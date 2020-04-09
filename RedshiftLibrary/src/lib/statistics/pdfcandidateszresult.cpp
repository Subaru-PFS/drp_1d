#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>
#include <RedshiftLibrary/statistics/pdfz.h>

#include <RedshiftLibrary/log/log.h>

#include <RedshiftLibrary/processflow/context.h>


using namespace NSEpic;
using namespace std;
#include <fstream>

CPdfCandidateszResult::CPdfCandidateszResult()
{
    optMethod = 0; //di
    //optMethod = 1; //gaussian fit
    dzDefault = 1e-3;
    Fullwidth = 6*dzDefault;//default value in case deltaz couldnt be calculted
    
}

CPdfCandidateszResult::~CPdfCandidateszResult()
{

}

void CPdfCandidateszResult::Resize(Int32 n)
{
    Redshifts.resize(n);
    ValSumProba.resize(n);
    Rank.resize(n);
    ExtremaIDs.resize(n);
    Deltaz.resize(n);
    //only for method 1
    GaussAmp.resize(n);
    GaussSigma.resize(n);
    GaussAmpErr.resize(n);
    GaussSigmaErr.resize(n);
}
/**
 * 1) Fix Deltaz problems: none is passed or none could be compute --> use default values
 * 2) Check if integration windows overlap, mostly for close candidates
 *      1) if no, keep deltaz value
 *      2) if overlapping, update the half-width of the left and right sides of the integration windows
 * Note: Output range includes the multiplication by (1+z).
 * Returns 0 if no overlapping; otherwise 1;
*/
Int32 CPdfCandidateszResult::SetIntegrationWindows( std::vector<Float64> redshifts,  std::vector<Float64> deltaz, std::vector<Float64>& range_right, std::vector<Float64>& range_left){ 
    Bool nodz = false;
    Int32 n = redshifts.size();
    if(!deltaz.size()){
        Log.LogInfo("    CPdfCandidateszResult::Compute pdf using Default full-window size, i.e., 6e-3" );      
        nodz = true;
    }
    std::vector<Float64> halfWidth;
    //check cases where deltaz couldnt be computed or wasnt set--> use default value, 
    for(Int32 i = 0; i< n; i++){
        if(deltaz[i] == -1 || nodz) 
            deltaz[i] = dzDefault;
        Deltaz[i] = deltaz[i]; 
        halfWidth.push_back(3*deltaz[i]*(1 + redshifts[i]));     
        //initialize range boundaries for each candidate
        range_right.push_back(redshifts[i] + halfWidth[i]);//higher boundary
        range_left.push_back(redshifts[i] - halfWidth[i]);//lower boundary
    };
    // sort zc values to facilitate comparison and keep track of initial order
    vector<pair<Float64,Int32 >> vp;
    vp.reserve(n);
    for (Int32 i = 0 ; i < n ; i++) {
        vp.push_back(make_pair(redshifts[i], i));
    }
    std::sort(vp.rbegin(), vp.rend());  //sort from the highest to the lowest!

    Int32 b = 0; 
    for(Int32 i = 0; i<n - 1; i++){ //i represents the higher candidate; go till n-1 since j increments i by one
        Int32 idx_h = vp[i].second; 
        Int32 j = i + 1; 
        Int32 idx_l = vp[j].second;
        Float64 overlap =  range_left[idx_h] - range_right[idx_l];
        if(overlap < 0 ){
            b = 1;
            Log.LogDebug("    CPdfCandidateszResult::Trimming: integration supports overlap for %f and %f", redshifts[idx_h], redshifts[idx_l] ); 
            range_left[idx_h] = ( std::max(redshifts[idx_l], range_left[idx_h]) +
                                  std::min(redshifts[idx_h], range_right[idx_l]) )/2;
            range_right[idx_l] = range_left[idx_h] - 1E-4;
         }
    }

    return b; //b is an indicator about overlapping presence
}
/**
 * @brief CPdfCandidateszResult::Compute
 */
Int32 CPdfCandidateszResult::Compute( std::vector<Float64> zc,  std::vector<Float64> Pdfz,  std::vector<Float64> PdfProbalog, std::vector<Float64> deltaz, std::vector<std::string> IDs)
{
    if(optMethod==0)
    {
        Log.LogInfo("    CPdfCandidateszResult::Compute pdf peaks info (method=direct integration)" );
    }else{
        Log.LogInfo("    CPdfCandidateszResult::Compute pdf peaks info (method=gaussian fitting)" );
    }
    Resize(zc.size());
 
    std::vector<Float64> range_right, range_left; 
    Int32 b = SetIntegrationWindows(zc, deltaz, range_right, range_left);
    //b == 0 --> no overlapping, b == 1 --> overlapping
    CPdfz pdfz;
    for(Int32 kc=0; kc<zc.size(); kc++)
    {
        Rank[kc] = -1;
        Redshifts[kc] = zc[kc];
        if(IDs.size()){
             ExtremaIDs[kc] = IDs[kc];
        }else{
            //generate IDs 
            ExtremaIDs[kc] = "Ext" + std::to_string(kc);
        }

        if(optMethod==0)
        {
            ValSumProba[kc] = pdfz.getCandidateSumTrapez( Pdfz, PdfProbalog, zc[kc], range_left[kc], range_right[kc]);
            GaussAmp[kc]=-1;
            GaussAmpErr[kc]=-1;
            GaussSigma[kc]=-1;
            GaussSigmaErr[kc]=-1;
        }else
        {
            //TODO: this requires further check ?...
            Int32 retGaussFit = pdfz.getCandidateRobustGaussFit( Pdfz, PdfProbalog, zc[kc], (range_left[kc] + range_right[kc]), GaussAmp[kc], GaussAmpErr[kc], GaussSigma[kc], GaussSigmaErr[kc]);
            if(retGaussFit==0)
            {
                ValSumProba[kc] = GaussAmp[kc]*GaussSigma[kc]*sqrt(2*M_PI);
            }else{
                ValSumProba[kc] = -1.;
            }
        }
    }

    SortByRank();
    return 0;
}


void CPdfCandidateszResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    stream  << "#fullwidth = " << "6* Deltaz" << std::endl;
    stream  << "#method = " << optMethod << std::endl;
    stream  << std::endl;

    stream  << "#" << store.GetSpectrumName() << "\t" << store.GetProcessingID() << "\t";
    stream  << std::endl;

    stream  << "#" << "rank" << "\t"  << "IDs" << "\t"<< "redshift" << "\t" << "intgProba"<< "\t" << "Rank_PDF" << "\t" <<"Deltaz";
    if(optMethod==1)
    {
        stream << "\t" << "gaussAmp" << "\t" << "gaussAmpErr" << "\t" << "gaussSigma" << "\t" << "gaussSigmaErr";
    }else{
        stream << "\t" << "gaussAmp_unused" << "\t" << "gaussAmpErr_unused" << "\t" << "gaussSigma_unused" << "\t" << "gaussSigmaErr_unused";
    }
    stream  << "\n";
    for(Int32 k=0; k<Redshifts.size(); k++)
    {
        stream << k << "\t"; 
        stream << ExtremaIDs[k] << "\t";
        stream << Redshifts[k] << "\t";
        stream << ValSumProba[k] << "\t";
        stream << Rank[k] << "\t";
        stream << Deltaz[k] << "\t";
        //only for method 1, but leave columns with -1 value ste in compute()
        stream << GaussAmp[k] << "\t";
        stream << GaussAmpErr[k] << "\t";
        stream << GaussSigma[k] << "\t";
        stream << GaussSigmaErr[k] << "\t";

        stream << "\n";
    }
    stream << std::endl;
}

void CPdfCandidateszResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    stream  << store.GetSpectrumName() << "\t" << store.GetProcessingID() << "\t";
    for(Int32 k=0; k<Redshifts.size(); k++)
    {
        stream << k << "\t";
        stream << ExtremaIDs[k] << "\t";
        stream << Redshifts[k] << "\t";
        stream << ValSumProba[k] << "\t";
        stream << Rank[k] << "\t";
        stream << Deltaz[k] << "\t";
        stream << GaussAmp[k] << "\t";
        stream << GaussSigma[k] << "\t"; 
    }
    stream << std::endl;
}

void CPdfCandidateszResult::SortByRank()
{
    for (Int32 i = 0 ; i < Redshifts.size() ; i++)
    {
        Rank[i] = i;
    }
    SortByValSumProbaInt(Rank);//update ranks based on valproba
    SortIDsByValSumProba(ExtremaIDs);//update ranks based on valproba
    SortByValSumProba(Redshifts);
    for (Int32 i = 0; i <Rank.size(); i++){
        if(Rank[i]!=i){
            Log.LogDebug("Zcand %f has his rank updated from %d to %d \n", Redshifts[i], i, Rank[i]);
        }
    }
    SortByValSumProba(ValSumProba);
    SortByValSumProba(Deltaz);
    if(optMethod==1)
    {
        SortByValSumProba(GaussAmp);
        SortByValSumProba(GaussAmpErr);
        SortByValSumProba(GaussSigma);
        SortByValSumProba(GaussSigmaErr);
    }  
}

void CPdfCandidateszResult::SortByValSumProba(TFloat64List& flist)
{
    //sort the valProbaSum and reorder flist accordingly
    TFloat64List sortedProba;
    TFloat64List sortedFlist;

    // This is a vector of {value,index} pairs
    vector<pair<Float64,Float64> > vp;
    vp.reserve(Redshifts.size());
    for (Int32 i = 0 ; i < Redshifts.size() ; i++) {
        vp.push_back(make_pair(ValSumProba[i], flist[i]));
    }
    std::sort(vp.rbegin(), vp.rend()); //sort reverse order
    for (Int32 i = 0 ; i < vp.size() ; i++) {
        sortedProba.push_back(vp[i].first);
        sortedFlist.push_back(vp[i].second);
    }

    for (Int32 i = 0 ; i < Redshifts.size() ; i++) {
        flist[i] = sortedFlist[i];
    }
}

void CPdfCandidateszResult::SortByValSumProbaInt(TInt32List& flist)
{
    //sort the valProbaSum and reorder flist accordingly
    TFloat64List sortedProba;
    TInt32List sortedFlist;

    // This is a vector of {value,index} pairs
    vector<pair<Float64,Int32> > vp;
    vp.reserve(Redshifts.size());
    for (Int32 i = 0 ; i < Redshifts.size() ; i++) {
        vp.push_back(make_pair(ValSumProba[i], flist[i]));
    }
    std::sort(vp.rbegin(), vp.rend()); //sort reverse order
    for (Int32 i = 0 ; i < vp.size() ; i++) {
        sortedProba.push_back(vp[i].first);
        sortedFlist.push_back(vp[i].second);
    }

    for (Int32 i = 0 ; i < Redshifts.size() ; i++) {
        flist[i] = sortedFlist[i];
    }
}

void CPdfCandidateszResult::SortIDsByValSumProba(std::vector<std::string>& flist)
{
    //sort the valProbaSum and reorder flist accordingly
    TFloat64List sortedProba;
    std::vector<std::string> sortedFlist;

    // This is a vector of {value,index} pairs
    vector<pair<Float64,std::string> > vp;
    vp.reserve(Redshifts.size());
    for (Int32 i = 0 ; i < Redshifts.size() ; i++) {
        vp.push_back(make_pair(ValSumProba[i], flist[i]));
    }
    std::sort(vp.rbegin(), vp.rend()); //sort reverse order
    for (Int32 i = 0 ; i < vp.size() ; i++) {
        sortedProba.push_back(vp[i].first);
        sortedFlist.push_back(vp[i].second);
    }

    for (Int32 i = 0 ; i < Redshifts.size() ; i++) {
        flist[i] = sortedFlist[i];
    }
}
