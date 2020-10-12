#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>
#include <RedshiftLibrary/statistics/pdfz.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/log/log.h>

#include <RedshiftLibrary/processflow/context.h>


using namespace NSEpic;
using namespace std;
#include <fstream>

CPdfCandidateszResult::CPdfCandidateszResult()
{
    optMethod = 0; //di
    //optMethod = 1; //gaussian fit
    dzDefault = 1e-3;//default value in case deltaz couldnt be calculted, should be instrument dependant (parameter ?)
}

CPdfCandidateszResult::~CPdfCandidateszResult()
{

}

void CPdfCandidateszResult::Init(TRedshiftList const & zc, const TRedshiftList & deltaz, const TStringList & IDs)
{

    Resize(zc.size());

    Redshifts = zc;
    Deltaz = deltaz;
    if(!IDs.empty())
        ExtremaIDs = IDs;
    else
    {
        for(Int32 kc=0; kc<zc.size(); kc++)
        {
            ExtremaIDs[kc] = "Ext" + std::to_string(kc);
        }

    }
    Rank = TInt32List(zc.size(), -1);
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
Int32 CPdfCandidateszResult::SetIntegrationWindows(const TRedshiftList & Pdfz, TFloat64RangeList & ranges)
{
    Bool nodz = false;
    Int32 n = Redshifts.size();
    if(Deltaz.empty()){
        Log.LogInfo("    CPdfCandidateszResult::Compute pdf using Default full-window size, i.e., 6e-3 *(1+z)" );
        nodz = true;
    }

    TRedshiftList halfWidth;

    ranges.resize(n);
    for(Int32 i = 0; i< n; i++){

        //check cases where deltaz couldnt be computed or wasnt set--> use default value,
        if(Deltaz[i] == -1 || nodz)
            Deltaz[i] = dzDefault*(1+Redshifts[i]);

        halfWidth.push_back(3*Deltaz[i]);

        //initialize range boundaries for each candidate]
        ranges[i].Set((Redshifts[i] - halfWidth[i]), Redshifts[i] + halfWidth[i]);
        ranges[i].IntersectWith(TFloat64Range(Pdfz));
    };
    // sort zc values to facilitate comparison and keep track of initial order
    vector<pair<Float64,Int32 >> vp;
    vp.reserve(n);
    for (Int32 i = 0 ; i < n ; i++) {
        vp.push_back(make_pair(Redshifts[i], i));
    }
    std::stable_sort(vp.rbegin(), vp.rend(), [](const pair<Float64,Int32 >& a, const pair<Float64,Int32 >& b) { return a.first < b.first; });
    Int32 b = 0; 
    for(Int32 i = 0; i<n - 1; i++){ //i represents the highest candidate; go till n-1 since j increments i by one
        Int32 idx_h = vp[i].second; 
        Int32 idx_l = vp[i+1].second;
        Redshift overlap = ranges[idx_h].GetBegin() - ranges[idx_l].GetEnd();
        if(overlap < 0){
            b = 1;
            //in the case of duplicates, trim completely the range of the second cand
            if(vp[i].first !=  vp[i+1].first){
                Log.LogDebug("    CPdfCandidateszResult::Trimming: integration supports overlap for %f and %f", Redshifts[idx_h], Redshifts[idx_l] );
                ranges[idx_h].SetBegin(( std::max(Redshifts[idx_l], ranges[idx_h].GetBegin()) +
                                     std::min(Redshifts[idx_h], ranges[idx_l].GetEnd()) )/2);
            }
            ranges[idx_l].SetEnd(ranges[idx_h].GetBegin() - 1E-4);
         }
    }

    return b; //b is an indicator about overlapping presence
}
/**
 * @brief CPdfCandidateszResult::Compute
 */
Int32 CPdfCandidateszResult::Compute(TRedshiftList const & zc , TRedshiftList const & Pdfz, TFloat64List const & PdfProbalog, const TRedshiftList & deltaz, const TStringList & IDs)
{
    if(optMethod==0)
    {
        Log.LogInfo("    CPdfCandidateszResult::Compute pdf peaks info (method=direct integration)" );
    }else{
        Log.LogInfo("    CPdfCandidateszResult::Compute pdf peaks info (method=gaussian fitting)" );
    }

    Init(zc, deltaz, IDs);

    TFloat64RangeList ranges;
    Int32 b = SetIntegrationWindows( Pdfz, ranges);
    //b == 0 --> no overlapping, b == 1 --> overlapping
    CPdfz pdfz;
    for(Int32 kc=0; kc<zc.size(); kc++)
    {
        if(optMethod==0)
        {
            ValSumProba[kc] = pdfz.getCandidateSumTrapez( Pdfz, PdfProbalog, zc[kc], ranges[kc]);
            GaussAmp[kc]=-1;
            GaussAmpErr[kc]=-1;
            GaussSigma[kc]=-1;
            GaussSigmaErr[kc]=-1;
        }else
        {
            //TODO: this requires further check ?...
            Int32 retGaussFit = pdfz.getCandidateRobustGaussFit( Pdfz, PdfProbalog, zc[kc], ranges[kc], GaussAmp[kc], GaussAmpErr[kc], GaussSigma[kc], GaussSigmaErr[kc]);
            if(retGaussFit==0)
            {
                ValSumProba[kc] = GaussAmp[kc]*GaussSigma[kc]*sqrt(2*M_PI);
            }else{
                ValSumProba[kc] = -1.;
            }
        }
    }

    for (Int32 i = 0 ; i < Redshifts.size() ; i++)
    {
        Rank[i] = i;
    }
    
    SortByValSumProbaInt(Rank);//update only ranks based on valproba

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
    for(Int32 i=0; i<Rank.size(); i++)
    {
        Int32 k = Rank[i]; //use final rank for the output order
        stream << i << "\t"; 
        stream << ExtremaIDs[k] << "\t";
        stream << Redshifts[k] << "\t";
        stream << ValSumProba[k] << "\t";
        stream << i << "\t";
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
    for(Int32 i=0; i<Rank.size(); i++)
    {
        Int32 k = Rank[i];
        stream << i << "\t";
        stream << ExtremaIDs[k] << "\t";
        stream << Redshifts[k] << "\t";
        stream << ValSumProba[k] << "\t";
        stream << i << "\t"; 
        stream << Deltaz[k] << "\t";
        stream << GaussAmp[k] << "\t";
        stream << GaussSigma[k] << "\t"; 
    }
    stream << std::endl;
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

Bool CPdfCandidateszResult::GetBestRedshiftsFromPdf(const CDataStore& store, 
                                                    TFloat64List Extrema,  
                                                    std::vector<TFloat64List> ExtremaExtendedRedshifts,
                                                    TFloat64List& candidates ) const
{
        std::string scope_res = "zPDF/logposterior.logMargP_Z_data";
        auto results_pdf =  store.GetGlobalResult( scope_res.c_str() );
        auto logzpdf1d = std::dynamic_pointer_cast<const CPdfMargZLogResult>( results_pdf.lock() );

        if(!logzpdf1d)
        {
            Log.LogError( "GetBestRedshiftFromPdf: no pdf results retrieved from scope: %s", scope_res.c_str());
            return false;
        }
        Float64 Fullwidth = 6e-3;//should be replaced with deltaz?
        Int32 method = 0; //0=maxpdf, 1=direct integration on peaks only
        for( Int32 i=0; i<Extrema.size(); i++ )
        {
            Float64 tmpIntgProba = -DBL_MAX;
            Float64 tmpRedshift = Extrema[i]; //by default
            //TODO: below commented code should replace executing code once the new finder is ready; find() arguments shd also be updated
            TPointList extremumList;
            TFloat64Range redshiftsRange = TFloat64Range( ExtremaExtendedRedshifts[i]);
            //call Find on each secondpass range and retrieve the best 10 peaks?
            /*CExtremum extremum(redshiftsRange, 5, 0.005, false);
            extremum.DeactivateSlidingWindow();
            extremum.Find(logzpdf1d->Redshifts, logzpdf1d->valProbaLog, extremumList);
            if(method == 0){
                candidates.push_back(extremumList[0].X);
                continue;
            }
            for(Int32 j = 0; j < extremumList.size(); j++){
                CPdfz pdfz;
                Float64 flux_integral = -1;
                flux_integral = pdfz.getCandidateSumTrapez( logzpdf1d->Redshifts, logzpdf1d->valProbaLog, extremumList[j].X, Fullwidth);
                if(flux_integral>tmpIntgProba){
                        tmpRedshift = extremumList[j].X;
                        tmpIntgProba = flux_integral;
                }
            }*/
            Float64 tmpProbaLog = -DBL_MAX;
            for(Int32 kval=0; kval<ExtremaExtendedRedshifts[i].size(); kval++)
            {
                Float64 zInCandidateRange = ExtremaExtendedRedshifts[i][kval];
                UInt32 solIdx = logzpdf1d->getIndex(zInCandidateRange);
                if(solIdx<0 || solIdx>=logzpdf1d->valProbaLog.size())
                {
                    Log.LogError( "GetBestRedshiftFromPdf: pdf proba value not found for extremumIndex = %d", i);
                    return false;
                }

                Float64 probaLog = logzpdf1d->valProbaLog[solIdx];
                
                if(method == 0){
                    if(probaLog>tmpProbaLog){
                        tmpRedshift = zInCandidateRange;
                        tmpProbaLog = probaLog;
                    }
                }    
                if(method == 1){//max integrated proba but only on peaks in this range
                    CPdfz pdfz;
                    Float64 flux_integral = -1;
                    Float64 prev, next;
                    if(solIdx == 0){
                        prev = logzpdf1d->valProbaLog[solIdx];
                    }else{
                        prev = logzpdf1d->valProbaLog[solIdx-1];
                    }

                    if(solIdx == logzpdf1d->valProbaLog.size() - 1){
                        next = logzpdf1d->valProbaLog[solIdx];
                    }else{
                        next = logzpdf1d->valProbaLog[solIdx+1];
                    }
                    if((probaLog > prev&& probaLog > next) ||
                        (solIdx == 0 && probaLog>next) ||
                        (solIdx == logzpdf1d->valProbaLog.size()-1 && probaLog > prev)){
                        //if current value is a peak
                        Float64 Fullwidth_zp1 = Fullwidth*(zInCandidateRange + 1);
                        TFloat64Range zrange(zInCandidateRange-Fullwidth_zp1/2, zInCandidateRange+Fullwidth_zp1/2);
                        flux_integral = pdfz.getCandidateSumTrapez( logzpdf1d->Redshifts, logzpdf1d->valProbaLog, zInCandidateRange, zrange);
                        if(flux_integral>tmpIntgProba){
                            tmpRedshift = zInCandidateRange;
                            tmpIntgProba = flux_integral;
                        }
                    }
                    else 
                    {
                        continue; //it doesnt work to compute here the pdfz
                    }
                } 
            }

            candidates.push_back(tmpRedshift);
        }

    return true;
}


  void CPdfCandidateszResult::getCandidateData(const int& rank,const std::string& name, Float64& v) const
  {
    if (name.compare("Redshift") == 0) v=Redshifts[rank];
    else if (name.compare("RedshiftError") == 0) v=Deltaz[rank];
    else if (name.compare("RedshiftProba") == 0) v=ValSumProba[rank];
    else Log.LogError("unknown candidate data %s",name.c_str());
  }

  void CPdfCandidateszResult::getCandidateData(const int& rank,const std::string& name, Int32& v) const
  {
    if (name.compare("Rank") == 0) v=Rank[rank];
    else Log.LogError("unknown candidate data %s",name.c_str());
  }

  void CPdfCandidateszResult::getCandidateData(const int& rank,const std::string& name, std::string& v) const{
    v=ExtremaIDs[rank];
  }

  void CPdfCandidateszResult::getCandidateData(const int& rank,const std::string& name, double **data, int *size) const{}

  void CPdfCandidateszResult::getData(const std::string& name, Int32& v) const{
    if (name.compare("NbCandidates") == 0) v= Rank.size();
  
  }
  void CPdfCandidateszResult::getData(const std::string& name, Float64& v) const{}
  void CPdfCandidateszResult::getData(const std::string& name, std::string& v) const{}
  void CPdfCandidateszResult::getData(const std::string& name, double **data, int *size) const
  {

  }
