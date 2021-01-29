#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>
#include <RedshiftLibrary/operator/pdfz.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/log/log.h>

#include <RedshiftLibrary/processflow/context.h>


using namespace NSEpic;
using namespace std;
#include <fstream>
#include <iostream>
#include <numeric>

CPdfCandidateszResult::CPdfCandidateszResult():
    optMethod(0),   // di
    //optMethod(1), // gaussian fit
    dzDefault(1e-3) // default value in case deltaz couldnt be calculted, should be instrument dependant (parameter ?)
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
 * Returns a list of identified very close candidates, at 2*1E-4
*/
TInt32List CPdfCandidateszResult::SetIntegrationWindows(const TRedshiftList & Pdfz, TFloat64RangeList & ranges)
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
    TInt32List b = {}; 
    for(Int32 i = 0; i<n - 1; i++){ //i represents the highest candidate; go till n-1 since j increments i by one
        Int32 idx_h = vp[i].second; 
        Int32 idx_l = vp[i+1].second;
        Redshift overlap = ranges[idx_h].GetBegin() - ranges[idx_l].GetEnd();
        if(overlap < 0){
            //in the case of duplicates, trim completely the range of the second cand
            if((vp[i].first -  vp[i+1].first)>2*1E-4){
                Log.LogDebug("    CPdfCandidateszResult::SetIntegrationWindows: integration supports overlap for %f and %f", Redshifts[idx_h], Redshifts[idx_l] );
                ranges[idx_h].SetBegin(( std::max(Redshifts[idx_l], ranges[idx_h].GetBegin()) +
                                     std::min(Redshifts[idx_h], ranges[idx_l].GetEnd()) )/2);
            }else{
                Log.LogInfo(" CPdfCandidateszResult::SetIntegrationWindows: very close candidates are identified %f and %f", vp[i].first,  vp[i+1].first);
                b.push_back(idx_l);
            }
            ranges[idx_l].SetEnd(ranges[idx_h].GetBegin() - 1E-4);
         }
    }

    //iterate over computed ranges and check that corresponding zcandidates belong to that range, otherwise throw error
    for(Int32 i = 0; i<n; i++){
        //if currend index belongs to the duplicates b vector, skip testing it and only test the others
        if(std::find(b.begin(), b.end(), i)!= b.end())
            continue;
        if( Redshifts[i]>= ranges[i].GetBegin() && 
            Redshifts[i]<= ranges[i].GetEnd()){
            continue;
        }else{
            Log.LogError("CPdfCandidateszResult::SetIntegrationWindows: Failed to identify a range including the candidate %f", Redshifts[i]);
            throw runtime_error("CPdfCandidateszResult::SetIntegrationWindows: Failed to identify a range including the candidate! Aborting");
        }
    } 
    return b; 
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
    TInt32List duplicates = SetIntegrationWindows( Pdfz, ranges);
    COperatorPdfz pdfz;
    for(Int32 kc=0; kc<zc.size(); kc++)
    {
        if(optMethod==0)
        {
            //check if current candidate belongs to the identified duplicates list
            //if yest, force its pdf value to 0 and avoid callling getCandidateSumTrapez
            if(std::find(duplicates.begin(), duplicates.end(),kc)!=duplicates.end())
                ValSumProba[kc] = 0;
            else
                ValSumProba[kc] = pdfz.getCandidateSumTrapez( Pdfz, PdfProbalog, zc[kc], ranges[kc]);
            GaussAmp[kc]=-1;
            GaussAmpErr[kc]=-1;
            GaussSigma[kc]=-1;
            GaussSigmaErr[kc]=-1;
        }else
        {
            //TODO: this requires further check ?...
            if(std::find(duplicates.begin(), duplicates.end(),kc)!=duplicates.end()){
                ValSumProba[kc] = 0;
                continue;
            }            
            Int32 retGaussFit = pdfz.getCandidateRobustGaussFit( Pdfz, PdfProbalog, zc[kc], ranges[kc], GaussAmp[kc], GaussAmpErr[kc], GaussSigma[kc], GaussSigmaErr[kc]);
            if(retGaussFit==0)
            {
                ValSumProba[kc] = GaussAmp[kc]*GaussSigma[kc]*sqrt(2*M_PI);
            }else{
                ValSumProba[kc] = -1.;
            }
        }
    }
    
    SortByValSumProbaInt();//update only ranks based on valproba

    return 0;
}


void CPdfCandidateszResult::Save( std::ostream& stream ) const
{
    stream  << "#fullwidth = " << "6* Deltaz" << std::endl;
    stream  << "#method = " << optMethod << std::endl;
    stream  << std::endl;

    //    stream  << "#" << store.GetSpectrumName() << "\t" << store.GetProcessingID() << "\t";
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

void CPdfCandidateszResult::SaveLine( std::ostream& stream ) const
{
  //    stream  << store.GetSpectrumName() << "\t" << store.GetProcessingID() << "\t";
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


void CPdfCandidateszResult::SortByValSumProbaInt()
{
    iota(Rank.begin(), Rank.end(), 0);
    
    const TFloat64List & v = ValSumProba;
    sort(Rank.begin(), Rank.end(), 
        [&v](Int32 i1, Int32 i2){return v[i1] > v[i2];});
}


  void CPdfCandidateszResult::getCandidateData(const int& rank,const std::string& name, Float64& v) const
  {
    if (name.compare("Redshift") == 0) v=Redshifts[Rank[rank]];
    else if (name.compare("RedshiftError") == 0) v=Deltaz[Rank[rank]];
    else if (name.compare("RedshiftProba") == 0) v=ValSumProba[Rank[rank]];
    else Log.LogError("unknown candidate data %s",name.c_str());
  }

  void CPdfCandidateszResult::getCandidateData(const int& rank,const std::string& name, Int32& v) const
  {
    if (name.compare("Rank") == 0) v=rank;
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
