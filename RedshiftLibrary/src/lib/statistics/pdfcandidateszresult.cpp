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
    Fullwidth = 6e-3;
}

CPdfCandidateszResult::~CPdfCandidateszResult()
{

}

void CPdfCandidateszResult::Resize(Int32 n)
{
    Redshifts.resize(n);
    ValSumProba.resize(n);
    Rank.resize(n);
    ExtremaIDs.resize(n) ;
    //only for method 1
    GaussAmp.resize(n);
    GaussSigma.resize(n);
    GaussAmpErr.resize(n);
    GaussSigmaErr.resize(n);
}

/**
 * @brief CPdfCandidateszResult::Compute
 */
Int32 CPdfCandidateszResult::Compute( std::vector<Float64> zc,  std::vector<Float64> Pdfz,  std::vector<Float64> PdfProbalog, std::vector<std::string> IDs)
{
    if(optMethod==0)
    {
        Log.LogInfo("    CPdfCandidateszResult::Compute pdf peaks info (method=direct integration)" );
    }else{
        Log.LogInfo("    CPdfCandidateszResult::Compute pdf peaks info (method=gaussian fitting)" );
    }
    Resize(zc.size());

    CPdfz pdfz;
    for(Int32 kc=0; kc<zc.size(); kc++)
    {
        Rank[kc] = -1;
        Redshifts[kc] = zc[kc];
        ExtremaIDs[kc] = IDs[kc];
        if(optMethod==0)
        {
            ValSumProba[kc] = pdfz.getCandidateSumTrapez( Pdfz, PdfProbalog, zc[kc], Fullwidth*(1+zc[kc]));
            GaussAmp[kc]=-1;
            GaussAmpErr[kc]=-1;
            GaussSigma[kc]=-1;
            GaussSigmaErr[kc]=-1;
        }else
        {
            Int32 retGaussFit = pdfz.getCandidateRobustGaussFit( Pdfz, PdfProbalog, zc[kc], Fullwidth*(1+zc[kc]), GaussAmp[kc], GaussAmpErr[kc], GaussSigma[kc], GaussSigmaErr[kc]);
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

void CPdfCandidateszResult::SetFullWidth(Float64 width)
{
    Fullwidth = width;
}

void CPdfCandidateszResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    stream  << "#fullwidth = " << Fullwidth << std::endl;
    stream  << "#method = " << optMethod << std::endl;
    stream  << std::endl;

    stream  << "#" << store.GetSpectrumName() << "\t" << store.GetProcessingID() << "\t";
    stream  << std::endl;

    stream  << "#" << "rank" << "\t"  << "IDs" << "\t"<< "redshift" << "\t" << "intgProba"<< "\t" << "Rank_PDF";
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
            Log.LogDebug("Zcand %f has his rank updated from %d to %f \n", Redshifts[i], i, Rank[i]);
        }
    }
    SortByValSumProba(ValSumProba);
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
