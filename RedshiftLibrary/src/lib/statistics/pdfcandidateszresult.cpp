#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>
#include <RedshiftLibrary/statistics/pdfz.h>

#include <RedshiftLibrary/log/log.h>

#include <RedshiftLibrary/processflow/context.h>


using namespace NSEpic;
using namespace std;
#include <fstream>

CPdfCandidateszResult::CPdfCandidateszResult()
{
    optMethod = 1; //gaussian fit
    Fullwidth = 1e-2;
}

CPdfCandidateszResult::~CPdfCandidateszResult()
{

}

void CPdfCandidateszResult::Resize(Int32 n)
{
    Redshifts.resize(n);
    ValSumProba.resize(n);
    Rank.resize(n);

    if(optMethod==1)
    {
        GaussAmp.resize(n);
        GaussSigma.resize(n);
        GaussAmpErr.resize(n);
        GaussSigmaErr.resize(n);
    }
}

/**
 * @brief CPdfCandidateszResult::Compute
 */
Int32 CPdfCandidateszResult::Compute( std::vector<Float64> zc,  std::vector<Float64> Pdfz,  std::vector<Float64> PdfProbalog )
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
        if(optMethod==0)
        {
            ValSumProba[kc] = pdfz.getCandidateSumTrapez( Pdfz, PdfProbalog, zc[kc], Fullwidth);
        }else
        {
            Int32 retGaussFit = pdfz.getCandidateRobustGaussFit( Pdfz, PdfProbalog, zc[kc], Fullwidth, GaussAmp[kc], GaussAmpErr[kc], GaussSigma[kc], GaussSigmaErr[kc]);
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

    stream  << "#" << "rank" << "\t" << "redshift" << "\t" << "intgProba";
    if(optMethod==1)
    {
        stream << "\t" << "gaussAmp" << "\t" << "gaussAmpErr" << "\t" << "gaussSigma" << "\t" << "gaussSigmaErr";
    }
    stream  << "\n";
    for(Int32 k=0; k<Redshifts.size(); k++)
    {
        stream << Rank[k] << "\t";
        stream << Redshifts[k] << "\t";
        stream << ValSumProba[k] << "\t";
        if(optMethod==1)
        {
            stream << GaussAmp[k] << "\t";
            stream << GaussAmpErr[k] << "\t";
            stream << GaussSigma[k] << "\t";
            stream << GaussSigmaErr[k] << "\t";
        }
        stream << "\n";
    }
    stream << std::endl;
}

void CPdfCandidateszResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    stream  << store.GetSpectrumName() << "\t" << store.GetProcessingID() << "\t";

    for(Int32 k=0; k<Redshifts.size(); k++)
    {
        stream << Rank[k] << "\t";
        stream << Redshifts[k] << "\t";
        stream << ValSumProba[k] << "\t";
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
    SortByValSumProba(Redshifts);
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
