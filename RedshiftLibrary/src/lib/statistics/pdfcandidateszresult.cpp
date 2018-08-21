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
    Redshifts.resize(zc.size());
    ValSumProba.resize(zc.size());
    if(optMethod==1)
    {
        GaussAmp.resize(zc.size());
        GaussSigma.resize(zc.size());
        GaussAmpErr.resize(zc.size());
        GaussSigmaErr.resize(zc.size());
    }
    CPdfz pdfz;
    for(Int32 kc=0; kc<zc.size(); kc++)
    {
        Redshifts[kc] = zc[kc];
        if(optMethod==0)
        {
            ValSumProba[kc] = pdfz.getCandidateSumTrapez( Pdfz, PdfProbalog, zc[kc], Fullwidth);
        }else
        {
            Int32 retGaussFit = pdfz.getCandidateGaussFit( Pdfz, PdfProbalog, zc[kc], Fullwidth, GaussAmp[kc], GaussAmpErr[kc], GaussSigma[kc], GaussSigmaErr[kc]);
            if(retGaussFit==0)
            {
                ValSumProba[kc] = GaussAmp[kc]*GaussSigma[kc]*sqrt(2*M_PI);
            }
        }
    }

    return 0;
}


void CPdfCandidateszResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    stream  << "#fullwidth = " << Fullwidth << std::endl;
    stream  << "#method = " << optMethod << std::endl;
    stream  << std::endl;

    stream  << "#" << store.GetSpectrumName() << "\t" << store.GetProcessingID() << "\t";
    stream  << std::endl;

    stream  << "#" << "redshift" << "\t" << "intgProba";
    if(optMethod==1)
    {
        stream << "\t" << "gaussAmp" << "\t" << "gaussAmpErr" << "\t" << "gaussSigma" << "\t" << "gaussSigmaErr";
    }
    stream  << "\n";
    for(Int32 k=0; k<Redshifts.size(); k++)
    {
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
        stream << Redshifts[k] << "\t";
        stream << ValSumProba[k] << "\t";
        stream << GaussAmp[k] << "\t";
        stream << GaussSigma[k] << "\t";
    }
    stream << std::endl;
}
