#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>
#include <RedshiftLibrary/statistics/pdfz.h>

#include <RedshiftLibrary/log/log.h>

#include <RedshiftLibrary/processflow/context.h>


using namespace NSEpic;
using namespace std;
#include <fstream>

CPdfCandidateszResult::CPdfCandidateszResult()
{
    Fullwidth = 1e-3;
}

CPdfCandidateszResult::~CPdfCandidateszResult()
{

}

/**
 * @brief CPdfCandidateszResult::Compute
 */
Int32 CPdfCandidateszResult::Compute( std::vector<Float64> zc,  std::vector<Float64> Pdfz,  std::vector<Float64> PdfProbalog )
{
    Redshifts.resize(zc.size());
    ValSumProba.resize(zc.size());
    CPdfz pdfz;
    for(Int32 kc=0; kc<zc.size(); kc++)
    {
        Redshifts[kc] = zc[kc];
        ValSumProba[kc] = pdfz.getCandidateSumTrapez( Pdfz, PdfProbalog, zc[kc], Fullwidth);
    }

    return 0;
}


Void CPdfCandidateszResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    stream  << "#fullwidth = " << Fullwidth << std::endl;

    stream  << store.GetSpectrumName() << "\t" << store.GetProcessingID() << "\t";

    for(Int32 k=0; k<Redshifts.size(); k++)
    {
        stream << Redshifts[k] << "\t";
        stream << ValSumProba[k] << "\t";
    }
    stream << std::endl;
}

Void CPdfCandidateszResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    stream  << store.GetSpectrumName() << "\t" << store.GetProcessingID() << "\t";

    for(Int32 k=0; k<Redshifts.size(); k++)
    {
        stream << Redshifts[k] << "\t";
        stream << ValSumProba[k] << "\t";
    }
    stream << std::endl;
}
