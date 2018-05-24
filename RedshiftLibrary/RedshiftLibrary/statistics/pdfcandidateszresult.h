#ifndef _REDSHIFT_PDF_CANDIDATESZRESULT_
#define _REDSHIFT_PDF_CANDIDATESZRESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>

namespace NSEpic
{


class CPdfCandidateszResult : public COperatorResult
{

public:

    CPdfCandidateszResult();
    virtual ~CPdfCandidateszResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;

    Int32 Compute(std::vector<Float64> zc , std::vector<Float64> Pdfz, std::vector<Float64> PdfProbalog);

    TFloat64List           		Redshifts;
    TFloat64List           		ValSumProba;
    Float64                     Fullwidth;
};


}

#endif
