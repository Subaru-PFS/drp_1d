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

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }

    Int32 Compute(std::vector<Float64> zc , std::vector<Float64> Pdfz, std::vector<Float64> PdfProbalog);

    TFloat64List           		Redshifts;
    Int32                       optMethod; //0: direct integration, 1:gaussian fit
    TFloat64List           		ValSumProba;
    Float64                     Fullwidth;
    //opt 1: direct integration
    //
    //opt 2: gaussian fit
    TFloat64List           		GaussAmp;
    TFloat64List           		GaussAmpErr;
    TFloat64List           		GaussSigma;
    TFloat64List           		GaussSigmaErr;
    //TFloat64List           		GaussSkewness;

};


}

#endif
