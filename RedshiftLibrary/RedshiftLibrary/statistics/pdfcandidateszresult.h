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

    void Resize(Int32 n);
    Int32 Compute(std::vector<Float64> zc , std::vector<Float64> Pdfz, std::vector<Float64> PdfProbalog, std::vector<std::string> IDs);
    void SetFullWidth(Float64 width);


    Int32                       optMethod; //0: direct integration, 1:gaussian fit
    Float64                     Fullwidth;
    std::vector<std::string> ExtremaIDs; //also sort ids

    TFloat64List           		Redshifts;
    TFloat64List           		ValSumProba;
    TFloat64List                Rank;
    //opt 1: direct integration
    //
    //opt 2: gaussian fit
    TFloat64List           		GaussAmp;
    TFloat64List           		GaussAmpErr;
    TFloat64List           		GaussSigma;
    TFloat64List           		GaussSigmaErr;
    //TFloat64List           		GaussSkewness; //todo !

private:
    void SortByRank();
    void SortByValSumProba(TFloat64List &flist);
    //to sort IDs
    void SortIDsByValSumProba(std::vector<std::string>& flist);

};


}

#endif
