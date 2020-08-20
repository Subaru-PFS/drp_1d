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

    Int32 Compute(TRedshiftList const & zc , TRedshiftList const & Pdfz, TFloat64List const & PdfProbalog,
                  const TRedshiftList & deltaz = TRedshiftList(), const TStringList & IDs= TStringList());

    void Init(TRedshiftList const & zc, const TRedshiftList & deltaz = TRedshiftList(), const TStringList & IDs= TStringList());

    Int32 SetIntegrationWindows(const TRedshiftList &Pdfz, TFloat64RangeList & ranges);

    Bool GetBestRedshiftsFromPdf(const CDataStore& store, 
                                TFloat64List Extrema,  
                                std::vector<TFloat64List> ExtremaExtendedRedshifts, 
                                TFloat64List& candidates) const;
  Float64 getDouble(std::string name,Int32 rank) const;
  std::string getString(std::string name,Int32 rank) const;
  Int32 getInt(std::string name,Int32 rank) const;
  Int32 getNbCandidates() const;
  
    Int32                       optMethod; //0: direct integration, 1:gaussian fit
    Float64                     dzDefault;
    std::vector<std::string> ExtremaIDs; //also sort ids

    TFloat64List           		Redshifts;
    TFloat64List           		ValSumProba;
    TInt32List                  Rank;
    TFloat64List                Deltaz;
    //opt 1: direct integration
    //
    //opt 2: gaussian fit
    TFloat64List           		GaussAmp;
    TFloat64List           		GaussAmpErr;
    TFloat64List           		GaussSigma;
    TFloat64List           		GaussSigmaErr;
    //TFloat64List           		GaussSkewness; //todo !

private:

    void SortByValSumProbaInt(TInt32List& flist);
};


}

#endif
