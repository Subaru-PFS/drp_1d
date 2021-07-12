#ifndef _REDSHIFT_OPERATOR_PDFZ_
#define _REDSHIFT_OPERATOR_PDFZ_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/operator/operator.h"

#include "RedshiftLibrary/operator/pdfMargZLogResult.h"
#include "RedshiftLibrary/statistics/pdfcandidateszresult.h"
#include "RedshiftLibrary/statistics/pdfcandidatesz.h"

#include <string>

namespace NSEpic
{

struct ChisquareArray {
    TFloat64List redshifts;
    std::vector<TFloat64List> chisquares;
    std::vector<TFloat64List> zpriors;
    TFloat64List modelpriors;
    Float64 cstLog=0.;
};   


/**
 * \ingroup Redshift
 * Pdfz
 */
class COperatorPdfz : public COperator
{

public:    
    COperatorPdfz(const std::string & opt_combine,
                  Float64 peakSeparation = 0.0, // no minimal seaparation
                  Float64 meritcut = 0.0,  // no cut
                  UInt32 maxCandidate=10, // max number of candidate at the end
                  const std::string & Id_prefix="EXT",
                  Bool allow_extrema_at_border=true,
                  UInt32 maxPeakCount_per_window=0,  // <=0 will be set to maxCandidate (default to one window)
                  const std::vector<TFloat64List> & candidatesRedshifts = std::vector<TFloat64List>(1) ,
                  const TStringList & candidatesIds = TStringList(1)
                );

    
  std::shared_ptr<CPdfCandidateszResult<TCandidateZ>> Compute(const ChisquareArray & chisquares, Bool integ=true);

    Int32 CombinePDF(const ChisquareArray & chisquares);



    Bool checkPdfSum();

    std::shared_ptr<CPdfMargZLogResult> m_postmargZResult;

    // static member function to do calculation on pdfs
    ///////////////////////////////////////////////////

    static Int32 getIndex( const std::vector<Float64> & redshifts, Float64 z );

    static Int32 ComputePdf(const TFloat64List &merits, const TFloat64List &redshifts, const Float64 cstLog, const TFloat64List &zPrior, TFloat64List &logPdf, Float64 &logEvidence);

    static Float64 getSumTrapez(const TRedshiftList &redshifts, const TFloat64List &valprobalog);
 
    static Float64 getSumRect(const TRedshiftList &redshifts, const TFloat64List &valprobalog);
     
    static Float64 logSumExpTrick(const TFloat64List & valproba, const TFloat64List & redshifts, Int32 sumMethod);


private:

    Int32 Marginalize(const ChisquareArray & chisquarearray);
    Int32 BestProba(const ChisquareArray & chisquarearray);
    Int32 BestChi2(const ChisquareArray & chisquarearray);

    TCandidateZbyID searchMaxPDFcandidates() const;

    const std::string m_opt_combine;
    TCandidateZRangebyID m_candidatesZRanges;
    const UInt32 m_maxPeakCount_per_window;
    const UInt32 m_maxCandidate;
    const Float64  m_peakSeparation;
    const Bool m_allow_extrema_at_border;
    const Float64  m_meritcut;
    const std::string m_Id_prefix;// =  "EXT"; // for "extrema"

};


}

#endif
