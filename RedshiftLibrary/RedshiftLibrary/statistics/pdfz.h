#ifndef _REDSHIFT_STATISTICS_PDFZ_
#define _REDSHIFT_STATISTICS_PDFZ_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>

#include <RedshiftLibrary/operator/pdfMargZLogResult.h>

#include <string>

namespace NSEpic
{

/**
 * \ingroup Redshift
 * Pdfz
 */
class CPdfz
{

public:

    CPdfz();
    ~CPdfz();

    Int32 Compute(const TFloat64List &merits, const TFloat64List &redshifts, const Float64 cstLog, const TFloat64List &zPrior, TFloat64List &logPdf, Float64 &logEvidence);
    TFloat64List GetConstantLogZPrior(UInt32 nredshifts);
    TFloat64List GetStrongLinePresenceLogZPrior(const TBoolList & linePresence, const Float64 penalization_factor);
    TFloat64List GetNLinesSNRAboveCutLogZPrior(const TInt32List & nlinesAboveSNR, const Float64 penalization_factor);
    TFloat64List GetEuclidNhaLogZPrior(const TFloat64List & redshifts, const Float64 aCoeff);
    TFloat64List CombineLogZPrior(const TFloat64List & logprior1, const TFloat64List & logprior2);

    Float64 getSumTrapez(const TRedshiftList &redshifts, const TFloat64List &valprobalog);
    Float64 getSumRect(const TRedshiftList &redshifts, const TFloat64List &valprobalog);
    Float64 getCandidateSumTrapez(const TRedshiftList &redshifts,
                                  const TFloat64List &valprobalog,
                                  const Redshift zcandidate,
                                  const TFloat64Range &zrange);//default: zwidth_left = zwidth_right

    Int32   getCandidateRobustGaussFit(const TRedshiftList &redshifts,
                                       const TFloat64List &valprobalog,
                                       const Float64 zcandidate,
                                       const TFloat64Range &zrange,
                                       Float64 &gaussAmp, Float64 &gaussAmpErr,
                                       Float64 &gaussSigma, Float64 &gaussSigmaErr);

    Int32   getPmis(const TRedshiftList &redshifts,
                    const TFloat64List &valprobalog,
                    const Float64 zbest,
                    TRedshiftList &zcandidates,
                    const Float64 zwidth,
                    Float64 &pmis);

    Int32 Marginalize(const TFloat64List &redshifts,
                      const std::vector<TFloat64List> &meritResults,
                      const std::vector<TFloat64List> &zPriors,
                      const Float64 cstLog,
                      std::shared_ptr<CPdfMargZLogResult> postmargZResult,
                      const TFloat64List &modelPriors=TFloat64List());
    Int32 BestProba(const TFloat64List &redshifts, const std::vector<TFloat64List> &meritResults, const std::vector<TFloat64List> &zPriors, const Float64 cstLog, std::shared_ptr<CPdfMargZLogResult> postmargZResult);
    Int32 BestChi2(const TFloat64List &redshifts, const std::vector<TFloat64List> &meritResults, const std::vector<TFloat64List> &zPriors, const Float64 cstLog, std::shared_ptr<CPdfMargZLogResult> postmargZResult);
    Int32 getIndex( std::vector<Float64> redshifts, Float64 z );
    Float64 logSumExpTrick(TFloat64List valproba, TFloat64List redshifts, Int32 sumMethod);
private:

    Int32   getCandidateGaussFit(const TRedshiftList &redshifts,
                                         const TFloat64List &valprobalog,
                                         const Float64 zcandidate,
                                         const TFloat64Range &zrange,
                                         Float64 &gaussAmp, Float64 &gaussAmpErr,
                                         Float64 &gaussSigma, Float64 &gaussSigmaErr);

};


}

#endif
