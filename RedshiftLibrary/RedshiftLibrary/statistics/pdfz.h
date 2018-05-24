#ifndef _REDSHIFT_RAY_PDFZ_
#define _REDSHIFT_RAY_PDFZ_

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

    Int32 Compute(TFloat64List merits, TFloat64List redshifts, Float64 cstLog, TFloat64List zPrior, TFloat64List &logPdf, Float64 &logEvidence);
    std::vector<Float64> GetConstantLogZPrior(UInt32 nredshifts);
    std::vector<Float64> GetStrongLinePresenceLogZPrior(std::vector<bool> linePresence);
    std::vector<Float64> GetEuclidNhaLogZPrior(std::vector<Float64> redshifts);
    std::vector<Float64> CombineLogZPrior(std::vector<Float64> logprior1, std::vector<Float64> logprior2);

    Float64 getSumTrapez(std::vector<Float64> redshifts, std::vector<Float64> valprobalog);
    Float64 getSumRect(std::vector<Float64> redshifts, std::vector<Float64> valprobalog);
    Float64 getCandidateSumTrapez(std::vector<Float64> redshifts, std::vector<Float64> valprobalog, Float64 zcandidate, Float64 zwidth);



    Int32 Marginalize(TFloat64List redshifts,
                      std::vector<TFloat64List> meritResults,
                      std::vector<TFloat64List> zPriors,
                      Float64 cstLog,
                      std::shared_ptr<CPdfMargZLogResult> postmargZResult,
                      std::vector<Float64> modelPriors=std::vector<Float64>());
    Int32 BestProba(TFloat64List redshifts, std::vector<TFloat64List> meritResults, std::vector<TFloat64List> zPriors, Float64 cstLog, std::shared_ptr<CPdfMargZLogResult> postmargZResult);
    Int32 BestChi2(TFloat64List redshifts, std::vector<TFloat64List> meritResults, std::vector<TFloat64List> zPriors, Float64 cstLog, std::shared_ptr<CPdfMargZLogResult> postmargZResult);


private:


};


}

#endif

