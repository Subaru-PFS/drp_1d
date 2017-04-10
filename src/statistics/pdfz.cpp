#include <epic/redshift/statistics/pdfz.h>

#include <epic/core/log/log.h>

using namespace NSEpic;
using namespace std;
#include <fstream>

#include <gsl/gsl_multifit.h>

CPdfz::CPdfz()
{

}

CPdfz::~CPdfz()
{

}

/**
 * @brief CPdfz::Compute
 * @param merits
 * @param redshifts
 * ...
 *
 * @return 0: success, 1:problem,
 */
Int32 CPdfz::Compute(TFloat64List merits, TFloat64List redshifts, Float64 cstLog, TFloat64List& logPdf)
{
    logPdf.clear();
    logPdf.resize(redshifts.size());
    Float64 logPrior = 0; //log(1.0)

    /* ------------------------------------------------------------------
    * NOTE (copied from dev_bayes branch by S. Jamal):
    * -------
    * The sum is realised in ( Z , TPL) 2D space
    * csteLOG = -N/2*LOG(2*pi) - LOG ( product ( sigma_RMSE) )
    *      =  -N/2*LOG(2*pi) - sum ( LOG ( sigma_RMSE) )
    *      csteLOG is independent from the (z, tpl) space; only dependant
    *      of input spectrum (rmse & nbr points for the overlapRate)
    *
    * LogEvidence  = LOG ( sum  ( Prior * Likelihood   )  )
    *      = LOG ( sum ( EXPONENTIAL ( logPrior + logLikelihood )  )  )
    *      = LOG ( sum ( EXPONENTIAL ( logPrior  - chi2/2 +  csteLOG )  )  )
    *      = LOG ( csteLOG * sum ( EXPONENTIAL  ( logPrior - chi2/2 )  )  )
    *      = csteLOG + LOG ( sum ( EXPONENTIAL ( logPrior  - chi2/2  )  )
    *
    * The values "logPrior  - chi2/2"   are very small.
    * Cannot comput EXP directly.
    * The LOG-SUM-EXP computational trick is used here:
    *              smallVALUES = logPrior  - chi2/2
    *              MAXI = max(smallVALUES)
    *              ModifiedEXPO = EXPONENTIAL  ( smallVALUES - MAXI  )
    *
    *              LogEvidence = csteLOG + MAXI +  LOG ( sum  ( ModifiedEXPO )
    *
    * ------------------------------------------------------------------  */

    //prepare logLikelihood and LogEvidence
    Float64 maxi = -DBL_MAX;
    std::vector<Float64> mchi2Sur2(redshifts.size(), 0.0);
    std::vector<Float64> smallVALUES(redshifts.size(), 0.0);
    for ( UInt32 k=0; k<redshifts.size(); k++)
    {
        mchi2Sur2[k] = -0.5*merits[k];
        smallVALUES[k] = mchi2Sur2[k] + logPrior;
        if(maxi<smallVALUES[k])
        {
            maxi = smallVALUES[k]; // maxi will be used to avoid underflows when summing exponential of small values
        }
    }

    Float64 sumModifiedExp = 0.0;
    for ( UInt32 k=0; k<redshifts.size(); k++)
    {
        Float64 modifiedEXPO = exp(smallVALUES[k]-maxi);
        sumModifiedExp += modifiedEXPO;
        //logEvidence;
    }
    Float64 logEvidence = cstLog + maxi + log(sumModifiedExp);



    for ( UInt32 k=0; k<redshifts.size(); k++)
    {
            logPdf[k] = logPrior + (mchi2Sur2[k] + cstLog) - logEvidence;
    }

    return 0;
}

