#include <RedshiftLibrary/statistics/pdfz.h>

#include <RedshiftLibrary/log/log.h>

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
 * @return 0: success, 1:problem, 2:dz not constant, 3 not enough z values
 */
Int32 CPdfz::Compute(TFloat64List merits, TFloat64List redshifts, Float64 cstLog, TFloat64List& logPdf)
{
    logPdf.clear();

    //check if there is more than 2 reshifts values
    if(redshifts.size()<=2)
    {
        return 2;
    }

    //check if the z step is constant. If not, pdf cannot be estimated by the current method.
    Float64 reldzThreshold = 0.05; //relative difference accepted
    bool constantdz = true;
    Float64 mindz = DBL_MAX;
    Float64 maxdz = -DBL_MAX;
    for ( UInt32 k=1; k<redshifts.size(); k++)
    {
        Float64 diff = redshifts[k]-redshifts[k-1];
        if(mindz > diff)
        {
            mindz = diff;
        }
        if(maxdz < diff)
        {
            maxdz = diff;
        }
    }
    Float64 zstep = (maxdz+mindz)/2.0;
    if(abs(maxdz-mindz)/zstep>reldzThreshold)
    {
        constantdz = false;
        return 2;
    }

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
    }
    Float64 logEvidence = cstLog + maxi + log(sumModifiedExp) + log(zstep);



    for ( UInt32 k=0; k<redshifts.size(); k++)
    {
            logPdf[k] = logPrior + (mchi2Sur2[k] + cstLog) - logEvidence;
    }

    return 0;
}

