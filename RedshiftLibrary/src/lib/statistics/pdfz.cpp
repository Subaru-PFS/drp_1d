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
 * @return 0: success, 1:problem, 2:dz not constant, 3 not enough z values, 4: zPrior not valid
 */
Int32 CPdfz::Compute(TFloat64List merits, TFloat64List redshifts, Float64 cstLog, TFloat64List logZPrior, TFloat64List& logPdf, Float64 &logEvidence)
{
    Bool verbose = false;
    logPdf.clear();

    //check if there is more than 2 redshifts values
    if(redshifts.size()<=2)
    {
        return 2;
    }

    //check that the zPrior is size-compatible
    if(logZPrior.size()!=redshifts.size())
    {
        return 4;
        //Float64 logPrior = log(1.0/redshifts.size()); //log(1.0)
    }
//    std::vector<Float64> logZPrior(zPrior.size(), 1.0);
//    for(UInt32 kz=0; kz<zPrior.size(); kz++)
//    {
//        logZPrior[kz] = log(zPrior[kz]);
//    }

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
        smallVALUES[k] = mchi2Sur2[k] + logZPrior[k];
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
    logEvidence = cstLog + maxi + log(sumModifiedExp) + log(zstep);

    if(verbose)
    {
        Log.LogInfo("chisquare2solve: Pdfz computation: using cstLog=%f", cstLog);
        Log.LogInfo("chisquare2solve: Pdfz computation: using logEvidence=%f", logEvidence);
        Log.LogInfo("chisquare2solve: Pdfz computation: using log(zstep)=%f",  log(zstep));
        //Log.LogInfo("chisquare2solve: Pdfz computation: using logPrior=%f",  logPrior); //logPrior can be variable with z
    }

    for ( UInt32 k=0; k<redshifts.size(); k++)
    {
            logPdf[k] = logZPrior[k] + (mchi2Sur2[k] + cstLog) - logEvidence;
    }

    return 0;
}

std::vector<Float64> CPdfz::GetConstantLogZPrior(UInt32 nredshifts)
{
    std::vector<Float64> zPrior(nredshifts, 1.0);
    for(UInt32 kz=0; kz<nredshifts; kz++)
    {
        zPrior[kz] = 1.0/nredshifts;
    }

    //switch to log
    std::vector<Float64> logzPrior(nredshifts, 0.0);
    for(UInt32 kz=0; kz<nredshifts; kz++)
    {
        logzPrior[kz] = log(zPrior[kz]);
    }

    return zPrior;
}

std::vector<Float64> CPdfz::GetStrongLinePresenceLogZPrior(std::vector<bool> linePresence)
{
    Float64 probaPresent = 1.0;
    Float64 probaAbsent = 1e-5;
    std::vector<Float64> zPrior(linePresence.size(), probaAbsent);
    Float64 sum = 0.0;
    for(UInt32 kz=0; kz<linePresence.size(); kz++)
    {
        if(linePresence[kz])
        {
            zPrior[kz] = probaPresent;
        }else{
            zPrior[kz] = probaAbsent;
        }
        sum += zPrior[kz];
    }

    if(sum>0)
    {
        for(UInt32 kz=0; kz<linePresence.size(); kz++)
        {
            zPrior[kz] /= sum;
        }
    }

    //switch to log
    std::vector<Float64> logzPrior(linePresence.size(), probaAbsent);
    for(UInt32 kz=0; kz<linePresence.size(); kz++)
    {
        logzPrior[kz] = log(zPrior[kz]);
    }

    return logzPrior;
}
