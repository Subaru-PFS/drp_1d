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
 * NB-2018-02-19 : this method works with IRREGULAR z grid. No need to have a regular grid z-pdf anymore
 * @return 0: success, 1:problem, 3 not enough z values, 4: zPrior not valid
 */
Int32 CPdfz::Compute(TFloat64List merits, TFloat64List redshifts, Float64 cstLog, TFloat64List logZPrior, TFloat64List& logPdf, Float64 &logEvidence)
{
    Bool verbose = false;
    Int32 sumMethod = 1; //0=rect, 1=trapez
    logPdf.clear();

    if(verbose)
    {
        Float64 meritmax = -DBL_MAX;
        Float64 meritmin = DBL_MAX;
        for ( UInt32 k=0; k<redshifts.size(); k++)
        {
            if(meritmax<merits[k])
            {
                meritmax = merits[k];
            }
            if(meritmin>merits[k])
            {
                meritmin = merits[k];
            }
        }
        Log.LogInfo("Pdfz: Pdfz computation: using merit min=%e", meritmin);
        Log.LogInfo("Pdfz: Pdfz computation: using merit max=%e", meritmax);
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
    Float64 zstep_mean = (maxdz+mindz)/2.0;
    if(abs(maxdz-mindz)/zstep_mean>reldzThreshold)
    {
        constantdz = false;
        return 2;
    }

    //deactivate logPrior just to see...
//    for ( UInt32 k=0; k<redshifts.size(); k++)
//    {
//        logZPrior[k] = 0;
//    }

    //check if there is at least 1 redshifts values
    if(redshifts.size()==1) //consider this as a success
    {
        logPdf.resize(redshifts.size());
        logPdf[0] = 1.0;
        logEvidence = 1.0;
        return 0;
    }else if(redshifts.size()<1)
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

    logPdf.resize(redshifts.size());

    /* ------------------------------------------------------------------
    * NOTE: this uses the LOG-SUM-EXP trick originally suggested by S. Jamal
    * ------------------------------------------------------------------  */

    //prepare logLikelihood and LogEvidence
    Float64 maxi = -DBL_MAX;
    std::vector<Float64> mchi2Sur2(redshifts.size(), 0.0);
    std::vector<Float64> smallVALUES(redshifts.size(), 0.0);
    for ( UInt32 k=0; k<redshifts.size(); k++)
    {
        mchi2Sur2[k] = -0.5*merits[k];
        Float64 zstep;
        // here, using rect. approx. (enough precision for maxi estimation)
        if(k==0)
        {
            zstep = (redshifts[k+1]-redshifts[k])*0.5;
        }
        else if(k==redshifts.size()-1)
        {
            zstep = (redshifts[k]-redshifts[k-1])*0.5;
        }else
        {
            zstep = (redshifts[k+1]+redshifts[k])*0.5 - (redshifts[k]+redshifts[k-1])*0.5;
        }
        smallVALUES[k] = mchi2Sur2[k] + logZPrior[k];
        if(maxi<smallVALUES[k] + log(zstep))
        {
            maxi = smallVALUES[k] + log(zstep); // maxi will be used to avoid underflows when summing exponential of small values
        }
    }


    Float64 sumModifiedExp = 0.0;
    if(sumMethod==0)
    {
        Log.LogDebug("Pdfz: Pdfz computation: summation method option = RECTANGLES");
        for ( UInt32 k=0; k<redshifts.size(); k++)
        {
            Float64 modifiedEXPO = exp(smallVALUES[k]-maxi);
            Float64 area = (modifiedEXPO);
            Float64 zstep;
            if(k==0)
            {
                zstep = (redshifts[k+1]-redshifts[k])*0.5;
            }
            else if(k==redshifts.size()-1)
            {
                zstep = (redshifts[k]-redshifts[k-1])*0.5;
            }else
            {
                zstep = (redshifts[k+1]+redshifts[k])*0.5 - (redshifts[k]+redshifts[k-1])*0.5;
            }
            area *= zstep;
            sumModifiedExp += area;
        }
    }else if(sumMethod==1)
    {
        Log.LogDebug("Pdfz: Pdfz computation: summation method option = TRAPEZOID");
        Float64 modifiedEXPO_previous = exp(smallVALUES[0]-maxi);
        for ( UInt32 k=1; k<redshifts.size(); k++)
        {
            Float64 modifiedEXPO = exp(smallVALUES[k]-maxi);
            Float64 trapezArea = (modifiedEXPO+modifiedEXPO_previous)/2.0;
            trapezArea *= (redshifts[k]-redshifts[k-1]);
            sumModifiedExp += trapezArea;
            modifiedEXPO_previous = modifiedEXPO;
        }
    }else{
        Log.LogError("Pdfz: Pdfz computation: unable to parse summation method option");
    }
    logEvidence = cstLog + maxi + log(sumModifiedExp);

    if(verbose)
    {
        Log.LogInfo("Pdfz: Pdfz computation: using cstLog=%e", cstLog);
        Log.LogInfo("Pdfz: Pdfz computation: using logEvidence=%e", logEvidence);
        //Log.LogInfo("Pdfz: Pdfz computation: using logPrior=%f",  logPrior); //logPrior can be variable with z
    }

    for ( UInt32 k=0; k<redshifts.size(); k++)
    {
            logPdf[k] = logZPrior[k] + (mchi2Sur2[k] + cstLog) - logEvidence;
    }

    return 0;
}

Float64 CPdfz::getSumTrapez(std::vector<Float64> redshifts, std::vector<Float64> valprobalog)
{
    Float64 sum=0.0;

    //prepare LogEvidence
    Float64 maxi = -DBL_MAX;
    std::vector<Float64> smallVALUES(redshifts.size(), 0.0);
    for ( UInt32 k=0; k<redshifts.size(); k++)
    {
        //find the smallest zstep in order to use most penalizing case fot the log-sum-exp trick
        Float64 zstepPrevious = -1.;
        if(k>0){
            zstepPrevious = (redshifts[k]-redshifts[k-1]);
        }else if(k<redshifts.size()-1){
            zstepPrevious = (redshifts[k+1]-redshifts[k]);
        }
        Float64 zstepNext = -1.;
        if(k<redshifts.size()-1){
            zstepNext = (redshifts[k+1]-redshifts[k]);
        }else if(k>0){
            zstepNext = (redshifts[k]-redshifts[k-1]);
        }
        Float64 zstepCurrent = min(zstepPrevious, zstepNext);
        smallVALUES[k] = valprobalog[k];
        if(maxi<smallVALUES[k] + log(zstepCurrent))
        {
            maxi = smallVALUES[k] + log(zstepCurrent); // maxi will be used to avoid underflows when summing exponential of small values
        }
    }

    Float64 sumModifiedExp = 0.0;
    Float64 modifiedEXPO_previous = exp(smallVALUES[0]-maxi);
    for ( UInt32 k=1; k<redshifts.size(); k++)
    {
        Float64 modifiedEXPO = exp(smallVALUES[k]-maxi);
        Float64 trapezArea = (modifiedEXPO+modifiedEXPO_previous)/2.0;
        trapezArea *= (redshifts[k]-redshifts[k-1]);
        sumModifiedExp += trapezArea;
        modifiedEXPO_previous = modifiedEXPO;
    }
    Float64 logSum = maxi + log(sumModifiedExp);

    sum = exp(logSum);

    return sum;
}

Float64 CPdfz::getSumRect(std::vector<Float64> redshifts, std::vector<Float64> valprobalog)
{
    Float64 sum=0.0;

    //prepare LogEvidence
    Float64 maxi = -DBL_MAX;
    std::vector<Float64> smallVALUES(redshifts.size(), 0.0);
    for ( UInt32 k=0; k<redshifts.size(); k++)
    {
        Float64 zstep;
        if(k==0)
        {
            zstep = (redshifts[k+1]-redshifts[k])*0.5;
        }
        else if(k==redshifts.size()-1)
        {
            zstep = (redshifts[k]-redshifts[k-1])*0.5;
        }else
        {
            zstep = (redshifts[k+1]+redshifts[k])*0.5 - (redshifts[k]+redshifts[k-1])*0.5;
        }
        smallVALUES[k] = valprobalog[k];
        if(maxi<smallVALUES[k] + log(zstep))
        {
            maxi = smallVALUES[k] + log(zstep); // maxi will be used to avoid underflows when summing exponential of small values
        }
    }

    Log.LogDebug("  pdfz: getSumRect - found maxi = %e", maxi);
    Float64 sumModifiedExp = 0.0;
    for ( UInt32 k=0; k<redshifts.size()-1; k++)
    {
        Float64 modifiedEXPO = exp(smallVALUES[k]-maxi);
        Float64 area = modifiedEXPO;
        Float64 zstep;
        if(k==0)
        {
            zstep = (redshifts[k+1]-redshifts[k])*0.5;
        }
        else if(k==redshifts.size()-1)
        {
            zstep = (redshifts[k]-redshifts[k-1])*0.5;
        }else
        {
            zstep = (redshifts[k+1]+redshifts[k])*0.5 - (redshifts[k]+redshifts[k-1])*0.5;
        }
        area *= zstep;
        sumModifiedExp += area;
    }
    Float64 logSum = maxi + log(sumModifiedExp);

    sum = exp(logSum);

    return sum;
}

/**
 * @brief CPdfz::getCandidateSumTrapez
 * @param redshifts
 * @param valprobalog
 * @param zcandidate
 * @param zwidth
 * @return -1 if error, else sum around the candidate
 */
Float64 CPdfz::getCandidateSumTrapez(std::vector<Float64> redshifts, std::vector<Float64> valprobalog, Float64 zcandidate, Float64 zwidth)
{
    //check that redshifts are sorted
    for ( UInt32 k=1; k<redshifts.size(); k++)
    {
        if(redshifts[k]<redshifts[k-1])
        {
            Log.LogError("    CPdfz::getCandidateSumTrapez - redshifts are not sorted for (at least) index={}", k);
            return -1.0;
        }
    }

    //find indexes kmin, kmax so that zmin and zmax are inside [ redshifts[kmin]:redshifts[kmax] ]
    Int32 kmin=0;
    Int32 kmax=redshifts.size()-1;
    Float64 halfzwidth = zwidth/2.0;
    for ( UInt32 k=0; k<redshifts.size(); k++)
    {
        if( redshifts[k] < (zcandidate-halfzwidth) )
        {
            kmin = k;
        }
    }
    Log.LogDebug("    CPdfz::getCandidateSumTrapez - kmin index=%d", kmin);

    for ( UInt32 k=redshifts.size()-1; k>0; k--)
    {
        if( redshifts[k] > (zcandidate+halfzwidth) )
        {
            kmax = k;
        }
    }
    Log.LogDebug("    CPdfz::getCandidateSumTrapez - kmax index=%d", kmax);

    //initialize the LOG-SUM-EXP trick
    Float64 maxi = -DBL_MAX;
    std::vector<Float64> smallVALUES(redshifts.size(), 0.0);
    for ( UInt32 k=kmin; k<=kmax; k++)
    {
        //estimate zstep for the log-sum-exp trick
        Float64 zstep;
        if(k==0)
        {
            zstep = (redshifts[k+1]-redshifts[k]);
        }
        else if(k==redshifts.size()-1)
        {
            zstep = (redshifts[k]-redshifts[k-1]);
        }else
        {
            zstep = (redshifts[k+1]+redshifts[k])*0.5 - (redshifts[k]+redshifts[k-1])*0.5;
        }
        smallVALUES[k] = valprobalog[k];
        if(maxi<smallVALUES[k] + log(zstep))
        {
            maxi = smallVALUES[k] + log(zstep); // maxi will be used to avoid underflows when summing exponential of small values
        }
    }

    //for now the sum is estimated between kmin and kmax.
    //todo: INTERPOLATE (linear) in order to start exactly at zmin and stop at zmax
    Float64 sum=0.0;
    Float64 sumModifiedExp = 0.0;
    Float64 modifiedEXPO_previous = exp(smallVALUES[kmin]-maxi);
    for ( UInt32 k=kmin+1; k<=kmax; k++)
    {
        Float64 modifiedEXPO = exp(smallVALUES[k]-maxi);
        Float64 trapezArea = (modifiedEXPO+modifiedEXPO_previous)/2.0;
        trapezArea *= (redshifts[k]-redshifts[k-1]);
        sumModifiedExp += trapezArea;
        modifiedEXPO_previous = modifiedEXPO;
    }
    Float64 logSum = maxi + log(sumModifiedExp);

    sum = exp(logSum);

    return sum;
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

    return logzPrior;
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

std::vector<Float64> CPdfz::GetEuclidNhaLogZPrior(std::vector<Float64> redshifts)
{
    std::vector<Float64> zPrior(redshifts.size(), 0.0);
    Float64 aCoeff = 1.0;
    for(UInt32 kz=0; kz<redshifts.size(); kz++)
    {
        zPrior[kz] = (1e-6)*( -3080.1*aCoeff*redshifts[kz]+15546.6 );
    }

    Float64 sum = 0.0;
    for(UInt32 kz=0; kz<redshifts.size(); kz++)
    {
        sum+=zPrior[kz];
    }
    if(sum>0)
    {
        for(UInt32 kz=0; kz<redshifts.size(); kz++)
        {
            zPrior[kz] /= sum;
        }
    }

    //switch to log
    std::vector<Float64> logzPrior(redshifts.size(), 0.0);
    for(UInt32 kz=0; kz<redshifts.size(); kz++)
    {
        logzPrior[kz] = log(zPrior[kz]);
    }

    return logzPrior;
}

/**
 * @brief CombineLogZPrior
 * returns a vector with log(prior1*prior2) for each z, with normalization so that sumPrior=1
 * @param logprior1
 * @param logprior2
 * @return
 */
std::vector<Float64> CPdfz::CombineLogZPrior(std::vector<Float64> logprior1, std::vector<Float64> logprior2)
{
    std::vector<Float64> logzPriorCombined;
    if(logprior1.size()!=logprior2.size())
    {
        return logzPriorCombined;
    }
    Int32 n=logprior1.size();

    Float64 maxi = DBL_MIN;
    for ( UInt32 k=0; k<n; k++)
    {
        Float64 val= logprior1[k]+logprior2[k];
        if(maxi<val)
        {
            maxi=val;
        }
    }
    std::vector<Float64> zPrior_modif(n, 0.0);
    for ( UInt32 k=0; k<n; k++)
    {
        zPrior_modif[k] = logprior1[k]+logprior2[k]-maxi;
    }

    Float64 sumExpModif=0.0;
    for ( UInt32 k=0; k<n; k++)
    {
        sumExpModif+=exp(zPrior_modif[k]);
    }

    Float64 logSum = log(exp(maxi))*sumExpModif;

    logzPriorCombined.resize(n);
    for(Int32 k=0; k<n; k++)
    {
        logzPriorCombined[k] =  logprior1[k]+logprior2[k]-logSum;
    }


    return logzPriorCombined;
}


Int32 CPdfz::Marginalize(TFloat64List redshifts, std::vector<TFloat64List> meritResults, std::vector<TFloat64List> zPriors, Float64 cstLog, std::shared_ptr<CPdfMargZLogResult> postmargZResult, std::vector<Float64> modelPriors)
{
    bool verbose = false;

    if(meritResults.size() != zPriors.size())
    {
        Log.LogError("Pdfz: Pdfz marginalize problem. merit.size (%d) != prior.size (%d)", meritResults.size(), zPriors.size());
        return -9;
    }
    if(meritResults.size()<1 || zPriors.size()<1 || redshifts.size()<1)
    {
        Log.LogError("Pdfz: Pdfz marginalize problem. merit.size (%d), prior.size (%d), or redshifts.size (%d) is zero !", meritResults.size(), zPriors.size(), redshifts.size());
        return -99;
    }

    //check merit curves. Maybe this should be assert stuff ?
    for(Int32 km=0; km<meritResults.size(); km++)
    {
        TFloat64List _merit = meritResults[km];
        Bool invalidFound = false;
        for(Int32 kz=0; kz<_merit.size(); kz++)
        {
            if(_merit[kz] != _merit[kz])
            {
                Log.LogError("    CPdfz::Marginalize - merit result #%d has at least one nan or invalid value at index=%d", km, kz);
                invalidFound=true;
            }
            if(invalidFound)
            {
                break;
            }
        }
    }


    Bool initPostMarg = false;
    std::vector<UInt32> nSum;

    Float64 MaxiLogEvidence = -DBL_MAX;
    TFloat64List LogEvidencesWPriorM;
    Float64 sumModifiedEvidences = 0;

    std::vector<Float64> logPriorModel;
    if( /*false &&*/ modelPriors.size()!=meritResults.size())
    {
        Float64 priorModelCst = 1.0/((Float64)meritResults.size());
        Log.LogInfo("Pdfz: Marginalize: no priors loaded, using constant priors (=%f)", priorModelCst);
        for(Int32 km=0; km<meritResults.size(); km++)
        {
            logPriorModel.push_back(log(priorModelCst));
        }
    }else
    {
        /*
        //override modelPriors with pypelid 10 knn templates priors
        logPriorModel.push_back(log(0.1490));
        logPriorModel.push_back(log(0.0794));
        logPriorModel.push_back(log(0.0744));
        logPriorModel.push_back(log(0.0836));
        logPriorModel.push_back(log(0.1089));
        logPriorModel.push_back(log(0.1124));
        logPriorModel.push_back(log(0.0786));
        logPriorModel.push_back(log(0.1563));
        logPriorModel.push_back(log(0.0509));
        logPriorModel.push_back(log(0.1060));
        //*/
        for(Int32 km=0; km<meritResults.size(); km++)
        {
            logPriorModel.push_back(log(modelPriors[km]));
        }

        //logging priors used
        Float64 sumPriors=0.;
        for(Int32 km=0; km<meritResults.size(); km++)
        {
            sumPriors += exp(logPriorModel[km]);
            Log.LogInfo("Pdfz: Marginalize: for model k=%d, using prior=%f", km, exp(logPriorModel[km]));
        }
        Log.LogInfo("Pdfz: Marginalize: sumPriors=%f", sumPriors);
        if(sumPriors>1.1 || sumPriors<0.9)
        {
            Log.LogError("Pdfz: sumPriors should be close to 1... !!!");
        }
    }

    for(Int32 km=0; km<meritResults.size(); km++)
    {
        //Todo: Check if the status is OK ?
        //meritResult->Status[i] == COperator::nStatus_OK

        CPdfz pdfz;
        TFloat64List logProba;
        Float64 logEvidence;
        Int32 retPdfz = pdfz.Compute(meritResults[km], redshifts, cstLog, zPriors[km], logProba, logEvidence);
        if(retPdfz!=0)
        {
            Log.LogError("Pdfz: Pdfz computation failed for result km=%d", km);
            return -1;
        }else{
//            if(verbose)
//            {
//                Log.LogInfo("Pdfz: Marginalize: for km=%d, logEvidence=%e", km, MaxiLogEvidence);
//            }
            Float64 logEvidenceWPriorM = logEvidence+logPriorModel[km];

            LogEvidencesWPriorM.push_back(logEvidenceWPriorM);
            if(MaxiLogEvidence<logEvidenceWPriorM){
                MaxiLogEvidence=logEvidenceWPriorM;
            }
        }
    }
    if(verbose)
    {
        Log.LogInfo("Pdfz: Marginalize: MaxiLogEvidence=%e", MaxiLogEvidence);
    }

    //Using computational trick to sum the evidences
    for(Int32 k=0; k<LogEvidencesWPriorM.size(); k++)
    {
        sumModifiedEvidences += exp(LogEvidencesWPriorM[k]-MaxiLogEvidence);
    }
    Float64 logSumEvidence = MaxiLogEvidence + log(sumModifiedEvidences);
    if(verbose)
    {
        Log.LogInfo("Pdfz: Marginalize: logSumEvidence=%e", logSumEvidence);
    }

    for(Int32 km=0; km<meritResults.size(); km++)
    {
        if(verbose)
        {
            Log.LogInfo("Pdfz: Marginalize: processing chi2-result km=%d", km);
        }

        //Todo: Check if the status is OK ?
        //meritResult->Status[i] == COperator::nStatus_OK

        CPdfz pdfz;
        TFloat64List logProba;
        Float64 logEvidence;
        Int32 retPdfz = pdfz.Compute(meritResults[km], redshifts, cstLog, zPriors[km], logProba, logEvidence);
        if(retPdfz!=0)
        {
            Log.LogError("Pdfz: Pdfz computation failed for result km=%d", km);
            return -1;
        }else{
            if(!initPostMarg)
            {
                nSum.resize(redshifts.size());
                postmargZResult->countTPL = redshifts.size(); // assumed 1 model per z
                postmargZResult->Redshifts.resize(redshifts.size());
                postmargZResult->valProbaLog.resize(redshifts.size());
                for ( UInt32 k=0; k<redshifts.size(); k++)
                {
                    postmargZResult->Redshifts[k] = redshifts[k] ;
                    postmargZResult->valProbaLog[k] = log(0.0);
                    nSum[k] = 0;
                }
                initPostMarg = true;
            }else
            {
                //check if the redshift bins are the same
                for ( UInt32 k=0; k<redshifts.size(); k++)
                {
                    if(postmargZResult->Redshifts[k] != redshifts[k])
                    {
                        Log.LogError("pdfz: Pdfz computation (z-bins comparison) failed for result km=%d", km);
                        break;
                    }
                }
            }

            postmargZResult->valEvidenceLog = logSumEvidence;
            for ( UInt32 k=0; k<redshifts.size(); k++)
            {
                if( true /*meritResult->Status[k]== COperator::nStatus_OK*/) //todo: check (temporarily considers status is always OK for linemodel tplshape)
                {
                    Float64 logValProba = postmargZResult->valProbaLog[k];
                    Float64 logValProbaAdd = logProba[k]+logPriorModel[km]+logEvidence-logSumEvidence;
                    Float64 maxP = logValProba;
                    if(maxP<logValProbaAdd)
                    {
                        maxP=logValProbaAdd;
                    }
                    Float64 valExp = exp(logValProba-maxP)+exp(logValProbaAdd-maxP);
                    postmargZResult->valProbaLog[k] = maxP+log(valExp);
                    nSum[k]++;
                }
            }
        }


    }

    //THIS DOES NOT ALLOW Marginalization with coverage<100% for ALL templates
    for ( UInt32 k=0; k<postmargZResult->Redshifts.size(); k++)
    {
        if(nSum[k]!=meritResults.size())
        {
            postmargZResult->valProbaLog[k] = NAN;
            Log.LogError("Pdfz: Pdfz computation failed. For z=%f, nSum=%d", postmargZResult->Redshifts[k], nSum[k]);
            Log.LogError("Pdfz: Pdfz computation failed. For z=%f, meritResults.size()=%d", postmargZResult->Redshifts[k], meritResults.size());
            Log.LogError("Pdfz: Pdfz computation failed. Not all templates have 100 percent coverage for all redshifts!");
        }
    }


    return 0;
}

//This mathematically does not correspond to any valid method for combining PDFs.
//TODO: problem while estimating best proba. is it best proba for each z ? In that case: what about sum_z P = 1 ?
//TODO: this methid should be replaced/modified to correspond to the MaxPDF technique.
Int32 CPdfz::BestProba(TFloat64List redshifts, std::vector<TFloat64List> meritResults, std::vector<TFloat64List> zPriors, Float64 cstLog, std::shared_ptr<CPdfMargZLogResult> postmargZResult)
{
    bool verbose = false;

    if(meritResults.size() != zPriors.size())
    {
        Log.LogError("Pdfz: Pdfz-bestproba problem. merit.size (%d) != prior.size (%d)", meritResults.size(), zPriors.size());
        return -9;
    }
    if(meritResults.size()<1 || zPriors.size()<1 || redshifts.size()<1)
    {
        Log.LogError("Pdfz: Pdfz-bestproba problem. merit.size (%d), prior.size (%d), or redshifts.size (%d) is zero !", meritResults.size(), zPriors.size(), redshifts.size());
        return -99;
    }

    Bool initPostMarg = false;

    for(Int32 km=0; km<meritResults.size(); km++)
    {
        if(verbose)
        {
            Log.LogInfo("Pdfz:-bestproba: processing chi2-result km=%d", km);
        }

        //Todo: Check if the status is OK ?
        //meritResult->Status[i] == COperator::nStatus_OK

        CPdfz pdfz;
        TFloat64List logProba;
        Float64 logEvidence;
        Int32 retPdfz = pdfz.Compute(meritResults[km], redshifts, cstLog, zPriors[km], logProba, logEvidence);
        if(retPdfz!=0)
        {
            Log.LogError("Pdfz: Pdfz-bestproba computation failed for result km=%d", km);
            return -1;
        }else{
            if(!initPostMarg)
            {
                postmargZResult->countTPL = redshifts.size(); // assumed 1 model per z
                postmargZResult->Redshifts.resize(redshifts.size());
                postmargZResult->valProbaLog.resize(redshifts.size());
                for ( UInt32 k=0; k<redshifts.size(); k++)
                {
                    postmargZResult->Redshifts[k] = redshifts[k] ;
                    postmargZResult->valProbaLog[k] = -DBL_MAX;
                }
                initPostMarg = true;
            }else
            {
                //check if the redshift bins are the same
                for ( UInt32 k=0; k<redshifts.size(); k++)
                {
                    if(postmargZResult->Redshifts[k] != redshifts[k])
                    {
                        Log.LogError("pdfz: Pdfz-bestproba computation (z-bins comparison) failed for result km=%d", km);
                        break;
                    }
                }
            }
            for ( UInt32 k=0; k<redshifts.size(); k++)
            {
                if( true /*meritResult->Status[k]== COperator::nStatus_OK*/) //todo: check (temporarily considers status is always OK for linemodel tplshape)
                {
                    if(logProba[k]>postmargZResult->valProbaLog[k])
                    {
                        postmargZResult->valProbaLog[k]=logProba[k];
                    }
                }
            }
        }
    }

    //normalize: sum_z P = 1
    //1. check if the z step is constant. If not, pdf cannot be estimated by the current method.
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

    //2. prepare LogEvidence
    Float64 maxi = -DBL_MAX;
    std::vector<Float64> smallVALUES(redshifts.size(), 0.0);
    for ( UInt32 k=0; k<redshifts.size(); k++)
    {
        smallVALUES[k] = postmargZResult->valProbaLog[k];
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
    Float64 logEvidence = maxi + log(sumModifiedExp) + log(zstep);

    if(verbose)
    {
        Log.LogInfo("Pdfz: Pdfz-bestproba computation: using logEvidence=%e", logEvidence);
        Log.LogInfo("Pdfz: Pdfz-bestproba computation: using log(zstep)=%e",  log(zstep));
    }

    for ( UInt32 k=0; k<redshifts.size(); k++)
    {
            postmargZResult->valProbaLog[k] = postmargZResult->valProbaLog[k] - logEvidence;
    }

    return 0;
}

/**
 * @brief CPdfz::BestChi2
 * Computes the combined pdf by taking the best chi2
 * WARNING: as long as the prior on the models are constant, it is equivalent to compute the bestchi2 and the MaxPDF. If this prior is not constant any more, the mas search has to be modified.
 * @param redshifts
 * @param meritResults
 * @param zPriors
 * @param cstLog
 * @param postmargZResult
 * @return
 */
Int32 CPdfz::BestChi2(TFloat64List redshifts, std::vector<TFloat64List> meritResults, std::vector<TFloat64List> zPriors, Float64 cstLog, std::shared_ptr<CPdfMargZLogResult> postmargZResult)
{
    bool verbose = false;

    if(meritResults.size() != zPriors.size())
    {
        Log.LogError("Pdfz: Pdfz-bestchi2 problem. merit.size (%d) != prior.size (%d)", meritResults.size(), zPriors.size());
        return -9;
    }
    if(meritResults.size()<1 || zPriors.size()<1 || redshifts.size()<1)
    {
        Log.LogError("Pdfz: Pdfz-bestchi2 problem. merit.size (%d), prior.size (%d), or redshifts.size (%d) is zero !", meritResults.size(), zPriors.size(), redshifts.size());
        return -99;
    }

    //build best chi2 vector
    std::vector<Float64> chi2Min(redshifts.size(), DBL_MAX);
    for(Int32 km=0; km<meritResults.size(); km++)
    {
        if(verbose)
        {
            Log.LogInfo("Pdfz:-bestchi2: processing chi2-result km=%d", km);
        }

        //Todo: Check if the status is OK ?
        //meritResult->Status[i] == COperator::nStatus_OK

        //Todo: use the priors for the min chi2 search ?
        for ( UInt32 k=0; k<redshifts.size(); k++)
        {
            if( true /*meritResult->Status[k]== COperator::nStatus_OK*/) //todo: check (temporarily considers status is always OK)
            {
                if(meritResults[km][k]<chi2Min[k])
                {
                    chi2Min[k] = meritResults[km][k];
                }
            }
        }

    }

    //estimate Posterior on the best chi2
    CPdfz pdfz;
    TFloat64List logProba;
    Float64 logEvidence;
    TFloat64List zprior;
    zprior = pdfz.GetConstantLogZPrior(redshifts.size());
    Int32 retPdfz = pdfz.Compute(chi2Min, redshifts, cstLog, zprior, logProba, logEvidence);
    if(retPdfz!=0)
    {
        Log.LogError("Pdfz: Pdfz-bestchi2: Pdfz computation failed");
        return -1;
    }else{
        postmargZResult->valEvidenceLog = logEvidence;
        postmargZResult->countTPL = redshifts.size(); // assumed 1 model per z
        postmargZResult->Redshifts.resize(redshifts.size());
        postmargZResult->valProbaLog.resize(redshifts.size());
        for ( UInt32 k=0; k<redshifts.size(); k++)
        {
            postmargZResult->Redshifts[k] = redshifts[k] ;
            postmargZResult->valProbaLog[k] = logProba[k];
        }
    }

    return 0;
}
