#include <RedshiftLibrary/operator/pdfz.h>
#include <RedshiftLibrary/statistics/zprior.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/extremum/extremum.h>

using namespace NSEpic;
using namespace std;
#include <fstream>
#include <algorithm>
#include <cfloat>

COperatorPdfz::COperatorPdfz(const std::string & opt_combine, 
                             Float64 peakSeparation,
                             Float64 meritcut,
                             UInt32 maxCandidate,
                             const std::string & Id_prefix,
                             UInt32 maxPeakCount_per_window,
                             const std::vector<TFloat64List> & candidatesRedshifts,
                             const TStringList & candidatesIds
                             ):
    m_opt_combine(opt_combine),
    m_peakSeparation(peakSeparation),
    m_meritcut(meritcut),
    m_maxPeakCount_per_window(maxPeakCount_per_window<=0? maxCandidate:maxPeakCount_per_window),
    m_maxCandidate(maxCandidate),
    m_Id_prefix(Id_prefix)
{
    TStringList::const_iterator begin1 = candidatesIds.begin();
    TStringList::const_iterator end1 = candidatesIds.end();
    std::vector<TFloat64List>::const_iterator begin2 = candidatesRedshifts.begin();
    std::vector<TFloat64List>::const_iterator end2 = candidatesRedshifts.end();
    TStringList::const_iterator i1;
    std::vector<TFloat64List>::const_iterator i2;
    for(i1=begin1, i2=begin2; (i1!=end1) && (i2!=end2); ++i1, ++i2)
    {
        m_candidatesZRanges[*i1] = TFloat64Range(*i2); 
        // with nocandidatesRedshifts given will be with Id="" and empty range value.
    }
}


/*
* Main Pdf operator entrance
*  combine pdf and search for candidates.
*/
std::shared_ptr<CPdfCandidateszResult> COperatorPdfz::Compute(const ChisquareArray & chisquarearray,
                                                              Bool integ)
{

    // build PDF from chisquares and priors
    CombinePDF(chisquarearray);

    // find candidates redshifts
    TCandidateZbyID zcandidates; // will be sorted by the id syntax inside each redhisftwindow
    searchMaxPDFcandidates(zcandidates);

    std::shared_ptr<CPdfCandidateszResult> CandidateszResult;
    if (integ){
        // compute pdf candidate properties (deltaz, integ, rank )
        CPdfCandidatesZ zcand_op =  CPdfCandidatesZ(zcandidates);
        CandidateszResult = zcand_op.Compute(m_postmargZResult->Redshifts, m_postmargZResult->valProbaLog);

        //eventually truncate candidates at maxcount
        size_t newsize = std::min(CandidateszResult->m_ranked_candidates.size(), size_t(m_maxCandidate));
        CandidateszResult->m_ranked_candidates.resize(newsize);
    }else{
        CandidateszResult = std::make_shared<CPdfCandidateszResult>();
            for (auto c:zcandidates){
                CandidateszResult->m_ranked_candidates.emplace_back(c);
            }
    }

    return CandidateszResult;
}



Int32 COperatorPdfz::CombinePDF(const ChisquareArray & chisquarearray)
{
    Log.LogInfo("Pdfz combination");

    Int32 retPdfz = -1;

    if(chisquarearray.chisquares.size()>0)
    {
        // initialize m_postmargZResult
        m_postmargZResult = std::make_shared<CPdfMargZLogResult>(chisquarearray.redshifts);

        if(m_opt_combine=="marg")
        {
            Log.LogInfo("    Pdfz combination - Marginalization");
            retPdfz = Marginalize( chisquarearray);
            
        }else if(m_opt_combine=="bestchi2")
        {
            Log.LogInfo("    Pdfz combination - BestChi2");
            retPdfz = BestChi2( chisquarearray);

        }else if(m_opt_combine=="bestproba"){
            Log.LogInfo("    Pdfz combination - BestProba");
            retPdfz = BestProba( chisquarearray);

        }else{
            Log.LogError("   Pdfz combination: Unable to parse pdf combination method option");
            throw runtime_error("Pdfz combination: Unable to parse pdf combination method option");

        }
    }else
    {
        Log.LogError("    Pdfz combination: Unable to find any chisquares prepared for combination. chisquares.size()=%d", chisquarearray.chisquares.size());
        throw runtime_error("Pdfz combination: Unable to find any chisquares prepared for combination");
    }

    if(retPdfz!=0)
    {
        Log.LogError("    Pdfz combination: Pdfz computation failed");
        throw runtime_error("Pdfz combination: computation failed");

    }

    //check pdf sum=1
    if(!checkPdfSum()){
        Log.LogError("%s: Pdfz normalization failed", __func__);
        throw std::runtime_error("Pdfz normalization failed");
    }

    return retPdfz;    
}

Bool COperatorPdfz::checkPdfSum()
{
    Bool ret = true;

    //check pdf sum=1
    Float64 sumRect = getSumRect(m_postmargZResult->Redshifts, m_postmargZResult->valProbaLog);
    Float64 sumTrapez = getSumTrapez(m_postmargZResult->Redshifts, m_postmargZResult->valProbaLog);
    Log.LogDetail("    COperatorPdfz: Pdfz normalization - sum rect. = %e", sumRect);
    Log.LogDetail("    COperatorPdfz: Pdfz normalization - sum trapz. = %e", sumTrapez);
    Bool pdfSumCheck = abs(sumRect-1.0)<1e-1 || abs(sumTrapez-1.0)<1e-1;
    if(!pdfSumCheck){
        Log.LogError("    COperatorPdfz: Pdfz normalization failed (rectsum = %f, trapzesum = %f)", sumRect, sumTrapez);
        ret = false;
    }

    return ret;
}


Bool COperatorPdfz::searchMaxPDFcandidates(TCandidateZbyID & candidates) const
{
    candidates.clear();

    for (const auto & cand : m_candidatesZRanges)
    {
        TPointList extremumList;
        const TFloat64Range & redshiftsRange = cand.second;
        std::string id = cand.first;

        if (id!="") id += "_";

        //call Find on each secondpass range and retrieve the best  peak
        Bool invertForMinSearch=false;
        CExtremum extremum_op = CExtremum( m_maxPeakCount_per_window, m_peakSeparation, m_meritcut, invertForMinSearch, redshiftsRange); // no peak separation, no cut
        Bool findok = extremum_op.Find(m_postmargZResult->Redshifts, m_postmargZResult->valProbaLog, extremumList);
        if (!findok)
        {   
            Log.LogError("COperatorPdfz: searchMaxPDFcandidates failed");
            throw runtime_error("COperatorPdfz: searchMaxPDFcandidates failed");
            return false;
        }
        Int32 i = 0 ;
        for (const auto & extremum : extremumList)
        {   
            std::string newid = id + m_Id_prefix + std::to_string(i);
            candidates[newid].Redshift = extremum.X;
            candidates[newid].ValProba = extremum.Y;    
            ++i;
        }
    }

    return true;
}




/**
 * Use the log sum exp trick to sum up small numbers while avoiding underflows
 *   * ------------------------------------------------------------------
     * NOTE: this uses the LOG-SUM-EXP trick originally suggested by S. Jamal
     * ------------------------------------------------------------------  *
*/
Float64 COperatorPdfz::logSumExpTrick(const TFloat64List & valproba, const TFloat64List & redshifts, Int32 sumMethod)
{

    Float64 logfactor = -DBL_MAX;
    if(redshifts.size()<2){
        Log.LogError("COperatorPdfz::logSumExpTrick Can't compute on a range of less than 2 points");
        throw runtime_error("COperatorPdfz::logSumExpTrick Can't compute on a range of less than 2 points");
    }
    
    for (UInt32 k = 0; k < redshifts.size(); k++)
    {
       
        Float64 zstep;
        // here, using rect. approx. (enough precision for maxi estimation)
        if (k == 0)
        {
            zstep = (redshifts[k + 1] - redshifts[k]) * 0.5;
        } else if (k == redshifts.size() - 1)
        {
            zstep = (redshifts[k] - redshifts[k - 1]) * 0.5;
        } else
        {
            zstep = (redshifts[k + 1] + redshifts[k]) * 0.5 -
                    (redshifts[k] + redshifts[k - 1]) * 0.5;
        }
        if (logfactor < valproba[k] + log(zstep))
        {
            logfactor = valproba[k] + log(zstep); // maxi will be used to avoid underflows when
                                                  // summing exponential of small values
        }
    }
    
    Log.LogDebug("Pdfz: Pdfz computation: using common factor value for log-sum-exp trick=%e", logfactor);

    Float64 sumModifiedExp = 0.0;
    if (sumMethod == 0)
    {
        Log.LogDebug("Pdfz: Pdfz computation: summation method option = RECTANGLES");
        for (UInt32 k = 0; k < redshifts.size(); k++)
        {
            Float64 modifiedEXPO = exp(valproba[k] - logfactor);
            Float64 area = (modifiedEXPO);
            Float64 zstep;
            if (k == 0)
            {
                zstep = (redshifts[k + 1] - redshifts[k]);
            } else if (k == redshifts.size() - 1)
            {
                zstep = (redshifts[k] - redshifts[k - 1]);
            } else
            {
                zstep = (redshifts[k + 1] - redshifts[k - 1]) * 0.5;
            }
            area *= zstep;
            sumModifiedExp += area;
        }
    } else if (sumMethod == 1)
    {
        Log.LogDebug("Pdfz: Pdfz computation: summation method option = TRAPEZOID");
        Float64 modifiedEXPO_previous = exp(valproba[0] - logfactor);
        for (UInt32 k = 1; k < redshifts.size(); k++)
        {
            Float64 modifiedEXPO = exp(valproba[k] - logfactor);
            Float64 trapezArea = (modifiedEXPO + modifiedEXPO_previous) / 2.0;
            trapezArea *= (redshifts[k] - redshifts[k - 1]);
            sumModifiedExp += trapezArea;
            modifiedEXPO_previous = modifiedEXPO;
        }
    } else
    {
        Log.LogError("Pdfz: Pdfz computation: unable to parse summation method option");
    }
    Float64 sum = logfactor + log(sumModifiedExp);

    return sum;
}

/**
 * @brief COperatorPdfz::Compute
 * @param merits
 * @param redshifts
 * ...
 * NB-2018-02-19 : this method works with IRREGULAR z grid. No need to have a
 * regular grid z-pdf anymore
 * @return 0: success, 1:problem, 3 not enough z values, 4: zPrior not valid
 */
Int32 COperatorPdfz::ComputePdf(const TFloat64List & merits, const TFloat64List & redshifts,
                                        const Float64 cstLog, const TFloat64List & logZPrior,
                                        TFloat64List &logPdf, Float64 &logEvidence)
{
    Bool verbose = true;
    Int32 sumMethod = 1; // 0=rect, 1=trapez
    logPdf.clear();

    if (verbose)
    {
        Float64 meritmax = -DBL_MAX;
        Float64 meritmin = DBL_MAX;
        for (UInt32 k = 0; k < redshifts.size(); k++)
        {
            if (meritmax < merits[k])
            {
                meritmax = merits[k];
            }
            if (meritmin > merits[k])
            {
                meritmin = merits[k];
            }
        }
        Log.LogDebug("Pdfz: Pdfz computation: using merit min=%e", meritmin);
        Log.LogDebug("Pdfz: Pdfz computation: using merit max=%e", meritmax);
    }

    //
    //    //check if the z step is constant. If not, pdf cannot be estimated by
    //    the current method. Float64 reldzThreshold = 0.05; //relative
    //    difference accepted bool constantdz = true; Float64 mindz = DBL_MAX;
    //    Float64 maxdz = -DBL_MAX;
    //    for ( UInt32 k=1; k<redshifts.size(); k++)
    //    {
    //        Float64 diff = redshifts[k]-redshifts[k-1];
    //        if(mindz > diff)
    //        {
    //            mindz = diff;
    //        }
    //        if(maxdz < diff)
    //        {
    //            maxdz = diff;
    //        }
    //    }
    //    Float64 zstep_mean = (maxdz+mindz)/2.0;
    //    if(abs(maxdz-mindz)/zstep_mean>reldzThreshold)
    //    {
    //        constantdz = false;
    //        return 2;
    //    }
    //

    // deactivate logPrior just to see...
    //    for ( UInt32 k=0; k<redshifts.size(); k++)
    //    {
    //        logZPrior[k] = 0;
    //    }

    // check if there is at least 1 redshifts values
    if (redshifts.size() == 1) // consider this as a success
    {
        logPdf.resize(redshifts.size());
        logPdf[0] = 1.0;
        logEvidence = 1.0;
        return 0;
    } else if (redshifts.size() < 1)
    {
        return 2;
    }

    // check that the zPrior is size-compatible
    if (logZPrior.size() != redshifts.size())
    {
        return 4;
        // Float64 logPrior = log(1.0/redshifts.size()); //log(1.0)
    }
    if (verbose)
    {
        Float64 logZPriorMax = -DBL_MAX;
        Float64 logZPriorMin = DBL_MAX;
        for (UInt32 k = 0; k < redshifts.size(); k++)
        {
            if (logZPriorMax < logZPrior[k])
            {
                logZPriorMax = logZPrior[k];
            }
            if (logZPriorMin > logZPrior[k])
            {
                logZPriorMin = logZPrior[k];
            }
        }
        Log.LogDebug("Pdfz: Pdfz computation: using logZPrior min=%e", logZPriorMax);
        Log.LogDebug("Pdfz: Pdfz computation: using logZPrior max=%e", logZPriorMin);
    }

    // renormalize zprior to 1
    Float64 logsumZPrior = logSumExpTrick(logZPrior, redshifts, sumMethod);

    logPdf.resize(redshifts.size());
    TFloat64List Xi2_2withPrior;
    for(Int32 i = 0; i < merits.size(); i++){
        Xi2_2withPrior.push_back(-0.5*merits[i] + logZPrior[i] - logsumZPrior);
    }
    // prepare logLikelihood and LogEvidence
    Float64 logsumexp = logSumExpTrick( Xi2_2withPrior, redshifts, sumMethod);
    logEvidence = cstLog + logsumexp;

    if (verbose)
    {
        Log.LogDebug("Pdfz: Pdfz computation: using cstLog=%e", cstLog);
        Log.LogDebug("Pdfz: Pdfz computation: using logEvidence=%e", logEvidence);
        // Log.LogDebug("Pdfz: Pdfz computation: using logPrior=%f", logPrior);
        // //logPrior can be variable with z
    }

    for (UInt32 k = 0; k < redshifts.size(); k++)
    {
        logPdf[k] = Xi2_2withPrior[k] + cstLog - logEvidence;
    }

    if (verbose)
    {
        Float64 pdfmax = -DBL_MAX;
        Float64 pdfmin = DBL_MAX;
        for (UInt32 k = 0; k < redshifts.size(); k++)
        {
            if (pdfmax < logPdf[k])
            {
                pdfmax = logPdf[k];
            }
            if (pdfmin > logPdf[k])
            {
                pdfmin = logPdf[k];
            }
        }
        Log.LogDebug("Pdfz: Pdfz computation: found pdf min=%e", pdfmin);
        Log.LogDebug("Pdfz: Pdfz computation: found pdf max=%e", pdfmax);
    }

    return 0;
}

Float64 COperatorPdfz::getSumTrapez(const TRedshiftList & redshifts,
                            const TFloat64List & valprobalog)
{
    Float64 sum = 0.0;
    if(redshifts.size()==0)
    {
        return sum;
    }
    if(redshifts.size()!=valprobalog.size()) //this should raise an exception ? or return some error values ?
    {
        return sum;
    }

    // prepare LogEvidence
    Int32 sumMethod = 1;
    Float64 logSum = logSumExpTrick( valprobalog, redshifts, sumMethod);
    sum = exp(logSum);

    return sum;
}

Float64 COperatorPdfz::getSumRect(const TRedshiftList & redshifts,
                          const TFloat64List & valprobalog)
{
    Float64 sum = 0.0;
    // prepare LogEvidence
    Int32 sumMethod = 0;
    Float64 logSum = logSumExpTrick( valprobalog, redshifts, sumMethod);
    sum = exp(logSum);
    return sum;
}

Int32 COperatorPdfz::getIndex( std::vector<Float64> redshifts, Float64 z )
{
    Int32 solutionIdx=-1;
    for ( UInt32 i2=0; i2<redshifts.size(); i2++)
    {
        if( redshifts[i2]==z )
        {
            solutionIdx = i2;
            break;
        }
    }
    return solutionIdx;
}





/*Int32 COperatorPdfz::Marginalize(const TFloat64List & redshifts,
                         const std::vector<TFloat64List> & meritResults,
                         const std::vector<TFloat64List> & zPriors, const Float64 cstLog,
                         std::shared_ptr<CPdfMargZLogResult> postmargZResult,
                         const TFloat64List &modelPriors)*/
Int32 COperatorPdfz::Marginalize(const ChisquareArray & chisquarearray)
{
    bool verbose = false;
    
    const TFloat64List & redshifts = chisquarearray.redshifts;
    const std::vector<TFloat64List> & meritResults = chisquarearray.chisquares;
    const std::vector<TFloat64List> & zPriors = chisquarearray.zpriors;
    const Float64 & cstLog  = chisquarearray.cstLog;
    const TFloat64List & modelPriors = chisquarearray.modelpriors;

    if (meritResults.size() != zPriors.size())
    {
        Log.LogError("Pdfz: Pdfz marginalize problem. merit.size (%d) != prior.size (%d)", meritResults.size(), zPriors.size());
        return -9;
    }
    if (meritResults.size() < 1 || zPriors.size() < 1 || redshifts.size() < 1)
    {
        Log.LogError("Pdfz: Pdfz marginalize problem. merit.size (%d), prior.size (%d), or redshifts.size (%d) is zero !", meritResults.size(), zPriors.size(), redshifts.size());
        return -99;
    }

    // check merit curves. Maybe this should be assert stuff ?
    for (Int32 km = 0; km < meritResults.size(); km++)
    {
        TFloat64List _merit = meritResults[km];
        Bool invalidFound = false;
        for (Int32 kz = 0; kz < _merit.size(); kz++)
        {
            if (_merit[kz] != _merit[kz])
            {
                Log.LogError("    COperatorPdfz::Marginalize - merit result #%d has at least one nan or invalid value at index=%d", km, kz);
                invalidFound = true;
                break;
            }
        }
    }

    std::vector<UInt32> nSum(redshifts.size(), 0);

    Float64 MaxiLogEvidence = -DBL_MAX;
    TFloat64List LogEvidencesWPriorM;
    Float64 sumModifiedEvidences = 0.0;

    std::vector<Float64> logPriorModel;
    if (/*false &&*/ modelPriors.size() != meritResults.size()){
    
        Float64 priorModelCst = 1.0 / Float64(meritResults.size());
        Log.LogInfo("Pdfz: Marginalize: no priors loaded, using constant priors (=%f)", priorModelCst);
        for (Int32 km = 0; km < meritResults.size(); km++)
        {
            logPriorModel.push_back(log(priorModelCst));
        }
    } else
    { //we need to check if modelPriors is a const vector passed from linemodelsolve.combinePDF

        for (Int32 km = 0; km < meritResults.size(); km++)
        {
            logPriorModel.push_back(log(modelPriors[km]));
        }

        // logging priors used
        Float64 sumPriors = 0.0;
        for (Int32 km = 0; km < meritResults.size(); km++)
        {
            sumPriors += exp(logPriorModel[km]);
            Log.LogDetail("Pdfz: Marginalize: for model k=%d, using prior=%f", km, exp(logPriorModel[km]));
        }
        Log.LogInfo("Pdfz: Marginalize: sumPriors=%f", sumPriors);
        if (sumPriors > 1.1 || sumPriors < 0.9)
        {
            Log.LogError("Pdfz: sumPriors should be close to 1... !!!");
            throw runtime_error("Pdfz: sumPriors should be close to 1");
        }
    }

    std::vector<TFloat64List> logProbaList;
    TFloat64List logEvidenceList;
    for (Int32 km = 0; km < meritResults.size(); km++)
    {
        // Todo: Check if the status is OK ?
        // meritResult->Status[i] == COperator::nStatus_OK

        logProbaList.emplace_back();
        TFloat64List & logProba = logProbaList.back();
        Float64 logEvidence;
        
        Int32 retPdfz = ComputePdf(meritResults[km], redshifts, cstLog, zPriors[km], logProba, logEvidence); //here we are passing cte priors over all Z;
        if (retPdfz != 0)
        {
            Log.LogError("Pdfz: Pdfz computation - compute logEvidence: failed for result km=%d", km);
            throw runtime_error("Pdfz: Pdfz computation - compute logEvidence failed");
        } else
        {
            logEvidenceList.push_back(logEvidence);

            Float64 logEvidenceWPriorM = logEvidence + logPriorModel[km];
            LogEvidencesWPriorM.push_back(logEvidenceWPriorM);
            if (MaxiLogEvidence < logEvidenceWPriorM)
            {
                MaxiLogEvidence = logEvidenceWPriorM;
            }
        }
    }
    if (verbose)
    {
        Log.LogDebug("Pdfz: Marginalize: MaxiLogEvidence=%e", MaxiLogEvidence);
    }

    Float64 & logSumEvidence = m_postmargZResult->valEvidenceLog;

    // Using computational trick to sum the evidences
    for (Int32 k = 0; k < LogEvidencesWPriorM.size(); k++)
    {
        sumModifiedEvidences += exp(LogEvidencesWPriorM[k] - MaxiLogEvidence);
    }
    logSumEvidence = MaxiLogEvidence + log(sumModifiedEvidences); //here is the marginalized evidence, used for classification
    if (verbose)
    {
        Log.LogDebug("Pdfz: Marginalize: logSumEvidence=%e", logSumEvidence);
    }

    for (Int32 km = 0; km < meritResults.size(); km++)
    {
        if (verbose)
        {
            Log.LogDebug("Pdfz: Marginalize: processing chi2-result km=%d", km);
        }

        // Todo: Check if the status is OK ?
        // meritResult->Status[i] == COperator::nStatus_OK
        TFloat64List logProba;
        Float64 logEvidence;
        logProba = logProbaList[km];
        logEvidence = logEvidenceList[km];

        // check if the redshift bins are the same
        for (UInt32 k = 0; k < redshifts.size(); k++)
        {
            if (m_postmargZResult->Redshifts[k] != redshifts[k])
            {
                Log.LogError("Pdfz: Pdfz computation (z-bins comparison) failed for result km=%d", km);
                throw runtime_error("Pdfz: Pdfz computation (z-bins comparison) failed");
            }
        }

        for (UInt32 k = 0; k < redshifts.size(); k++)
        {
            if (true /*meritResult->Status[k]== COperator::nStatus_OK*/) // todo: check (temporarily considers status is always OK for linemodel tplshape)
            {
                Float64 & logValProba = m_postmargZResult->valProbaLog[k];
                Float64 logValProbaAdd = logProba[k] + logPriorModel[km] + logEvidence - logSumEvidence;
                Float64 maxP = logValProba;
                if (maxP < logValProbaAdd)
                {
                    maxP = logValProbaAdd;
                }
                Float64 valExp = exp(logValProba - maxP) + exp(logValProbaAdd - maxP);
                logValProba = maxP + log(valExp);
                nSum[k]++;
            }
        }

    }

    // THIS DOES NOT ALLOW Marginalization with coverage<100% for ALL templates
    for (UInt32 k = 0; k < m_postmargZResult->Redshifts.size(); k++)
    {
        if (nSum[k] != meritResults.size())
        {
            m_postmargZResult->valProbaLog[k] = NAN;
            Log.LogError("Pdfz: Pdfz computation failed. For z=%f, nSum=%d", m_postmargZResult->Redshifts[k], nSum[k]);
            Log.LogError("Pdfz: Pdfz computation failed. For z=%f, meritResults.size()=%d", m_postmargZResult->Redshifts[k], meritResults.size());
            Log.LogError("Pdfz: Pdfz computation failed. Not all templates have 100 percent coverage for all redshifts!");
            throw runtime_error("Pdfz: Pdfz computation failed. Not all templates have 100 percent coverage for all redshifts");
        }
    }

    return 0;
}

// This mathematically does not correspond to any valid method for combining
// PDFs.
// TODO: problem while estimating best proba. is it best proba for each z ? In
// that case: what about sum_z P = 1 ?
// TODO: this methid should be replaced/modified to correspond to the MaxPDF
// technique.
/*Int32 COperatorPdfz::BestProba(const TFloat64List & redshifts,
                       const std::vector<TFloat64List> & meritResults,
                       const std::vector<TFloat64List> & zPriors, const Float64 cstLog,
                       std::shared_ptr<CPdfMargZLogResult> m_postmargZResult)*/
Int32 COperatorPdfz::BestProba(const ChisquareArray & chisquarearray)
{
    Log.LogError("Pdfz: Pdfz-bestproba computation ! This method is currently not working !! It will produce bad results as is....");
    bool verbose = false;

    const TFloat64List & redshifts = chisquarearray.redshifts;
    const std::vector<TFloat64List> & meritResults = chisquarearray.chisquares;
    const std::vector<TFloat64List> & zPriors = chisquarearray.zpriors;
    const Float64 & cstLog  = chisquarearray.cstLog;

    if (meritResults.size() != zPriors.size())
    {
        Log.LogError("Pdfz: Pdfz-bestproba problem. merit.size (%d) != prior.size (%d)",
                     meritResults.size(), zPriors.size());
        throw std::runtime_error("Pdfz: Pdfz-bestproba problem, merit.size != prior.size" );
    }
    if (meritResults.size() < 1 || zPriors.size() < 1 || redshifts.size() < 1)
    {
        Log.LogError("Pdfz: Pdfz-bestproba problem. merit.size (%d), prior.size (%d), or redshifts.size (%d) is zero !",
                     meritResults.size(), zPriors.size(), redshifts.size());
        throw std::runtime_error("Pdfz: Pdfz-bestproba problem, merit.size, prior.size or redshifts.size is zero !");
    }

    Bool initPostMarg = false;  

    for (Int32 km = 0; km < meritResults.size(); km++)
    {
        if (verbose)
        {
            Log.LogDebug("Pdfz:-bestproba: processing chi2-result km=%d", km);
        }

        // Todo: Check if the status is OK ?
        // meritResult->Status[i] == COperator::nStatus_OK

        TFloat64List logProba;
        Float64 logEvidence;
        Int32 retPdfz = ComputePdf(meritResults[km], redshifts, cstLog, zPriors[km], logProba, logEvidence);
        if (retPdfz != 0)
        {
            Log.LogError("Pdfz: Pdfz-bestproba computation failed for result km=%d", km);
            throw std::runtime_error("Pdfz: Pdfz-bestproba: Pdfz computation failed");
        } else
        {
            // check if the redshift bins are the same
            for (UInt32 k = 0; k < redshifts.size(); k++)
            {
                if (m_postmargZResult->Redshifts[k] != redshifts[k])
                {
                    Log.LogError("Pdfz: Pdfz-bestproba, computation (z-bins comparison) failed for result km=%d", km);
                    throw std::runtime_error("Pdfz: Pdfz-bestproba, computation (z-bins comparison) failed");
                }
            }
            for (UInt32 k = 0; k < redshifts.size(); k++)
            {
                if (true /*meritResult->Status[k]== COperator::nStatus_OK*/) // todo: check (temporarily considers status is always OK for linemodel tplshape)
                {
                    if (logProba[k] > m_postmargZResult->valProbaLog[k])
                    {
                        m_postmargZResult->valProbaLog[k] = logProba[k];
                    }
                }
            }
        }
    }

    // normalize: sum_z P = 1
    // 1. check if the z step is constant. If not, pdf cannot be estimated by
    // the current method.
    Float64 reldzThreshold = 0.05; // relative difference accepted
    Float64 mindz = DBL_MAX;
    Float64 maxdz = -DBL_MAX;
    for (UInt32 k = 1; k < redshifts.size(); k++)
    {
        Float64 diff = redshifts[k] - redshifts[k - 1];
        if (mindz > diff)
        {
            mindz = diff;
        }
        if (maxdz < diff)
        {
            maxdz = diff;
        }
    }
    Float64 zstep = (maxdz + mindz) / 2.0;
    if (abs(maxdz - mindz) / zstep > reldzThreshold)
    {
        return 2;
    }

    // 2. prepare LogEvidence
    Float64 maxi = -DBL_MAX;
    std::vector<Float64> smallVALUES(redshifts.size(), 0.0);
    for (UInt32 k = 0; k < redshifts.size(); k++)
    {
        smallVALUES[k] = m_postmargZResult->valProbaLog[k];
        if (maxi < smallVALUES[k])
        {
            maxi = smallVALUES[k]; // maxi will be used to avoid underflows when
                                   // summing exponential of small values
        }
    }

    Float64 sumModifiedExp = 0.0;
    for (UInt32 k = 0; k < redshifts.size(); k++)
    {
        Float64 modifiedEXPO = exp(smallVALUES[k] - maxi);
        sumModifiedExp += modifiedEXPO;
    }
    Float64 logEvidence = maxi + log(sumModifiedExp) + log(zstep);

    if (verbose)
    {
        Log.LogDebug("Pdfz: Pdfz-bestproba computation: using logEvidence=%e", logEvidence);
        Log.LogDebug("Pdfz: Pdfz-bestproba computation: using log(zstep)=%e", log(zstep));
    }

    for (UInt32 k = 0; k < redshifts.size(); k++)
    {
        m_postmargZResult->valProbaLog[k] = m_postmargZResult->valProbaLog[k] - logEvidence;
    }

    return 0;
}

/**
 * @brief COperatorPdfz::BestChi2
 * Computes the combined pdf by taking the best chi2
 * WARNING: as long as the prior on the models are constant, it is equivalent to
 * compute the bestchi2 and the MaxPDF. If this prior is not constant any more,
 * the mas search has to be modified.
 * @param redshifts
 * @param meritResults
 * @param zPriors
 * @param cstLog
 * @param postmargZResult
 * @return
 */
/* Int32 COperatorPdfz::BestChi2(const TFloat64List & redshifts,
                      const std::vector<TFloat64List> & meritResults,
                      const std::vector<TFloat64List> & zPriors, const Float64 cstLog,
                      std::shared_ptr<CPdfMargZLogResult> postmargZResult) */ 
Int32 COperatorPdfz::BestChi2(const ChisquareArray & chisquarearray)
{
    bool verbose = false;

    const TFloat64List & redshifts = chisquarearray.redshifts;
    const std::vector<TFloat64List> & meritResults = chisquarearray.chisquares;
    const std::vector<TFloat64List> & zPriors = chisquarearray.zpriors;
    const Float64 & cstLog  = chisquarearray.cstLog;

    if (meritResults.size() != zPriors.size())
    {
        Log.LogError("Pdfz: Pdfz-bestchi2 problem, merit.size (%d) != prior.size (%d)",
                     meritResults.size(), zPriors.size());
        throw std::runtime_error("Pdfz: Pdfz-bestchi2 problem, merit.size != prior.size" );
    }
    if (meritResults.size() < 1 || zPriors.size() < 1 || redshifts.size() < 1)
    {
        Log.LogError("Pdfz: Pdfz-bestchi2 problem, merit.size (%d), prior.size (%d), or redshifts.size (%d) is zero !",
                     meritResults.size(), zPriors.size(), redshifts.size());
        throw std::runtime_error("Pdfz: Pdfz-bestchi2 problem, merit.size, prior.size or redshifts.size is zero !");
    }

    // build best chi2 vector
    TFloat64List chi2Min(redshifts.size(), DBL_MAX);
    for (Int32 km = 0; km < meritResults.size(); km++)
    {
        if (verbose)
        {
            Log.LogDebug("Pdfz:-bestchi2: processing chi2-result km=%d", km);
        }

        // Todo: Check if the status is OK ?
        // meritResult->Status[i] == COperator::nStatus_OK

        // Todo: use the priors for the min chi2 search ?
        for (UInt32 k = 0; k < redshifts.size(); k++)
        {
            if (true /*meritResult->Status[k]== COperator::nStatus_OK*/) // todo:
                                                                         // check
                                                                         // (temporarily
                                                                         // considers
                                                                         // status
                                                                         // is
                                                                         // always
                                                                         // OK)
            {
                if (meritResults[km][k] < chi2Min[k])
                {
                    chi2Min[k] = meritResults[km][k];
                }
            }
        }
    }

    // estimate Posterior on the best chi2
    CZPrior zpriorhelper;
    TFloat64List logProba;
    Float64 logEvidence;
    TFloat64List zprior;
    zprior = zpriorhelper.GetConstantLogZPrior(redshifts.size());

    Int32 retPdfz = ComputePdf(chi2Min, redshifts, cstLog, zprior, logProba, logEvidence);
    if (retPdfz != 0)
    {
        Log.LogError("Pdfz: Pdfz-bestchi2: Pdfz computation failed");
        throw std::runtime_error("Pdfz: Pdfz-bestchi2: Pdfz computation failed");
    } else
    {
        m_postmargZResult->valEvidenceLog = logEvidence;
        m_postmargZResult->valProbaLog = std::move(logProba);
    }

    return 0;
}
