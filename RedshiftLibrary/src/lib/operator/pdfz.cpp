// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/statistics/zprior.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/extremum/extremum.h"

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
                             bool allow_extrema_at_border,
                             UInt32 maxPeakCount_per_window,
                             const std::vector<TFloat64List> & candidatesRedshifts,
                             const TStringList & candidatesIds
                             ):
    m_opt_combine(opt_combine),
    m_peakSeparation(peakSeparation),
    m_meritcut(meritcut),
    m_allow_extrema_at_border(allow_extrema_at_border),
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
std::shared_ptr<PdfCandidatesZResult> COperatorPdfz::Compute(const ChisquareArray & chisquarearray,
                                                              bool integ)
{
    Log.LogInfo("%s: Pdfz computation", __func__);    

    // build PDF from chisquares and priors
    CombinePDF(chisquarearray);

    // find candidates redshifts
    TCandidateZbyID zcandidates = searchMaxPDFcandidates(); // will be sorted by the id syntax inside each redhisftwindow

    std::shared_ptr<PdfCandidatesZResult> CandidateszResult;
    if (integ){
        // compute pdf candidate properties (deltaz, integ, rank )
        CPdfCandidatesZ zcand_op =  CPdfCandidatesZ(zcandidates);
        CandidateszResult = zcand_op.Compute(m_postmargZResult->Redshifts, m_postmargZResult->valProbaLog);

        //eventually truncate candidates at maxcount
        size_t newsize = std::min(CandidateszResult->m_ranked_candidates.size(), size_t(m_maxCandidate));
        CandidateszResult->m_ranked_candidates.resize(newsize);
    }else{
        CandidateszResult = std::make_shared<PdfCandidatesZResult>();
            for (auto c:zcandidates){
                CandidateszResult->m_ranked_candidates.emplace_back(c);
            }
    }

    return CandidateszResult;
}



void COperatorPdfz::CombinePDF(const ChisquareArray & chisquarearray)
{
    Log.LogInfo("COperatorPdfz::CombinePDF: Pdfz combination");

    if(chisquarearray.chisquares.size()>0)
    {
        // initialize m_postmargZResult
        m_postmargZResult = std::make_shared<CPdfMargZLogResult>(chisquarearray.redshifts);

        if(m_opt_combine=="marg")
        {
            Log.LogInfo("COperatorPdfz::CombinePDF: Marginalization");
            Marginalize( chisquarearray);
            
        }else if(m_opt_combine=="bestchi2")
        {
            Log.LogInfo("COperatorPdfz::CombinePDF: BestChi2");
            BestChi2( chisquarearray);

        }else if(m_opt_combine=="bestproba"){
            Log.LogInfo("COperatorPdfz::CombinePDF: BestProba");
            BestProba( chisquarearray);

        }else{
            throw GlobalException(INTERNAL_ERROR,"COperatorPdfz::CombinePDF: Unable to parse pdf combination method option");

        }
    }else
    {
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"COperatorPdfz::CombinePDF: Unable to find any chisquares prepared for combination. chisquares.size()="<< chisquarearray.chisquares.size());
    }

    //check pdf is ok
    isPdfValid(); //will throw an error if not
}

bool COperatorPdfz::checkPdfSum() const
{
    bool ret = true;

    //check pdf sum=1
    Float64 sumTrapez = getSumTrapez(m_postmargZResult->Redshifts, m_postmargZResult->valProbaLog);
    Log.LogDetail("COperatorPdfz::checkPdfSum: Pdfz normalization - sum trapz. = %e", sumTrapez);
    if(abs(sumTrapez-1.0)>1e-1){
        Log.LogError("COperatorPdfz::checkPdfSum: Pdfz normalization failed (trapzesum = %f)", sumTrapez);
        ret = false;
    }

    return ret;
}


TCandidateZbyID COperatorPdfz::searchMaxPDFcandidates() const
{
    TCandidateZbyID candidates;
    
    for (const auto & cand : m_candidatesZRanges)
    {
        TPointList extremumList;
        const TFloat64Range & redshiftsRange = cand.second;
        std::string id = cand.first;

        //call Find on each secondpass range and retrieve the best  peak
        bool invertForMinSearch=false;
        CExtremum extremum_op = CExtremum( m_maxPeakCount_per_window, m_peakSeparation, m_meritcut, 
                                           invertForMinSearch, m_allow_extrema_at_border, redshiftsRange);
        bool findok = extremum_op.Find(m_postmargZResult->Redshifts, m_postmargZResult->valProbaLog, extremumList);
        if (!findok){
            if (m_candidatesZRanges.size() >1){ // we are in 2nd pass (several redshift ranges)
                Log.LogInfo("COperatorPdfz::searchMaxPDFcandidates: Second-pass fitting degenerates the first-pass results of candidate:%s in range [%f , %f]\n", 
                                cand.first.c_str(), redshiftsRange.GetBegin(), redshiftsRange.GetEnd());
                Log.LogInfo(" Flag - Eliminating a second-pass candidate");
                continue; 
            }
            else{
                throw GlobalException(INTERNAL_ERROR,"COperatorPdfz: searchMaxPDFcandidates failed");
            }       
        }
        Int32 i = 0 ;
        const std::string Id_prefix = (id=="" ? m_Id_prefix : id + "_" + m_Id_prefix);
        for (const auto & extremum : extremumList) // ranked by ValProba (from CExtremum::Find)
        {   
            std::string newid = Id_prefix + std::to_string(i++);
            candidates[newid].Redshift = extremum.X;
            candidates[newid].ValProba = extremum.Y;
            if (id!="") // this extrema has a parent id 
                candidates[newid].ParentId = id;   
        }
    }

    if (candidates.empty()){
        throw GlobalException(INTERNAL_ERROR,"COperatorPdfz: searchMaxPDFcandidates failed");
    }

    return candidates;
}



const std::shared_ptr<const CPdfMargZLogResult> COperatorPdfz::compressFirstpassPDF(Int32 ratio)
{
    std::shared_ptr<CPdfMargZLogResult> compressed_postmargZResult = std::make_shared<CPdfMargZLogResult>();
    Int32 s = (Int32)m_postmargZResult->Redshifts.size()/ratio;
    compressed_postmargZResult->Redshifts.reserve(s);
    compressed_postmargZResult->valProbaLog.reserve(s);
    for(Int32 i = 0;i<m_postmargZResult->Redshifts.size(); i+=ratio)
    {
        compressed_postmargZResult->Redshifts.push_back(m_postmargZResult->Redshifts[i]);
        compressed_postmargZResult->valProbaLog.push_back(m_postmargZResult->valProbaLog[i]);
    }
    return compressed_postmargZResult;
}

/**
 * Use the log sum exp trick to sum up small numbers while avoiding underflows
 *   * ----------------------------------------------------------------------
     * NOTE: this uses the LOG-SUM-EXP trick originally suggested by S. Jamal
     * ----------------------------------------------------------------------
*/
Float64 COperatorPdfz::logSumExpTrick(const TFloat64List & valproba, const TFloat64List & redshifts)
{

    Float64 logfactor = -DBL_MAX;
    if(redshifts.size()<2){
        throw GlobalException(INTERNAL_ERROR,"COperatorPdfz::logSumExpTrick Can't compute on a range of less than 2 points");
    }
    
    for (UInt32 k = 0; k < redshifts.size(); k++)
    {
       
        Float64 zstep;
        if (k == 0)
        {
            zstep = (redshifts[k + 1] - redshifts[k]) * 0.5;
        } else if (k == redshifts.size() - 1)
        {
            zstep = (redshifts[k] - redshifts[k - 1]) * 0.5;
        } else
        {
            zstep = (redshifts[k + 1] - redshifts[k - 1]) * 0.5;
        }
        if (logfactor < valproba[k] + log(zstep))
        {
            logfactor = valproba[k] + log(zstep); // maxi will be used to avoid underflows when
                                                  // summing exponential of small values
        }
    }
    
    Log.LogDebug("COperatorPdfz::logSumExpTrick: using common factor value for log-sum-exp trick=%e", logfactor);

    Float64 sumModifiedExp = 0.0;
        Float64 modifiedEXPO_previous = exp(valproba[0] - logfactor);
    for (UInt32 k = 1; k < redshifts.size(); k++)
    {
        Float64 modifiedEXPO = exp(valproba[k] - logfactor);
        Float64 trapezArea = (modifiedEXPO + modifiedEXPO_previous) / 2.0;
        trapezArea *= (redshifts[k] - redshifts[k - 1]);
        sumModifiedExp += trapezArea;
        modifiedEXPO_previous = modifiedEXPO;
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
void COperatorPdfz::ComputePdf(const TFloat64List & merits, const TFloat64List & redshifts,
                                        const Float64 cstLog, const TFloat64List & logZPrior,
                                        TFloat64List &logPdf, Float64 &logEvidence)
{
    bool verbose = true;
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
        Log.LogDebug("COperatorPdfz::ComputePdf: using merit min=%e", meritmin);
        Log.LogDebug("COperatorPdfz::ComputePdf: using merit max=%e", meritmax);
    }

    // check if there is at least 1 redshift value
    if (redshifts.size() == 1) // consider this as a success
    {
        logPdf.resize(redshifts.size());
        logPdf[0] = 1.0;
        logEvidence = 1.0;
        return ;
    } else if (redshifts.size() < 1)
    {
        throw GlobalException(INTERNAL_ERROR,"COperatorPdfz::ComputePdf, redshifts is empty");
    }

    // check that the zPrior is size-compatible
    if (logZPrior.size() != redshifts.size())
    {
        throw GlobalException(INTERNAL_ERROR,"COperatorPdfz::ComputePdf, redshifts and logZPrior have different sizes");
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
        Log.LogDebug("COperatorPdfz::ComputePdf: using logZPrior min=%e", logZPriorMax);
        Log.LogDebug("COperatorPdfz::ComputePdf: using logZPrior max=%e", logZPriorMin);
    }

    // renormalize zprior to 1
    Float64 logsumZPrior = logSumExpTrick(logZPrior, redshifts);

    logPdf.resize(redshifts.size());
    TFloat64List Xi2_2withPrior;
    for(Int32 i = 0; i < merits.size(); i++){
        Xi2_2withPrior.push_back(-0.5*merits[i] + logZPrior[i] - logsumZPrior);
    }
    // prepare logLikelihood and LogEvidence
    Float64 logsumexp = logSumExpTrick( Xi2_2withPrior, redshifts);
    logEvidence = cstLog + logsumexp;

    if (verbose)
    {
        Log.LogDebug("COperatorPdfz::ComputePdf: using cstLog=%e", cstLog);
        Log.LogDebug("COperatorPdfz::ComputePdf: using logEvidence=%e", logEvidence);
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
        Log.LogDebug("COperatorPdfz::ComputePdf: found pdf min=%e", pdfmin);
        Log.LogDebug("COperatorPdfz::ComputePdf: found pdf max=%e", pdfmax);
    }
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
    Float64 logSum = logSumExpTrick( valprobalog, redshifts);
    sum = exp(logSum);

    return sum;
}

Int32 COperatorPdfz::getIndex( const std::vector<Float64> & redshifts, Float64 z )
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

void COperatorPdfz::ComputeEvidenceAll(const TFloat64List &LogEvidencesWPriorM, Float64 &MaxiLogEvidence)
{
    MaxiLogEvidence = *std::max_element(LogEvidencesWPriorM.cbegin(), LogEvidencesWPriorM.cend());

    Float64 & logSumEvidence = m_postmargZResult->valMargEvidenceLog;

    // Using computational trick to sum the evidences
    Float64 sumModifiedEvidences = 0.0;
    for (auto & logEv:LogEvidencesWPriorM)
        sumModifiedEvidences += exp(logEv - MaxiLogEvidence);
    logSumEvidence = MaxiLogEvidence + log(sumModifiedEvidences); //here is the marginalized evidence, used for classification
}

void COperatorPdfz::ComputeAllPdfs( const ChisquareArray & chisquarearray, 
                                    std::vector<TFloat64List> &logProbaList,
                                    TFloat64List &LogEvidencesWPriorM,
                                    TFloat64List &logPriorModel,
                                    Float64 &MaxiLogEvidence)
{
    const TFloat64List & redshifts = chisquarearray.redshifts;
    const std::vector<TFloat64List> & meritResults = chisquarearray.chisquares;
    const std::vector<TFloat64List> & zPriors = chisquarearray.zpriors;
    const Float64 & cstLog  = chisquarearray.cstLog;
    const TFloat64List & modelPriors = chisquarearray.modelpriors;

    if (meritResults.size() != zPriors.size())
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"COperatorPdfz::ComputeAllPdfs: merit.size ("<<meritResults.size()<<") != prior.size ("<<zPriors.size()<<")");
    
    if (meritResults.size() < 1 || zPriors.size() < 1 || redshifts.size() < 1)
        throw GlobalException(INTERNAL_ERROR,Formatter()<<"COperatorPdfz::ComputeAllPdfs: merit.size("<<meritResults.size()<<"), prior.size("<<zPriors.size()<<") or redshifts.size("<<redshifts.size()<<") is zero !");

    // check merit curves. Maybe this should be assert stuff ?
    for (const TFloat64List & _merit : meritResults)
        for (const Float64 & m : _merit)
            if (m != m) // test NAN value
                throw GlobalException(INTERNAL_ERROR,"COperatorPdfz::ComputeAllPdfs - merit result has at least one nan value");

    if (modelPriors.empty()){
        const Float64 priorModelCst = 1.0 / Float64(meritResults.size());
        Log.LogInfo("COperatorPdfz::ComputeAllPdfs: no priors loaded, using constant priors (=%f)", priorModelCst);
        logPriorModel = TFloat64List(meritResults.size(), log(priorModelCst));

    } else {
        if (modelPriors.size() != meritResults.size())
            throw GlobalException(INTERNAL_ERROR, "COperatorPdfz::ComputeAllPdfs: modelPriors has wrong size");

        logPriorModel.resize(meritResults.size());
        std::transform(modelPriors.cbegin(), modelPriors.cend(), logPriorModel.begin(), 
            [](Float64 v){
                return log(v);
                }
            );

        // logging priors used
        Float64 sumPriors = 0.0;
        for (Int32 km = 0; km < meritResults.size(); km++)
        {
            sumPriors += exp(logPriorModel[km]);
            Log.LogDetail("COperatorPdfz::ComputeAllPdfs: for model k=%d, using prior=%f", km, exp(logPriorModel[km]));
        }

        Log.LogInfo("COperatorPdfz::ComputeAllPdfs: sumPriors=%f", sumPriors);
        if (sumPriors > 1.1 || sumPriors < 0.9)
            throw GlobalException(INTERNAL_ERROR,"Pdfz::ComputeAllPdfs: sumPriors should be close to 1... !!!");

    }

    TFloat64List logEvidenceList(meritResults.size());
    MaxiLogEvidence = -DBL_MAX;
    logProbaList.resize(meritResults.size());
    LogEvidencesWPriorM.resize(meritResults.size());
    for (Int32 km = 0; km < meritResults.size(); km++)
        ComputePdf(meritResults[km], redshifts, cstLog, zPriors[km], logProbaList[km], logEvidenceList[km] ); 
    
    std::transform( logEvidenceList.cbegin(), logEvidenceList.cend(), 
                    logPriorModel.cbegin(), LogEvidencesWPriorM.begin(), std::plus<Float64>());
    
    ComputeEvidenceAll(LogEvidencesWPriorM, MaxiLogEvidence);
}

void COperatorPdfz::Marginalize(const ChisquareArray & chisquarearray)
{
    bool verbose = false;
    
    const auto nmodel = chisquarearray.chisquares.size();
    const TFloat64List & redshifts = chisquarearray.redshifts;
    const auto zsize = redshifts.size();

    std::vector<TFloat64List> logProbaList;
    TFloat64List LogEvidencesWPriorM;
    TFloat64List logPriorModel;
    Float64 MaxiLogEvidence;
    ComputeAllPdfs(chisquarearray, logProbaList, LogEvidencesWPriorM, logPriorModel, MaxiLogEvidence);

    if (verbose)
        Log.LogDebug("COperatorPdfz::Marginalize: MaxiLogEvidence=%e", MaxiLogEvidence);

    m_postmargZResult->valEvidenceLog = m_postmargZResult->valMargEvidenceLog;

    if (verbose)
        Log.LogDebug("COperatorPdfz::Marginalize: logSumEvidence=%e", m_postmargZResult->valMargEvidenceLog);

    // marginalize: ie sum all PDFS
    std::vector<UInt32> nSum(zsize, 0);
    for (Int32 km = 0; km < nmodel; km++)
    {
        if (verbose)
            Log.LogDebug("COperatorPdfz::Marginalize: processing chi2-result km=%d", km);

        // Todo: Check if the status is OK ?
        // meritResult->Status[i] == COperator::nStatus_OK
        const TFloat64List &logProba = logProbaList[km];
        const Float64 logWeight = LogEvidencesWPriorM[km] - m_postmargZResult->valMargEvidenceLog;

        // check if the redshift bins are the same
        if (m_postmargZResult->Redshifts != redshifts)
            throw GlobalException(INTERNAL_ERROR, "COperatorPdfz::Marginalize z-bins comparison failed");

        for (UInt32 k = 0; k < zsize; k++)
        {
            Float64 & logValProba = m_postmargZResult->valProbaLog[k];
            const Float64 logValProbaAdd = logProba[k] + logWeight;
            const Float64 maxP = std::max(logValProbaAdd,logValProba);
            const Float64 valExp = exp(logValProba - maxP) + exp(logValProbaAdd - maxP);
            logValProba = maxP + log(valExp);
            nSum[k]++;
        }

    }

    // THIS DOES NOT ALLOW Marginalization with coverage<100% for ALL templates
    for (UInt32 k = 0; k < zsize; k++)
    {
        if (nSum[k] != nmodel)
        {
            m_postmargZResult->valProbaLog[k] = NAN;
            Log.LogError("Pdfz: Pdfz computation failed. For z=%f, nSum=%d", m_postmargZResult->Redshifts[k], nSum[k]);
            Log.LogError("Pdfz: Pdfz computation failed. For z=%f, meritResults.size()=%d", m_postmargZResult->Redshifts[k], nmodel);
            throw GlobalException(INTERNAL_ERROR,"COperatorPdfz::Marginalize computation failed. Not all templates have 100 percent coverage for all redshifts!");
        }
    }
}

// This mathematically does not correspond to any valid method for combining
// PDFs.
// TODO: problem while estimating best proba. is it best proba for each z ? In
// that case: what about sum_z P = 1 ?
// TODO: this method should be replaced/modified to correspond to the MaxPDF
// technique.
void COperatorPdfz::BestProba(const ChisquareArray & chisquarearray)
{
    Log.LogError("Pdfz: Pdfz-bestproba computation ! This method is currently not working !! It will produce bad results as is....");
    bool verbose = false;

    const TFloat64List & redshifts = chisquarearray.redshifts;
    const std::vector<TFloat64List> & meritResults = chisquarearray.chisquares;
    const std::vector<TFloat64List> & zPriors = chisquarearray.zpriors;
    const Float64 & cstLog  = chisquarearray.cstLog;

    if (meritResults.size() != zPriors.size())
    {
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"COperatorPdfz: Pdfz-bestproba problem, merit.size ("<<meritResults.size()<<") != prior.size ("<<zPriors.size()<<")");
    }
    if (meritResults.size() < 1 || zPriors.size() < 1 || redshifts.size() < 1)
    {
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"COperatorPdfz: Pdfz-bestproba problem, merit.size("<<meritResults.size()<<"), prior.size("<<zPriors.size()<<") or redshifts.size("<<redshifts.size()<<") is zero !");
    }

    for (Int32 km = 0; km < meritResults.size(); km++)
    {
        if (verbose)
        {
            Log.LogDebug("COperatorPdfz::BestProba: processing chi2-result km=%d", km);
        }

        // Todo: Check if the status is OK ?
        // meritResult->Status[i] == COperator::nStatus_OK

        TFloat64List logProba;
        Float64 logEvidence;
        ComputePdf(meritResults[km], redshifts, cstLog, zPriors[km], logProba, logEvidence);

        // check if the redshift bins are the same
        for (UInt32 k = 0; k < redshifts.size(); k++)
        {
            if (m_postmargZResult->Redshifts[k] != redshifts[k])
            {
		  throw GlobalException(INTERNAL_ERROR,Formatter()<<"Pdfz: Pdfz-bestproba, computation (z-bins comparison) failed for result km="<< km);
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
        throw GlobalException(INTERNAL_ERROR,"COperatorPdfz::BestProba: zstep is not constant, cannot normalize");
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
        Log.LogDebug("COperatorPdfz::BestProba: using logEvidence=%e", logEvidence);
        Log.LogDebug("COperatorPdfz::BestProba: using log(zstep)=%e", log(zstep));
    }

    for (UInt32 k = 0; k < redshifts.size(); k++)
    {
        m_postmargZResult->valProbaLog[k] = m_postmargZResult->valProbaLog[k] - logEvidence;
    }
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
void COperatorPdfz::BestChi2(const ChisquareArray & chisquarearray)
{
    bool verbose = false;

    const auto nmodel = chisquarearray.chisquares.size();
    const TFloat64List & redshifts = chisquarearray.redshifts;
    const auto zsize = redshifts.size();
    const std::vector<TFloat64List> & meritResults = chisquarearray.chisquares;

    std::vector<TFloat64List> logProbaList;
    TFloat64List LogEvidencesWPriorM;
    TFloat64List logPriorModel;
    Float64 MaxiLogEvidence;
    ComputeAllPdfs(chisquarearray, logProbaList, LogEvidencesWPriorM, logPriorModel, MaxiLogEvidence);

/*    // build min chi2 vector
    // build instead maxlikelihood vector including priors: ie pdf x evidence
    TFloat64List maxLikelihood(zsize, -DBL_MAX);
    TFloat64List likelihood(zsize);
    for (Int32 km = 0; km < nmodel; km++)
    {
        if (verbose)
            Log.LogDebug("COperatorPdfz::BestChi2: processing chi2-result km=%d", km);

        // rescale the pdf by the evidence
        const Float64 &logev = LogEvidencesWPriorM[km];
        std::transform(logProbaList[km].cbegin(), logProbaList[km].cend(), likelihood.begin(), 
            [logev](Float64 logpdfval){return logpdfval+logev; });

        for (UInt32 k = 0; k < zsize; k++)
            if (likelihood[k]>maxLikelihood[k]) maxLikelihood[k] = likelihood[k];
    }
*/

    // build best chi2 vector
    TFloat64List chi2Min(zsize, DBL_MAX);
    for (Int32 km = 0; km < nmodel; km++)
        // Todo: use the priors for the min chi2 search ?
        for (UInt32 k = 0; k < zsize; k++)
            if (meritResults[km][k] < chi2Min[k])
                chi2Min[k] = meritResults[km][k];
                
    // estimate Posterior on the best chi2
    CZPrior zpriorhelper;
    TFloat64List zprior = zpriorhelper.GetConstantLogZPrior(zsize);

    ComputePdf(chi2Min, redshifts, chisquarearray.cstLog, zprior, m_postmargZResult->valProbaLog, m_postmargZResult->valEvidenceLog);



}

/**
 * @brief isPdfValid
 * @return
 */
void COperatorPdfz::isPdfValid() const
{
    if(!m_postmargZResult)
    {
      throw GlobalException(INTERNAL_ERROR,"COperatorPdfz::isPdfValid: PDF ptr is null");
    }

    if(m_postmargZResult->Redshifts.size()<2)
    {
      throw GlobalException(INTERNAL_ERROR,"COperatorPdfz::isPdfValid: PDF has size less than 2");
    }

    //is it completely flat ?
    Float64 minVal=DBL_MAX;
    Float64 maxVal=-DBL_MAX;
    for(Int32 k=0; k<m_postmargZResult->valProbaLog.size(); k++)
    {
        if(m_postmargZResult->valProbaLog[k]<minVal)
        {
            minVal = m_postmargZResult->valProbaLog[k];
        }
        if(m_postmargZResult->valProbaLog[k]>maxVal)
        {
            maxVal = m_postmargZResult->valProbaLog[k];
        }
    }
    if(minVal==maxVal){
        throw GlobalException(INTERNAL_ERROR,"COperatorPdfz::isPdfValid: PDF is flat !");
    }

    //is pdf any value nan ?
    for(Int32 k=0; k<m_postmargZResult->valProbaLog.size(); k++)
    {
        if(m_postmargZResult->valProbaLog[k] != m_postmargZResult->valProbaLog[k])
        {
            throw GlobalException(INTERNAL_ERROR,"COperatorPdfz::isPdfValid: PDF has nan or invalid values !");
        }
    }

    // is sum equal to 1
    if (!checkPdfSum()){
        throw GlobalException(INTERNAL_ERROR,"COperatorPdfz::isPdfValid: Pdfz normalization failed");
    };
}
