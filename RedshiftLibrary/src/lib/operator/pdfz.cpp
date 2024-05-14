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
#include <algorithm>
#include <cfloat>
#include <fstream>
#include <numeric>

#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/linemodel/linemodelextremaresult.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/statistics/zprior.h"

using namespace NSEpic;
using namespace std;

COperatorPdfz::COperatorPdfz(const std::string &opt_combine,
                             Float64 peakSeparation, Float64 meritcut,
                             Int32 maxCandidate, bool redshiftLogSampling,
                             const std::string &Id_prefix,
                             bool allow_extrema_at_border,
                             Int32 maxPeakCount_per_window)
    : m_opt_combine(opt_combine), m_peakSeparation(peakSeparation),
      m_meritcut(meritcut), m_allow_extrema_at_border(allow_extrema_at_border),
      m_maxPeakCount_per_window(maxPeakCount_per_window <= 0
                                    ? maxCandidate
                                    : maxPeakCount_per_window),
      m_maxCandidate(maxCandidate), m_redshiftLogSampling(redshiftLogSampling),
      m_Id_prefix(Id_prefix) {}

/*
 * Main Pdf operator entrance
 *  combine pdf and search for candidates.
 */
std::shared_ptr<PdfCandidatesZResult>
COperatorPdfz::Compute(const ChisquareArray &chisquarearray) {
  Log.LogInfo(Formatter() << __func__ << ": Pdfz computation");

  if (chisquarearray.zgridParams.size() !=
      chisquarearray.parentCandidates.size())
    THROWG(INTERNAL_ERROR, "candidates zgridParams and parentCandidates do not "
                           "have the same size");

  const bool isFirstPass = chisquarearray.parentCandidates.empty();

  // Set m_parentcandidates & m_candidateszrange depending on 1st / 2nd pass
  if (isFirstPass) {
    m_parentCandidates = {{"", nullptr}};
    m_candidatesZRanges = {{"", {}}};
  } else {
    // second pass
    m_parentCandidates =
        TCandidateZbyID(chisquarearray.parentCandidates.cbegin(),
                        chisquarearray.parentCandidates.cend());
    for (Int32 i = 0; i < chisquarearray.zgridParams.size(); ++i)
      m_candidatesZRanges[chisquarearray.parentCandidates[i].first] =
          TFloat64Range(chisquarearray.zgridParams[i].zmin,
                        chisquarearray.zgridParams[i].zmax);
  }

  // create ClogZPdfResult
  createPdfResult(chisquarearray);

  // build PDF from chisquares and priors
  CombinePDF(chisquarearray);

  // find candidates redshifts
  TCandidateZbyID zcandidates =
      searchMaxPDFcandidates(); // will be sorted by the id syntax inside each
                                // redhisftwindow

  CPdfCandidatesZ zcand_op = CPdfCandidatesZ(zcandidates);

  std::shared_ptr<PdfCandidatesZResult> CandidateszResult;

  // compute pdf candidate properties (deltaz)
  CandidateszResult = zcand_op.Compute(m_postmargZResult->redshifts,
                                       m_postmargZResult->valProbaLog);

  // eventually truncate candidates at maxcount
  size_t newsize = std::min(CandidateszResult->m_ranked_candidates.size(),
                            size_t(m_maxCandidate));
  CandidateszResult->m_ranked_candidates.resize(newsize);

  // Sets rank
  for (int r = 0; r < CandidateszResult->size(); r++)
    CandidateszResult->m_ranked_candidates[r].second->Rank = r;

  // Checks window size
  TFloat64Range window_range;
  for (int r = 0; r < CandidateszResult->size(); r++) {
    auto rth_ranked_candidate = CandidateszResult->m_ranked_candidates[r];
    TFloat64Range integration_range{
        rth_ranked_candidate.second->ValSumProbaZmin,
        rth_ranked_candidate.second->ValSumProbaZmax};

    window_range =
        isFirstPass
            ? TFloat64Range{m_postmargZResult->redshifts.front(),
                            m_postmargZResult->redshifts.back()}
            : m_candidatesZRanges.at(rth_ranked_candidate.second->ParentId);
    checkWindowSize(integration_range, window_range);
  }

  return CandidateszResult;
}

void COperatorPdfz::checkWindowSize(const TFloat64Range &integration_range,
                                    const TFloat64Range &window_range) {

  if (integration_range.GetBegin() < window_range.GetBegin() ||
      integration_range.GetEnd() > window_range.GetEnd()) {
    Flag.warning(WarningCode::PDF_INTEGRATION_WINDOW_TOO_SMALL,
                 Formatter()
                     << "COperatorPdfz::" << __func__
                     << " Window is too small : \n"
                     << "integration range: " << integration_range.GetBegin()
                     << ", " << integration_range.GetEnd()
                     << ", window range : [" << window_range.GetBegin() << ", "
                     << window_range.GetEnd());
  }
}

void COperatorPdfz::createPdfResult(const ChisquareArray &chisquarearray) {

  // initialize m_postmargZResult
  TZGridListParams zparams(chisquarearray.zgridParams.size() + 1);
  zparams[0] = CZGridParam(TFloat64Range(chisquarearray.redshifts),
                           chisquarearray.zstep);
  std::copy(chisquarearray.zgridParams.cbegin(),
            chisquarearray.zgridParams.cend(), zparams.begin() + 1);
  m_postmargZResult =
      std::make_shared<CLogZPdfResult>(zparams, m_redshiftLogSampling);

  // check if the redshift grids are the same
  if (m_postmargZResult->redshifts != chisquarearray.redshifts)
    THROWG(INTERNAL_ERROR, "z-grid comparison failed");

  m_postmargZResult->valProbaLog =
      TFloat64List(chisquarearray.redshifts.size(), -DBL_MAX);
}

void COperatorPdfz::CombinePDF(const ChisquareArray &chisquarearray) {
  Log.LogInfo("COperatorPdfz::CombinePDF: Pdfz combination");
  if (!chisquarearray.chisquares.size()) {
    THROWG(INTERNAL_ERROR, Formatter() << "chisquarearray is empty");
  }

  if (m_opt_combine == "marg") {
    Log.LogInfo("COperatorPdfz::CombinePDF: Marginalization");
    Marginalize(chisquarearray);

  } else if (m_opt_combine == "bestChi2") {
    Log.LogInfo("COperatorPdfz::CombinePDF: BestChi2");
    BestChi2(chisquarearray);

  } else if (m_opt_combine == "bestproba") {
    Log.LogInfo("COperatorPdfz::CombinePDF: BestProba");
    BestProba(chisquarearray);

  } else {
    THROWG(INTERNAL_ERROR,
           Formatter() << "Invalid pdfcombination method: " << m_opt_combine);
  }

  // check pdf is ok
  if (!m_postmargZResult)
    THROWG(INTERNAL_ERROR, " PDF ptr is null");

  m_postmargZResult->isPdfValid(); // will throw an error if not
}

TCandidateZbyID COperatorPdfz::searchMaxPDFcandidates() const {
  TCandidateZbyID candidates;
  const TFloat64List zgrid = m_postmargZResult->redshifts;
  for (const auto &cand : m_candidatesZRanges) {
    const TFloat64Range &redshiftsRange = cand.second;
    std::string id = cand.first;

    // call Find on each secondpass range and retrieve the best  peak
    bool invertForMinSearch = false;
    CExtremum extremum_op = CExtremum(
        m_maxPeakCount_per_window, m_peakSeparation, m_meritcut,
        invertForMinSearch, m_allow_extrema_at_border, redshiftsRange);
    bool isFirstPass = m_candidatesZRanges.size() == 1;
    TPointList extremumList =
        extremum_op.Find(zgrid, m_postmargZResult->valProbaLog, isFirstPass);
    if (!isFirstPass && extremumList.empty()) {
      // we are in 2nd pass (several redshift ranges) (error in Find if first
      // pass with not findok)
      Log.LogInfo(Formatter()
                  << "COperatorPdfz::searchMaxPDFcandidates: Second-pass "
                     "fitting degenerates the first-pass results of "
                     "candidate:"
                  << cand.first << " in range [" << redshiftsRange.GetBegin()
                  << " , " << redshiftsRange.GetEnd() << "]\n");
      Log.LogInfo(" Flag - Eliminating a second-pass candidate");
    }
    Int32 i = 0;
    const std::string Id_prefix =
        (id == "" ? m_Id_prefix : id + "_" + m_Id_prefix);

    // loop ranked by ValProba (from CExtremum::Find)
    for (const auto &extremum : extremumList) {
      std::string newid = Id_prefix + std::to_string(i++);
      candidates[newid] = std::make_shared<TCandidateZ>();
      candidates[newid]->Redshift = extremum.X;
      candidates[newid]->ValProba = extremum.Y;
      if (id != "") { // this extrema has a parent id
        candidates[newid]->ParentId = id;
        candidates[newid]->ParentObject = m_parentCandidates.at(id);
      }
    }
  }

  if (candidates.empty()) {
    THROWG(INTERNAL_ERROR, " Failed to identify pdf candidates");
  }

  return candidates;
}

/**
 * Use the log sum exp trick to sum up small numbers while avoiding underflows
 *   * ----------------------------------------------------------------------
 * NOTE: this uses the LOG-SUM-EXP trick originally suggested by S. Jamal
 * ----------------------------------------------------------------------
 */
Float64 COperatorPdfz::logSumExpTrick(const TFloat64List &valproba,
                                      const TFloat64List &redshifts) {

  Float64 logfactor = -DBL_MAX;
  if (redshifts.size() < 2) {
    THROWG(INTERNAL_ERROR, "Can't compute on a "
                           "range of less than 2 points");
  }

  for (Int32 k = 0; k < redshifts.size(); k++) {

    Float64 zstep;
    if (k == 0) {
      zstep = (redshifts[k + 1] - redshifts[k]) * 0.5;
    } else if (k == redshifts.size() - 1) {
      zstep = (redshifts[k] - redshifts[k - 1]) * 0.5;
    } else {
      zstep = (redshifts[k + 1] - redshifts[k - 1]) * 0.5;
    }
    if (logfactor < valproba[k] + log(zstep)) {
      logfactor =
          valproba[k] + log(zstep); // maxi will be used to avoid underflows
                                    // when summing exponential of small values
    }
  }

  Float64 sumModifiedExp = 0.0;
  Float64 modifiedEXPO_previous = exp(valproba[0] - logfactor);
  for (Int32 k = 1; k < redshifts.size(); k++) {
    Float64 modifiedEXPO = exp(valproba[k] - logfactor);
    Float64 trapezArea = (modifiedEXPO + modifiedEXPO_previous) / 2.0;
    trapezArea *= (redshifts[k] - redshifts[k - 1]);
    sumModifiedExp += trapezArea;
    modifiedEXPO_previous = modifiedEXPO;
  }
  if (sumModifiedExp == 0.0)
    THROWG(INTERNAL_ERROR, "sumModifiedExp cannot be null");
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
void COperatorPdfz::ComputePdf(const TFloat64List &merits,
                               const TFloat64List &redshifts,
                               const Float64 cstLog,
                               const TFloat64List &logZPrior,
                               TFloat64List &logPdf, Float64 &logEvidence) {
  logPdf.clear();
  if (redshifts.size() < 1)
    THROWG(INTERNAL_ERROR, "redshifts is empty");
  // check if there is at least 1 redshift value
  if (redshifts.size() == 1) // consider this as a success
  {
    logPdf.resize(redshifts.size());
    logPdf[0] = 1.0;
    logEvidence = 1.0;
    return;
  }

  // check that the zPrior is size-compatible
  if (logZPrior.size() != redshifts.size())
    THROWG(INTERNAL_ERROR, "Size do not match between redshifts and logZPrior");

  // renormalize zprior to 1
  Float64 logsumZPrior = logSumExpTrick(logZPrior, redshifts);

  logPdf.resize(redshifts.size());
  TFloat64List Xi2_2withPrior;
  for (Int32 i = 0; i < merits.size(); i++)
    Xi2_2withPrior.push_back(-0.5 * merits[i] + logZPrior[i] - logsumZPrior);

  // prepare logLikelihood and LogEvidence
  Float64 logsumexp = logSumExpTrick(Xi2_2withPrior, redshifts);
  logEvidence = cstLog + logsumexp;

  for (Int32 k = 0; k < redshifts.size(); k++)
    logPdf[k] = Xi2_2withPrior[k] + cstLog - logEvidence;
}

Int32 COperatorPdfz::getIndex(const TFloat64List &redshifts, Float64 z) {
  Int32 solutionIdx = -1;
  for (Int32 i2 = 0; i2 < redshifts.size(); i2++)
    if (redshifts[i2] == z) {
      solutionIdx = i2;
      break;
    }
  return solutionIdx;
}

void COperatorPdfz::ComputeEvidenceAll(const TFloat64List &LogEvidencesWPriorM,
                                       Float64 &MaxiLogEvidence) {
  MaxiLogEvidence = *std::max_element(LogEvidencesWPriorM.cbegin(),
                                      LogEvidencesWPriorM.cend());

  Float64 &logSumEvidence = m_postmargZResult->valMargEvidenceLog;

  // Using computational trick to sum the evidences
  Float64 sumModifiedEvidences = 0.0;
  for (auto &logEv : LogEvidencesWPriorM)
    sumModifiedEvidences += exp(logEv - MaxiLogEvidence);
  logSumEvidence =
      MaxiLogEvidence +
      log(sumModifiedEvidences); // here is the marginalized evidence, used for
                                 // classification
}

void COperatorPdfz::validateChisquareArray(
    const ChisquareArray &chisquarearray) const {

  const TFloat64List &redshifts = chisquarearray.redshifts;
  const std::vector<TFloat64List> &meritResults = chisquarearray.chisquares;
  const std::vector<TFloat64List> &zPriors = chisquarearray.zpriors;

  if (meritResults.size() != zPriors.size())
    THROWG(INTERNAL_ERROR, Formatter()
                               << " merit.size (" << meritResults.size()
                               << ") != prior.size (" << zPriors.size() << ")");

  if (meritResults.size() < 1 || zPriors.size() < 1 || redshifts.size() < 1)
    THROWG(INTERNAL_ERROR, Formatter() << "merit.size(" << meritResults.size()
                                       << "), prior.size(" << zPriors.size()
                                       << ") or redshifts.size("
                                       << redshifts.size() << ") is zero");

  // check merit curves. Maybe this should be assert stuff ?
  for (const TFloat64List &_merit : meritResults)
    for (const Float64 &m : _merit)
      if (m != m) // test NAN value
        THROWG(INTERNAL_ERROR, "merit result "
                               "has at least one nan value");
  return;
}

void COperatorPdfz::ComputeAllPdfs(const ChisquareArray &chisquarearray,
                                   std::vector<TFloat64List> &logProbaList,
                                   TFloat64List &LogEvidencesWPriorM,
                                   TFloat64List &logPriorModel,
                                   Float64 &MaxiLogEvidence) {
  const TFloat64List &redshifts = chisquarearray.redshifts;
  const std::vector<TFloat64List> &meritResults = chisquarearray.chisquares;
  const std::vector<TFloat64List> &zPriors = chisquarearray.zpriors;
  const Float64 &cstLog = chisquarearray.cstLog;
  const TFloat64List &modelPriors = chisquarearray.modelpriors;

  validateChisquareArray(chisquarearray);

  if (modelPriors.empty()) {
    const Float64 priorModelCst = 1.0 / Float64(meritResults.size());
    Log.LogInfo(Formatter()
                << "COperatorPdfz::ComputeAllPdfs: no priors loaded, using "
                   "constant priors (="
                << priorModelCst << ")");
    logPriorModel = TFloat64List(meritResults.size(), log(priorModelCst));
  } else {
    if (modelPriors.size() != meritResults.size())
      THROWG(INTERNAL_ERROR, "modelPriors has wrong size");

    logPriorModel.resize(meritResults.size());
    std::transform(modelPriors.cbegin(), modelPriors.cend(),
                   logPriorModel.begin(), [](Float64 v) { return log(v); });

    // summing priors used
    Float64 sumPriors = accumulate(
        logPriorModel.begin(), logPriorModel.end(), 0.0,
        [](Float64 sum, Float64 logprior) { return sum + exp(logprior); });

    Log.LogInfo(Formatter()
                << "COperatorPdfz::ComputeAllPdfs: sumPriors=" << sumPriors);
    if (sumPriors > 1.1 || sumPriors < 0.9)
      THROWG(INTERNAL_ERROR, "sumPriors should be between ]0.9;1.1[");
  }

  TFloat64List logEvidenceList(meritResults.size());
  MaxiLogEvidence = -DBL_MAX;
  logProbaList.resize(meritResults.size());
  LogEvidencesWPriorM.resize(meritResults.size());
  for (Int32 km = 0; km < meritResults.size(); km++)
    ComputePdf(meritResults[km], redshifts, cstLog, zPriors[km],
               logProbaList[km], logEvidenceList[km]);

  std::transform(logEvidenceList.cbegin(), logEvidenceList.cend(),
                 logPriorModel.cbegin(), LogEvidencesWPriorM.begin(),
                 std::plus<Float64>());

  ComputeEvidenceAll(LogEvidencesWPriorM, MaxiLogEvidence);
}

void COperatorPdfz::Marginalize(const ChisquareArray &chisquarearray) {

  const auto nmodel = chisquarearray.chisquares.size();
  const TFloat64List &redshifts = chisquarearray.redshifts;
  const auto zsize = redshifts.size();

  std::vector<TFloat64List> logProbaList;
  TFloat64List LogEvidencesWPriorM;
  TFloat64List logPriorModel;
  Float64 MaxiLogEvidence;
  ComputeAllPdfs(chisquarearray, logProbaList, LogEvidencesWPriorM,
                 logPriorModel, MaxiLogEvidence);

  Log.LogDebug(Formatter() << "COperatorPdfz::Marginalize: MaxiLogEvidence="
                           << MaxiLogEvidence);

  m_postmargZResult->valEvidenceLog = m_postmargZResult->valMargEvidenceLog;

  Log.LogDebug(Formatter() << "COperatorPdfz::Marginalize: logSumEvidence="
                           << m_postmargZResult->valMargEvidenceLog);

  // marginalize: ie sum all PDFS
  TInt32List nSum(zsize, 0);
  for (Int32 km = 0; km < nmodel; km++) {
    const TFloat64List &logProba = logProbaList[km];
    const Float64 logWeight =
        LogEvidencesWPriorM[km] - m_postmargZResult->valMargEvidenceLog;

    for (Int32 k = 0; k < zsize; k++) {
      Float64 &logValProba = m_postmargZResult->valProbaLog[k];
      const Float64 logValProbaAdd = logProba[k] + logWeight;
      const Float64 maxP = std::max(logValProbaAdd, logValProba);
      const Float64 valExp =
          exp(logValProba - maxP) + exp(logValProbaAdd - maxP);
      logValProba = maxP + log(valExp);
      nSum[k]++;
    }
  }

  // THIS DOES NOT ALLOW Marginalization with coverage<100% for ALL templates
  for (Int32 k = 0; k < zsize; k++) {
    if (nSum[k] != nmodel)
      THROWG(INTERNAL_ERROR, "Computation failed. Not all templates "
                             "have 100 percent coverage for all redshifts");
  }
}

// This mathematically does not correspond to any valid method for combining
// PDFs.
// TODO: problem while estimating best proba. is it best proba for each z ? In
// that case: what about sum_z P = 1 ?
// TODO: this method should be replaced/modified to correspond to the MaxPDF
// technique.
void COperatorPdfz::BestProba(const ChisquareArray &chisquarearray) {

  validateChisquareArray(chisquarearray);

  const TFloat64List &redshifts = chisquarearray.redshifts;
  const std::vector<TFloat64List> &meritResults = chisquarearray.chisquares;
  const std::vector<TFloat64List> &zPriors = chisquarearray.zpriors;
  const Float64 &cstLog = chisquarearray.cstLog;

  for (Int32 km = 0; km < meritResults.size(); km++) {
    Log.LogDebug(Formatter()
                 << "COperatorPdfz::BestProba: processing chi2-result km="
                 << km);

    TFloat64List logProba;
    Float64 logEvidence;
    ComputePdf(meritResults[km], redshifts, cstLog, zPriors[km], logProba,
               logEvidence);

    // check if the redshift bins are the same
    if (m_postmargZResult->redshifts != redshifts)
      THROWG(INTERNAL_ERROR, "z-bins comparison failed");

    for (Int32 k = 0; k < redshifts.size(); k++)
      if (true)
        m_postmargZResult->valProbaLog[k] =
            std::max(logProba[k], m_postmargZResult->valProbaLog[k]);
  }

  // normalize: sum_z P = 1
  // 1. check if the z step is constant. If not, pdf cannot be estimated by
  // the current method.
  Float64 reldzThreshold = 0.05; // relative difference accepted
  Float64 mindz = DBL_MAX;
  Float64 maxdz = -DBL_MAX;
  for (Int32 k = 1; k < redshifts.size(); k++) {
    Float64 diff = redshifts[k] - redshifts[k - 1];
    mindz = std::min(mindz, diff);
    maxdz = std::max(maxdz, diff);
  }
  Float64 zstep = (maxdz + mindz) / 2.0;
  if (abs(maxdz - mindz) / zstep > reldzThreshold)
    THROWG(INTERNAL_ERROR, "Impossible to normalize: zstep is not constant");

  // 2. prepare LogEvidence
  Float64 maxi = -DBL_MAX;
  TFloat64List smallVALUES(redshifts.size(), 0.0);
  for (Int32 k = 0; k < redshifts.size(); k++) {
    smallVALUES[k] = m_postmargZResult->valProbaLog[k];
    maxi = std::max(maxi,
                    smallVALUES[k]); // maxi will be used to avoid underflows
                                     // when summing exponential of small values
  }

  Float64 sumModifiedExp = accumulate(
      smallVALUES.begin(), smallVALUES.end(), 0.0,
      [maxi](Float64 sum, Float64 value) { return sum + exp(value - maxi); });

  Float64 logEvidence = maxi + log(sumModifiedExp) + log(zstep);

  Log.LogDebug(Formatter() << "COperatorPdfz::BestProba: using logEvidence="
                           << logEvidence);
  Log.LogDebug(Formatter() << "COperatorPdfz::BestProba: using log(zstep)="
                           << log(zstep));

  std::transform(m_postmargZResult->valProbaLog.cbegin(),
                 m_postmargZResult->valProbaLog.cend(),
                 m_postmargZResult->valProbaLog.begin(),
                 [logEvidence](Float64 v) { return v - logEvidence; });
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
void COperatorPdfz::BestChi2(const ChisquareArray &chisquarearray) {

  const auto nmodel = chisquarearray.chisquares.size();
  const TFloat64List &redshifts = chisquarearray.redshifts;
  const auto zsize = redshifts.size();
  const std::vector<TFloat64List> &meritResults = chisquarearray.chisquares;

  std::vector<TFloat64List> logProbaList;
  TFloat64List LogEvidencesWPriorM;
  TFloat64List logPriorModel;
  Float64 MaxiLogEvidence;
  ComputeAllPdfs(chisquarearray, logProbaList, LogEvidencesWPriorM,
                 logPriorModel, MaxiLogEvidence);

  /*    // build min chi2 vector
      // build instead maxlikelihood vector including priors: ie pdf x evidence
      TFloat64List maxLikelihood(zsize, -DBL_MAX);
      TFloat64List likelihood(zsize);
      for (Int32 km = 0; km < nmodel; km++)
      {
          if (verbose)

          // rescale the pdf by the evidence
          const Float64 &logev = LogEvidencesWPriorM[km];
          std::transform(logProbaList[km].cbegin(), logProbaList[km].cend(),
     likelihood.begin(), [logev](Float64 logpdfval){return logpdfval+logev; });

          for (Int32 k = 0; k < zsize; k++)
              if (likelihood[k]>maxLikelihood[k]) maxLikelihood[k] =
     likelihood[k];
      }
  */

  // build best chi2 vector
  TFloat64List chi2Min(zsize, DBL_MAX);
  for (Int32 km = 0; km < nmodel; km++)
    // Todo: use the priors for the min chi2 search ?
    for (Int32 k = 0; k < zsize; k++)
      if (meritResults[km][k] < chi2Min[k])
        chi2Min[k] = meritResults[km][k];

  // estimate Posterior on the best chi2
  CZPrior zpriorhelper;
  TFloat64List zprior = zpriorhelper.GetConstantLogZPrior(zsize);

  ComputePdf(chi2Min, redshifts, chisquarearray.cstLog, zprior,
             m_postmargZResult->valProbaLog, m_postmargZResult->valEvidenceLog);
}
