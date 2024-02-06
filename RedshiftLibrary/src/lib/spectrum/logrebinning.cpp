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
#include "RedshiftLibrary/spectrum/logrebinning.h"
#include "RedshiftLibrary/log/log.h"

namespace bfs = boost::filesystem;

using namespace NSEpic;
using namespace std;

CSpectrumLogRebinning::CSpectrumLogRebinning(CInputContext &inputContext)
    : m_inputContext(inputContext) {
  m_logGridStep = m_inputContext.getLogGridStep();
  std::shared_ptr<CSpectrum> spc;
  if (inputContext.GetSpectrum()->GetSpectralAxis().IsLogSampled()) {
    spc =
        m_inputContext
            .GetRebinnedSpectrum(); // retrieve the corrected rebinned spectrum
  } else {
    spc = m_inputContext.GetSpectrum();
  }
  setupRebinning(*spc, *(m_inputContext.getLambdaRange()));
}

/**
 * Brief: Get loglambdastep and update the zrange accordingly
 * Below code relies on the fact that both loglambda grid and the
 log(Redshift+1) grid follows the same arithmetic progession with a common step
 * if spectrum is already rebinned, then it imposes the rebinning and the
 creation of zGrid
 * Otherwise, it's the input zrange that decides on the rebinning param.
 * Note : If we want the log step to be log(1.0+redshiftstep) then spreadOverlog
 should be modified to use 1+redshiftstep for the common ratio (or construct the
 grid using arithmetic log progression).
*/
void CSpectrumLogRebinning::setupRebinning(CSpectrum &spectrum,
                                           const TFloat64Range &lambdaRange) {
  if (spectrum.GetSpectralAxis().IsLogSampled(m_logGridStep)) {
    // compute reference lambda range
    // (the effective lambda range of log-sampled spectrum when initial spectrum
    // overlaps lambdaRange assuming all spectra are aligned (this to avoid
    // rebining the template several times)
    TLambdaRange lambda_range_spc = spectrum.GetLambdaRange();
    Float64 loglambda_start_spc = log(lambda_range_spc.GetBegin());
    Float64 loglambda_end_spc = log(lambda_range_spc.GetEnd());
    Float64 loglambda_start_ref = log(lambdaRange.GetBegin());
    Float64 loglambda_end_ref = log(lambdaRange.GetEnd());
    loglambda_start_ref =
        loglambda_start_spc +
        ceil((loglambda_start_ref - loglambda_start_spc) / m_logGridStep) *
            m_logGridStep;
    loglambda_end_ref =
        loglambda_end_spc +
        floor((loglambda_end_ref - loglambda_end_spc) / m_logGridStep) *
            m_logGridStep;
    m_lambdaRange_ref = TFloat64Range(
        exp(loglambda_start_ref),
        exp(loglambda_end_ref)); // this is the effective lambda range, to be
                                 // used in inferTemplateRebinningSetup
  } else if (!spectrum.GetSpectralAxis().IsLogSampled()) {
    // compute reference lambda range (the effective lambda range of rebinned
    // spectrum when initial spectrum overlaps lambdaRange)
    Int32 loglambda_count_ref;
    Float64 loglambda_start_ref = log(lambdaRange.GetBegin());
    Float64 loglambda_end_ref = log(lambdaRange.GetEnd());
    loglambda_count_ref = Int32(
        floor((loglambda_end_ref - loglambda_start_ref) /
              m_logGridStep)); // we should sample inside range hence floor
    loglambda_end_ref =
        loglambda_start_ref + loglambda_count_ref * m_logGridStep;
    m_lambdaRange_ref = TFloat64Range(
        lambdaRange.GetBegin(),
        exp(loglambda_end_ref)); // this is the effective lambda range, to be
                                 // used in inferTemplateRebinningSetup
  } else {
    THROWG(INTERNAL_ERROR,
           Formatter() << "Log-sampled spectrum has wrong logGridStep : "
                       << m_logGridStep);
  }
  Log.LogDetail("  Log-Rebin: logGridStep = %f", m_logGridStep);
  return;
}

/**
 *  Rebin the spectrum with the calculated logGridStep if spectrum not already
 * rebinned: step1: construct the spectralAxis step2: do the rebin
 */
std::shared_ptr<CSpectrum> CSpectrumLogRebinning::loglambdaRebinSpectrum(
    CSpectrum const &spectrum, std::string const &errorRebinMethod) const {
  TFloat64Range lambdaRange_spc;
  Int32 loglambda_count_spc;
  if (spectrum.GetSpectralAxis().IsLogSampled()) {
    THROWG(INTERNAL_ERROR, Formatter() << "spectrum is already log-sampled");
  } else {
    // compute rebinned spectrum lambda range lambdaRange_spc (ie clamp on
    // reference grid) to be passed to computeTargetLogSpectralAxis in
    // loglambdaRebinSpectrum
    spectrum.GetSpectralAxis().ClampLambdaRange(m_lambdaRange_ref,
                                                lambdaRange_spc);
    Float64 loglambda_start_spc = log(lambdaRange_spc.GetBegin());
    Float64 loglambda_end_spc = log(lambdaRange_spc.GetEnd());
    Float64 loglambda_start_ref = log(m_lambdaRange_ref.GetBegin());
    Float64 loglambda_end_ref = log(m_lambdaRange_ref.GetEnd());
    if (lambdaRange_spc.GetBegin() > m_lambdaRange_ref.GetBegin()) {
      loglambda_start_spc =
          loglambda_start_ref +
          ceil((loglambda_start_spc - loglambda_start_ref) / m_logGridStep) *
              m_logGridStep; // ceil to be bigger or equal the first sample
      lambdaRange_spc.SetBegin(exp(loglambda_start_spc));
    }
    if (lambdaRange_spc.GetEnd() < m_lambdaRange_ref.GetEnd()) {
      loglambda_end_spc =
          loglambda_end_ref -
          ceil((loglambda_end_ref - loglambda_end_spc) / m_logGridStep) *
              m_logGridStep; // ceil to be less or equal the last sample
      lambdaRange_spc.SetEnd(exp(loglambda_end_spc));
    }
    Float64 count_ = (loglambda_end_spc - loglambda_start_spc) /
                     m_logGridStep; // should integer at numerical precision...
    loglambda_count_spc = round(count_) + 1;

    if (loglambda_count_spc < 2) {
      THROWG(INTERNAL_ERROR, Formatter()
                                 << "Cannot rebin spectrum of a grid of size:  "
                                 << loglambda_count_spc << "<2");
    }
  }

  // prepare return rebinned vector
  auto spectrumRebinedLog = make_shared<CSpectrum>(spectrum.GetName());
  CMask mskRebinedLog;

  const CSpectrumSpectralAxis targetSpectralAxis =
      computeTargetLogSpectralAxis(lambdaRange_spc, loglambda_count_spc);

  TFloat64Range spcLbdaRange(targetSpectralAxis[0] - 0.5 * m_logGridStep,
                             targetSpectralAxis[loglambda_count_spc - 1] +
                                 0.5 * m_logGridStep);

  // rebin the spectrum
  spectrum.setRebinInterpMethod(m_rebinMethod);
  spectrum.Rebin(spcLbdaRange, targetSpectralAxis, *spectrumRebinedLog,
                 mskRebinedLog, errorRebinMethod);

  spectrumRebinedLog->GetSpectralAxis().IsLogSampled(
      m_logGridStep); // double make sure that sampling is well done

  Log.LogDetail("  Log-regular lambda resampling FINISHED");

  return spectrumRebinedLog;
}

/*
    Aims at computing the template lambda range once for all template catalog
*/
Int32 CSpectrumLogRebinning::inferTemplateRebinningSetup(
    const TFloat64Range &zrange, TFloat64Range &lambdaRange_tpl) const {
  Float64 loglbdamin =
      log(m_lambdaRange_ref.GetBegin() / (1.0 + zrange.GetEnd()));
  Float64 loglbdamax =
      log(m_lambdaRange_ref.GetEnd() / (1.0 + zrange.GetBegin()));
  Int32 _round = std::round((loglbdamax - loglbdamin) / m_logGridStep) + 1;
  Float64 _neat =
      (loglbdamax - loglbdamin) / m_logGridStep +
      1; // we expect to get an int value with no need to any rounding
  if (std::abs(_round - _neat) > 1E-8) {
    THROWG(INTERNAL_ERROR, "Problem in logrebinning setup");
  }
  Int32 loglambda_count_tpl =
      std::round((loglbdamax - loglbdamin) / m_logGridStep) + 1;

  Float64 tgt_loglbdamax = loglbdamax;
  Float64 tgt_loglbdamin =
      loglbdamax - (loglambda_count_tpl - 1) * m_logGridStep;
  lambdaRange_tpl = TFloat64Range(exp(tgt_loglbdamin), exp(tgt_loglbdamax));
  Log.LogDetail("  Operator-TemplateFittingLog: Log-Rebin: tpl raw "
                "loglbdamin=%f : raw loglbdamax=%f",
                loglbdamin, loglbdamax);
  Log.LogDetail(
      "  Operator-TemplateFittingLog: zmin_new = %f, tpl->lbdamax = %f",
      zrange.GetBegin(), exp(loglbdamax));
  Log.LogDetail(
      "  Operator-TemplateFittingLog: zmax_new = %f, tpl->lbdamin = %f",
      zrange.GetEnd(), exp(loglbdamin));
  return loglambda_count_tpl;
}
/**
 *
 * Brief: Log Rebin the template spectral axis
 * Important: below function considered that we have already constructed the new
 *log spectralAxis of our spectrum Step1: align the rebined spectrum max lambda
 *at min redshift. Step2: get grid count for target template and update size
 *accordingly Step3: construct target loglambda axis for the template and check
 *borders Step4: rebin the template --> rebinned flux is saved in
 *templateRebinedLog
 **/
std::shared_ptr<CTemplate> CSpectrumLogRebinning::loglambdaRebinTemplate(
    std::shared_ptr<const CTemplate> tpl, TFloat64Range &lambdaRange_tpl,
    const Int32 loglambda_count_tpl) const {
  Log.LogInfo("  Operator-TemplateFittingLog: Log-regular lambda resampling "
              "START for template %s",
              tpl->GetName().c_str());
  // check template coverage is enough for zrange and spectrum coverage
  bool overlapFull = true;
  if (lambdaRange_tpl.GetBegin() < tpl->GetSpectralAxis()[0])
    overlapFull = false;
  if (lambdaRange_tpl.GetEnd() >
      tpl->GetSpectralAxis()[tpl->GetSampleCount() - 1])
    overlapFull = false;
  if (!overlapFull) {
    THROWG(INTERNAL_ERROR,
           Formatter() << "overlap found to be lower than 1.0 for template "
                       << tpl->GetName().c_str());
  }

  const CSpectrumSpectralAxis targetSpectralAxis =
      computeTargetLogSpectralAxis(lambdaRange_tpl, loglambda_count_tpl);

  if (targetSpectralAxis[loglambda_count_tpl - 1] < targetSpectralAxis[0])
    THROWG(INTERNAL_ERROR,
           " Last element of the target spectral axis is not valid. Template "
           "count is not well computed due to exp/conversions");

  auto templateRebinedLog =
      make_shared<CTemplate>(tpl->GetName(), tpl->GetCategory());
  templateRebinedLog->m_ismCorrectionCalzetti = tpl->m_ismCorrectionCalzetti;
  templateRebinedLog->m_igmCorrectionMeiksin = tpl->m_igmCorrectionMeiksin;
  CMask mskRebinedLog;

  TFloat64Range tplLbdaRange(targetSpectralAxis[0] - 0.5 * m_logGridStep,
                             targetSpectralAxis[loglambda_count_tpl - 1] +
                                 0.5 * m_logGridStep);

  tpl->setRebinInterpMethod(m_rebinMethod);
  tpl->Rebin(tplLbdaRange, targetSpectralAxis, *templateRebinedLog,
             mskRebinedLog);

  templateRebinedLog->GetSpectralAxis().IsLogSampled(m_logGridStep);

  return templateRebinedLog;
}

CSpectrumSpectralAxis CSpectrumLogRebinning::computeTargetLogSpectralAxis(
    const TFloat64Range &lambdarange,
    Int32 count) const { // spreadoverlog expects m_Begin to be non-log value
  TFloat64List axis = lambdarange.SpreadOverLogEpsilon(m_logGridStep);
  if (axis.size() != count) {
    THROWG(INTERNAL_ERROR,
           "computed axis does not have the expected samples number");
  }
  CSpectrumSpectralAxis targetSpectralAxis(std::move(axis));
  return targetSpectralAxis;
}

bool CSpectrumLogRebinning::checkTemplateAlignment(
    const std::shared_ptr<const CTemplate> &tpl,
    const TFloat64Range &lambdaRange_tpl) const {
  const TAxisSampleList &w = tpl->GetSpectralAxis().GetSamplesVector();
  const Float64 &lstart = lambdaRange_tpl.GetBegin();
  Int32 idx = CIndexing<Float64>::getCloserIndex(w, lstart);
  return (lstart - w[idx]) / w[idx] <= 2E-7;
}

TFloat64Range CSpectrumLogRebinning::logRebinTemplateCatalog(
    const std::string &spectrumModel) const {
  TFloat64Range redshiftRange =
      m_inputContext.GetParameterStore()->Get<TFloat64Range>(spectrumModel +
                                                             ".redshiftRange");
  Int32 SSratio = 1;
  if (m_inputContext.GetParameterStore()->Get<std::string>(
          spectrumModel + ".redshiftSolver.method") == "lineModelSolve" &&
      m_inputContext.GetParameterStore()->Has<Int32>(
          spectrumModel + ".redshiftSolver.lineModelSolve.lineModel.firstPass."
                          "largeGridStepRatio"))
    SSratio = m_inputContext.GetParameterStore()->Get<Int32>(
        spectrumModel + ".redshiftSolver.lineModelSolve.lineModel.firstPass."
                        "largeGridStepRatio");

  // compute the effective zrange of the new redshift grid
  // set the min to the initial min
  // set the max to an interger number of log(z+1) steps
  Float64 zmin_new = redshiftRange.GetBegin(),
          zmax_new = redshiftRange.GetEnd();
  {
    Float64 log_zmin_new_p1 = log(zmin_new + 1.);
    Float64 log_zmax_new_p1 = log(zmax_new + 1.);
    Int32 nb_z = Int32(
        ceil((log_zmax_new_p1 - log_zmin_new_p1) / m_logGridStep / SSratio));
    zmax_new = exp(log_zmin_new_p1 + nb_z * m_logGridStep * SSratio) - 1.;
  }
  TFloat64Range zrange(zmin_new,
                       zmax_new); // updating zrange based on the new logstep
  TFloat64Range lambdaRange_tpl;
  Int32 loglambda_count_tpl =
      inferTemplateRebinningSetup(zrange, lambdaRange_tpl);
  // rebin templates using previously identified parameters,
  //  rebin only if rebinning parameters are different from previously used ones
  std::shared_ptr<CTemplateCatalog> tplcat =
      m_inputContext.GetTemplateCatalog();
  const TStringList categoryList = tplcat->GetCategoryList();
  tplcat->m_logsampling = true;

  // check existence of already  & correctly logsampled templates
  const Int32 ntpl = tplcat->GetTemplateCount(spectrumModel);
  if (ntpl > 0) {
    for (Int32 i = 0; i < ntpl; i++) {
      std::shared_ptr<const CTemplate> tpl =
          tplcat->GetTemplate(spectrumModel, i);
      bool needrebinning = isRebinningNeeded(tpl, lambdaRange_tpl);
      if (!needrebinning)
        continue;

      {
        Log.LogDetail(
            " CInputContext::RebinInputs: need to rebin again the template: %s",
            tpl->GetName().c_str());
        tplcat->m_logsampling = false;
        std::shared_ptr<const CTemplate> input_tpl = tplcat->GetTemplateByName(
            TStringList{spectrumModel}, tpl->GetName());
        tplcat->m_logsampling = true;
        tplcat->SetTemplate(
            loglambdaRebinTemplate(input_tpl, lambdaRange_tpl,
                                   loglambda_count_tpl),
            i); // assigin the tpl pointer to a new rebined template
      }
    }
  } else {
    // no rebined templates in the category: rebin all templates
    tplcat->m_logsampling = false;
    const TTemplateConstRefList TplList =
        std::const_pointer_cast<const CTemplateCatalog>(tplcat)
            ->GetTemplateList(TStringList{spectrumModel});
    tplcat->m_logsampling = true;
    for (auto tpl : TplList) {
      std::shared_ptr<CTemplate> rebinnedTpl =
          loglambdaRebinTemplate(tpl, lambdaRange_tpl, loglambda_count_tpl);
      tplcat->Add(rebinnedTpl);
    }
  }
  tplcat->m_logsampling = false;

  return zrange;
}

bool CSpectrumLogRebinning::isRebinningNeeded(
    const std::shared_ptr<const CTemplate> &tpl,
    const TFloat64Range &lambdaRange_tpl) const {
  bool needrebinning = false;
  if (!tpl->GetSpectralAxis().IsLogSampled(m_logGridStep))
    needrebinning = true;
  if (lambdaRange_tpl.GetBegin() < tpl->GetSpectralAxis()[0])
    needrebinning = true;
  if (lambdaRange_tpl.GetEnd() >
      tpl->GetSpectralAxis()[tpl->GetSampleCount() - 1])
    needrebinning = true;
  if (!checkTemplateAlignment(tpl, lambdaRange_tpl))
    needrebinning = true;
  return needrebinning;
};
