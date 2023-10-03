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
#include "RedshiftLibrary/linemodel/templatesortho.h"
#include "RedshiftLibrary/linemodel/linemodelfitting.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/spectrum/logrebinning.h"

using namespace NSEpic;
void CTemplatesOrthogonalization::Orthogonalize(
    CInputContext &inputContext, const std::string category,
    std::shared_ptr<const CLSF> lsf) {
  // retrieve params from InputContext
  m_LSF = lsf;
  m_enableOrtho =
      inputContext.GetParameterStore()->EnableTemplateOrthogonalization(
          category);
  CAutoScope autoscope_category(Context.m_ScopeStack, category);
  CAutoScope autoscope_method(Context.m_ScopeStack, "LineModelSolve");

  // retrieve templateCatalog
  std::shared_ptr<CTemplateCatalog> tplCatalog =
      inputContext.GetTemplateCatalog();
  bool currentsampling = tplCatalog->m_logsampling;
  TBoolList samplingList{0, 1};

  bool needOrthogonalization = prepareTplCatForOrthogonalization(
      tplCatalog, inputContext, category, samplingList);

  if (!needOrthogonalization) {
    tplCatalog->m_logsampling = currentsampling;
    tplCatalog->m_orthogonal = 0;
    return;
  }

  for (bool sampling : samplingList) {
    tplCatalog->m_logsampling = sampling;
    // orthogonalize all templates
    tplCatalog->m_orthogonal = 0;
    const TTemplateConstRefList TplList =
        std::const_pointer_cast<const CTemplateCatalog>(tplCatalog)
            ->GetTemplateList(TStringList{category});
    tplCatalog->m_orthogonal = 1;
    bool partial = (sampling && tplCatalog->GetTemplateCount(category))
                       ? true
                       : false; // need to replace only some of the ortho
    for (Int32 i = 0; i < TplList.size(); i++) {
      const CTemplate tpl = *TplList[i];
      if (partial && tplCatalog->GetTemplate(category, i))
        break; // ortho template  is already there
      Log.LogDetail(Formatter()
                    << "    TplOrthogonalization: now processing tpl"
                    << tpl.GetName());
      std::shared_ptr<CTemplate> _orthoTpl = OrthogonalizeTemplate(tpl);
      if (partial)
        tplCatalog->SetTemplate(_orthoTpl, i);
      else
        tplCatalog->Add(_orthoTpl);
    }
  }
  tplCatalog->m_logsampling = currentsampling;
  tplCatalog->m_orthogonal = 0;
  return;
}

bool CTemplatesOrthogonalization::prepareTplCatForOrthogonalization(
    std::shared_ptr<CTemplateCatalog> &tplCatalog, CInputContext &inputContext,
    const std::string category, TBoolList &samplingList) const {

  // check if category has never been orthogonalized
  samplingList = {0, 1};
  bool first_time_ortho = !tplCatalog->GetTemplateCount(category, true, false);

  // check if LSF has changed, if yes reorthog all
  bool differentLSF = false;
  Float64 lambda = (inputContext.getLambdaRange()->GetBegin() +
                    inputContext.getLambdaRange()->GetEnd()) /
                   2;
  if (tplCatalog->m_ortho_LSFWidth !=
      m_LSF->GetWidth(lambda)) // true also if m_ortho_LSFWidth is NAN
  {
    differentLSF = true;
    tplCatalog->m_ortho_LSFWidth =
        m_LSF->GetWidth(lambda); // thus initialize if m_ortho_LSFWidth is NAN
  }

  if (first_time_ortho)
    return true; // no need to continue analysis since no more required setup,
                 // other than lsfwidth

  if (differentLSF) {
    tplCatalog->ClearTemplateList(category, 1,
                                  0); // clear orthog templates - non-rebinned
    tplCatalog->ClearTemplateList(category, 1,
                                  1); // clear orthog templates - rebinned
    return true; // no need to continue analysis since we cleared all
                 // templates
  }

  // check if log-sampled templates have changed and need ortho
  bool need_ortho_logsampling =
      hasLogRebinnedTemplatesChanged(tplCatalog, category);

  if (need_ortho_logsampling) // need ortho recompute,
    samplingList = {1};

  return need_ortho_logsampling;
}

bool CTemplatesOrthogonalization::hasLogRebinnedTemplatesChanged(
    std::shared_ptr<CTemplateCatalog> &tplCatalog, std::string category) const {

  tplCatalog->m_logsampling = 1;
  tplCatalog->m_orthogonal = 0; // orig log

  bool need_ortho_logsampling =
      tplCatalog->GetTemplateCount(category) > 0; // orig log list is not empty

  if (!need_ortho_logsampling)
    return false;
  // log sampling is there
  // check if logsampling has changed and thus need orthogonalization
  // has it changed for all or some templates ? (if changed, the orthog
  // template will be cleared)
  tplCatalog->m_orthogonal = 1;
  Int32 tpl_log_ortho_count = tplCatalog->GetTemplateCount(category);
  if (tpl_log_ortho_count > 0) // log-ortho list still there
  {
    // check if ortho missing for some templates only
    Int32 nonNullCount = tplCatalog->GetNonNullTemplateCount(category);
    if (nonNullCount == tpl_log_ortho_count)
      need_ortho_logsampling = false; // no need to recompute ortho
  }
  return need_ortho_logsampling;
}
/**
 * @brief getOrthogonalTemplate
 * Computes as follows:
 * - process linemodel on the spectrum
 * - creates the output template from the subtraction of the input spectrum
 * and the fitted linemodel
 * - add the newly created template to the tplCatalogOrtho member
 * @return
 */
std::shared_ptr<CTemplate> CTemplatesOrthogonalization::OrthogonalizeTemplate(
    const CTemplate &inputTemplate) {

  std::shared_ptr<CTemplate> tplOrtho =
      std::make_shared<CTemplate>(inputTemplate);

  if (!m_enableOrtho)
    return tplOrtho;

  std::string opt_continuumcomponent = "fromspectrum";
  Float64 opt_continuum_neg_threshold =
      -INFINITY; // not relevant in the "fromspectrum" case
  Float64 opt_continuum_nullamp_threshold =
      0.; // not relevant in the "fromspectrum" case;
  tplOrtho->SetLSF(m_LSF);

  std::string saveContinuumEstimationMethod =
      tplOrtho->GetContinuumEstimationMethod();
  tplOrtho->SetContinuumEstimationMethod("zero");

  // Compute linemodel on the template
  TLambdaRange lambdaRange = inputTemplate.GetLambdaRange();

  std::shared_ptr<COperatorTemplateFitting> TFOperator;
  CLineModelFitting model(tplOrtho, lambdaRange, TFOperator);

  Float64 redshift = 0.0;
  Float64 contreest_iterations = 0;
  bool enableLogging = true;
  CLineModelSolution modelSolution(Context.getCLineMap());
  CTplModelSolution continuumModelSolution;

  model.fit(redshift, modelSolution, continuumModelSolution,
            contreest_iterations, enableLogging);

  // Restore the continuum estimation method
  tplOrtho->SetContinuumEstimationMethod(saveContinuumEstimationMethod);

  // get mtm
  Float64 mtm = model.EstimateMTransposeM();

  // Subtract the fitted model from the original template
  model.getSpectrumModel().refreshModel();
  CSpectrum modelSpc = model.getSpectrumModel().GetModelSpectrum();
  /*//debug:
  FILE* f = fopen( "templatesortho_fittedmodel_dbg.txt", "w+" );
  for( Int32 t=0;t<modelSpc.GetSampleCount();t++)
  {
      fprintf( f, "%f %e\n", modelSpc.GetSpectralAxis()[t],
  modelSpc.GetFluxAxis()[t]);
  }
  fclose( f );
  //*/

  const CSpectrumFluxAxis &modelFluxAxis = modelSpc.GetFluxAxis();
  CSpectrumFluxAxis continuumOrthoFluxAxis = std::move(tplOrtho->GetFluxAxis());
  for (Int32 i = 0; i < continuumOrthoFluxAxis.GetSamplesCount(); i++) {
    continuumOrthoFluxAxis[i] -= modelFluxAxis[i];
  }
  tplOrtho->SetFluxAxis(std::move(continuumOrthoFluxAxis));

  /*//debug:
  FILE* f2 = fopen( "templatesortho_orthotemplate_dbg.txt", "w+" );
  for( Int32 t=0;t<modelSpc.GetSampleCount();t++)
  {
      fprintf( f2, "%f %e\n", modelSpc.GetSpectralAxis()[t],
  continuumOrthoFluxAxis[t]);
  }
  fclose( f2 );
  //*/

  return tplOrtho;
}
