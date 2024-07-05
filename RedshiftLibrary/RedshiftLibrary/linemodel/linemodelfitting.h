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
#ifndef _REDSHIFT_LINEMODEL_FITTING_
#define _REDSHIFT_LINEMODEL_FITTING_

#include <memory>
#include <unordered_set>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/line/regulament.h"
#include "RedshiftLibrary/linemodel/abstractfitter.h"
#include "RedshiftLibrary/linemodel/continuummanager.h"
#include "RedshiftLibrary/linemodel/element.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/linemodel/obsiterator.h"
#include "RedshiftLibrary/linemodel/spectrummodel.h"
#include "RedshiftLibrary/operator/continuumfitting.h"
#include "RedshiftLibrary/operator/linemodelresult.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

namespace NSEpic {

class CLineRatioManager;

class CLineModelFitting {

public:
  CLineModelFitting(
      const std::shared_ptr<COperatorContinuumFitting>
          &continuumFittingOperator,
      ElementComposition element_composition = ElementComposition::Default);
  CLineModelFitting(
      const std::shared_ptr<const CSpectrum> &template_,
      const TLambdaRange &lambdaRange,
      const std::shared_ptr<COperatorContinuumFitting>
          &continuumFittingOperator); // only used for template
                                      // orthogonalization, TODO use only one of
                                      // the future subclasses ? at least
                                      // inherit from clinemodelfitting

  void initParameters();
  void initMembers(const std::shared_ptr<COperatorContinuumFitting>
                       &continuumFittingOperator,
                   ElementComposition element_composition);

  void LogCatalogInfos();

  void setRedshift(Float64 redshift, bool reinterpolatedContinuum = false);
  void setContinuumComponent(std::string component);
  const std::string &getContinuumComponent() const {
    return m_continuumManager->getContinuumComponent();
  };

  bool initDtd();
  Float64 EstimateDTransposeD(const std::string &spcComponent) const;
  Float64 EstimateMTransposeM() const;
  Float64 EstimateLikelihoodCstLog() const;
  Float64 getDTransposeD();
  Float64 getLikelihood_cstLog();

  void SetVelocityEmission(Float64 vel);
  void SetVelocityAbsorption(Float64 vel);
  void setVelocityAbsorptionByGroup(Float64 vel, const TInt32List &inds);
  void setVelocityEmissionByGroup(Float64 vel, const TInt32List &inds);

  Float64 GetVelocityEmission() const;
  Float64 GetVelocityAbsorption() const;

  Float64 fit(Float64 redshift, CLineModelSolution &modelSolution,
              CContinuumModelSolution &continuumModelSolution,
              Int32 contreest_iterations = 0, bool enableLogging = 0);
  TFloat64Range &getDTDLambdaRange() { return m_dTransposeDLambdaRange; };

  void SetFittingMethod(const std::string &fitMethod,
                        bool enableAmplitudeOffsets = false,
                        bool enableLambdaOffsetsFit = false);
  void setLineRatioType(const std::string &lineratio);
  void SetAbsLinesLimit(Float64 limit);

  Float64 GetRedshift() const;

  CMask getOutsideLinesMask() const;
  Float64 getOutsideLinesSTD(Int32 which) const;

  Int32 computeSpcNSamples() const;

  Float64 getLeastSquareContinuumMeritFast() const;
  Float64 getLeastSquareContinuumMerit() const;

  Float64 getScaleMargCorrection(Int32 Eltidx = undefIdx) const {
    return m_ElementsVector->getScaleMargCorrection(Eltidx);
  }

  std::unordered_set<std::string> getLinesAboveSNR(Float64 snrcut = 3.5) const {
    // TODO temp basic impl
    m_spectraIndex.reset();
    return getSpectrumModel().getLinesAboveSNR(getLambdaRange(), snrcut);
  }

  std::pair<Float64, Float64> getCumulSNRStrongEL() const;
  std::pair<Float64, Float64> getSNROnRange(TInt32Range idxRange) const;

  void LoadModelSolution(const CLineModelSolution &modelSolution);
  CLineModelSolution GetModelSolution(Int32 opt_level = 0);

  Float64 getModelFluxVal(Int32 idx) const;
  void logParameters();

  Int32 setPassMode(Int32 iPass);
  Int32 GetPassNumber() const;

  void prepareAndLoadContinuum(Int32 icontfitting, Float64 redshift);
  void computeSpectrumFluxWithoutContinuum();
  bool isContinuumComponentTplfitxx() const {
    return m_continuumManager->isContinuumComponentTplfitxx();
  }

  bool isContinuumComponentPowerLaw() const {
    return m_continuumManager->isContinuumComponentPowerLaw();
  }

  Int32 GetModelNonZeroElementsNDdl() const {
    return m_ElementsVector->GetModelNonZeroElementsNDdl();
  }

  const CSpectrum &getSpectrum() const {
    return *((*m_inputSpcs).at(m_spectraIndex.get()));
  }
  shared_ptr<const CSpectrum> getSpectrumPtr() {
    return (*m_inputSpcs).at(m_spectraIndex.get());
  }

  CSpectrumModel &getSpectrumModel() { return m_models->getSpectrumModel(); }

  const CSpectrumModel &getSpectrumModel() const {
    return m_models->getSpectrumModel();
  }

  CLineModelElementList &getElementList() {
    return m_ElementsVector->getElementList();
  }
  const CLineModelElementList &getElementList() const {
    return m_ElementsVector->getElementList();
  }
  std::vector<TLineModelElementParam_ptr> &getElementParam() {
    return m_ElementsVector->getElementParam();
  }

  const TLambdaRange &getLambdaRange() const {
    return *(m_lambdaRanges.at(m_spectraIndex.get()));
  }

  const std::string &getFittingMethod() const { return m_fittingmethod; }

  std::shared_ptr<CContinuumManager> getContinuumManager() {
    return m_continuumManager;
  }
  std::shared_ptr<const CContinuumModelSolution> getContinuumFitValues() const {
    return m_continuumFitValues;
  }

  // we keep that getters temporarily, waiting for refactoring linemodel extrema
  // result
  Int32 getTplratio_count() const;
  TFloat64List getTplratio_priors() const;

  std::string const &getLineRatioType() const { return m_lineRatioType; }

  void loadFitContinuumParameters(Int32 icontinuum, Float64 redshift);

  std::shared_ptr<CAbstractFitter> m_fitter;
  std::shared_ptr<CLineRatioManager> m_lineRatioManager;

  const CLineMap m_RestLineList;

  Int32 m_pass = 1;
  bool m_enableAmplitudeOffsets;
  bool m_enableLbdaOffsets;

  Float64 m_LambdaOffsetMin = -400.0;
  Float64 m_LambdaOffsetMax = 400.0;
  Float64 m_LambdaOffsetStep = 25.0;

  // Multi obs combination/aggregation methods on elements Lists

  std::pair<Float64, Float64>
  GetMeanContinuumUnderLine(Int32 eltIdx, Int32 line_index, Float64 redshift);

  std::pair<Float64, Float64>
  getFluxDirectIntegration(const TInt32List &eIdx_list,
                           const TInt32List &subeIdx_list,
                           bool substract_abslinesmodel) const;

  CSpectraGlobalIndex &getSpectraIndex() { return m_spectraIndex; }
  void refreshAllModels();

private:
  void SetLSF();

  void applyPolynomCoeffs(Int32 eIdx, const TPolynomCoeffs &polynom_coeffs);

  std::shared_ptr<CContinuumModelSolution> m_continuumFitValues;
  std::shared_ptr<CContinuumManager> m_continuumManager;

  CSpcModelVectorPtr m_models;

  Float64 m_dTransposeD; // the cached dtd (maximum chisquare value)
  TFloat64Range m_dTransposeDLambdaRange; // the lambdaRange used to computed
                                          // cached dTransposeD values
  Float64 m_likelihood_cstLog; // constant term for the Likelihood calculation

  // Float64* m_unscaleContinuumFluxAxisDerivZ;

  Float64 m_nominalWidthDefault;

  std::string m_fittingmethod;

  std::string m_lineRatioType;

  bool m_forcedisableMultipleContinuumfit = false;

  CTLambdaRangePtrVector m_lambdaRanges;
  CCSpectrumVectorPtr m_inputSpcs;
  std::shared_ptr<CLMEltListVector> m_ElementsVector;

  bool m_opt_firstpass_forcedisableMultipleContinuumfit = true;
  Int32 m_opt_fitcontinuum_maxN;
  std::string m_opt_firstpass_fittingmethod = "hybrid";
  std::string m_opt_secondpass_fittingmethod = "hybrid";

  //  bool m_opt_enable_improveBalmerFit = false;

  bool m_useloglambdasampling = false;

  mutable CSpectraGlobalIndex m_spectraIndex;
};

} // namespace NSEpic

#endif // _REDSHIFT_LINEMODEL_ELEMENTLIST_
