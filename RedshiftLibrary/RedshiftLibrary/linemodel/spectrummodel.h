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

#ifndef __REDSHIFT_LM_SPECTRUM_MODEL__
#define __REDSHIFT_LM_SPECTRUM_MODEL__

#include <unordered_set>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/linemodel/continuumcomponent.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/operator/continuumfitting.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

namespace NSEpic {

class CLineModelSolution;
class CContinuumModelSolution;
class CTemplate;
class COperatorTemplateFittingBase;
class CSpectrumModel {
public:
  CSpectrumModel(
      const std::shared_ptr<CLineModelElementList> &elements,
      const std::shared_ptr<const CSpectrum> &spc,
      const CLineMap &m_RestLineList,
      const std::shared_ptr<CContinuumModelSolution> &continuumModelSolution,
      const std::shared_ptr<COperatorContinuumFitting>
          &continuumFittingOperator,
      Int32 spcIndex);

  void reinitModel() { m_SpectrumModel.SetFluxAxis(m_ContinuumFluxAxis); };
  void refreshModel(CLine::EType lineTypeFilter = CLine::EType::nType_All);
  void reinitModelUnderElements(const TInt32List &filterEltsIdx, Int32 lineIdx);
  void refreshModelUnderElements(const TInt32List &filterEltsIdx,
                                 Int32 lineIdx = undefIdx);

  CSpectrumFluxAxis
  getModel(const TInt32List &eIdx_list,
           CLine::EType lineTypeFilter = CLine::EType::nType_All) const;
  void setContinuumToInputSpc();
  void setContinuumComponent(TContinuumComponent const &component);
  void EstimateSpectrumContinuum(Float64 opt_enhance_lines);

  const CSpectrum &GetModelSpectrum() const;
  const CSpectrumFluxAxis &GetModelContinuum() const;

  CSpectrum GetObservedSpectrumWithLinesRemoved(
      CLine::EType lineTypeFilter = CLine::EType::nType_All);
  Float64 GetWeightingAnyLineCenterProximity(Int32 sampleIndex,
                                             const TInt32List &EltsIdx) const;
  std::pair<TInt32Range, TFloat64List>
  GetLineRangeAndProfile(Int32 eIdx, Int32 line_id, Float64 redshift) const;

  std::tuple<Float64, Float64, Float64>
  GetContinuumWeightedSumInRange(TInt32Range const &indexRange,
                                 TFloat64List const &weights,
                                 TPolynomCoeffs const &polynomCoeffs) const;

  std::tuple<Float64, Float64, Float64>
  getContinuumSquaredResidualInRange(TInt32Range const &indexRange);

  Float64 getMaxContinuumUnderElement(Int32 eIdx) const;

  std::pair<Float64, Float64>
  getModelSquaredResidualUnderElements(TInt32List const &EltsIdx,
                                       bool with_continuum,
                                       bool with_weight = true) const;

  std::pair<Float64, Float64> getFluxDirectIntegration(
      const TInt32List &eIdx_list, const TInt32List &subeIdx_list,
      bool substract_abslinesmodel, const TFloat64Range &lambdaRange) const;
  std::unordered_set<std::string>
  getLinesAboveSNR(const TFloat64Range &lambdaRange,
                   Float64 snrcut = 3.5) const;

  bool m_enableAmplitudeOffsets = false;
  Float64 m_Redshift = 0.;
  // new methods
  Int32 m__count = 0;
  std::shared_ptr<COperatorContinuumFitting> m_continuumFittingOperator;

  void initModelWithContinuum();
  void setContinuumFromTplFit(Float64 alpha, Float64 tplAmp,
                              const TFloat64List &polyCoeffs);
  const CSpectrumFluxAxis &getSpcFluxAxis() const { return m_SpcFluxAxis; }
  const CSpectrumFluxAxis &getContinuumFluxAxis() const {
    return m_ContinuumFluxAxis;
  }
  const CSpectrumFluxAxis &getSpcFluxAxisNoContinuum() const {
    return m_spcFluxAxisNoContinuum;
  }

  Int32 ApplyContinuumPowerLawOnGrid(
      std::shared_ptr<CContinuumModelSolution> const &continuum);

  Int32 ApplyContinuumTplOnGrid(const std::shared_ptr<const CTemplate> &tpl,
                                Float64 zcontinuum);
  void initObserveGridContinuumFlux(Int32 size);
  const TPhotVal &getPhotValues() const { return m_photValues; };

private:
  CSpectrumFluxAxis
  getContinuumUnderLines(const TInt32RangeList &indexRangeList,
                         const TInt32List &eIdx_list,
                         bool substract_abslinesmodel) const;
  std::shared_ptr<const CSpectrum> m_inputSpc; // model
  const CLineMap &m_RestLineList;
  std::shared_ptr<CContinuumModelSolution> m_fitContinuum;

  CSpectrum m_SpectrumModel; // model
  std::shared_ptr<CLineModelElementList> m_Elements;
  CSpectrumFluxAxis m_ContinuumFluxAxis;
  CSpectrumFluxAxis m_SpcFluxAxis;
  CSpectrumFluxAxis
      m_spcFluxAxisNoContinuum; // observed spectrum for line fitting

  Int32 m_spcIndex = 0;

  TAxisSampleList
      m_observeGridContinuumFlux; // the continuum spectre without the
  // amplitude coeff; m_ContinuumFLux = amp *
  // m_observeGridContinuumFlux
  TPhotVal m_photValues;
};

class CSpcModelVector {
public:
  CSpcModelVector(const CSpectraGlobalIndex &spcIndex)
      : m_spectraIndex(spcIndex) {}

  void push_back(const CSpectrumModel &model) { m_models.push_back(model); }

  CSpectrumModel &getSpectrumModel() {
    m_spectraIndex.assertIsValid();
    return m_models.at(m_spectraIndex.get());
  }

  const CSpectrumModel &getSpectrumModel() const {
    m_spectraIndex.assertIsValid();
    return m_models.at(m_spectraIndex.get());
  }

  void refreshAllModels() {
    for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
      getSpectrumModel().refreshModel();
    }
  }

  void reinitAllModels() {
    for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
      getSpectrumModel().reinitModel();
    }
  }

  void refreshAllModelsUnderElements(const TInt32List &filterEltsIdx,
                                     Int32 line_index = undefIdx) {
    for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
      getSpectrumModel().refreshModelUnderElements(filterEltsIdx, line_index);
    }
  }

  Float64 getModelResidualRmsUnderElements(TInt32List const &EltsIdx,
                                           bool with_continuum,
                                           bool with_weight = true) {
    Float64 fit_allObs = 0;
    Float64 sumErr_allObs = 0;
    std::size_t nb_nan = 0;
    for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
      auto [fit, sumErr] =
          getSpectrumModel().getModelSquaredResidualUnderElements(
              EltsIdx, with_continuum, with_weight);
      if (fit == 0.0)
        continue;
      if (std::isnan(fit)) {
        nb_nan++;
        continue;
      }
      fit_allObs += fit;
      sumErr_allObs += sumErr;
    }
    if (nb_nan == m_models.size())
      return NAN;
    return sumErr_allObs != 0.0 ? sqrt(fit_allObs / sumErr_allObs) : NAN;
  }

  Float64 getMaxContinuumUnderElement(Int32 eIdx) const {
    Float64 max_all = -INFINITY;
    for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
      auto const &model = getSpectrumModel();
      Float64 const max = model.getMaxContinuumUnderElement(eIdx);
      max_all = std::max(max, max_all);
    }
    return max_all;
  }

  void setEnableAmplitudeOffsets(bool enableAmplitudeOffsets) {
    for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
      getSpectrumModel().m_enableAmplitudeOffsets = enableAmplitudeOffsets;
    }
  }

private:
  std::vector<CSpectrumModel> m_models;
  mutable CSpectraGlobalIndex m_spectraIndex;
};
using CSpcModelVectorPtr = std::shared_ptr<CSpcModelVector>;
} // namespace NSEpic
#endif
