#ifndef __REDSHIFT_LM_SPECTRUM_MODEL__
#define __REDSHIFT_LM_SPECTRUM_MODEL__

#include <unordered_set>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

namespace NSEpic {

class CLineModelSolution;
class CTplModelSolution;
class CTemplate;
class COperatorTemplateFittingBase;
class CSpectrumModel {
public:
  CSpectrumModel(
      const std::shared_ptr<CLineModelElementList> &elements,
      const std::shared_ptr<const CSpectrum> &spc,
      const CLineMap &m_RestLineList,
      const std::shared_ptr<CTplModelSolution> &tfv,
      const std::shared_ptr<COperatorTemplateFittingBase> &TFOperator,
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
  void setContinuumComponent(const std::string &component);
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

  std::pair<Float64, Int32>
  getContinuumSquaredResidualInRange(TInt32Range const &indexRange);

  std::pair<Float64, Float64>
  getModelSquaredResidualUnderElements(TInt32List const &EltsIdx,
                                       bool with_continuum) const;

  std::pair<Float64, Float64> getFluxDirectIntegration(
      const TInt32List &eIdx_list, const TInt32List &subeIdx_list,
      bool substract_abslinesmodel, const TFloat64Range &lambdaRange) const;
  std::unordered_set<std::string>
  getLinesAboveSNR(const TFloat64Range &lambdaRange,
                   Float64 snrcut = 3.5) const;

  void integrateFluxes_usingTrapez(const CSpectrumFluxAxis &continuumFlux,
                                   const TInt32RangeList &indexRangeList,
                                   Float64 &sumFlux, Float64 &sumErr) const;

  bool m_enableAmplitudeOffsets = false;
  Float64 m_Redshift = 0.;
  // new methods
  Int32 m__count = 0;
  std::shared_ptr<COperatorTemplateFittingBase> m_templateFittingOperator;

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

  Int32 ApplyContinuumOnGrid(const std::shared_ptr<const CTemplate> &tpl,
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
  std::shared_ptr<CTplModelSolution> m_fitContinuum;

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

using CSpcModelVectorPtr = std::shared_ptr<std::vector<CSpectrumModel>>;
} // namespace NSEpic
#endif
