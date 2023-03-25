#ifndef __REDSHIFT_LM_SPECTRUM_MODEL__
#define __REDSHIFT_LM_SPECTRUM_MODEL__

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

#include <unordered_set>

namespace NSEpic {

class CLineModelSolution;
class CTplModelSolution;
class CTemplate;
class COperatorTemplateFitting;
class CSpectrumModel {
public:
  CSpectrumModel(CLineModelElementList &elements,
                 std::shared_ptr<const CSpectrum> spc,
                 const TLineVector &m_RestLineList,
                 std::shared_ptr<CTplModelSolution> tfv);

  void reinitModel() { m_SpectrumModel.SetFluxAxis(m_ContinuumFluxAxis); };
  void refreshModel(Int32 lineTypeFilter = -1);
  void reinitModelUnderElements(const TInt32List &filterEltsIdx, Int32 lineIdx);
  void refreshModelUnderElements(const TInt32List &filterEltsIdx,
                                 Int32 lineIdx = -1);

  CSpectrumFluxAxis getModel(const TInt32List &eIdx_list,
                             Int32 lineTypeFilter = -1) const;
  void setContinuumToInputSpc();
  void setContinuumComponent(const std::string &component);
  void EstimateSpectrumContinuum(Float64 opt_enhance_lines);

  const CSpectrum &GetModelSpectrum() const;
  const CSpectrumFluxAxis &GetModelContinuum() const;

  const CSpectrum &
  GetObservedSpectrumWithLinesRemoved(Int32 lineTypeFilter = -1);
  Float64 GetWeightingAnyLineCenterProximity(Int32 sampleIndex,
                                             const TInt32List &EltsIdx) const;

  Float64 GetContinuumError(Int32 eIdx, Int32 subeIdx);
  Float64 getModelErrorUnderElement(Int32 eltId,
                                    const CSpectrumFluxAxis &fluxRef) const;
  void getFluxDirectIntegration(const TInt32List &eIdx_list,
                                const TInt32List &subeIdx_list,
                                bool substract_abslinesmodel, Float64 &fluxdi,
                                Float64 &snrdi,
                                const TFloat64Range &lambdaRange) const;
  std::unordered_set<std::string>
  getLinesAboveSNR(const TFloat64Range &lambdaRange,
                   Float64 snrcut = 3.5) const;

  void
  integrateFluxes_usingTrapez(const CSpectrumFluxAxis &continuumFlux,
                              const std::vector<CRange<Int32>> &indexRangeList,
                              Float64 &sumFlux, Float64 &sumErr) const;

  bool m_enableAmplitudeOffsets = false;
  Float64 m_Redshift = 0.;
  // new methods
  Int32 m__count = 0;
  std::shared_ptr<COperatorTemplateFitting> m_templateFittingOperator;
  Int32 m_spcIndex = 0;
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

private:
  CSpectrumFluxAxis
  getContinuumUnderLines(const TInt32RangeList &indexRangeList,
                         const TInt32List &eIdx_list,
                         bool substract_abslinesmodel) const;
  std::shared_ptr<const CSpectrum> m_inputSpc; // model
  const TLineVector &m_RestLineList;
  std::shared_ptr<CTplModelSolution> m_fitContinuum;

  CSpectrum m_SpectrumModel; // model
  CLineModelElementList &m_Elements;
  CSpectrumFluxAxis m_ContinuumFluxAxis;
  CSpectrum m_spcCorrectedUnderLines;
  CSpectrumFluxAxis m_SpcFluxAxis;
  CSpectrumFluxAxis
      m_spcFluxAxisNoContinuum; // observed spectrum for line fitting

  TAxisSampleList
      m_observeGridContinuumFlux; // the continuum spectre without the
  // amplitude coeff; m_ContinuumFLux = amp *
  // m_observeGridContinuumFlux
};

} // namespace NSEpic
#endif
