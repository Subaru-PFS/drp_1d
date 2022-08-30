#ifndef __REDSHIFT_LM_SPECTRUM_MODEL__
#define __REDSHIFT_LM_SPECTRUM_MODEL__

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

namespace NSEpic {

class CLineModelSolution;
class CContinuumModelSolution;
class CSpectrumModel {
public:
  CSpectrumModel(CLineModelElementList &elements,
                 std::shared_ptr<const CSpectrum> spc);

  void reinitModel();
  void refreshModel(Int32 lineTypeFilter = -1);
  void reinitModelUnderElements(const TInt32List &filterEltsIdx, Int32 lineIdx);
  void refreshModelInitAllGrid();
  void refreshModelUnderElements(const TInt32List &filterEltsIdx,
                                 Int32 lineIdx = -1);

  // methods modifying/initializing continuum (m_ContinuumFluxAxis)
  // void prepareContinuum(std::string continuumComponent);
  void PrepareContinuum();
  void SetContinuumComponent(std::string component);
  void EstimateSpectrumContinuum(Float64 opt_enhance_lines);

  const CSpectrum &GetModelSpectrum() const;
  CSpectrum GetSpectrumModelContinuum() const;
  const CSpectrumFluxAxis &GetModelContinuum() const;

  const CSpectrum &
  GetObservedSpectrumWithLinesRemoved(Int32 lineTypeFilter = -1);
  Float64 GetWeightingAnyLineCenterProximity(Int32 sampleIndex,
                                             const TInt32List &EltsIdx) const;

  Float64 GetContinuumError(Int32 eIdx, Int32 subeIdx);
  Float64 getModelErrorUnderElement(Int32 eltId) const;

  bool m_enableAmplitudeOffsets = false;
  Float64 m_Redshift = 0.;

  // new methods
  void setContinuum(const CSpectrumFluxAxis &continuum);
  void setRedshift(Float64 redshift) { m_Redshift = redshift; }
  void initModelWithContinuum();
  void substractContToFlux();
  void buildFromParameters(const CLineModelSolution &lm_solution,
                           const CContinuumModelSolution &c_solution);
  void setContinuumFromTplFit(Float64 alpha, Float64 tplAmp,
                              const TFloat64List &polyCoeffs,
                              const TAxisSampleList &observeGridContinuumFlux);
  const CSpectrumFluxAxis &getSpcFluxAxis() const { return m_SpcFluxAxis; }
  const CSpectrumFluxAxis &getContinuumFluxAxis() const {
    return m_ContinuumFluxAxis;
  }
  const CSpectrumFluxAxis &getSpcFluxAxisNoContinuum() const {
    return m_spcFluxAxisNoContinuum;
  }

private:
  std::shared_ptr<const CSpectrum> m_inputSpc; // model

  CSpectrum m_SpectrumModel; // model
  CLineModelElementList &m_Elements;
  CSpectrumFluxAxis m_ContinuumFluxAxis;
  CSpectrumFluxAxis m_modelFluxAxis;
  CSpectrum m_spcCorrectedUnderLines;
  CSpectrumFluxAxis m_SpcFluxAxis; // observed spectrum //this is DUBIOUS ! why
                                   // isnt'it a const ref ? why initializing
                                   // this way ? why  const auto &SpcFluxAxis =
                                   // m_SpcFluxAxis; and not directly a const
                                   // ref ?
  CSpectrumFluxAxis
      m_spcFluxAxisNoContinuum; // observed spectrum for line fitting

  std::string m_ContinuumComponent;
  bool isContinuumComponentTplfitxx() const {
    return m_ContinuumComponent == "tplfit" ||
           m_ContinuumComponent == "tplfitauto";
  }
};

} // namespace NSEpic
#endif
