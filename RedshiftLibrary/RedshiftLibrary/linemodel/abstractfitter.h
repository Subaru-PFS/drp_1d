#ifndef _ABSTRACT_FITTER_H
#define _ABSTRACT_FITTER_H

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/linemodel/spectrummodel.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

namespace NSEpic

{
class CAbstractFitter {
public:
  CAbstractFitter(CLineModelElementList &elements,
                  std::shared_ptr<const CSpectrum> inputSpectrum,
                  std::shared_ptr<const TLambdaRange> lambdaRange,
                  std::shared_ptr<CSpectrumModel> spectrumModel);

  virtual void fit(Float64 redshift) = 0;

  void enableAmplitudeOffsets() { m_enableAmplitudeOffsets = true; }

  Int32 m_cont_reestim_iterations = 0;

protected:
  CLineModelElementList &m_Elements;
  const CSpectrum &m_inputSpc;
  const CLineCatalog::TLineVector &m_RestLineList;
  const TFloat64Range &m_lambdaRange;
  std::shared_ptr<CSpectrumModel> m_model;

  // hard coded options
  bool m_enableAmplitudeOffsets = false;
  bool m_enableLambdaOffsetsFit = true;

  Int32 m_AmplitudeOffsetsDegree = 2;
  Float64 m_LambdaOffsetMin = -400.0;
  Float64 m_LambdaOffsetMax = 400.0;
  Float64 m_LambdaOffsetStep = 25.0;

  Int32 fitAmplitudesLinSolveAndLambdaOffset(
      TInt32List EltsIdx, const CSpectrumSpectralAxis &spectralAxis,
      const CSpectrumFluxAxis &fluxAxis,
      const CSpectrumFluxAxis &continuumfluxAxis,
      std::vector<Float64> &ampsfitted, std::vector<Float64> &errorsfitted,
      bool enableOffsetFitting, Float64 redshift);

  Int32 fitAmplitudesLinSolve(const TInt32List &EltsIdx,
                              const CSpectrumSpectralAxis &spectralAxis,
                              const CSpectrumFluxAxis &fluxAxis,
                              const CSpectrumFluxAxis &continuumfluxAxis,
                              std::vector<Float64> &ampsfitted,
                              std::vector<Float64> &errorsfitted,
                              Float64 redshift);
};
} // namespace NSEpic

#endif
