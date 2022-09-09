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
#ifndef _ABSTRACT_FITTER_H
#define _ABSTRACT_FITTER_H

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/linemodel/spectrummodel.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

namespace NSEpic

{
class CContinuumManager;
class CAbstractFitter {
public:
  CAbstractFitter(CLineModelElementList &elements,
                  std::shared_ptr<const CSpectrum> inputSpectrum,
                  std::shared_ptr<const TLambdaRange> lambdaRange,
                  std::shared_ptr<CSpectrumModel> spectrumModel,
                  const CLineCatalog::TLineVector &restLineList);

  virtual void fit(Float64 redshift) = 0;

  void enableAmplitudeOffsets() { m_enableAmplitudeOffsets = true; }

  static std::unique_ptr<CAbstractFitter>
  makeFitter(std::string fittingMethod, CLineModelElementList &elements,
             std::shared_ptr<const CSpectrum> inputSpectrum,
             std::shared_ptr<const TLambdaRange> lambdaRange,
             std::shared_ptr<CSpectrumModel> spectrumModel,
             const CLineCatalog::TLineVector &restLineList,
             std::shared_ptr<CContinuumManager> continuumManager);
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
      std::vector<Float64> &ampsfitted, std::vector<Float64> &errorsfitted,
      bool enableOffsetFitting, Float64 redshift);

  Int32 fitAmplitudesLinSolve(const TInt32List &EltsIdx,
                              const CSpectrumSpectralAxis &spectralAxis,
                              const CSpectrumFluxAxis &fluxAxis,
                              const CSpectrumFluxAxis &continuumfluxAxis,
                              std::vector<Float64> &ampsfitted,
                              std::vector<Float64> &errorsfitted,
                              Float64 redshift);

private:
  void setOffset(const TInt32List &EltsIdx, Int32 offsetCount) const;
};
} // namespace NSEpic

#endif
