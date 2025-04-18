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
#include "RedshiftLibrary/linemodel/obsiterator.h"
#include "RedshiftLibrary/linemodel/spectrummodel.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

namespace NSEpic

{
class CContinuumManager;
class CAbstractFitter {
public:
  CAbstractFitter(const std::shared_ptr<CLMEltListVector> &elementsVector,
                  const CCSpectrumVectorPtr &inputSpcs,
                  const CTLambdaRangePtrVector &lambdaRanges,
                  const CSpcModelVectorPtr &spectrumModels,
                  const CLineMap &restLineList,
                  const CSpectraGlobalIndex &spcIndex,
                  bool enableAmplitudeOffsets = false,
                  bool enableLambdaOffsetsFit = false);

  void fit(Float64 redshift);

  virtual void resetSupport(Float64 redshift);

  void enableAmplitudeOffsets() { m_enableAmplitudeOffsets = true; }
  void enableLambdaOffsets() { m_enableLambdaOffsetsFit = true; }

  static std::shared_ptr<CAbstractFitter> makeFitter(
      std::string fittingMethod,
      const std::shared_ptr<CLMEltListVector> &elementsVector,
      const CCSpectrumVectorPtr &inputSpcs,
      const CTLambdaRangePtrVector &lambdaRanges,
      const CSpcModelVectorPtr &spectrumModels, const CLineMap &restLineList,
      std::shared_ptr<CContinuumManager> continuumManager,
      const CSpectraGlobalIndex &spcIndex, bool enableAmplitudeOffsets = false,
      bool enableLambdaOffsetsFit = false);

  void logParameters();

  TAsymParams fitAsymParameters(Float64 redshift, Int32 idxLyaE,
                                const Int32 &idxLineLyaE);

  Int32 fitAsymIGMCorrection(
      Float64 redshift,
      std::vector<std::pair<Int32, TInt32List>> const &idxLines);

  Float64 getModelResidualRmsUnderElements(TInt32List const &EltsIdx,
                                           bool with_continuum,
                                           bool with_weight = true) {
    return m_models->getModelResidualRmsUnderElements(EltsIdx, with_continuum,
                                                      with_weight);
  }
  Int32 m_cont_reestim_iterations = 0;

protected:
  virtual void doFit(Float64 redshift) = 0;

  void initFit(Float64 redshift);

  void resetLambdaOffsets() { m_ElementsVector->resetLambdaOffsets(); };

  void fitLyaProfile(Float64 redshift);

  void fitAmplitude(Int32 eltIndex, Float64 redshift, Int32 lineIdx = undefIdx);

  virtual void fitAmplitudeAndLambdaOffset(Int32 eltIndex, Float64 redshift,
                                           Int32 lineIdx = undefIdx,
                                           bool enableOffsetFitting = true);

  Float64 getLeastSquareMeritFast(Int32 eltIdx = undefIdx) const;

  void setLambdaOffset(const TInt32List &EltsIdx, Int32 offsetCount);

  bool HasLambdaOffsetFitting(TInt32List EltsIdx,
                              bool enableOffsetFitting) const;
  Int32 GetLambdaOffsetSteps(bool atLeastOneOffsetToFit) const;

  CSpectrumModel &getModel() { return m_models->getSpectrumModel(); }
  const CSpectrumModel &getModel() const {
    m_spectraIndex.assertIsValid();
    return m_models->getSpectrumModel();
  }
  const CSpectrum &getSpectrum() {
    m_spectraIndex.assertIsValid();
    return *((*m_inputSpcs).at(m_spectraIndex.get()));
  }
  const TLambdaRange &getLambdaRange() {
    m_spectraIndex.assertIsValid();
    return *(m_lambdaRanges.at(m_spectraIndex.get()));
  }
  CLineModelElementList &getElementList() {
    m_spectraIndex.assertIsValid();
    return m_ElementsVector->getElementList();
  }
  const CLineModelElementList &getElementList() const {
    m_spectraIndex.assertIsValid();
    return m_ElementsVector->getElementList();
  }

  std::vector<TLineModelElementParam_ptr> &getElementParam() {
    return m_ElementsVector->getElementParam();
  }
  const std::vector<TLineModelElementParam_ptr> &getElementParam() const {
    return m_ElementsVector->getElementParam();
  }

  Int32 m_nbElements;

  CCSpectrumVectorPtr m_inputSpcs;
  const CLineMap &m_RestLineList;
  CTLambdaRangePtrVector m_lambdaRanges;
  CSpcModelVectorPtr m_models;
  std::shared_ptr<CLMEltListVector> m_ElementsVector;

  CSpectraGlobalIndex m_spectraIndex;
  // hard coded options
  bool m_enableAmplitudeOffsets = false;
  bool m_enableLambdaOffsetsFit = false;

  Float64 m_LambdaOffsetMin = -MAX_LAMBDA_OFFSET;
  Float64 m_LambdaOffsetMax = MAX_LAMBDA_OFFSET;
  Float64 m_LambdaOffsetStep = LAMBDA_OFFSET_STEP;
  Float64 m_enlarge_line_supports = 0.0;

  Float64 m_absLinesLimit = 1.0;

  Float64 m_opt_lya_fit_asym_min = NAN;
  Float64 m_opt_lya_fit_asym_max = NAN;
  Float64 m_opt_lya_fit_asym_step = NAN;
  Float64 m_opt_lya_fit_width_min = NAN;
  Float64 m_opt_lya_fit_width_max = NAN;
  Float64 m_opt_lya_fit_width_step = NAN;
  Float64 m_opt_lya_fit_delta_min = NAN;
  Float64 m_opt_lya_fit_delta_max = NAN;
  Float64 m_opt_lya_fit_delta_step = NAN;
};
} // namespace NSEpic

#endif
