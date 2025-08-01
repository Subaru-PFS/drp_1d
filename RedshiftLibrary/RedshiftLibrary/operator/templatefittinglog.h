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
#ifndef _REDSHIFT_OPERATOR_CHISQUARELOGLAMBDA_
#define _REDSHIFT_OPERATOR_CHISQUARELOGLAMBDA_

#include <fftw3.h>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/operator/templatefittingBase.h"
#include "RedshiftLibrary/operator/templatefittingresult.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"
#include "RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h"
#include "RedshiftLibrary/spectrum/fullspectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/statistics/priorhelper.h"

namespace templateFittingLog_test {
class EstimateXtY_test;
};

namespace NSEpic {

enum class EPrecomputedFFT {
  none,
  spcFluxOverErr2,
  spcFlux2OverErr2,
  spcOneOverErr2,
  spcMask,
  tplMask
};

class FFTPlans {
public:
  FFTPlans(Int32 n) { initFFT(n); };
  ~FFTPlans() { freeFFTPlans(); };
  FFTPlans(FFTPlans const &other) = delete;
  FFTPlans(FFTPlans &&other);
  FFTPlans &operator=(FFTPlans const &other) = delete;
  FFTPlans &operator=(FFTPlans &&other);
  void storeFFT(EPrecomputedFFT precomputed, fftw_complex *fft);
  void copyFFT(fftw_complex *src, fftw_complex *dest);
  bool isPrecomputed(EPrecomputedFFT precomputed);
  void allocatePrecomputedFFT(EPrecomputedFFT precomputed);

  // buffers for fft computation
  Int32 nPaddedSamples = 0;
  Float64 *inX = nullptr;
  fftw_complex *outX = nullptr;
  fftw_plan pX = nullptr;
  Float64 *inY = nullptr;
  fftw_complex *outY = nullptr;
  fftw_plan pY = nullptr;
  fftw_complex *outXY = nullptr;
  Float64 *inXY = nullptr;
  fftw_plan pBackward = nullptr;
  std::map<EPrecomputedFFT, fftw_complex *> precomputedFFT;

private:
  void initFFT(Int32 n);
  void freeFFTPlans();
  void freeFFTPrecomputedBuffers();
};

class COperatorTemplateFittingLog : public COperatorTemplateFittingBase {

public:
  COperatorTemplateFittingLog(const TFloat64List &redshifts);

  void SetRedshifts(const TFloat64List &redshifts) override;
  void CheckRedshifts();

  std::shared_ptr<CTemplateFittingResult> Compute(
      const CTemplate &tpl, Float64 overlapThreshold, std::string opt_interp,
      bool opt_extinction, bool opt_dustFitting,
      Float64 opt_continuum_null_amp_threshold = 0.,
      const CPriorHelper::TPriorZEList &logprior = CPriorHelper::TPriorZEList(),
      Int32 FitEbmvIdx = allIdx, Int32 FitMeiksinIdx = allIdx,
      TInt32Range zIdxRangeToCompute = TInt32Range(undefIdx, undefIdx),
      std::shared_ptr<CTemplateFittingResult> const &result = nullptr) override;

  inline bool IsFFTProcessing() override { return true; };

  // made public for unit-testing
  TInt32Range FindTplSpectralIndex(const TFloat64Range &redshiftrange) const;
  TInt32Range FindTplSpectralIndex(const CSpectrumSpectralAxis &spcSpectralAxis,
                                   const CSpectrumSpectralAxis &tplSpectralAxis,
                                   const TFloat64Range &redshiftrange) const;

private:
  friend class templateFittingLog_test::EstimateXtY_test;

  // shadow CoperatorContinuumFitting::m_spectra
  std::vector<std::shared_ptr<const CFullSpectrum>> m_spectra;

  Float64 m_logstep;
  Int32 m_ssRatio;

  // hardcoded config: FIT_RANGEZ
  Int32 exportIGMIdx = 5;
  Int32 exportISMIdx = -1;

  void FitAllz(std::shared_ptr<CTemplateFittingResult> result,
               const TInt32List &MeiksinList = TInt32List(1, 0),
               const TInt32List &EbmvList = TInt32List(1, 0),
               const CPriorHelper::TPriorZEList &logpriorze =
                   CPriorHelper::TPriorZEList(),
               CMask const &lineMask = {});

  void FitRangez(const TFloat64List &inv_err2, const TInt32Range &range,
                 const std::shared_ptr<CTemplateFittingResult> &result,
                 const TInt32List &MeiksinList, const TInt32List &EbmvList,
                 const Float64 &dtd, CMask const &lineMask);

  TInt32RangeList FindZRanges(const TFloat64List &redshifts);

  CMask maskTemplate();

  void updateGlobalResults(
      const std::shared_ptr<CTemplateFittingResult> &result,
      const std::shared_ptr<const CTemplateFittingResult> &subResult,
      Int32 resultIdx);

  void
  applyPrior(const std::shared_ptr<CTemplateFittingResult> &result,
             const std::shared_ptr<const CTemplateFittingResult> &subResult,
             Int32 resultIdx, const CPriorHelper::TPriorZEList &logpriorze,
             Float64 dtd);

  void computeFitQuality(const std::shared_ptr<CTemplateFittingResult> &result,
                         Int32 resultIdx, Int32 subResultSize,
                         Int32 firstTplIdx, CMask const &lineMask);

  void EstimateXtY(const TFloat64List &X, const TFloat64List &Y,
                   TFloat64List &XtY, FFTPlans &fftPlans,
                   EPrecomputedFFT fftX = EPrecomputedFFT::none,
                   EPrecomputedFFT fftY = EPrecomputedFFT::none);

  Int32 EstimateMtMFast(const TFloat64List &X, const TFloat64List &Y,
                        Int32 nShifts, TFloat64List &XtY);

  bool m_enableISM = true;
  bool m_enableIGM = true;
};
} // namespace NSEpic

#endif
