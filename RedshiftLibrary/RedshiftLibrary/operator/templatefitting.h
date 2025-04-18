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
#ifndef _REDSHIFT_OPERATOR_TEMPLATE_FITTING_
#define _REDSHIFT_OPERATOR_TEMPLATE_FITTING_

#include <cfloat>
#include <climits>
#include <numeric>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/operator/templatefittingBase.h"
#include "RedshiftLibrary/operator/templatefittingresult.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"
#include "RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/statistics/priorhelper.h"

namespace templateFitting_test {
class fitQuality_test;
}

namespace NSEpic {

struct TCrossProductResult {
  Float64 sumCross = 0.;
  Float64 sumT = 0.;
  Float64 sumS = 0.;
  Float64 sumCross_phot = 0.0;
  Float64 sumT_phot = 0.0;
  Float64 sumS_phot = 0.0;

  TCrossProductResult &operator+=(const TCrossProductResult &other) {
    sumT += other.sumT;
    sumS += other.sumS;
    sumCross += other.sumCross;
    sumCross_phot += other.sumCross_phot;
    sumT_phot += other.sumT_phot;
    sumS_phot += other.sumS_phot;
    return *this;
  }
};
struct TFittingResult {
  Float64 chiSquare = INFINITY;
  Float64 chiSquare_phot =
      0.0; // should be initialized at zero to be summed unconditionnaly
  Float64 ampl = NAN;
  Float64 ampl_err = NAN;
  Float64 ampl_sigma = NAN;
  Float64 logprior = 0.;
  TCrossProductResult cross_result;
};

struct TFittingIsmIgmResult : TFittingResult {
  TFittingIsmIgmResult(Int32 EbmvListSize, Int32 MeiksinListSize,
                       Int32 spcsize = 1)
      : overlapFraction(spcsize, NAN),
        ChiSquareInterm(EbmvListSize, TFloat64List(MeiksinListSize, DBL_MAX)),
        IsmCalzettiIdxInterm(EbmvListSize, undefIdx),
        IgmMeiksinIdxInterm(MeiksinListSize, undefIdx) {}

  TFloat64List overlapFraction;
  Float64 reducedChiSquare = INFINITY;
  Float64 pValue = 0;
  Float64 ebmvCoef = NAN;
  Int32 meiksinIdx = undefIdx;
  std::vector<TFloat64List> ChiSquareInterm;
  TInt32List IsmCalzettiIdxInterm;
  TInt32List IgmMeiksinIdxInterm;
};

class COperatorTemplateFitting : public COperatorTemplateFittingBase {

public:
  COperatorTemplateFitting(const TFloat64List &redshifts)
      : COperatorTemplateFittingBase(redshifts), m_kStart(m_spectra.size()),
        m_kEnd(m_spectra.size()){

        };
  virtual ~COperatorTemplateFitting() = default;

  std::shared_ptr<CTemplateFittingResult> Compute(
      const std::shared_ptr<const CTemplate> &tpl, Float64 overlapThreshold,
      std::string opt_interp, bool opt_extinction, bool opt_dustFitting,
      Float64 opt_continuum_null_amp_threshold = 0.,
      const CPriorHelper::TPriorZEList &logprior = CPriorHelper::TPriorZEList(),
      Int32 FitEbmvIdx = allIdx, Int32 FitMeiksinIdx = allIdx,
      TInt32Range zIdxRangeToCompute = TInt32Range(undefIdx, undefIdx),
      std::shared_ptr<CTemplateFittingResult> const &result = nullptr) override;

protected:
  friend class templateFitting_test::fitQuality_test;
  TFittingIsmIgmResult BasicFit(const std::shared_ptr<const CTemplate> &tpl,
                                Float64 redshift, Float64 overlapThreshold,
                                bool opt_extinction, bool opt_dustFitting,
                                const CPriorHelper::TPriorEList &logpriore,
                                const TInt32List &MeiksinList,
                                const TInt32List &EbmvList);

  virtual std::pair<TList<CMask>, Int32>
  getMaskListAndNSamples(Float64 redshift) const;

  virtual void init_fast_igm_processing(Int32 EbmvListSize);

  virtual bool igmIsInRange(const TFloat64RangeList &ranges) const;

  virtual TCrossProductResult ComputeCrossProducts(Int32 kM, Int32 kEbmv_,
                                                   Float64 redshift,
                                                   CMask const &mask,
                                                   Int32 spcIndex = 0);

  virtual void
  ComputeAmplitudeAndChi2(TFittingResult &fitres,
                          const CPriorHelper::SPriorTZE &logpriorTZ) const;

  bool m_option_igmFastProcessing;
  TInt32List m_kStart, m_kEnd;

  std::vector<TFloat64List> m_sumCross_outsideIGM;
  std::vector<TFloat64List> m_sumT_outsideIGM;
  std::vector<TFloat64List> m_sumS_outsideIGM;
};

} // namespace NSEpic

#endif
