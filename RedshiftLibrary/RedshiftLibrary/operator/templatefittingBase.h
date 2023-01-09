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
#ifndef _REDSHIFT_OPERATOR_TEMPLATE_FITTING_BASE_
#define _REDSHIFT_OPERATOR_TEMPLATE_FITTING_BASE_

#include "RedshiftLibrary/operator/operator.h"

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/processflow/result.h"
#include "RedshiftLibrary/spectrum/maskBuilder.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/statistics/priorhelper.h"

#include <vector>

namespace NSEpic {

class CSpectrum;
class COperatorResult;
class CModelSpectrumResult;

/**
 * \ingroup Redshift
 */
class COperatorTemplateFittingBase : public COperator {

public:
  COperatorTemplateFittingBase(const TFloat64List &redshifts = TFloat64List());

  virtual ~COperatorTemplateFittingBase() = default;

  virtual void SetRedshifts(TFloat64List redshifts) {
    m_redshifts = std::move(redshifts);
  };

  virtual std::shared_ptr<COperatorResult> Compute(
      const std::shared_ptr<const CTemplate> &tpl, Float64 overlapThreshold,
      std::string opt_interp, bool opt_extinction, bool opt_dustFitting,
      Float64 opt_continuum_null_amp_threshold = 0.,
      const CPriorHelper::TPriorZEList &logprior = CPriorHelper::TPriorZEList(),
      Int32 FitEbmvIdx = undefIdx, Int32 FitMeiksinIdx = undefIdx) = 0;

  std::shared_ptr<CModelSpectrumResult>
  ComputeSpectrumModel(const std::shared_ptr<const CTemplate> &tpl,
                       Float64 redshift, Float64 EbmvCoeff, Int32 meiksinIdx,
                       Float64 amplitude, const Float64 overlapThreshold,
                       Int32 index = 0);

  inline virtual bool IsFFTProcessing() { return false; };

  static Float64 GetIGMStartingRedshiftValue(const Float64 spcLbda0);
  void setMaskBuilder(std::shared_ptr<CMaskBuilder> maskBuilder) {
    m_maskBuilder = maskBuilder;
  }

protected:
  virtual void RebinTemplate(const std::shared_ptr<const CTemplate> &tpl,
                             Float64 redshift, TFloat64Range &currentRange,
                             Float64 &overlaprate,
                             const Float64 overlapThreshold,
                             Int32 spcIndex = 0);

  void RebinTemplateFAS(const std::shared_ptr<const CTemplate> &tpl,
                        Float64 redshift, TFloat64Range &currentRange,
                        Float64 &overlaprate, const Float64 overlapThreshold);
  // Likelihood
  virtual Float64 EstimateLikelihoodCstLog() const;
  void applyPositiveAndNonNullConstraint(Float64 amp_sigma,
                                         Float64 &ampl) const;

  std::vector<std::shared_ptr<const CSpectrum>> m_spectra;
  std::vector<std::shared_ptr<const TFloat64Range>> m_lambdaRanges;
  TFloat64List m_redshifts;

  std::vector<CTemplate> m_templateRebined_bf;
  std::vector<CSpectrumSpectralAxis> m_spcSpectralAxis_restframe;
  std::vector<CMask> m_mskRebined_bf;
  std::shared_ptr<CMaskBuilder> m_maskBuilder;
  Float64 m_continuum_null_amp_threshold; // in SNR
};

} // namespace NSEpic

#endif
