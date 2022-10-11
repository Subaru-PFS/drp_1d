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
  COperatorTemplateFittingBase(const CSpectrum &spectrum,
                               const TFloat64Range &lambdaRange,
                               const TFloat64List &redshifts = TFloat64List());
  COperatorTemplateFittingBase(const CSpectrum &&spectrum,
                               const TFloat64Range &lambdaRange,
                               const TFloat64List &redshifts) = delete;

  COperatorTemplateFittingBase(COperatorTemplateFittingBase const &other) =
      delete; // ref member inside
  COperatorTemplateFittingBase &
  operator=(COperatorTemplateFittingBase const &other) = delete;
  virtual ~COperatorTemplateFittingBase() = default;

  virtual void SetRedshifts(TFloat64List redshifts) {
    m_redshifts = std::move(redshifts);
  };

  virtual std::shared_ptr<COperatorResult> Compute(
      const std::shared_ptr<const CTemplate> &tpl, Float64 overlapThreshold,
      const std::vector<CMask> &additional_spcMasks, std::string opt_interp,
      Int32 opt_extinction, Int32 opt_dustFitting,
      const CPriorHelper::TPriorZEList &logprior = CPriorHelper::TPriorZEList(),
      bool keepigmism = false, Float64 FitEbmvCoeff = -1.,
      Int32 FitMeiksinIdx = -1) = 0;

  std::shared_ptr<CModelSpectrumResult>
  ComputeSpectrumModel(const std::shared_ptr<const CTemplate> &tpl,
                       Float64 redshift, Float64 EbmvCoeff, Int32 meiksinIdx,
                       Float64 amplitude, const Float64 overlapThreshold);

  inline virtual bool IsFFTProcessing() { return false; };

  static Float64 GetIGMStartingRedshiftValue(const Float64 spcLbda0);

protected:
  virtual void RebinTemplate(const std::shared_ptr<const CTemplate> &tpl,
                             Float64 redshift, TFloat64Range &currentRange,
                             Float64 &overlaprate,
                             const Float64 overlapThreshold);
  // Likelihood
  virtual Float64
  EstimateLikelihoodCstLog(const CSpectrum &spectrum,
                           const TFloat64Range &lambdaRange) const;

  const CSpectrum &m_spectrum;
  TFloat64Range m_lambdaRange;
  TFloat64List m_redshifts;

  CTemplate m_templateRebined_bf;
  CSpectrumSpectralAxis m_spcSpectralAxis_restframe;
  CMask m_mskRebined_bf;
};

} // namespace NSEpic

#endif
