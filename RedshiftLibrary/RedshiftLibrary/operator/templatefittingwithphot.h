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
#ifndef _REDSHIFT_OPERATOR_TEMPLATE_FITTING_WITHPHOT
#define _REDSHIFT_OPERATOR_TEMPLATE_FITTING_WITHPHOT

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/operator/templatefitting.h"
#include "RedshiftLibrary/photometry/photometricband.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/statistics/priorhelper.h"

#include <map>
#include <memory>
#include <string>
namespace NSEpic {

class COperatorTemplateFittingPhot : public COperatorTemplateFitting {

public:
  explicit COperatorTemplateFittingPhot(
      const CSpectrum &spectrum, const TFloat64Range &lambdaRange,
      const std::shared_ptr<const CPhotBandCatalog> &photbandcat,
      const Float64 weight = 1.0,
      const TFloat64List &redshifts = TFloat64List());

private:
  void checkInputPhotometry() const;

  void RebinTemplate(const std::shared_ptr<const CTemplate> &tpl,
                     Float64 redshift, const std::string &opt_interp,
                     TFloat64Range &currentRange, Float64 &overlaprate,
                     const Float64 overlapThreshold) override;

  void RebinTemplateOnPhotBand(const std::shared_ptr<const CTemplate> &tpl,
                               Float64 redshift, const std::string &opt_interp);

  void InitIsmIgmConfig(Float64 redshift,
                        const std::shared_ptr<CSpectrumFluxCorrectionCalzetti>
                            &ismCorrectionCalzetti,
                        const std::shared_ptr<CSpectrumFluxCorrectionMeiksin>
                            &igmCorrectionMeiksin,
                        Int32 EbmvListSize) override;

  bool
  CheckLyaIsInCurrentRange(const TFloat64Range &currentRange) const override;

  bool ApplyMeiksinCoeff(Int32 meiksinIdx) override;
  bool ApplyDustCoeff(Int32 kEbmv) override;

  TFittingResult ComputeLeastSquare(Int32 kM, Int32 kEbmv,
                                    const CPriorHelper::SPriorTZE &logprior,
                                    const CMask &spcMaskAdditional) override;

  void ComputePhotCrossProducts(Int32 kM, Int32 kEbmv_,
                                TFittingResult &fitResult,
                                Float64 &sumCross_phot, Float64 &sumT_phot,
                                Float64 &sumS_phot);

  Float64
  EstimateLikelihoodCstLog(const CSpectrum &spectrum,
                           const TFloat64Range &lambdaRange) const override;

  std::map<std::string, CSpectrumSpectralAxis> m_photSpectralAxis_restframe;
  std::map<std::string, CTemplate> m_templateRebined_phot;

  Float64 m_weight;
  std::shared_ptr<const CPhotBandCatalog> m_photBandCat;
  TStringList m_sortedBandNames;

  TFloat64List m_sumCross_outsideIGM_phot;
  TFloat64List m_sumT_outsideIGM_phot;
  TFloat64List m_sumS_outsideIGM_phot;
};

} // namespace NSEpic

#endif
