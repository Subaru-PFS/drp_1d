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

#include <map>
#include <memory>
#include <string>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/operator/templatefitting.h"
#include "RedshiftLibrary/photometry/photometricband.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/statistics/priorhelper.h"

namespace NSEpic {
#include <string>

class COperatorTemplateFittingPhot : public COperatorTemplateFitting {

public:
  explicit COperatorTemplateFittingPhot(
      const std::shared_ptr<const CPhotBandCatalog> &photbandcat,
      const TFloat64List &redshifts, const Float64 weight);

private:
  void checkInputPhotometry() const;

  void RebinTemplate(const std::shared_ptr<const CTemplate> &tpl,
                     Float64 redshift, TFloat64Range &currentRange,
                     Float64 &overlapFraction, const Float64 overlapThreshold,
                     Int32 spcIndex = 0) override;

  void RebinTemplateOnPhotBand(const std::shared_ptr<const CTemplate> &tpl,
                               Float64 redshift);

  void InitIsmIgmConfig(Float64 redshift, Int32 kstart, Int32 kend,
                        Int32 spcIndex) override;

  void init_fast_igm_processing(Int32 EbmvListSize) override;

  bool igmIsInRange(const TFloat64RangeList &ranges) const override;

  bool ApplyMeiksinCoeff(Int32 meiksinIdx, Int32 spcIndex = 0) override;
  bool ApplyDustCoeff(Int32 kEbmv, Int32 spcIndex = 0) override;

  std::pair<TList<CMask>, Int32>
  getMaskListAndNSamples(Float64 redshift) const override;

  void ComputePhotCrossProducts(Int32 kM, Int32 kEbmv_,
                                TCrossProductResult &fitResult);

  TCrossProductResult ComputeCrossProducts(Int32 kM, Int32 kEbmv_,
                                           Float64 redshift, CMask const &mask,
                                           Int32 spcIndex = 0) override;

  void ComputeAmplitudeAndChi2(
      TFittingResult &fitres,
      const CPriorHelper::SPriorTZE &logpriorTZ) const override;

  Float64 EstimateLikelihoodCstLog() const override;
  TPhotVal getIntegratedFluxes(Float64 ampl = 1.0) const override;
  std::map<std::string, CSpectrumSpectralAxis> m_photSpectralAxis_restframe;
  std::map<std::string, CTemplate> m_templateRebined_phot;

  Float64 m_weight;
  std::shared_ptr<const CPhotBandCatalog> m_photBandCat;
  TStringList m_sortedBandNames;

  TFloat64List m_sumCross_outsideIGM_phot;
  TFloat64List m_sumT_outsideIGM_phot;
  TFloat64List m_sumS_outsideIGM_phot;
  TStringList::const_iterator m_BandIgmEnd;
};

} // namespace NSEpic

#endif
