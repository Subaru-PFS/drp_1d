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

#ifndef _REDSHIFT_TPLRATIO_MANAGER_
#define _REDSHIFT_TPLRATIO_MANAGER_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/linemodel/lineratiomanager.h"
#include "RedshiftLibrary/statistics/priorhelper.h"

namespace NSEpic {

class CLineCatalogsTplRatio;

class CTplratioManager : public CLineRatioManager {
public:
  CTplratioManager(const std::shared_ptr<CLMEltListVector> &elementsVector,
                   const CSpcModelVectorPtr &models,
                   const CCSpectrumVectorPtr &inputSpcs,
                   const CTLambdaRangePtrVector &lambdaRanges,
                   std::shared_ptr<CContinuumManager> continuumManager,
                   const CLineMap &restLineList,
                   const CSpectraGlobalIndex &spcIndex);
  CTplratioManager() = delete;
  virtual ~CTplratioManager() = default;
  CTplratioManager(CTplratioManager const &other) = default;
  CTplratioManager &operator=(CTplratioManager const &other) = default;

  CTplratioManager(CTplratioManager &&other) = default;
  CTplratioManager &operator=(CTplratioManager &&other) = default;

  int prepareFit(Float64 redshift) override;
  bool init(Float64 redshift, Int32 itratio) override;

  virtual std::pair<Float64, Float64> computeMerit(Int32 itratio) override;
  void resetToBestRatio(Float64 redshift) override;
  void setPassMode(Int32 iPass) override;
  void saveResults(Int32 itratio) override { m_savedIdxFitted = itratio; };
  Int32 getTplratio_count() const override;
  TFloat64List getTplratio_priors() const override;

  void logParameters() override;
  const std::string &getTplratio_bestTplName() const;
  Float64 getTplratio_bestTplIsmCoeff() const;
  Float64 getTplratio_bestAmplitudeEm() const;
  Float64 getTplratio_bestAmplitudeAbs() const;
  Float64 getTplratio_bestAmplitudeUncertaintyEm() const;
  Float64 getTplratio_bestAmplitudeUncertaintyAbs() const;
  Float64 getTplratio_bestDtmEm() const;
  Float64 getTplratio_bestDtmAbs() const;
  Float64 getTplratio_bestMtmEm() const;
  Float64 getTplratio_bestMtmAbs() const;
  const TFloat64List &GetChisquareTplratio() const;
  TFloat64List GetPriorLinesTplratio() const;
  const TFloat64List &GetScaleMargTplratio() const;
  const TBoolList &GetStrongELPresentTplratio() const;
  const TBoolList &getHaELPresentTplratio() const;
  const TInt32List &GetNLinesAboveSNRTplratio() const;

  void setTplratioModel(Int32 itplratio, Float64 redshift,
                        bool enableSetVelocity = false);
  void SetLeastSquareFastEstimationEnabled(Int32 enabled);

  void SetForcedisableTplratioISMfit(bool opt);
  void duplicateTplratioResult(Int32 ifitting);
  void updateTplratioResults(Int32 ifitting, Float64 _merit,
                             Float64 _meritprior);
  Float64 computelogLinePriorMerit(
      Int32 itratio,
      const std::vector<CPriorHelper::SPriorTZE> &logPriorDataTplRatio);

  bool m_opt_firstpass_forcedisableTplratioISMfit = true;

protected:
  void initMerit(Int32 ntplratio);
  void SetTplratio_PriorHelper();
  void initTplratioCatalogs(Int32 opt_tplratio_ismFit);
  void SetNominalAmplitudes(Int32 iCatalog);

  Float64 GetIsmCoeff(Int32 idx) const;

  std::vector<std::vector<TFloat64List>> m_LineCatalogCorrespondingNominalAmp;
  Int32 m_savedIdxFitted = undefIdx;
  TFloat64List m_MeritTplratio;
  TFloat64List m_PriorMeritTplratio;
  std::vector<CPriorHelper::SPriorTZE> m_logPriorDataTplRatio;

  std::vector<TFloat64List> m_FittedAmpTplratio;
  std::vector<TFloat64List> m_FittedErrorTplratio;
  std::vector<TFloat64List> m_MtmTplratio;
  std::vector<TFloat64List> m_DtmTplratio;
  std::vector<TFloat64List> m_LyaAsymCoeffTplratio;
  std::vector<TFloat64List> m_LyaWidthCoeffTplratio;
  std::vector<TFloat64List> m_LyaDeltaCoeffTplratio;
  std::vector<TInt32List> m_LyaIgmIdxTplratio;
  std::vector<TFloat64List> m_LinesLogPriorTplratio;
  Int32 m_tplratioLeastSquareFast =
      0; // for rigidity=tplratio: switch to use fast least square estimation

  Int32 m_EmEltIdx = undefIdx;
  Int32 m_AbsEltIdx = undefIdx;

  TFloat64List m_ScaleMargCorrTplratio;
  TBoolList m_StrongELPresentTplratio;
  TBoolList m_StrongHalphaELPresentTplratio;
  TInt32List m_NLinesAboveSNRTplratio;

  std::shared_ptr<const CLineCatalogsTplRatio> m_CatalogTplRatio;
  std::shared_ptr<CPriorHelper> m_tplratio_priorhelper;

  bool m_forcedisableTplratioISMfit = false;
  Float64 m_NSigmaSupport;
  Float64 m_opt_haprior = -1.;
  bool m_opt_dust_calzetti;

private:
  void fillHalphaArray(Int32 idx);
};

} // namespace NSEpic

#endif
