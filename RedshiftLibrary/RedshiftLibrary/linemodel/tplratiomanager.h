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
#include "RedshiftLibrary/linemodel/rigiditymanager.h"
#include "RedshiftLibrary/statistics/priorhelper.h"

namespace NSEpic {

class CLineCatalogsTplRatio;

class CTplratioManager : public CRigidityManager {
public:
  CTplratioManager(CLineModelElementList &elements,
                   std::shared_ptr<CSpectrumModel> model,
                   std::shared_ptr<const CSpectrum> inputSpc,
                   std::shared_ptr<const TFloat64Range> lambdaRange,
                   std::shared_ptr<CContinuumManager> continuumManager,
                   const CLineCatalog::TLineVector &restLineList);

  int prepareFit(Float64 redshift) override;
  void init(Float64 redshift, Int32 itratio) override;

  Float64 computeMerit(Int32 itratio) override;
  void finish(Float64 redshift) override;
  void setPassMode(Int32 iPass) override;
  void saveResults(Int32 itratio) override;

  void logParameters();
  const std::string &getTplratio_bestTplName() const;
  Float64 getTplratio_bestTplIsmCoeff() const;
  Float64 getTplratio_bestAmplitude() const;
  Float64 getTplratio_bestDtm() const;
  Float64 getTplratio_bestMtm() const;
  Int32 getTplratio_count() const;
  TFloat64List getTplratio_priors();
  const TFloat64List &GetChisquareTplratio() const;
  TFloat64List GetPriorLinesTplratio() const;
  const TFloat64List &GetScaleMargTplratio() const;
  const TBoolList &GetStrongELPresentTplratio() const;
  const TBoolList &getHaELPresentTplratio() const;
  const TInt32List &GetNLinesAboveSNRTplratio() const;

  bool setTplratioModel(Int32 itplratio, bool enableSetVelocity = false);
  bool setTplratioAmplitude(const TFloat64List &ampsElts,
                            const TFloat64List &errorsElts);
  void SetLeastSquareFastEstimationEnabled(Int32 enabled);

  void SetForcedisableTplratioISMfit(bool opt);
  void duplicateTplratioResult(Int32 ifitting, TFloat64List &bestTplratioMerit,
                               TFloat64List &bestTplratioMeritPrior);
  void updateTplratioResults(Int32 ifitting, Float64 _merit,
                             Float64 _meritprior,
                             TFloat64List &bestTplratioMerit,
                             TFloat64List &bestTplratioMeritPrior);
  Float64 computelogLinePriorMerit(
      Int32 itratio,
      const std::vector<CPriorHelper::SPriorTZE> &logPriorDataTplRatio);

  bool m_opt_firstpass_forcedisableTplratioISMfit = true;

protected:
  Int32 getLineIndexInCatalog(
      Int32 iElts, Int32 idxLine,
      const CLineCatalog::TLineVector &catalog) const override;

  void initMerit(Int32 ntplratio);
  void SetTplratio_PriorHelper();
  void initTplratioCatalogs(Int32 opt_tplratio_ismFit);
  bool SetMultilineNominalAmplitudesFast(Int32 iCatalog);

  Float64 GetIsmCoeff(Int32 idx) const;

  std::vector<std::vector<TFloat64List>>
      m_LineCatalogLinesCorrespondingNominalAmp;
  Int32 m_savedIdxFitted = -1; // for rigidity=tplratio
  TFloat64List m_bestTplratioMerit;
  TFloat64List m_bestTplratioMeritPrior;
  std::vector<CPriorHelper::SPriorTZE> m_logPriorDataTplRatio;

  TFloat64List m_ChisquareTplratio;
  std::vector<TFloat64List> m_FittedAmpTplratio;
  std::vector<TFloat64List> m_FittedErrorTplratio;
  std::vector<TFloat64List> m_MtmTplratio;
  std::vector<TFloat64List> m_DtmTplratio;
  std::vector<TFloat64List> m_LyaAsymCoeffTplratio;
  std::vector<TFloat64List> m_LyaWidthCoeffTplratio;
  std::vector<TFloat64List> m_LyaDeltaCoeffTplratio;
  std::vector<TInt32List> m_LyaIgmIdxTplratio;
  std::vector<TFloat64List> m_LinesLogPriorTplratio;

  std::string m_tplratioBestTplName = "undefined";
  Float64 m_tplratioBestTplIsmCoeff = NAN;
  Float64 m_tplratioBestTplAmplitude = NAN;
  Float64 m_tplratioBestTplDtm = NAN;
  Float64 m_tplratioBestTplMtm = NAN;
  Int32 m_tplratioLeastSquareFast =
      0; // for rigidity=tplratio: switch to use fast least square estimation

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
};

} // namespace NSEpic

#endif
