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

#ifndef _REDSHIFT_LINE_RATIO_MANAGER_
#define _REDSHIFT_LINE_RATIO_MANAGER_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/line/catalog.h"

namespace NSEpic {

class CLineModelElementList;
class CSpectrumModel;
class CSpectrum;
class CContinuumManager;
class CLineModelSolution;
class CTplModelSolution;
class CAbstractFitter;
class CLineRatioManager {
public:
  CLineRatioManager(CLineModelElementList &elements,
                    std::shared_ptr<CSpectrumModel> model,
                    std::shared_ptr<const CSpectrum> inputSpc,
                    std::shared_ptr<const TFloat64Range> lambdaRange,
                    std::shared_ptr<CContinuumManager> continuumManager,
                    const CLineCatalog::TLineVector &restLineList);
  CLineRatioManager() = delete;
  virtual ~CLineRatioManager() {}
  CLineRatioManager(CLineRatioManager const &other) = default;
  CLineRatioManager &operator=(CLineRatioManager const &other) = default;

  CLineRatioManager(CLineRatioManager &&other) = default;
  CLineRatioManager &operator=(CLineRatioManager &&other) = default;

  virtual int prepareFit(Float64 redshift) { return 1; }
  virtual bool init(Float64 redshift, Int32 itratio = -1);
  virtual Float64 computeMerit(Int32 itratio) = 0;
  virtual void finish(Float64 redshift) {}
  virtual void saveResults(Int32 itratio) {}
  virtual void setPassMode(Int32 iPass);
  virtual TFloat64List getTplratio_priors() {
    static TFloat64List dumb;
    return dumb;
  }
  virtual Int32 getTplratio_count() const { return 0; }

  virtual void logParameters();

  void SetLeastSquareFastEstimationEnabled(Int32 enabled) {
  } // TODO, called in computeFirstPass in the general case but only active when
    // line_ratio_type = tplratio -> should be reviewed

  void setFitter(std::shared_ptr<CAbstractFitter> fitter) { m_fitter = fitter; }
  static std::shared_ptr<CLineRatioManager>
  makeLineRatioManager(const std::string &lineRatioType,
                       CLineModelElementList &elements,
                       std::shared_ptr<CSpectrumModel> model,
                       std::shared_ptr<const CSpectrum> inputSpc,
                       std::shared_ptr<const TFloat64Range> lambdaRange,
                       std::shared_ptr<CContinuumManager> continuumManager,
                       const CLineCatalog::TLineVector &restLineList,
                       std::shared_ptr<CAbstractFitter> fitter);

protected:
  void setLyaProfile(Float64 redshift,
                     const CLineCatalog::TLineVector &lineList);
  void setAsymProfile(Int32 idxLyaE, Int32 idxLineLyaE, Float64 redshift,
                      const CLineCatalog::TLineVector &lineList);
  virtual Int32
  getLineIndexInCatalog(Int32 iElts, Int32 idxLine,
                        const CLineCatalog::TLineVector &catalog) const;

  void setSymIgmProfile(Int32 iElts, const TInt32List &idxLineIGM,
                        Float64 redshift);

  Float64 getLeastSquareMerit() const;

  CLineModelElementList &m_Elements;
  std::shared_ptr<CSpectrumModel> m_model;
  std::shared_ptr<const CSpectrum> m_inputSpc;
  std::shared_ptr<CContinuumManager> m_continuumManager;
  std::shared_ptr<CAbstractFitter> m_fitter;
  const CLineCatalog::TLineVector &m_RestLineList;
  std::shared_ptr<const TFloat64Range> m_lambdaRange;

  bool m_forceDisableLyaFitting = false;
  bool m_forceLyaFitting = false;

  bool m_opt_lya_forcefit = false;
  bool m_opt_lya_forcedisablefit = false;
};

} // namespace NSEpic

#endif
