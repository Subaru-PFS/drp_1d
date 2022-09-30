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
#ifndef _REDSHIFT_LINEMODEL_FITTING_
#define _REDSHIFT_LINEMODEL_FITTING_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/line/regulament.h"
#include "RedshiftLibrary/linemodel/abstractfitter.h"
#include "RedshiftLibrary/linemodel/continuummanager.h"
#include "RedshiftLibrary/linemodel/element.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/linemodel/spectrummodel.h"
#include "RedshiftLibrary/operator/linemodelresult.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <boost/shared_ptr.hpp>
#include <memory>
#include <unordered_set>

namespace NSEpic {

class CLineRatioManager;

class CLineModelFitting {

public:
  CLineModelFitting();
  CLineModelFitting(
      const std::shared_ptr<const CSpectrum> &template_,
      const TLambdaRange
          &lambdaRange); // only used for template orthogonalization, TODO use
                         // only one of the future subclasses ? at least inherit
                         // from clinemodelfitting

  void initParameters();
  void initMembers();
  void LoadCatalog(const CLineCatalog::TLineVector &restLineList);
  void LoadCatalogOneMultiline(const CLineCatalog::TLineVector &restLineList);
  void
  LoadCatalogTwoMultilinesAE(const CLineCatalog::TLineVector &restLineList);

  void LogCatalogInfos();

  void setRedshift(Float64 redshift, bool reinterpolatedContinuum);
  void setContinuumComponent(std::string component);
  const std::string &getContinuumComponent() const {
    return m_continuumManager->getContinuumComponent();
  };

  bool initDtd();
  Float64 EstimateDTransposeD(const std::string &spcComponent) const;
  Float64 EstimateMTransposeM() const;
  Float64 EstimateLikelihoodCstLog() const;
  Float64 getDTransposeD();
  Float64 getLikelihood_cstLog();

  Int32 GetNElements() const;

  void SetVelocityEmission(Float64 vel);
  void SetVelocityAbsorption(Float64 vel);
  void SetVelocityEmissionOneElement(Float64 vel, Int32 idxElt);
  void SetVelocityAbsorptionOneElement(Float64 vel, Int32 idxElt);

  void setVelocity(Float64 vel, Int32 lineType);
  void setVelocity(Float64 vel, Int32 idxElt, Int32 lineType);

  Float64 GetVelocityEmission() const;
  Float64 GetVelocityAbsorption() const;
  Int32 ApplyVelocityBound(Float64 inf, Float64 sup);

  bool initModelAtZ(Float64 redshift,
                    const CSpectrumSpectralAxis &spectralAxis);

  Float64 fit(Float64 redshift, CLineModelSolution &modelSolution,
              CContinuumModelSolution &continuumModelSolution,
              Int32 contreest_iterations = 0, bool enableLogging = 0);
  TFloat64Range &getLambdaRange() { return m_dTransposeDLambdaRange; };

  void SetFittingMethod(const std::string &fitMethod);
  void setLineRatioType(const std::string &lineratio);
  void SetSecondpassContinuumFitPrms();

  void SetAbsLinesLimit(Float64 limit);

  Float64 GetRedshift() const;

  CMask getOutsideLinesMask() const;
  Float64 getOutsideLinesSTD(Int32 which) const;

  Int32 getSpcNSamples() const;

  Float64 getLeastSquareContinuumMeritFast() const;
  Float64 getLeastSquareContinuumMerit() const;
  Float64 getLeastSquareMeritUnderElements() const;
  Float64 getScaleMargCorrection(Int32 idxLine = -1) const {
    return m_Elements.getScaleMargCorrection(idxLine);
  }

  Float64 getStrongerMultipleELAmpCoeff() const;

  std::unordered_set<std::string> getLinesAboveSNR(Float64 snrcut = 3.5) const {
    return m_model->getLinesAboveSNR(*m_lambdaRange, snrcut);
  }

  Float64 getCumulSNRStrongEL() const;
  Float64 getCumulSNROnRange(TInt32Range idxRange) const;

  void LoadModelSolution(const CLineModelSolution &modelSolution);
  CLineModelSolution GetModelSolution(Int32 opt_level = 0) const;

  Float64 getModelFluxVal(Int32 idx) const;
  void logParameters();
  CLineModelElementList m_Elements;
  std::unique_ptr<CAbstractFitter> m_fitter;
  std::shared_ptr<CLineRatioManager> m_lineRatioManager;
  std::shared_ptr<const CSpectrum> m_inputSpc;
  const CLineCatalog::TLineVector m_RestLineList;

  Int32 setPassMode(Int32 iPass);
  Int32 GetPassNumber() const;

  void prepareAndLoadContinuum(Int32 icontfitting, Float64 redshift);
  void computeSpectrumFluxWithoutContinuum();
  bool isContinuumComponentTplfitxx() const {
    return m_continuumManager->isContinuumComponentTplfitxx();
  }

  std::shared_ptr<CSpectrumModel> getSpectrumModel() {
    return m_model;
  } // not const because of tplortho

  std::shared_ptr<const CSpectrumModel> getConstSpectrumModel() {
    return m_model;
  }

  // we keep that getters temporarily, waiting for refactoring linemodel extrema
  // result
  Int32 getTplratio_count() const;
  TFloat64List getTplratio_priors();

  std::string getLineRatioType() { return m_lineRatioType; }

  Int32 m_pass = 1;
  bool m_enableAmplitudeOffsets;

  Float64 m_LambdaOffsetMin = -400.0;
  Float64 m_LambdaOffsetMax = 400.0;
  Float64 m_LambdaOffsetStep = 25.0;
  bool m_enableLambdaOffsetsFit;

  std::shared_ptr<CContinuumManager> m_continuumManager;

private:
  void fitAmplitudesSimplex();

  void SetLSF();

  TInt32List findLineIdxInCatalog(const CLineCatalog::TLineVector &restLineList,
                                  const std::string &strTag, Int32 type) const;

  void applyPolynomCoeffs(Int32 eIdx, const TPolynomCoeffs &polynom_coeffs);
  void addDoubleLine(const CLine &r1, const CLine &r2, Int32 index1,
                     Int32 index2, Float64 nominalWidth, Float64 a1,
                     Float64 a2);

  void applyRules(bool enableLogs = false);

  std::shared_ptr<CSpectrumModel> m_model;

  Float64 m_dTransposeD; // the cached dtd (maximum chisquare value)
  TFloat64Range m_dTransposeDLambdaRange; // the lambdaRange used to computed
                                          // cached dTransposeD values
  Float64 m_likelihood_cstLog; // constant term for the Likelihood calculation

  // Float64* m_unscaleContinuumFluxAxisDerivZ;

  std::string m_LineWidthType;

  Float64 m_velocityEmission;
  Float64 m_velocityAbsorption;
  Float64 m_velocityEmissionInit;
  Float64 m_velocityAbsorptionInit;

  Float64 m_nominalWidthDefaultEmission;
  Float64 m_nominalWidthDefaultAbsorption;

  std::string m_fittingmethod;

  std::string m_lineRatioType;

  Int32 m_secondpass_fitContinuum_dustfit;
  Int32 m_secondpass_fitContinuum_igm;

  bool m_forcedisableMultipleContinuumfit = false;

  std::shared_ptr<const TFloat64Range> m_lambdaRange;

  linetags ltags;

  bool m_opt_firstpass_forcedisableMultipleContinuumfit = true;

  std::string m_opt_firstpass_fittingmethod = "hybrid";
  std::string m_opt_secondpass_fittingmethod = "hybrid";
  bool m_ignoreLinesSupport = false;

  //  bool m_opt_enable_improveBalmerFit = false;

  bool m_useloglambdasampling = false;
};

} // namespace NSEpic

#endif // _REDSHIFT_LINEMODEL_ELEMENTLIST_
