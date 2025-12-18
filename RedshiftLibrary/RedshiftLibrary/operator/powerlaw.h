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
#ifndef _REDSHIFT_OPERATOR_POWER_LAW_
#define _REDSHIFT_OPERATOR_POWER_LAW_

#include "RedshiftLibrary/common/curve3d.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/linemodel/continuummodelsolution.h"
#include "RedshiftLibrary/operator/continuumfitting.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/operator/pass.h"
#include "RedshiftLibrary/processflow/result.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"
#include <Eigen/Dense>

class PowerLaw_fixture;

namespace powerLawOperator_test {
class init;
class basicfit_powerlaw;
class basicfit_without_extinction;
class basicfit_without_extinction_only_one_coef;
class basicfit_double_without_extinction;
class basicfit_double_default;
class basicfit_double_with_var;
class basicfit_multiobs;
class basicfit_negative;
class basicfit_default;
class simple_powerlaw;
} // namespace powerLawOperator_test

namespace NSEpic {

struct TPowerLawCoefs {
  // a * x^b
  Float64 a = NAN;
  Float64 b = NAN;
  Float64 stda = NAN;
  Float64 stdb = NAN;
};

typedef std::vector<std::vector<TPowerLawCoefs>> T2DPowerLawCoefs;
typedef std::pair<TPowerLawCoefs, TPowerLawCoefs> TPowerLawCoefsPair;
typedef std::vector<std::vector<TPowerLawCoefsPair>> T2DPowerLawCoefsPair;
typedef std::pair<T2DPowerLawCoefs, T2DPowerLawCoefs> TPair2DPowerLawCoefs;

struct TPowerCoefLimits {
  Float64 min = -INFINITY;
  Float64 max = INFINITY;
};
typedef std::pair<TPowerCoefLimits, TPowerCoefLimits> TPowerCoefsPairLimits;

struct TPowerLawCalcStorage {
  // 1 <-> power law first part / 2 <-> second part

  // General values
  Float64 xc;

  // n pixels
  Int32 N1 = 0;
  Int32 N2 = 0;
  // sum pixels weights
  Float64 n1 = 0;
  Float64 n2 = 0;
  // sum ln xi * wi
  Float64 sx1 = 0;
  Float64 sx2 = 0;
  // sum (ln xi)^2 * wi
  Float64 sxx1 = 0;
  Float64 sxx2 = 0;
  // sum ln yi * wi
  Float64 sy1 = 0;
  Float64 sy2 = 0;
  // sum ln xi * ln yi * wi
  Float64 sxy1 = 0;
  Float64 sxy2 = 0;

  // sum wi*(lnxc - lnxi)**2 (on second part only)
  Float64 sx2mc2 = 0;

  Eigen::Matrix3d m; // The M^TN^-1M matrix

  // Values of the M^TN^-1Y vector
  Eigen::Vector3d v;
};

// For one z
struct TPowerLawResult : TContinuumResult {
  Float64 chiSquare = INFINITY;
  TPowerLawCoefsPair coefs;
};

struct TChi2Result {
  Int32 igmIdx;
  Int32 ismIdx;
  Float64 chi2;
};

class COperatorPowerLaw : public COperatorContinuumFitting,
                          public COperatorPass {

public:
  COperatorPowerLaw(const TFloat64List &redshifts = TFloat64List(),
                    Float64 lambdaCut = POWER_LOW_WAVELENGTH_CUT);

  COperatorPowerLaw(COperatorPowerLaw const &other) = default;
  COperatorPowerLaw &operator=(COperatorPowerLaw const &other) = default;

  COperatorPowerLaw(COperatorPowerLaw &&other) = default;
  COperatorPowerLaw &operator=(COperatorPowerLaw &&other) = default;
  ~COperatorPowerLaw() = default;

  std::shared_ptr<const COperatorResult>
  Compute(bool opt_extinction, bool opt_dustFitting, Float64 nullFluxThreshold,
          Int32 FitEbmvIdx, Int32 FitMeiksinIdx);
  CModelSpectrumResult
  ComputeSpectrumModel(const CContinuumModelSolution &continuum,
                       Int32 spcIndex);
  bool checkCoefsOrNull(TPowerLawCoefs &coefs) const;
  bool checkCoefsOrNull(TPowerLawCoefsPair &coefs) const;
  TPowerLawCoefs DEFAULT_COEFS = {NAN, NAN, INFINITY, INFINITY};
  TPowerLawCoefsPair DEFAULT_COEFS_PAIR = {DEFAULT_COEFS, DEFAULT_COEFS};
  TPowerLawCoefs NULL_COEFS = {0, 0, INFINITY, INFINITY};
  TPowerLawCoefsPair NULL_COEFS_PAIR = {NULL_COEFS, NULL_COEFS};

protected:
  friend ::PowerLaw_fixture;
  friend powerLawOperator_test::basicfit_powerlaw;
  friend powerLawOperator_test::basicfit_without_extinction;
  friend powerLawOperator_test::basicfit_without_extinction_only_one_coef;
  friend powerLawOperator_test::init;
  friend powerLawOperator_test::basicfit_double_without_extinction;
  friend powerLawOperator_test::basicfit_double_default;
  friend powerLawOperator_test::basicfit_double_with_var;
  friend powerLawOperator_test::basicfit_multiobs;
  friend powerLawOperator_test::basicfit_negative;
  friend powerLawOperator_test::basicfit_default;
  friend powerLawOperator_test::simple_powerlaw;

  TPowerLawResult BasicFit(Float64 redshift, bool opt_extinction,
                           bool opt_dustFitting, Float64 nullFluxThreshold);

private:
  // igm ism curves
  std::shared_ptr<const CSpectrumFluxCorrectionMeiksin> m_igmCorrectionMeiksin;
  std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
      m_ismCorrectionCalzetti;

  // Indexes of igm / ism elements to take into account
  TList<Int32> m_igmIdxList;
  TList<Int32> m_ismIdxList;
  Int32 m_nIgmCurves;
  Int32 m_nIsmCurves;

  TList<Int32> m_nPixels;
  Float64 m_lambdaCut;
  Int32 m_nSpectra;
  std::vector<CSpectrumSpectralAxis> m_spcSpectralAxis_restframe;
  TPowerCoefsPairLimits m_powerCoefsLimits;

  void initIgmIsm(bool opt_extinction, bool opt_dustFitting, Int32 FitEbmvIdx,
                  Int32 FitMeiksinIdx);
  void addTooFewSamplesWarning(Int32 N, Int32 igmIdx, Int32 ismIdx,
                               const char *funcName) const;
  TPowerLawCoefsPair computeConstantLawCoefs(TFloat64List const &flux,
                                             TFloat64List const &error) const;
  TPowerLawCoefsPair computeFullPowerLawCoefs(Int32 N1, Int32 N2,
                                              TCurve const &lnCurve);
  TAxisSampleList lnLambda(TAxisSampleList const &lambda) const;
  T2DPowerLawCoefsPair powerLawCoefs3D(T3DCurve const &lnCurves);
  TBoolList computeSNRCompliantPixels(TFloat64List const &spectrumFlux,
                                      TFloat64List const &spectrumFluxError,
                                      Float64 nullFluxThreshold) const;
  T3DCurve computeLnCurve(T3DCurve const &emittedCurve) const;
  T2DList<Float64> computeChi2(T3DCurve const &curve3D,
                               T2DPowerLawCoefsPair const &coefs);
  TChi2Result findMinChi2OnIgmIsm(T3DCurve const &curve,
                                  T2DPowerLawCoefsPair const &coefs);
  Float64 computeDoublePowerLaw(TPowerLawCoefsPair const &coefs,
                                Float64 lambda) const;
  Float64 computePowerLaw(TPowerLawCoefs const &coefs, Float64 lambda) const;
  TPowerLawCoefs compute2PassSimplePowerLawCoefs(TCurve const &lnCurves);
  TPowerLawCoefsPair compute2PassDoublePowerLawCoefs(TCurve const &lnCurves);
  TPowerLawCoefsPair computeDoublePowerLawCoefs(
      TCurve const &lnCurve,
      std::optional<TPowerLawCoefsPair> const &coefsFirstEstim = std::nullopt);
  TCurve initializeFluxCurve(Float64 redshift, Float64 nullFluxThreshold);
  T3DList<Float64>
  computeIsmIgmCorrections(Float64 redshift,
                           CSpectrumSpectralAxis const &spectrumLambdaRest,
                           bool opt_extinction, bool opt_dustFitting) const;
  TList<Float64>
  computeIsmIgmCorrection(Float64 redshift,
                          CSpectrumSpectralAxis const &spectrumLambdaRest,
                          Int32 igmIdx, Float64 ismCoef) const;
  T3DCurve computeEmittedCurve(Float64 redshift, bool opt_extinction,
                               bool opt_dustFitting, TCurve &&fluxCurve);
  TPowerLawCoefs computeSimplePowerLawCoefs(
      TCurve const &lnCurve,
      std::optional<TPowerLawCoefs> const &coefsFirstEstim = std::nullopt);
  TFloat64List computeModelFlux(const CSpectrumSpectralAxis &lambdaRestAxis,
                                const Float64 redshift, const Int32 meiksinIdx,
                                const Float64 ebmvCoef,
                                const TPowerLawCoefsPair &coefs) const;
  Float64 limitCoef(Float64 coef, TPowerCoefLimits limits) const;
  void limitCoefs(TPowerLawCoefsPair &coefs);

  void updatePowerLawCalcStorage(
      TCurve const &lnCurve,
      std::optional<TPowerLawCoefsPair> const &coefsFirstEstim);
  void updatePowerLawCalcStorageForSimple(
      TCurve const &lnCurve,
      std::optional<TPowerLawCoefs> const &coefsFirstEstim);
  TPowerLawCoefsPair computeDoublePowerLawCoefs_b2_fixed(Float64 b2);
  TPowerLawCoefsPair computeDoublePowerLawCoefs_b1_fixed(Float64 b1);
  TPowerLawCoefsPair computeDoublePowerLawCoefs_b1_b2_fixed(Float64 b1,
                                                            Float64 b2);
  TPowerLawCoefs computeSimplePowerLawCoefs_b_fixed(Float64 b) const;
  std::pair<Float64, Float64> computea2(const Float64 a1, const Float64 b1,
                                        const Float64 b2, Float64 varA1,
                                        Float64 varb1, Float64 varb2,
                                        Float64 covb1b2, Float64 covA1b1,
                                        Float64 covA1b2) const;
  std::pair<Float64, Float64> computea1(const Float64 a2, const Float64 b1,
                                        const Float64 b2, Float64 varA2,
                                        Float64 varb1, Float64 varb2,
                                        Float64 covb1b2, Float64 covA2b1,
                                        Float64 covA2b2) const;
  Float64 stdExpA(Float64 a, Float64 varA) const;
  Float64 computeVarA2(Float64 varA1, Float64 varb1, Float64 varb2,
                       Float64 covb1b2, Float64 covA1b1, Float64 covA1b2) const;
  Float64 computeVarA1(Float64 varA2, Float64 varb1, Float64 varb2,
                       Float64 covb1b2, Float64 covA2b1, Float64 covA2b2) const;
  TPowerLawCalcStorage m_powerLawCalcStorage;
};

inline Float64 COperatorPowerLaw::computePowerLaw(TPowerLawCoefs const &coefs,
                                                  Float64 lambda) const {
  return coefs.a * std::pow(lambda, coefs.b);
}

inline Float64
COperatorPowerLaw::computeDoublePowerLaw(TPowerLawCoefsPair const &fullCoefs,
                                         Float64 lambda) const {
  return computePowerLaw(
      lambda < m_lambdaCut ? fullCoefs.first : fullCoefs.second, lambda);
}

} // namespace NSEpic

#endif