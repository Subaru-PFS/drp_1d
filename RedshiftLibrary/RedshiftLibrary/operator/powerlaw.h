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

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/operator/powerlawbase.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"

class PowerLaw_fixture;

namespace powerLawOperator_test {
class init;
class basicfit_powerlaw;
class basicfit_without_extinction;
class basicfit_without_extinction_only_one_coef;
class basicfit_simple_without_extinction;
class basicfit_simple_weighted_without_extinction;
class basicfit_double_without_extinction;
} // namespace powerLawOperator_test

namespace NSEpic {

struct TPowerLawCoefs {
  // a * x^b
  Float64 a = NAN;
  Float64 b = NAN;
};

typedef std::vector<std::vector<TPowerLawCoefs>> T2DPowerLawCoefs;
typedef std::pair<TPowerLawCoefs, TPowerLawCoefs> TPowerLawCoefsPair;
typedef std::vector<std::vector<TPowerLawCoefsPair>> T2DPowerLawCoefsPair;
typedef std::pair<T2DPowerLawCoefs, T2DPowerLawCoefs> TPair2DPowerLawCoefs;

struct TPowerLawResult {
  Float64 redshift = NAN;
  Float64 chiSquare = INFINITY;
  Int32 igmIdx = -1;
  Int32 ismIdx = -1;
  TPowerLawCoefsPair coefs;
  Float64 lambdaCut;
  TFloat64List lambdaRest;
  TFloat64List emittedFlux;
  TFloat64List fluxError;
  TBoolList pixelsToUse;
};

struct TCurve {
  // TODO here replace by axis ?
  TAxisSampleList lambda;
  TFloat64List flux;
  TBoolList pixelsToUse;
  TFloat64List fluxError;
};

struct T3DCurve {
  // CSpectrumSpectralAxis lambda;
  // T3DFloatList lambda;
  TFloat64List lambda;
  T3DFloatList flux;
  T3DBoolList pixelsToUse;
  T3DFloatList fluxError;
};

struct TChi2Result {
  Int32 igmIdx;
  Int32 ismIdx;
  Float64 chi2;
};

class COperatorPowerLaw : public COperatorPowerLawBase {

public:
  COperatorPowerLaw(const TFloat64List &redshifts = TFloat64List(),
                    Float64 lambdaCut = POWER_LOW_WAVELENGTH_CUT);

  COperatorPowerLaw(COperatorPowerLaw const &other) = default;
  COperatorPowerLaw &operator=(COperatorPowerLaw const &other) = default;

  COperatorPowerLaw(COperatorPowerLaw &&other) = default;
  COperatorPowerLaw &operator=(COperatorPowerLaw &&other) = default;
  ~COperatorPowerLaw() = default;

  std::shared_ptr<COperatorResult> Compute();

protected:
  friend ::PowerLaw_fixture;
  friend powerLawOperator_test::basicfit_powerlaw;
  friend powerLawOperator_test::basicfit_without_extinction;
  friend powerLawOperator_test::basicfit_without_extinction_only_one_coef;
  friend powerLawOperator_test::init;
  friend powerLawOperator_test::basicfit_simple_without_extinction;
  friend powerLawOperator_test::basicfit_simple_weighted_without_extinction;
  friend powerLawOperator_test::basicfit_double_without_extinction;

  TPowerLawResult BasicFit(Float64 redshift, bool opt_extinction,
                           bool opt_dustFitting, Float64 nullFluxThreshold,
                           Int32 nLogSamplesMin, std::string method);
  TAxisSampleList lnLambda(TAxisSampleList const &lambda);
  T2DPowerLawCoefsPair powerLawCoefs3D(T3DCurve const &lnCurves,
                                       TInt32Range const pixelsRange,
                                       Int32 nLogSamplesMin,
                                       std::string method) const;
  TBoolList computeSNRCompliantPixels(TFloat64List const &spectrumFlux,
                                      TFloat64List const &spectrumFluxError,
                                      Float64 nullFluxThreshold) const;
  T3DCurve computeLnCurves(T3DCurve const &emittedCurve);
  T2DFloatList computeChi2(T3DCurve const &curve3D,
                           T2DPowerLawCoefsPair const &coefs);
  TChi2Result findMinChi2OnIgmIsm(T3DCurve const &curve,
                                  T2DPowerLawCoefsPair const &coefs);
  Float64 theoreticalFluxAtLambda(TPowerLawCoefsPair coefs, Float64 lambda);
  Float64 computePowerLaw(TPowerLawCoefs coefs, Float64 lambda);
  TPowerLawCoefsPair computeFullPowerLawCoefs(TCurve const &lnCurve,
                                              TInt32Range const pixelsRange,
                                              Int32 nLogSamplesMin) const;
  T3DCurve initializeFluxCurve(Float64 redshift,
                               Float64 nullFluxThreshold) const;
  T3DFloatList
  computeIsmIgmCorrections(Float64 redshift,
                           CSpectrumSpectralAxis const &spectrumLambdaRest,
                           bool opt_extinction, bool opt_dustFitting) const;
  T3DCurve computeEmittedCurve(Float64 redshift, Float64 nullFluxThreshold,
                               bool opt_extinction, bool opt_dustFitting) const;
  TPowerLawCoefs computeSimplePowerLawCoefs(TCurve const &lnCurve,
                                            TInt32Range const pixelsRange,
                                            Int32 nLogSamplesMin,
                                            bool applyWeight) const;
  void checkInputPowerLawCoefs(TCurve const &lnCurve,
                               TInt32Range const pixelsRange) const;

private:
  // igm ism curves
  std::shared_ptr<CSpectrumFluxCorrectionCalzetti> m_ismCorrectionCalzetti;
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> m_igmCorrectionMeiksin;
  Int32 m_nIsmCurves;
  Int32 m_nIgmCurves;
  Int32 m_nPixels;
  Int32 m_firstPixelIdxInRange;
  Int32 m_lastPixelIdxInRange;

  Float64 m_lambdaCut;
};
} // namespace NSEpic

#endif