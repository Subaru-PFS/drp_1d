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
#ifndef _REDSHIFT_SPECTRUM_FLUXCORRECTIONMEIKSIN_
#define _REDSHIFT_SPECTRUM_FLUXCORRECTIONMEIKSIN_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/range.h"
namespace fluxcorrectionmeiksin_test { // boost test suite
// all boost_auto_test_case that use private method
class overall_test;
class convolveByLSF_test;
} // namespace fluxcorrectionmeiksin_test
namespace NSEpic {
class CLSF;
typedef struct MeiksinCorrection {
  MeiksinCorrection(TFloat64List _lbda, std::vector<TFloat64List> _fluxcorr)
      : lbda(_lbda), fluxcorr(_fluxcorr){};

  MeiksinCorrection() = default;
  TFloat64List lbda;                  // wavelength
  std::vector<TFloat64List> fluxcorr; // 7 flux correction lists
} MeiksinCorrection;
/**
 * \ingroup Redshift
 */
class CSpectrumFluxCorrectionMeiksin {

public:
  CSpectrumFluxCorrectionMeiksin(
      std::vector<MeiksinCorrection> meiksinCorrectionCurves,
      TFloat64List zbins);

  void convolveByLSF(const std::shared_ptr<const CLSF> &lsf,
                     const TFloat64Range &lambdaRange);
  bool isConvolved() { return m_convolved; };
  Int32 getIdxCount() const {
    return m_rawCorrections[0].fluxcorr.size();
  }; // harcoded value from the number of cols in the ascii files
  Int32 getRedshiftIndex(Float64 z) const;
  Int32 getWaveIndex(Float64 w) const;

  std::pair<TFloat64List, TFloat64List>
  getWaveAndCorrectionVector(const TFloat64Range &wrange, Float64 redshift,
                             Int32 meiksinIdx) const;

  Float64 getCorrection(Int32 zIdx, Int32 meiksinIdx, Int32 lbdaIdx) const {
    return m_corrections[zIdx].fluxcorr[meiksinIdx].at(lbdaIdx);
  };
  Float64 getCorrectionDerivLbdaRest(Int32 zIdx, Int32 meiksinIdx,
                                     Int32 lbdaIdx) const;

  const TFloat64List &getRedshiftBins() const { return m_zbins; };

  Float64 getLambdaMin() const { return m_LambdaMin; };
  Float64 getLambdaMax() const { return m_LambdaMax; };
  Float64 getCorrection(Float64 redshift, Int32 meiksinIdx,
                        Float64 lambda) const;
  std::tuple<Float64, Float64>
  getCorrectionAndDerivLbdaRest(Float64 redshift, Int32 meiksinIdx,
                                Float64 lambda) const;

private:
  friend class fluxcorrectionmeiksin_test::overall_test;
  friend class fluxcorrectionmeiksin_test::convolveByLSF_test;

  TFloat64List
  ConvolveByLSFOneCurve(const TFloat64List &arr, const TFloat64List &lambdas,
                        const TFloat64List &fineLambdas,
                        const TFloat64Range &zbin,
                        const std::shared_ptr<const CLSF> &lsf) const;

  TInt32Range getWaveRangeIndices(const TFloat64Range &wrange, bool raw) const;
  TFloat64List getWaveVector(const TFloat64Range &wrange, bool raw) const;
  TFloat64List getWaveVector(const TInt32Range &wrange, bool raw) const;

  TFloat64List m_zbins;
  std::vector<MeiksinCorrection> m_rawCorrections;
  std::vector<MeiksinCorrection> m_corrections;
  const Float64 m_LambdaMin = NAN;
  const Float64 m_LambdaMax = NAN;
  Int32 m_LambdaSize = -1;
  Int32 m_fineLambdaSize = -1;
  TFloat64Range m_convolRange;
  bool m_convolved = false;
  Float64 m_finegridstep = IGM_RAW_STEP / IGM_OVERSAMPLING;
};

} // namespace NSEpic

#endif
