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
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/spectrum/LSF.h"

namespace fluxcorrectionmeiksin_test { // boost test suite
// all boost_auto_test_case that use private method
class correction_multiply_test;
class correction_multiply_test_CteResolution;
class correction_multiply_test_CteResolution25_4;
class correction_multiply_test_CteResolution25_4_incontext;
class correction_test;
} // namespace fluxcorrectionmeiksin_test
namespace NSEpic {
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
    return 7;
  }; // harcoded value from the number of cols in the ascii files
  Int32 getRedshiftIndex(Float64 z) const;
  Float64 getCorrection(Int32 zIdx, Int32 meiksinIdx, Int32 lbdaIdx) const {
    return m_corrections[zIdx].fluxcorr[meiksinIdx][lbdaIdx];
  };
  const TFloat64List &getRedshiftBins() const { return m_zbins; };
  Float64 getLambdaMin() const { return m_LambdaMin; };
  Float64 getLambdaMax() const { return m_LambdaMax; };

private:
  friend class fluxcorrectionmeiksin_test::correction_multiply_test;
  friend class fluxcorrectionmeiksin_test::
      correction_multiply_test_CteResolution;
  friend class fluxcorrectionmeiksin_test::
      correction_multiply_test_CteResolution25_4;
  friend class fluxcorrectionmeiksin_test::
      correction_multiply_test_CteResolution25_4_incontext;
  friend class fluxcorrectionmeiksin_test::correction_test;

  TFloat64List applyAdaptativeKernel(const TFloat64List &arr,
                                     const Float64 z_center,
                                     const std::shared_ptr<const CLSF> &lsf,
                                     const TFloat64List &lambdas);

  TFloat64List m_zbins;
  std::vector<MeiksinCorrection> m_rawCorrections;
  std::vector<MeiksinCorrection> m_corrections;
  Float64 m_LambdaMin;
  Float64 m_LambdaMax;
  TFloat64Range m_convolRange;
  bool m_convolved = false;
};

} // namespace NSEpic

#endif
