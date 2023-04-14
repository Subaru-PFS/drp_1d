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
#ifndef _REDSHIFT_LINE_PROFILE_SYMIGM_
#define _REDSHIFT_LINE_PROFILE_SYMIGM_
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/line/lineprofile.h"
#include "RedshiftLibrary/line/lineprofileSYM.h"
#include "RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h"
#include <cmath>
#include <string>

namespace lineProfile_test {
class lineprofileSYMIGM_test;
} // namespace lineProfile_test

namespace NSEpic {
/**
 * \ingroup Redshift
 */
class CLineProfileSYMIGM : public CLineProfileSYM {
public:
  CLineProfileSYMIGM(const std::shared_ptr<CSpectrumFluxCorrectionMeiksin>
                         &igmcorrectionMeiksin,
                     const Float64 nsigmasupport = N_SIGMA_SUPPORT);

  Float64 GetLineProfileVal(Float64 x, Float64 x0,
                            Float64 sigma) const override;
  Float64 GetLineFlux(Float64 x0, Float64 sigma,
                      Float64 A = 1.0) const override;
  Float64 GetLineProfileDerivZ(Float64 x, Float64 x0, Float64 redshift,
                               Float64 sigma) const override;
  Float64 GetLineProfileDerivSigma(Float64 x, Float64 x0,
                                   Float64 sigma) const override;
  TSymIgmParams GetSymIgmParams() const override;
  bool isSymIgm() const override { return true; };
  bool isSymIgmFit() const override { return m_igmFit; };

  void SetSymIgmParams(const TSymIgmParams &params) override;
  void SetSymIgmFit(bool val = true) override { m_igmFit = val; };
  void SetSymIgmFixed() override { m_igmFit = false; };

  void resetParams() override;
  Int32 getIGMIdxCount() const override;

private:
  friend class lineProfile_test::lineprofileSYMIGM_test;
  CLineProfile *CloneImplementation() const override {
    return new CLineProfileSYMIGM(*this);
  }
  void CheckMeiksinInit() const;
  Float64 getIGMCorrection(Float64 x) const;

  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> m_igmCorrectionMeiksin;
  Float64 m_redshift = NAN;
  Int32 m_igmidx = -1;
  bool m_igmFit = false;
};
} // namespace NSEpic
#endif
