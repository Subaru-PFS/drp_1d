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
#include "RedshiftLibrary/line/lineprofileASYMFIT.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/log/log.h"
using namespace NSEpic;
using namespace std;

CLineProfileASYMFIT::CLineProfileASYMFIT(Float64 nsigmasupport,
                                         const TAsymParams &params,
                                         const std::string &centeringMethod)
    : CLineProfileASYM(ASYMFIT, nsigmasupport, params, centeringMethod) {}

CLineProfileASYMFIT::CLineProfileASYMFIT(const CLineProfileASYM &other)
    : CLineProfileASYMFIT(other.GetNSigmaSupport(), other.GetAsymParams(),
                          other.GetCenteringMethod()) {}

std::unique_ptr<CLineProfile> CLineProfileASYMFIT::cloneToASYM() const {
  return std::unique_ptr<CLineProfile>(new CLineProfileASYM(*this));
}

void CLineProfileASYMFIT::SetAsymParams(const TAsymParams &params) {
  if (std::isnan(params.sigma) || std::isnan(params.alpha) ||
      std::isnan(params.delta)) {
    Flag.warning(Flag.ASYMFIT_NAN_PARAMS,
                 Formatter() << "CLineProfileASYMFIT::" << __func__
                             << " AsymFit params are NaN");
  }
  m_asym_sigma_coeff = params.sigma;
  m_asym_alpha = params.alpha;
  m_asym_delta = params.delta;
}

void CLineProfileASYMFIT::resetParams() {
  m_asym_sigma_coeff = 2.;
  m_asym_alpha = 0.;
  m_asym_delta = 0.;
}