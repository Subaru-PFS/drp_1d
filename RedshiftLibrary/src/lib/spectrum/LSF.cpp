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
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/spectrum/LSFConstantResolution.h"
#include "RedshiftLibrary/spectrum/LSFConstantWidth.h"
#include "RedshiftLibrary/spectrum/LSF_NISPSIM_2016.h"
#include "RedshiftLibrary/spectrum/LSF_NISPVSSPSF_201707.h"
using namespace NSEpic;
using namespace std;

CLSF::CLSF(TLSFType name, CLineProfile_ptr &&profile)
    : m_name(name), m_profile(std::move(profile)) {}

Float64 CLSF::GetProfileVal(Float64 lambda, Float64 lambda0) const {
  return m_profile->GetLineProfileVal(lambda, lambda0, GetWidth(lambda0));
}

Float64 CLSF::GetProfileVal(Float64 lambda, Float64 lambda0,
                            Float64 sigma0) const {
  return m_profile->GetLineProfileVal(lambda, lambda0, sigma0);
}

const CLineProfile_ptr &CLSF::GetProfile() const { return m_profile; }

TFloat64List CLSF::getNormalizedProfileVector(TFloat64List lambda,
                                              Float64 lambda0) const {

  if (!IsValid())
    THROWG(ErrorCode::INTERNAL_ERROR, "lsf is not valid");

  Float64 sigma = GetWidth(lambda0, true);

  Float64 norm = 0.0;
  TFloat64List kernel;
  kernel.reserve(lambda.size());
  for (auto l : lambda) {
    kernel.push_back(GetProfileVal(l, lambda0, sigma));
    norm += kernel.back();
  }
  // normalizing  values
  Float64 inv_norm = 1.0 / norm;
  for (auto &k : kernel)
    k *= inv_norm;

  return kernel;
}
