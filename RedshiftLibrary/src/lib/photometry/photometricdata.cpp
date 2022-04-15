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
#include "RedshiftLibrary/photometry/photometricdata.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/log/log.h"

#include <stdexcept>

using namespace NSEpic;

CPhotometricData::CPhotometricData(const TStringList &name,
                                   const TFloat64List &flux,
                                   const TFloat64List &fluxerr) {
  // check all sizes are the same
  if (name.size() != flux.size() || name.size() != fluxerr.size()) {
    Log.LogError("CPhotometryData::CPhotometryData name, flux or fluxerr have "
                 "not the same size");
    throw GlobalException(ErrorCode::INTERNAL_ERROR,
                          "CPhotometryData::CPhotometryData name, flux or "
                          "fluxerr have not the same size");
  }

  for (Int32 i = 0; i < name.size(); i++) {
    m_flux[name[i]] = flux[i];
    m_fluxErr[name[i]] = fluxerr[i];
    m_oneOverErr2[name[i]] = 1. / (fluxerr[i] * fluxerr[i]);
    m_fluxOverErr2[name[i]] = flux[i] * flux[i] * m_oneOverErr2[name[i]];
  }
}

TStringList CPhotometricData::GetNameList() const {
  TStringList names;
  for (auto const &flux : m_flux) {
    names.push_back(flux.first);
  }
  return (names);
}
