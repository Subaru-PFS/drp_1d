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
#include "RedshiftLibrary/line/lineprofileLOR.h"
#include "RedshiftLibrary/common/exception.h"

using namespace NSEpic;
using namespace std;

CLineProfileLOR::CLineProfileLOR(Float64 nsigmasupport)
    : CLineProfile(nsigmasupport, LOR) {}

Float64 CLineProfileLOR::GetLineProfileVal(Float64 x, Float64 x0,
                                           Float64 sigma) const {
  Float64 xc = x - x0;
  Float64 val = 0.0;
  Float64 xsurc;

  xsurc = xc / sigma;
  val = 1.0 / (1 + xsurc * xsurc);
  return val;
}

Float64 CLineProfileLOR::GetNSigmaSupport() const {
  return m_nsigmasupport * 2.0;
}

Float64 CLineProfileLOR::GetLineFlux(Float64 x0, Float64 sigma,
                                     Float64 A) const {
  return A * sigma * M_PI;
}

Float64 CLineProfileLOR::GetLineProfileDerivX0(Float64 x, Float64 lambda0,
                                               Float64 sigma) const {
  THROWG(ErrorCode::INTERNAL_ERROR,
         "Deriv for Z not IMPLEMENTED for profile LOR");
}

Float64 CLineProfileLOR::GetLineProfileDerivSigma(Float64 x, Float64 x0,
                                                  Float64 sigma) const {
  THROWG(ErrorCode::INTERNAL_ERROR, "No derivate sigma for LOR profile");
}
