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
#include "RedshiftLibrary/spectrum/LSFConstantResolution.h"
#include "RedshiftLibrary/line/lineprofileSYM.h"
#include "RedshiftLibrary/log/log.h"

using namespace NSEpic;
using namespace std;

CLSFGaussianConstantResolution::CLSFGaussianConstantResolution(
    const Float64 resolution)
    : CLSF(GaussianConstantResolution,
           std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM())),
      m_Resolution(resolution) {
  if (!IsValid())
    THROWG(ErrorCode::INVALID_LSF,
           Formatter() << "invalid LSF, resolution=" << m_Resolution);
}

Float64 CLSFGaussianConstantResolution::GetWidth(Float64 lambda,
                                                 bool cliplambda) const {

  if (cliplambda)
    lambda = getSpectralRange().Clamp(lambda);

  Float64 defaultSigma = lambda / m_Resolution * INSTRUMENT_RESOLUTION_FACTOR;
  return defaultSigma;
}

bool CLSFGaussianConstantResolution::IsValid() const {
  return (m_Resolution > 1.);
}

Float64 CLSFGaussianConstantResolution::computeResolution(Float64 lambda,
                                                          Float64 width) {
  return lambda / width * INSTRUMENT_RESOLUTION_FACTOR;
}
