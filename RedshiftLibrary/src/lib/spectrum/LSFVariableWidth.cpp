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
#include "RedshiftLibrary/spectrum/LSFVariableWidth.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/line/lineprofileSYM.h"
#include "RedshiftLibrary/log/log.h"
using namespace NSEpic;
using namespace std;

// used to create mapping between width and lambda values
// Instead of this
CLSFGaussianVariableWidth::CLSFGaussianVariableWidth(
    const std::shared_ptr<const TLSFGaussianVarWidthArgs> &args)
    : CLSF(GaussianVariableWidth,
           std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM())),
      m_width(args->width), m_spcAxis(args->lambdas) {
  IsValid();
}

Float64 CLSFGaussianVariableWidth::GetWidth(Float64 lambda,
                                            bool cliplambda) const {

  if (cliplambda)
    lambda = getSpectralRange().Clamp(lambda);

  if (!checkAvailability(lambda))
    THROWG(INTERNAL_ERROR, " lambda outside spectralAxis range");

  Int32 idx = undefIdx;
  TFloat64Index::getClosestLowerIndex(m_spcAxis.GetSamplesVector(), lambda,
                                      idx);

  if (m_spcAxis[idx] == lambda)
    return m_width[idx];

  // interpolation
  Float64 t = (lambda - m_spcAxis[idx]) / (m_spcAxis[idx + 1] - m_spcAxis[idx]);
  return m_width[idx] * (1 - t) + m_width[idx + 1] * t;
}

bool CLSFGaussianVariableWidth::IsValid() const {
  if (!m_width.size()) {
    THROWG(BAD_COUNTMATCH, "Width array cannot be null ");
  }
  if (m_spcAxis.GetSamplesCount() != m_width.size()) {
    THROWG(BAD_COUNTMATCH, "Sizes do not match between Spectral axis "
                           "and width axis");
  }
  for (Float64 w : m_width)
    if (w <= 0.)
      return false;
  return true;
}
