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
#ifndef _REDSHIFT_GAUSSIANFIT_GAUSSIANFIT_
#define _REDSHIFT_GAUSSIANFIT_GAUSSIANFIT_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace NSEpic {

class CSpectrum;

/**
 * \ingroup Redshift
 * Single gaussian equation fit.
 **/
class CGaussianFit {

public:
  enum EStatus {
    nStatus_Success = (1 << 0),

    // Results corresponding to an error
    nStatus_IllegalInput = (1 << 1),
    nStatus_InvalidStudyRange = (1 << 2),
    nStatus_IterationHasNotConverged = (1 << 3),

    // Results corresponding to a success
    nStatus_FailToReachTolerance = (1 << 4)
  };

  CGaussianFit();

  EStatus Compute(const CSpectrum &s, const TInt32Range &studyRange);

  void GetResults(Float64 &amplitude, Float64 &position, Float64 &width) const;
  void GetResultsPolyCoeff0(Float64 &coeff0) const;
  void GetResultsError(Float64 &amplitude, Float64 &position,
                       Float64 &width) const;

private:
  static int GaussF(const gsl_vector *param, void *data, gsl_vector *f);
  static int GaussDF(const gsl_vector *param, void *data, gsl_matrix *J);
  static int GaussFDF(const gsl_vector *param, void *data, gsl_vector *f,
                      gsl_matrix *J);

  void ComputeFirstGuess(const CSpectrum &spectrum,
                         const TInt32Range &studyRange, Int32 polyOrder,
                         Float64 &peakValue, Float64 &peakPos,
                         Float64 &gaussAmp);
  EStatus getReturnCode(int status) const;
  struct SUserData {
    const CSpectrum *spectrum;
    const TInt32Range *studyRange;
    Int32 polyOrder;
  };

  Float64 m_AbsTol;
  Float64 m_Amplitude;
  Float64 m_AmplitudeErr;
  Float64 m_C;
  Float64 m_CErr;
  Float64 m_Mu;
  Float64 m_MuErr;
  Float64 m_RelTol;
  Float64 m_coeff0;
  Int32 m_PolyOrder;
};

} // namespace NSEpic

#endif
