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
#ifndef _REDSHIFT_OPERATOR_CONTINUUM_FITTING_BASE_
#define _REDSHIFT_OPERATOR_CONTINUUM_FITTING_BASE_

#include <vector>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/operator/operator.h"
#include "RedshiftLibrary/spectrum/maskBuilder.h"

namespace NSEpic {

class CSpectrum;
class COperatorResult;
class CModelSpectrumResult;

struct TFitQuality {
  Float64 reducedChiSquare = INFINITY;
  Float64 pValue = 0;
  Float64 meanResiduals = INFINITY;
  Float64 stdResiduals = INFINITY;
  Float64 skewnessResiduals = INFINITY;
  Float64 kurtosisResiduals = INFINITY;
  Float64 ksResiduals = INFINITY;
  Float64 ksStdResiduals = INFINITY;
  Float64 ksStdMeanResiduals = INFINITY;
  Float64 andersonResiduals = INFINITY;
  Int32 nPixels = 0;
};

struct TContinuumResult {
  Float64 ebmvCoef = NAN;
  Int32 meiksinIdx = undefIdx;
  TFitQuality fitQuality;
};

/**
 * \ingroup Redshift
 */
class COperatorContinuumFitting : public COperator {

public:
  COperatorContinuumFitting();
  virtual ~COperatorContinuumFitting() = default;
  COperatorContinuumFitting(const COperatorContinuumFitting &) = default;
  COperatorContinuumFitting(COperatorContinuumFitting &&) = default;
  COperatorContinuumFitting &
  operator=(const COperatorContinuumFitting &) = default;
  COperatorContinuumFitting &operator=(COperatorContinuumFitting &&) = default;

  virtual bool IsFFTProcessing() { return false; };
  void setMaskBuilder(const std::shared_ptr<CMaskBuilder> &maskBuilder) {
    m_maskBuilder = maskBuilder;
  }

protected:
  std::shared_ptr<CMaskBuilder> m_maskBuilder;
  std::vector<std::shared_ptr<const CSpectrum>> m_spectra;
  std::vector<std::shared_ptr<const TFloat64Range>> m_lambdaRanges;
  TInt32List m_kStart, m_kEnd;

  void checkTemplateOverlap(const Float64 overlapFraction,
                            const Float64 overlapThreshold);
  virtual Float64 EstimateLikelihoodCstLog() const;
  Int32 m_nSamplesMinForContinuumFit;
};
} // namespace NSEpic

#endif
