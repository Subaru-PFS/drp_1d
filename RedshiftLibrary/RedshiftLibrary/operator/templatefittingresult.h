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
#ifndef _REDSHIFT_OPERATOR_TEMPLATEFITTINGRESULT_
#define _REDSHIFT_OPERATOR_TEMPLATEFITTINGRESULT_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/operator/operator.h"
#include "RedshiftLibrary/processflow/result.h"

namespace NSEpic {

class TFittingIsmIgmResult;
class CTemplateFittingResult : public COperatorResult {

public:
  CTemplateFittingResult(Int32 n);
  CTemplateFittingResult(Int32 n, Int32 EbmvListSize, Int32 MeiksinListSize);
  void set_at_redshift(Int32 i, TFittingIsmIgmResult val);

  virtual ~CTemplateFittingResult() = default;
  CTemplateFittingResult(const CTemplateFittingResult &) = default;
  CTemplateFittingResult(CTemplateFittingResult &&) = default;
  CTemplateFittingResult &operator=(const CTemplateFittingResult &) = default;
  CTemplateFittingResult &operator=(CTemplateFittingResult &&) = default;

  Float64 SNRCalculation(Float64 dtm, Float64 mtm);

  TFloat64List Redshifts;

  // best fit results
  TFloat64List ChiSquare;
  TFloat64List ChiSquarePhot;
  TFloat64List FitAmplitude;
  TFloat64List FitAmplitudeError;
  TFloat64List FitAmplitudeSigma;
  TFloat64List FitEbmvCoeff;
  TInt32List FitMeiksinIdx;
  TFloat64List FitDtM;
  TFloat64List FitMtM;
  TFloat64List LogPrior;
  TFloat64List SNR;

  // intermediate chisquare results
  std::vector<std::vector<TFloat64List>>
      ChiSquareIntermediate; // chi2 for each intermediate results (for each
                             // config [z][Calzetti][Meiksin])
  std::vector<std::vector<TFloat64List>>
      IsmEbmvCoeffIntermediate; // calzetti dust coeff for each intermediate
                                // result (for each config
                                // [z][Calzetti][Meiksin])
  std::vector<std::vector<TInt32List>>
      IgmMeiksinIdxIntermediate; // meiksin idx for each intermediate result
                                 // (for each config [z][Calzetti][Meiksin])
  // TODO: std::vector<std::vector<TFloat64List>> LogPriorIntermediate
  // TODO: std::vector<std::vector<TFloat64List>> AmpIntermediate //is needed
  // for correct prior use in marg. mode tplmodel method

  Float64 CstLog = NAN;
  std::vector<TFloat64List> Overlap; // overlap rate by redshift by spectra
};

} // namespace NSEpic

#endif
