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

class TExtremaResult : public TCandidateZ {

public:
  TExtremaResult() = default;
  TExtremaResult(const TExtremaResult &res) = default;
  TExtremaResult(TExtremaResult &&res) = default;
  TExtremaResult &operator=(const TExtremaResult &res) = default;
  TExtremaResult &operator=(TExtremaResult &&res) = default;
  virtual ~TExtremaResult() = default;
  TExtremaResult(const TCandidateZ &candz) : TCandidateZ(candz) {
    m_type = "TExtremaResult";
  }
  TExtremaResult(const CContinuumModelSolution &cms)
      : FittedTplName(cms.tplName), FittedTplAmplitude(cms.tplAmplitude),
        FittedTplAmplitudeError(cms.tplAmplitudeError),
        FittedTplMerit(cms.tplMerit), FittedTplMeritPhot(cms.tplMeritPhot),
        FittedTplEbmvCoeff(cms.tplEbmvCoeff),
        FittedTplMeiksinIdx(cms.tplMeiksinIdx), FittedTplDtm(cms.tplDtm),
        FittedTplMtm(cms.tplMtm), FittedTplLogPrior(cms.tplLogPrior) {
    m_type = "TExtremaResult";
  }

  // Name of the best template fitted for continuum
  std::string FittedTplName = "";
  // Amplitude for the best template fitted for continuum
  Float64 FittedTplAmplitude = NAN;
  // Amplitude error for the best template fitted for continuum
  Float64 FittedTplAmplitudeError = NAN;
  // Chisquare for the best template fitted for continuum
  Float64 FittedTplMerit = NAN;
  // extra chisquare term due to photometry (set to 0.0 if not enabled)
  Float64 FittedTplMeritPhot = NAN;
  // Calzetti ebmvcoeff for the best template fitted for continuum
  Float64 FittedTplEbmvCoeff = NAN;
  // Meiksin igm index for the best template fitted for continuum
  Int32 FittedTplMeiksinIdx = -1;
  // DTM for the best template fitted for continuum
  Float64 FittedTplDtm = NAN;
  // MTM for the best template fitted for continuum
  Float64 FittedTplMtm = NAN;
  // log prior for the best template fitted for continuum
  Float64 FittedTplLogPrior = NAN;
  Float64 FittedTplSNR = NAN;
};
