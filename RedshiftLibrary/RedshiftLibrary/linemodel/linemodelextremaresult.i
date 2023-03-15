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

class TLineModelResult : public TExtremaResult {
public:
  TLineModelResult(const TCandidateZ &candz) : TExtremaResult(candz) {
    m_type = "TLineModelResult";
  }

  void updateFromContinuumModelSolution(std::shared_ptr<const CTplModelSolution> cms);

  void updateFromLineModelSolution(const CLineModelSolution &cms);

  void
  updateTplRatioFromModel(const std::shared_ptr<const CTplratioManager> &lmel);

  void updateFromModel(const std::shared_ptr<const CLineModelFitting> &lmel,
                       const std::shared_ptr<const CLineModelResult> &lmresult,
                       bool estimateLeastSquareFast, int indx);

  Float64 MeritContinuum; // best continum  chi2
  Float64 Merit;          // fullmodel best chi2


  Float64 CorrScaleMarg; // extrema scale marg. correction
  Int32 NDof;            // non zero elements in the lambdarange
  Float64 snrHa;
  Float64 lfHa;
  Float64 snrOII;
  Float64 lfOII;

  Float64 NLinesOverThreshold;

  Float64 StrongELSNR;
  std::unordered_set<std::string> StrongELSNRAboveCut;
  Float64 bic; // bayesian information criterion for each extrema

  Float64 OutsideLinesSTDFlux; // STD measured on the spectrum continuum
                               // substracted outside lines
  Float64
      OutsideLinesSTDError; // STD measured on the error spectrum outside lines

  // line width
  Float64 Elv;            // emission line width
  Float64 Alv;            // absorption line width


  // template ratio
  std::string FittedTplratioName = "undefined"; // Name of the best template
                                                // fitted for tplcorr/tplratio
  Float64 FittedTplratioAmplitudeEm = NAN; // amp of the best template fitted for
                                         // tplcorr/tplratio
  Float64 FittedTplratioAmplitudeAbs = NAN;

  Float64 FittedTplratioIsmCoeff = NAN; // IsmCoeff/EBMV of the best template
                                        // fitted for tplcorr/tplratio
};
