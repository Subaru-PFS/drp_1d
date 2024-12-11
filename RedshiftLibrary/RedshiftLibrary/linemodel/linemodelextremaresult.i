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

  void updateFromContinuumModelSolution(
      std::shared_ptr<const CContinuumModelSolution> cms);
  void updateFromContinuumModelSolution(
      const CContinuumModelSolution &cms);

  void updateFromLineModelSolution(const CLineModelSolution &cms);

  void
  updateTplRatioFromModel(const std::shared_ptr<const CTplratioManager> &ratioMgr);

  void updateFromModel(const std::shared_ptr<const CLineModelFitting> &lmel,
                       const std::shared_ptr<const CLineModelResult> &lmresult,
                       bool estimateLeastSquareFast, int indx);

  Float64 MeritContinuum = NAN; // best continum  chi2
  Float64 Merit = NAN;          // fullmodel best chi2

  Float64 CorrScaleMarg = NAN; // extrema scale marg. correction
  Int32 NDof = undefIdx;            // non zero elements in the lambdarange
  Float64 snrHa = NAN;
  Float64 lfHa = NAN;
  Float64 snrHa_DI = NAN;
  Float64 lfHa_DI = NAN;
  Float64 snrOII = NAN;
  Float64 lfOII = NAN;
  Float64 snrOII_DI = NAN;
  Float64 lfOII_DI = NAN;
  Float64 LyaWidthCoeff = NAN;
  Float64 LyaAlpha = NAN;
  Float64 LyaDelta = NAN;
  Int32 LyaIgm = undefIdx;

  Float64 NLinesOverThreshold = NAN;

  Float64 ELSNR = NAN;
  Float64 StrongELSNR = NAN;
  std::unordered_set<std::string> StrongELSNRAboveCut;
  Float64 bic = NAN; // bayesian information criterion for each extrema

  Float64 OutsideLinesResidualRMS = NAN; // STD measured on the spectrum continuum
                               // substracted outside lines
  Float64
      OutsideLinesInputStDevRMS = NAN; // STD measured on the error spectrum outside lines

  // line width
  Float64 Elv = NAN; // emission line width
  Float64 Alv = NAN; // absorption line width

  CMask OutsideLinesMask;
  // template ratio
  std::string FittedTplratioName = undefStr; // Name of the best template
                                             // fitted for tplcorr/tplratio
  Float64 FittedTplratioAmplitudeEm = NAN;   // amp of the best template fitted
                                             // for tplcorr/tplratio
  Float64 FittedTplratioAmplitudeUncertaintyEm = NAN;
  Float64 FittedTplratioSNREm = NAN;
  Float64 FittedTplratioAmplitudeAbs = NAN;
  Float64 FittedTplratioAmplitudeUncertaintyAbs = NAN;
  Float64 FittedTplratioSNRAbs = NAN;

  Float64 FittedTplratioIsmCoeff = NAN; // IsmCoeff/EBMV of the best template
                                        // fitted for tplcorr/tplratio
};
