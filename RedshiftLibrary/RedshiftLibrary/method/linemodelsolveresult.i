class CLineModelSolveResult : public CPdfSolveResult {

public:
  CLineModelSolveResult(
      const std::shared_ptr<const TCandidateZ> &BestExtremumResult,
      const std::string &opt_pdfcombination, Float64 evidence,
      Float64 continuumEvidence = NAN);
  Float64 getContinuumEvidence() const override;
  bool getSwitchedToFromSpectrum() const override;
  void setSwitchedToFromSpectrum(const bool switched);



  std::string tplratioName = undefStr;
  std::string tplContinuumName = undefStr;
  Float64 sigma = NAN;
  Float64 snrHa = NAN;
  Float64 lfHa = NAN;
  Float64 snrOII = NAN;
  Float64 lfOII = NAN;
  Float64 LyaWidthCoeff = NAN;
  Float64 LyaAlpha = NAN;
  Float64 LyaDelta = NAN;
  Int32 LyaIgm = undefIdx;
  bool m_switchedToFromSpectrum = false;
  Float64 m_continuumEvidence =
      NAN; // For cases where switches in fromSpectrum, keeps evidence of
           // continuum for usage in classification
  Float64 minContinuumReducedChi2 = 0.;
};
