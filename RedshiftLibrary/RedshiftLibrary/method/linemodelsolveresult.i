class CLineModelSolveResult : public CPdfSolveResult {

public:
  CLineModelSolveResult(
      const std::shared_ptr<const TCandidateZ> &BestExtremumResult,
      const std::string &opt_pdfcombination, Float64 evidence,
      Float64 continuumEvidence = NAN);
  Float64 getContinuumEvidence() const override;
  bool getSwitchedToFromSpectrum() const override;
  void setSwitchedToFromSpectrum(const bool switched);
  bool m_switchedToFromSpectrum = false;
  Float64 m_continuumEvidence =
      NAN; // For cases where switches in fromSpectrum, keeps evidence of
           // continuum for usage in classification
  Float64 minContinuumReducedChi2 = 0;
  Float64 maxFitAmplitudeSigma = 0;
  Float64 maxPValue = 0;
};
