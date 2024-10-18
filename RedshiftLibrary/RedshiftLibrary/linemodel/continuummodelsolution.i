struct CContinuumModelSolution {

  // template continuum
  std::string name = undefStr;
  Float64 tplAmplitude = NAN;
  Float64 tplAmplitudeError = NAN;
  Float64 tplAmplitudeSigma = NAN;

  Float64 tplMeritPhot = NAN;
  Float64 tplDtM = NAN;
  Float64 tplMtM = NAN;
  Float64 tplLogPrior = 0.;

  // polynom
  TFloat64List pCoeffs;

  // Power law
  Float64 a1 = NAN;
  Float64 a2 = NAN;
  Float64 b1 = NAN;
  Float64 b2 = NAN;

  Float64 a1std = NAN;
  Float64 a2std = NAN;
  Float64 b1std = NAN;
  Float64 b2std = NAN;

  // Common
  Float64 merit = NAN;
  Float64 reducedChi2 = NAN;
  Float64 SNR = NAN;
  Int32 meiksinIdx = undefIdx;
  Float64 ebmvCoef = NAN;
  Float64 redshift = NAN;

};
