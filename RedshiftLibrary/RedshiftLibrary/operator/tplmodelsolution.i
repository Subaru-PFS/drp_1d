struct CTplModelSolution {

  // template continuum
  std::string tplName = undefStr;
  Float64 tplAmplitude = NAN;
  Float64 tplAmplitudeError = NAN;
  Float64 tplAmplitudeSigma = NAN;
  Float64 tplEbmvCoeff = NAN;
  Int32 tplMeiksinIdx = undefIdx;
  Float64 tplRedshift = NAN;

  Float64 tplMerit = NAN;
  Float64 tplMeritPhot = NAN;
  Float64 tplDtM = NAN;
  Float64 tplMtM = NAN;
  Float64 tplLogPrior = 0.;

  // polynom
  TFloat64List pCoeffs;
};
