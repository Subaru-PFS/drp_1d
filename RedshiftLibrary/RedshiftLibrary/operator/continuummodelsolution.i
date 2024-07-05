struct CContinuumModelSolution {

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
  Float64 tplSNR = NAN;

  // polynom
  TFloat64List pCoeffs;

  // Power law
  Float64 a1 = NAN;
  Float64 a2 = NAN;
  Float64 b1 = NAN;
  Float64 b2 = NAN;
  // TODO see which values to put for coefs uncertainties, if possible
  Int32 MeiksinIdx = undefIdx;
  Float64 EbmvCoeff = NAN;
  Float64 redshift = NAN;
  Float64 SNR = NAN;

  // TODO put in common chi2 and tplMerit (with rename) ?
};
