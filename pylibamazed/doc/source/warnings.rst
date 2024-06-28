Warnings
========

Possible warning codes list 

**AIR_VACUUM_CONVERSION_IGNORED**
   An air to vacuum conversion has been defined in the parameters. 
   However, it is specified in the input spectrum that the data are in vacuum.
   The air to vacuum conversion is therefore ignored.

**PDF_PEAK_NOT_FOUND**
   For a given redshift candidate, in the second pass window, no pdf peak has
   been found.

**ESTIMATED_STD_FAR_FROM_INPUT**
   Residual std (obtained from rms) (spectrum - model) very different
   (ratio > 1.5 & 1/ 15 <-> 50 + or 50%-) from input spectrum std

**LINEMATCHING_REACHED_ENDLOOP**
   In linematching, when trying to match lines on the spectrum, the maximum
   number of iterations has been reached without converging.

**FORCED_IGNORELINESUPPORT_TO_FALSE**
   In line model, for continuum fitting, in the parameters fft processing is
   active and ignore line support too (lineModel.continuumFit.ignoreLineSupport).
   However masking lines to fit continuum is in this case impossible.
   To fix this, ignorelinesupport is forced to False i.e.full spectrum with lines
   is used for template fitting.

**FORCED_CONTINUUM_COMPONENT_TO_FROMSPECTRUM**
   In line model, while fitting the continuum, a negative continuum amplitude
   has been found.
   To fix this, continuumComponent is forced to fromSpectrum.

**AIR_VACUUM_REACHED_MAX_ITERATIONS**
  Wavelengths convertion from air to vacuum is made by iterating an inverse
  translation, until a minimal precision is reached or until the max number
  of iterations is reached.
  The marning means that max precision was not reached at the end ot the iterations.

**ASYMFIT_NAN_PARAMS**
   At least one asymetric fitting param is NaN. All parameters of asymetric line
   profile are therefore set to NaN.

**DELTAZ_COMPUTATION_FAILED**
   Some error appeared when trying to calculate Deltaz on a redshift candidate
   probability.

**INVALID_FOLDER_PATH**
   Absolute path to piors calibration file does not exist. No priors used

**FORCED_CONTINUUM_TO_NOCONTINUUM**
   The calculated continuum amplitude is weaker than the min continuum amplitude
   defined in parameters (`continuumFit.nullThreshold`).
   Therefore, continuum is set to noContinuum.

**FORCED_CONTINUUM_REESTIMATION_TO_NO**
  Continuum reestimation was set to "onlyextrema" but it is not possible since
  second pass is skipped. It is automatically set to no continuum reestimation.

**LESS_OBSERVED_SAMPLES_THAN_AMPLITUDES_TO_FIT**
  Deficient rank for lines fitting: when trying to fit lines with line model,
  the number of parameters to find is greater than the number of observed samples.

**LBFGSPP_ERROR**
   There was an error while fitting whith lbfgspp. We then fix line position and
   width and use the initial svd guess for the amplitude.

**PDF_INTEGRATION_WINDOW_TOO_SMALL**
   PDF integration range goes beyond redshift window size. Integration range is
   therefore cropped.
   If this warning appears many times, consider increasing in the parameter
   second pass `halfWindowSize` value.

**UNUSED_PARAMETER**
   A parameter has been defined in the parameters file but it is not used. 
   NB: Some parameters are usefull only under some conditions.

**SPECTRUM_WAVELENGTH_TIGHTER_THAN_PARAM**
   Parameters lambda range goes beyond spectrum range.

**MULTI_OBS_ARBITRARY_LSF**
  Happens in the case of multi obs, for which only one lsf is support.
  Only the lsf of the first observation is taken into account and applied
  to all observations.