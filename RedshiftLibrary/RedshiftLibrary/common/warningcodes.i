enum class WarningCode {
  WARNING_NONE = 0,
  AIR_VACUUM_CONVERSION_IGNORED,                // 1
  PDF_PEAK_NOT_FOUND,                                  // 2
  ESTIMATED_STD_FAR_FROM_INPUT,                 // 3
  LINEMATCHING_REACHED_ENDLOOP,                 // 4
  FORCED_IGNORELINESUPPORT_TO_FALSE,            // 5
  FORCED_CONTINUUM_COMPONENT_TO_FROMSPECTRUM,   // 6
  AIR_VACUUM_REACHED_MAX_ITERATIONS,            // 7   --
  ASYMFIT_NAN_PARAMS,                           // 8
  DELTAZ_COMPUTATION_FAILED,                    // 9
  INVALID_FOLDER_PATH,                          // 10
  FORCED_CONTINUUM_TO_NOCONTINUUM,              // 11
  FORCED_CONTINUUM_REESTIMATION_TO_NO,          // 12
  LESS_OBSERVED_SAMPLES_THAN_AMPLITUDES_TO_FIT,      // 13
  LBFGSPP_ERROR,                                // 14
  PDF_INTEGRATION_WINDOW_TOO_SMALL,                           
  // Python
  UNUSED_PARAMETER,                             // 15
  SPECTRUM_WAVELENGTH_TIGHTER_THAN_PARAM,       // 16 
  MULTI_OBS_ARBITRARY_LSF                       // 17
};
