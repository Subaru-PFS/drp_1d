enum class WarningCode {
  AIR_VACUUM_CONVERSION_IGNORED = 0,            // 0 (1st bit, ie bit number 0 or 2**0)
  PDF_PEAK_NOT_FOUND,                           // 1 (2nd bit, ie bit number 1 or 2**1)
  ESTIMATED_STD_FAR_FROM_INPUT,                 // 2 (etc...)
  LINEMATCHING_REACHED_ENDLOOP,                 // 3
  FORCED_IGNORELINESUPPORT_TO_FALSE,            // 4
  FORCED_CONTINUUM_COMPONENT_TO_FROMSPECTRUM,   // 5
  AIR_VACUUM_REACHED_MAX_ITERATIONS,            // 6
  ASYMFIT_NAN_PARAMS,                           // 7
  DELTAZ_COMPUTATION_FAILED,                    // 8
  INVALID_FOLDER_PATH,                          // 9
  FORCED_CONTINUUM_TO_NOCONTINUUM,              // 10
  FORCED_CONTINUUM_REESTIMATION_TO_NO,          // 11
  LESS_OBSERVED_SAMPLES_THAN_AMPLITUDES_TO_FIT, // 12
  LBFGSPP_ERROR,                                // 13
  PDF_INTEGRATION_WINDOW_TOO_SMALL,             // 14  
  FORCED_POWERLAW_TO_ZERO,                      // 15
            
  // Python
  UNUSED_PARAMETER,                             // 16
  SPECTRUM_WAVELENGTH_TIGHTER_THAN_PARAM,       // 17 
  MULTI_OBS_ARBITRARY_LSF                       // 18
};
