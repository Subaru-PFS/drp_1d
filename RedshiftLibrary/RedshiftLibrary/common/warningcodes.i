enum class WarningCode {
  WARNING_NONE = 0,
  AIR_VACCUM_CONVERSION_IGNORED,             // 1
  CRANGE_VALUE_OUTSIDERANGE,                 // 2
  CRANGE_VECTBORDERS_OUTSIDERANGE,           // 3
  CRANGE_NO_INTERSECTION,                    // 4
  FINDER_NO_PEAKS,                           // 5
  STDESTIMATION_NO_MATCHING,                 // 6
  STDESTIMATION_FAILED,                      // 7
  MULTIROLL_STRTAG_NOTFOUND,                 // 8
  LINEMATCHING_REACHED_ENDLOOP,              // 9
  IGNORELINESSUPPORT_DISABLED_FFT,           // 10
  FORCE_FROMSPECTRUM_NEG_CONTINUUMAMP,       // 11
  INVALID_MERIT_VALUES,                      // 12
  AIR_VACCUM_REACHED_MAX_ITERATIONS,         // 13
  ASYMFIT_NAN_PARAMS,                        // 14
  DELTAZ_COMPUTATION_FAILED,                 // 15
  INVALID_FOLDER_PATH,                       // 16
  TPL_NAME_EMPTY,                            // 17
  FORCE_NOCONTINUUM_WEAK_CONTINUUMAMP,       // 18
  TEMPLATEFITTINGLOG_NO_MASK,                // 19
  CORRECT_SPECTRUM_NOMINFLUX,                // 20
  DEACTIVATE_CONTREESTIMATION_SKIPSECONDPASS, // 21
  LINEARFIT_RANK_DEFICIENT,                   // 22
  LBFGSPP_ERROR,                              // 23
  WINDOW_TOO_SMALL,                           // 24
  UNUSED_PARAMETER,                           // 24
};
