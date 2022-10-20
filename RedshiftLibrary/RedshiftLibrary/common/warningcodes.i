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
    FORCE_LOGSAMPLING_FFT,                     // 10
    IGNORELINESSUPPORT_DISABLED_FFT,           // 11
    FORCE_FROMSPECTRUM_NEG_CONTINUUMAMP,       // 12
    INVALID_MERIT_VALUES,                      // 13
    AIR_VACCUM_REACHED_MAX_ITERATIONS,         // 14
    ASYMFIT_NAN_PARAMS,                        // 15
    DELTAZ_COMPUTATION_FAILED,                 // 16
    INVALID_FOLDER_PATH,                       // 17
    TPL_NAME_EMPTY,                            // 18
    FORCE_NOCONTINUUM_WEAK_CONTINUUMAMP,       // 20
    TEMPLATEFITTINGLOG_NO_MASK,                // 21
    CORRECT_SPECTRUM_NOMINFLUX,                // 22
    DEACTIVATE_CONTREESTIMATION_SKIPSECONDPASS // 23
  };
