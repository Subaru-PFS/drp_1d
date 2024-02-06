enum class ErrorCode {
  INTERNAL_ERROR = 1,
  EXTERNAL_LIB_ERROR,
  INVALID_SPECTRUM_WAVELENGTH,
  INVALID_SPECTRUM_FLUX,
  INVALID_NOISE,
  INVALID_WAVELENGTH_RANGE,
  NEGATIVE_CONTINUUMFIT,
  BAD_CONTINUUMFIT,
  NULL_AMPLITUDES,
  PEAK_NOT_FOUND_PDF,
  MAX_AT_BORDER_PDF,
  MISSING_PARAMETER,
  BAD_PARAMETER_VALUE,
  UNKNOWN_ATTRIBUTE,
  BAD_LINECATALOG,
  BAD_LOGSAMPLEDSPECTRUM,
  BAD_COUNTMATCH,
  BAD_TEMPLATECATALOG,
  INVALID_SPECTRUM,
  OVERLAPFRACTION_NOTACCEPTABLE,
  DZ_NOT_COMPUTABLE,
  INCOHERENT_INPUTPARAMETERS,
  BAD_CALZETTICORR,
  // below are api errorCodes
  SPECTRUM_NOT_LOADED,
  LSF_NOT_LOADED,
  UNALLOWED_DUPLICATES,
  UNSORTED_ARRAY,
  INVALID_DIRECTORY,
  INVALID_FILEPATH,
  INVALID_PARAMETER,
  MISSING_CONFIG_OPTION,
  BAD_FILEFORMAT,
  INCOHERENT_CONFIG_OPTIONS,
  ATTRIBUTE_NOT_SUPPORTED,
  INCOMPATIBLE_PDF_MODELSHAPES,
  UNKNOWN_RESULT_TYPE,
  RELIABILITY_NEEDS_TENSORFLOW,
  OUTPUT_READER_ERROR,
  PYTHON_API_ERROR,
  INVALID_NAME,
  INVALID_FILTER_INSTRUCTION,
  INVALID_FILTER_KEY,
  NO_CLASSIFICATION,
  INVALID_PARAMETER_FILE,
    DUPLICATED_LINES,
  STAGE_NOT_RUN_BECAUSE_OF_PREVIOUS_FAILURE
};
