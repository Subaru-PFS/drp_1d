enum class ErrorCode {
  INTERNAL_ERROR = 0,
  EXTERNAL_LIB_ERROR,
  INVALID_SPECTRA_FLUX,
  INVALID_NOISE,
  SMALL_WAVELENGTH_RANGE,
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
  OVERLAPRATE_NOTACCEPTABLE,
  DZ_NOT_COMPUTABLE,
  INCOHERENT_INPUTPARAMETERS,
  BAD_CALZETTICORR,
  // below are api errorCodes
  SPECTRUM_NOT_LOADED,
  MULTILSF_NOT_HANDLED,
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
  PYTHON_API_ERROR
};
