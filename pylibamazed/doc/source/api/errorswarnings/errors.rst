Errors
========

``ErrorCode`` is an enum of the different possible error codes.

Possible error codes list (incomplete)

**INVALID_MERIT_VALUES**
  Looping on all templates, could not find one which resulted in a chi2 < INFINITY

**LSF_NOT_LOADED**
If ``lsfType`` is ``fromSpectrumData`` in input parameters, LSF must be loaded before calling ``get_spectrum()`` on a ``AbstractSpectrumReader`` object.

**SPECTRUM_CORRECTION_ERROR**
  When correcting input spectrum values (parameter `autoCorrectInput`` set to true), we try for the missing values, to :
  - use the lowest flux abs value
  - use the highest noise value
  However, we were unable to find a min flux value or a max noise value.

**SPECTRUM_NOT_LOADED**
When reading the input spectra, before launching ``get_spectrum()`` on a ``AbstractSpectrumReader`` object, you must first load the spectrum data, for instance using ``load_all()``.
In the multiobs case, the spectrum items must have the same dimension than and observations. E.g if there are 3 observations, wave, flux etc must have been loaded for the 3 observatuions.

**UNALLOWED_DUPLICATES**
For a given observation, all wavelengths must be unique. This error is raised when the input spectrum contains duplicate wavelengths.

**INVALID_DIRECTORY**
A calibration directory indicated in the config or parameters files cannot be accessed.

**INVALID_FILEPATH**
Raised for any file which cannot be accessed.

**INVALID_PARAMETER**

**MISSING_CONFIG_OPTION**

**BAD_FILEFORMAT**

**INCOHERENT_CONFIG_OPTIONS**

**ATTRIBUTE_NOT_SUPPORTED**

**INCOMPATIBLE_PDF_MODELSHAPES**

**RELIABILITY_NEEDS_TENSORFLOW**

**OUTPUT_READER_ERROR**

**PYTHON_API_ERROR**

**INVALID_NAME**

**INVALID_FILTER_INSTRUCTION**

**INVALID_FILTER_KEY**

**NO_CLASSIFICATION**

**INVALID_PARAMETER_FILE**

**DUPLICATED_LINES**

**STAGE_NOT_RUN_BECAUSE_OF_PREVIOUS_FAILURE**

**LINE_RATIO_UNKNOWN_LINE**