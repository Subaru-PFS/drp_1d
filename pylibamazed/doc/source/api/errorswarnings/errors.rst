Errors
========

``ErrorCode`` is an enum of the different possible error codes.
All errors beginning with ``IE_`` or named ``INTERNAL_ERROR`` are internal errors and should not appear. If you are confronted to an internal error, don't hesitate to contact the pylibamazed development team.

Possible error codes list and their meaning:

**ATTRIBUTE_NOT_SUPPORTED**
  There is an unsupported attribute in the config file.

**BAD_CALZETTICORR**
  Calzetti correction curve is invalid.

**BAD_CONTINUUMFIT**
  The fitted continuum is considered as "bad" i.e. its value is higher than the specified ``badChi2Threshold``. If ``continuumComponent`` is not set to an auto mode (``tplFitAuto`` / ``powerLawFitAuto``) which allows to switch to ``fromSPectrum`` and filter out the continuum, the fit is not made.

**BAD_FILEFORMAT**
  Generic error for all file format errors.

**BAD_LINE_FORCE**
  Line force is unknown. It must be "S" or "W".

**BAD_LINE_TYPE**
  Line type is unknown. It must be "A" or "E".

**BAD_TEMPLATECATALOG**
  An empty template catalalog has been found for some spectrum model.

**DUPLICATED_LINES**
  Line catalog contains several descriptions of the same line.

**EXTERNAL_LIB_ERROR**
  An error likely caused by an external library. Please contact amazed development team.

**FLAT_ZPDF**
  Pdf is flat.

**IMPORT_ERROR**
  An import error occurred, likely due to a missing dependency to ``scikit-learn`` (reliability stage only).

**INCOHERENT_CONFIG_OPTIONS**
  In configruration file, some elements are conctradictory. It could be:

  - A spectrum model is present in ``linemeascatalog`` but missing from ``linemeas_catalog_columns``
  - For a given spectrum model, linemeas is activated "alone" i.e. without redshift solver, but there is not ``linemeascatalog`` specified in config
  - For a given spectrum model, linemeas is activated "piped" i.e. with redshift solver, but there is ``linemeascatalog`` specified in config

**INCOMPATIBLE_PDF_MODELSHAPES**
  When reliability is activated, redshift range and step from parameters and reliability models  must be identitical.

**INSUFFICIENT_LSF_COVERAGE**
  The LSF is not defined for the whole wavelength range of the spectrum.

**INTERNAL_ERROR / IE_<...>**
  Should not appear. Contact pylibamazed development team.

**INVALID_DIRECTORY**
  A calibration directory indicated in the config or parameters files cannot be accessed.

**INVALID_FILEPATH**
  Raised for any file which cannot be accessed.

**INVALID_LSF**
  There is an error in LSF.

**INVALID_MERIT_VALUES**
  Looping on all templates, could not find one which resulted in a chi2 < INFINITY

**INVALID_NOISE**
  Input spectrum noise is empty or has invalid values.

**INVALID_PARAMETER_FILE**
  The input parameter file is not valid. See error description for more details.

**INVALID_SPECTRUM_FLUX**
  Input spectrum flux is empty, has invalid values, or contains only null values.

**INVALID_SPECTRUM**
  Spectrum canot be validated Causes can be :

  -  axis (flux, error etc) are not of the same size than spectral axis
  - invalid values encountered in spectral axis
  - spectral axis does not contain enough samples (can be after clamping to specified ``lambdaRange``)

**INVALID_SPECTRUM_WAVELENGTH**
  Input spectrum has less than 2 pixels once clamped i.e. reduced to the input lambda range (``lambdaRange`` in parameters file).

**INVALID_WAVELENGTH_RANGE**
  Lambda range is empty (for example when begin and end values are the same), from the input parameters file or once clamped to spectral axis range.

**LESS_OBSERVED_SAMPLES_THAN_AMPLITUDES_TO_FIT**
  When using svd fitter, if there are more elements to fit (number of lines + 1 - for continuum) than samples, the fit cannot be made.

**LINE_CATALOG_ERROR**
  There is an error in the line catalog. See detailed message for more information.

**LINE_NOT_FOUND**
  Could not find the requested line in the line catalog.

**LINE_RATIO_UNKNOWN_LINE**
  A line is missing from line ratio catalog.

**LINEMEAS_CATALOG_ERROR**
  There is an error in the linemeas catalog: some sources are missing.

**LSF_NOT_LOADED**
  If ``lsfType`` is ``fromSpectrumData`` in input parameters, LSF must be loaded before calling ``get_spectrum()`` on a ``AbstractSpectrumReader`` object.

**MISSING_CONFIG_OPTION**
  Some keys are missing from configuration file.

**MISSING_PARAMETER**
  A required parameter is missing from input parameters file.

**MISSING_PHOTOMETRIC_DATA**
  Photometry has ben activated in parameters, but it could not be found in spectrum or a band is missing.

**MISSING_PHOTOMETRIC_TRANSMISSION**
  Photometry has ben activated in parameters, but its transmission could not be found in spectrum or a band is missing.

**NEGATIVE_CONTINUUMFIT**
  A negative continnuum amplitude has been fitted on the input spectrum. If ``continuumComponent`` is not set to an auto mode (``tplFitAuto`` / ``powerLawFitAuto``) which allows to switch to ``fromSPectrum`` and filter out the continuum, the fit cannot be done.

**NO_CLASSIFICATION**
  Classification was not run or failed. It can happen if all redshift solvers failed, or il the classification itself failed.

**NULL_MODEL**
  All amplitudes (continuum & lines) are null at all redshifts.

**OUTPUT_READER_ERROR**
  There was an error when retrieving results from resultStore: the attribute level, name or method is unknown.

**PDF_PEAK_NOT_FOUND**
  A peak could not be identified in the pdf. It can happen in first pass if there is absolutely no peak in the pdf, or in the second pass if there is no peak around any of the first pass candidates.

**PDF_NORMALIZATION_FAILED**
  There was an error during pdf normalization.

**PYTHON_API_ERROR**

**RELIABILITY_NEEDS_TENSORFLOW**
  When reliability is activated, tensorflow must be installed.

**SPECTRUM_CORRECTION_ERROR**
  When correcting input spectrum values (parameter `autoCorrectInput`` set to true), we try for the missing values, to :
  
  - use the lowest flux abs value
  - use the highest noise value

  However, we were unable to find a min flux value or a max noise value.

**SPECTRUM_NOT_LOADED**
  When reading the input spectra, before launching ``get_spectrum()`` on a ``AbstractSpectrumReader`` object, you must first load the spectrum data, for instance using ``load_all()``.
  In the multiobs case, the spectrum items must have the same dimension than and observations. E.g if there are 3 observations, wave, flux etc must have been loaded for the 3 observatuions.

**STAGE_NOT_RUN_BECAUSE_OF_PREVIOUS_FAILURE**
  A stage was not run because a previous stage failed. For example, if the redshift solver failed, the classification stage will not be run.

**TEMPLATE_OVERLAP_TOO_SMALL**
  The overlap between spectrum and template is too small for the z range

**UNALLOWED_DUPLICATES**
  For a given observation, all wavelengths must be unique. This error is raised when the input spectrum contains duplicate wavelengths.

**UNKNOWN_ATTRIBUTE**
  Access to an unkown result has been requested.
