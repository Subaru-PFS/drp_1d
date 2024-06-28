Errors
========

Possible error codes list (incomplete)

**INVALID_MERIT_VALUES**
  Looping on all templates, could not find one which resulted in a chi2 < INFINITY

**TPL_NAME_EMPTY**
  The template name of a continuum model solution is empty (redshift line model solver).

**SPECTRUM_CORRECTION_ERROR**
  When correcting input spectrum values (parameter `autoCorrectInput`` set to true), we try for the missing values, to :
  - use the lowest flux abs value
  - use the highest noise value
  However, we were unable to find a min flux value or a max noise value.
