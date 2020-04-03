# 0.8.1 (2020-04-03)

## Bug Fixes

* Fix 5610 : Make output files: redshift.csv and candidateresults.csv consistent
* Fix 5618 : Generating IDs for candidates obtained with methods other than linemodel
* Fix 5620 : Add missing column header in candidates.csv
* Fix 5633 : Throw error if the number of identified candidates is greater than a constant set to 20 candidates
* Fix 5701 : Correction of unit tests in spectralaxis_test.cpp
* Fix 5710 : Fix output files consistency (between redshift.csv and candidateresults.csv)

# 0.8.0 (2020-02-27)

## New Features

* Add support for IEEE 754 raw floating-point and special values
* Add GSL error handling

## API Changes
None

## Bug Fixes

* Fix ranking redshift candidate
* Fix chi2 computation in chisquareloglamba
* Fix clang compiler detection for OSX plateform
* Fix build sequence for python build

## Other Changes and Additions

* Update minimum thirdparties version (`cmake` > 3.0, `boost` > 1.53, `cfistio` > 3.36, `gsl` > 2.5, `fftw` > 3.3.8)
* Python module name change to `pylibamazed`
* Thirdparty script is now compatible with python 2.7


# 0.6.1 (2019-04-11)

## Bug Fixes

* added missing CParameterStore::FromString and wrapper

# 0.6.0 (2019-04-11)

## New Features

* Added the new option `tplfitauto` to `continuumcomponent`
  parameter. With this option, the code selects the best method
  (`linemodel` or `fullmodel`) for redshift determination depending on
  the signal to noise ratio (SNR) of all template fitting. If the best
  SNR is greater than 50 then the `fullmodel` is selected. Otherwise
  the `linemodel` is selected.

* Added the new option `svdlcp2` to
  `secondpasslcfittingmethod`. Enable second degree polynomial on line
  support on second pass (correct continuum removal residue)

* Implemented new fit options for LyAlpha line: `lyafit.asymfitmin`,
  `lyafit.asymfitmax`, `lyafit.asymfitstep`, `lyafit.widthfitmin`,
  `lyafit.widthfitmax`, `lyafit.widthfitstep`, `lyafit.deltafitmin`,
  `lyafit.deltafitmax`, `lyafit.deltafitstep`

## API Changes
None

## Bug Fixes
None

## Other Changes and Additions
None
