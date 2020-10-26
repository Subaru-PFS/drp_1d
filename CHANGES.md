# Changelog

## 0.14.0 (2020-10-09)
### Added
* [5929] HDF5 output
* [5903] Spectrum model for chisquare2 method
### Changed
* [5971] Continuum templates redshifted using NGP interpolation from a fine grid (computed with cubic spline interpolation) in linemodel method
* [5637] Detect null amplitudes for every z values (flat PDF) and raise an error
* [5539] Option to build static library in cmake command
* [5900] Delta_z computing
### Deprecated
* [5929] CSV output
### Removed
* [5539] `TEST_COVERAGE` and `GPROF` options from cmake command
### Fixed
* [5839] Log messages overflow

## 0.12.0 (2020-08-11)
### Added
* Add the probability of the source to be a star / galaxy / qso
### Changed
* [5898] tellar solver output to comply with galaxy solver
* [5834] Candidate filtering when the solver provide a unique candidate
### Fixed
* [5907] Coredump caused by `skipsecondpass : yes` parameter
 
## 0.10.1 (2020-07-02)
### Fixed
* [5858] Fix memory bug
* [5859] Fix Tplcombination operator

## 0.10.0 (2020-06-23)
### Changed
* [5619] Peak surch
* [5770] Compute reliability with merit
* [5417] Detect negative continuum amplitude 
* [5431] Delta z computation based on the marginalized pdf rather than the secondpass Xi2 results
* [5800] Refactoring on method name in API
### Fixed
* [5772] Pdf integration (especially when Zcand close to Zmin or Zmax)
* [5659] Wrong Lya profile
* [5739] Redshift candidates rank
* [5757] Redshift candidates rank
* [5652] Redshift measurement on close candidates (overlap of integration support)
* [5778] Avoid negative SNR for Ha and OII lines
### Removed
* [5627] Extrema continuum computation for linemodel method

## 0.8.1 (2020-04-03)
### Fixed
* [5610] Make output files: redshift.csv and candidateresults.csv consistent
* [5618] Generating IDs for candidates obtained with methods other than linemodel
* [5620] Add missing column header in candidates.csv
* [5633] Throw error if the number of identified candidates is greater than a constant set to 20 candidates
* [5701] Correction of unit tests in spectralaxis_test.cpp
* [5710] Fix output files consistency (between redshift.csv and candidateresults.csv)

# 0.8.0 (2020-02-27)
## Added
* Add support for IEEE 754 raw floating-point and special values
* Add GSL error handling
## Changed
* Update minimum thirdparties version (`cmake` > 3.0, `boost` > 1.53, `cfistio` > 3.36, `gsl` > 2.5, `fftw` > 3.3.8)
* Python module name change to `pylibamazed`
* Thirdparty script is now compatible with python 2.7
## Fixed
* Fix ranking redshift candidate
* Fix chi2 computation in chisquareloglamba
* Fix clang compiler detection for OSX plateform
* Fix build sequence for python build

## 0.6.1 (2019-04-11)
### Fixed
* added missing CParameterStore::FromString and wrapper

## 0.6.0 (2019-04-11)
### Added
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