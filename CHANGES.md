# Changelog

## 0.32.0 (2022-05-31)
### Added
* [7093] : Centralize numerical constants in `defaults.h`
* [7221] : Extend unit test coverage
* [6681] : Add optional extended results output
### Changed
* [7123] : `CProcessFlowContext` is singleton
* [7173] : Substract polynomial from continuum for flux `directIntegration`
* [7110] : Warning flags are stored by object category in HDF5 file
* [7133] : Unify all error handling in the library (C++ and python)
* [7127] : Move Meiksin file loading from C++ to python
* [6882] : All C++ code reformatted using clang-format tool
* [5645] : Apply IGM on lines below 1216A
### Removed
* [7029] : Clean unused files
### Fixed
* [7209] : Use disctinct line catalogs for line measurement and z estimation
* [7309] : Fix wrong polynomial index and raise exception if invalid index
* [7149] : Fix IGM above last IGM z bin in `fftprocessing`
* [7295] : Memory management and profiling
* [7327] : Correct heap overflow bugs

## 0.30.0 (2022-03-11)
### Added
* [6754] : Add generic object categories
* [7027] : Save new output results for linemeas method
* [7056] : Save first pass PDF
* [6716] : Add account for polynom contribution to continuum flux
* [6683] : Add linemodel first pass data to hdf5 result file
* [6724] : Add linemeas for QSO
* [5638] : Add processing warning flags
### Changed
* [6499] : Optimizing processing time, when using hAlpha priors
* [6965] : Change processflow context API
* [5811] : Change invalid values for output results (mainly from -1 to NaN)
* [6856] : Change line catalog and line ratio catalog file formats
* [6963] : Change reliability model
### Fixed
* [6982] : Fix photometry flux unit (compute photometric flux in density/Hz)
* [5936] : Fix heap-buffer-overflow in fluxes computation
* [7019] : Fix best model solution for Linemeas, by selecting the minimum Chi2 value
* [6893] : Fix memory leak in raydetection
* [6503] : Fix star/qso/galaxy classification

## 0.28.2 (2022-01-21)
### Fixed
* [6986] : Fix photometric data loading

## 0.28.1 (2022-01-13)
### Fixed
* [6841] : Fix LSF type in API
* [6970] : Fix python Unit Test

## 0.28.0 (2021-12-14)
### Added
* [6783] : Add `AbstractSpectrumReader`, abstract class for writing spectrum reader adaptated for user instrument file.
* [6783] : Add `CalibrationLibrary`class for loading calibration files.
* [6783] : Add basic multi-observation processing. Several spectra can be provided for one source. Wavelength ranges of spectra can overlaped or be disjoined
* [6616] : Add option to use photometric data in the chisquare for TemplatefittingSolve & LineModelSolve 
### Changed
* None
### Deprecated
* None
### Removed
* None
### Fixed
* [6785] : Take into account the input lambdaRange for cropping input log sampled spectra
* [6872] : Fix conditions to orthogonalize all templates from both QSO and galaxy objects
* [6873] : Limit LSF access to observable lambda range

## 0.26.0 (2021-10-25)
### Added
* [5994] : Add option to logrebin input spectra for linemodel processing
* [5775] : Add Air to Vacuum conversion in `spectralAxis`
* [6489] : Activate fftprocessing for qso
* [6687] : Add `photometricdata` and `photometricband` classes
### Changed
* [6701] : Replace `yes` / `no` with `True` / `False` booleans for amazed params
* [6379] : Replace runtime_exception with throw GlobalException
* [6773] : Reactivate the template continuum removal feature
### Deprecated
* None
### Removed
* [6686] : Remove `CDataStore` and `zqual` reliability classes
* [6684] : Remove saving `CPdfCandidateszResult` from `ResultStore` & unnecessary operator result `CModelFittingResult`
### Fixed
* [6695] : Use all ISM coefficients when fitting ISM
* [6765] : Fix crash with empty spectrum
* Miscellaneous bug fixes and code quality improvement

## 0.24.0 (2021-09-01)
### Added
* [6520] : Add a new class for LSF with wavelength dependent width
* [5905] : Add line measurement feature
### Changed
* [6123] : Set default meiksinIdx to -1 when IGM does not apply
* [6659] : Set required minimum swif version to 4.0
### Deprecated
* None
### Removed
* None
### Fixed
* [6657] : Fix SNR of best model fitted in `templatecombination` method
* [6627] : Fix lambda offsets activation
* Miscellaneous bug fixes and code quality improvement

## 0.22.0 (2021-05-07)
### Added
* [5644] : Adding IGM convolution with LSF
* [6469] : Add coverage option in CMakeLists
* [6329] : Activate ISM/IGM on template combination method (PCA)
* [6086] : Save spectrum model and corresponding parameters for template combination method
* [6483] : Add qso linemodel, move `linecatalog` parameter inside `<qso|galaxy>.linemodelsolve.linemodel.linemodelcatalog`
### Changed
* [6343] : `ResultStore` getters changed in python API classes
* [6481] : Disable reading `templatefittingsolve.extinction` parameter for stars
### Deprecated
* None
### Removed
* None
### Fixed
* [6535] : Make functional both template combination method and operator
* [6512] : Fix bug on `tplfitauto` in second pass

## 0.20.1 (2021-05-19)
### Fixed
* [6494] : Rebin only galaxy templates if enablestar is true

## 0.20.0 (2021-05-07)
### Added
* [6061] : Add QSO section in HDF5 output 			
* [6400] : Add new parameter: `linemodel.continuumfit.negativethreshold`, referring to the maximal acceptable fitted amplitude of the continuum (unit: sigma of fitted amplitude, put the minus sign for a negative threshold)
* [5989] : Adapt the library to take log-sampled spectra as inputs for all methods (`linemodelsolve`, `templatefittingsolve`, `tplcombinationsolve`). Replace parameter `linemodelsolve.linemodel.firstpass.largegridstep` with `linemodelsolve.linemodel.firstpass.largegridstepratio`
* [6426] : Add new parameters ebmv.start, ebmv.step and ebmv.count to set the sampling of dust extinction
### Changed
* [6204] : Change CProcessFlowContext method definition
* [6382] : Change version number strategy
* [6425] : Change parameter name from `linemodelsolve.linemodel.continuumfit.method` to `linemodelsolve.linemodel.continuumfit.fftprocessing` (`yes`/`no`)
* [6265] : Merge templatefittinglogsolve method into templatefittingsolve method with a new parameter: `fftprocessing` (`yes`/`no`)
* [6387] : Create dedicated method for computing reliability
* [5933] : Change `gslcBlas` third party dependency to `openBlas` third party dependency
### Deprecated
* None
### Removed
* [6433] : Remove csv ouput file writing
### Fixed
* [6356] : Fix critical bug in templateFittingLog method
* [6392] : Fix bad logger level when storing 1st pass results
* [6434] : Restore TemplateCombination Method
* [6474] : Fix spectral axis of output spectrum-model (templatefittingsolve only): shifting axis from restframe to selected redshift

## 0.18.1 (2021-03-19)
### Fixed
* [6361] : Fix bug with linemodelsolve.linemodel.continuumcomponent=fromspectrum parameter

## 0.18.0 (2021-02-26)
### Added
* [6120] : Add Line Spread Function	
* [6167] : Add `redshiftrange` and `redshiftstep` parameters for stellar solver to tune the redshift computationnal grid
* [6133] : Add precompute continuum fit in the 2nd pass of linemodel method
* [6290] : Add firstpass data to HDF5 outputs

### Changed
* [6261] : Rename chisquare2 method into templatefitting
* [5815] : New design of PDF computation and processing
* [6039] : New loading method of template catalogs
* [6163] : New exceptions management 

### Deprecated

### Removed
* [6146] : Remove deprecated computeContinuumStat parameter
* [6167] : Remove deprecated CMethods and COperators

### Fixed
* [6018] : Fix negative evidence values
* [6328] : Correct bug in star models

## 0.16.0 (2020-12-15)
### Added
* [5997] : Add model parameters into HDF5 output for chisquare2 method
* [6057] : Add ExtremaRedshiftSeparation parameter for candidate finder
* [6057] : Add halfwindowsize parameter (size of the redshift half-range) for the 2nd pass in linemodel method
* [6027] : Add continuumfit.method in linemodel, chisquare2 and chisquarelog, methods
* [6125] : Add pylibamazed python API
## Changed
* [5932] : Disable default AVX and SSE build option for third parties
* [6043] : Upgrade minimum Boost version to 1.57
* [6051] : Avoid strict overlap of 2nd pass redshfit windows in linemodel method
* [5976] : Enable method selecting for QSO solver
## Removed
* [5648] : Remove systematic continuum estimation
* [5753] : Remove unused csv output files
* [6129] : Remove extinction parameter (IGM) for stellar solver
## Fixed
* [6130] : Fix core segmentation error in the 2nd of the linemodel method when called with only one redshift
* [6127] : Fix slightly randomly changing pdf for hight redshift galaxies on different compilers
* [6131] : Fix wrong redshift values in HDF5 file
* [6129] : Fix memory leak from one spectra to the other

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