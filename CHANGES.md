# Changelog

## 1.6.0 (2025-01-23)
## Added
* [#6715] : Add velocity uncertainty and lambda offset uncertainty
* [#9300] : Template fitting solve : implements the two pass coarse/ fine zgrid fitting
* [#9299] : Add python filename,method,line number in APIexceptions, fix some exception handling
* [#9365] : [api] Add optional class AbstractExternalStorage
* [#9360] : Filters should be defined by observation
* [#8308] : Enable linemeas after template fitting
## Changed
* [#9245] : Switched to pyproject.toml installation and build procedure (compatible with python 3.12)
* [#9255] : Refactor AbstractSpectrumReader
* [#9311] : Makes referencing lib easily optional (if not using json scheme checker)
* [#9033] : Add a minimal number of samples required to process spectra
* [#9366] : Upgrading numpy : small code adaptations for upgrading numpy to 2.0.2
## Removed
* [#9218] : Removes parameters v1 support. Only v2 parameters format are now supported.
* [#9310] : Remove linemodel firstpass b
## Fixed
* [#9295] : Parameters converter : Fixes v1 to v2 parameters conversion
* [#9271] : Fix when few samples in one of the power law
* [#9294] : Handle negative or null spectrum in power Law fitting
* [#9185] : Fix linemeas SNR
* [#9374] : Fix templateFittingSolve with continuum filtered out


## 1.4.0 (2024-11-07)
## Added
* [9058] : Linemeas only for classified object
* [8313] : Power law fit for QSO continuum 
* [5803] : Evaluate the goodness of fit for the continuum, then raise an error if below a given threshold or switch to a median removal
## Fixed
* [9093] : Parameters checker : change default columns name to check in filters
* [9184] : In template ratio modes, all the line flux uncertainties and SNR with Gaussian fit are disabled (set to NAN), and (re-)enabled with direct integration 
* [9223] : Mispelled error code
* [9221] : Fix the estimation of the continuum under each lines
* [7047] : Better track all the conditions when the lines cannot be fitted (abs lines with null continuum, only visible lines with null ratio, ...)  
* [9243] : Should check if the classified type has a lineMeasSolver listed in its stages

## 1.2.0 (2024-06-26)
## Added
* [8312] : Change all the structure of the code to fit jointly several spectra of the same source to a unique model
* [8117] : Implement amplitude uncertainty in lbfgsbfitter
* [8641] : Add new error codes
* [8235] : Set max number of candidates to retain after 1st pass to a dedicated parameter `lineModel.firstPass.extremaCount`
* [8471] : Add new warning output `InitWarningFlags` related to warning raised during spectrum preprocessing stage
## Changed
* [8219] : Changes input parameters file structure. Keeps compatibility with old parameters file
* [8536] : Improve the estimation of the continuum flux level and uncertainty under a line by averaging the fitted continuum weighted by the line profile and using the residual for the uncertainty. This has an impact on the estimated flux level and uncertainty of absorption lines since they depend on the estimated continuum level
* [8470] : Upgrade boost version to 1.74
* [8435] : Change Lya section in parameters 
* [8588] : Clean warnings (enhance naming, convert some warnings into errors, adds documentation)
* [8596] : Adds timestamp in logs
* [8813] : Relax line doublet and triplet constrain in linemeas
* [8734] : Context renamed in ProcessFlow
* [8867] : Rename output "Evidence" to "LogEvidence"
* [8005] : Do not switch to nocontinuum or fromspectrum in second pass continuum refit
* [8528] : Change amazed output access API to retrieve fitted line fluxes, linemeas fluxes and PDF
## Removed
* [8541] : Remove line individual SNR for template ratio fitting
* [6459] : Remove ism extinction application on amplitude of absorption lines
## Fixed
* [8634] : Check spectra is in lambda range before validating its samples
* [8476] : Fix parameter checker
* [8669] : In linemeas, better deals with line offsets to ensure the lines remains in the wavelength range where the LSF is available
* [8664] : Throw amazed exception when catching exceptions in lbfgsbpp
* [8748] : Fix continum template sorting by merit for line model
* [8778] : Fix flux estimation for Lya profile (gaussian+IGM)
* [8441] : Refactor process flow, error and warning handling
* [8812] : Fix IGM fitting and photometry fitting
* [8914] : Fix border effect for next spectrum affecting templates when continuum fitting fails (negative amplitude)

## 1.0.0 (2024-01-17)
## Added
* The API is now checking the syntax of the json input parameter file, raising an error if the syntax is wrong or if some required parameters are missing, issueing a warning for each unused parameter
* The API is now handling one lsf per observation (but the current multi-observation in interleave mode is only able to handle one lsf: the first of the list is chosen)
* In full multi-observation mode, Amazed provides now one model per observation
* "not and" operator have been added to filtering feature
* A warning is raised when lambda range is clamped (it occurs when user lambda range is not compliant with the lambda grid of the input spectrum)
* String vectors type are now available in Amazed output
* The flux uncertainty with direct integration method is now available in Amazed output
* The cumulative SNR for all emission lines is now available in Amazed output
## Changed
* Change line ID as line number in linecatalog, ensuring uniqueness of line name, wavelength (rounded to 2 decimals) and type (E or A)
* Change version of LBFGSpp to v0.3.0
* Changes the way of dealing with redshifts candidates with overlapping integration ranges. In case of overlap greater than 30%, only the best candidates (i.e. highest proba) is selected. In other case, the overlap is equally splitted and candidates are both kept
## Fixed
* The linemeas method is no more executed when redshift solver has failed
* Better protect bad calibration directory content from segfault, and raise informative exceptions
* The output of fitted Line ratio main amplitude of absorption lines is fixed
* The output of fitted Lya asym parameters is fixed
* The output of fitted photometric model is fixed
* The computation of cumulative SNR of strong emmission lines is fixed

## 0.44.0 (2023-09-06)
### Added
* [8033] : Read additional columns from input spectrum data and filter input spectrum based on column values
* [8112] : Add Lya model to amazed output
* [8065] : Add classification warning to amazed output
* [8058] : Add Continuum SNR to amazed output
* [7529] : Add a warning when 2nd pass window size is too small
* [8223] : Raise an error when classification is not run
### Changes
* [8051] : Enhance error message description when spectrum validation fails
* [7706] : Updates rule "strong higher than weak"
* [8048] : Change parameters access. All accesses are made through object. There no dictionary access any more
* [8222] : Change error enumeration
## Removed
* [8120] : End of support for python 3.6
## Fixed
* [8070] : Fix wrong signatures of load_XXX methods in ASCIISpectrumReader
* [8143] : Fix values accessor in AmazedSpectrumWriter
* [8171] : Fix overlapping lines
* [8199] : Fix result specification filtering

## 0.42.0 (2023-04-20)
### Added
* [7875] : Add 38 new solver results in the extended results list
* [7967] : Reliability is now able to handle three classes
* [7961] : Add eigen 3.3 and lbfgspp as thirdparties (requiered for the new line measurement solver based on Gaussian Fit)
* [7892] : Specific dataset result can be excluded from HDF5 output
* [7213] : Template line-ratio amplitude are now available in the results
* [7759] : Implement multi observation for template fitting method solver
* [7378] : Implement new line measurement solver based on full Gaussian fit
### Changes
* [7078] : change resolution conversion to sigma to 1/2.35
### Fixed
* [7983] : Fix `CLineModelElement::getSupport/SubElt()`  : last sample return is skipped by caller
* [7977] : Fix polynomial under lines handling (linemeas)
* [7919] : Fix skipsecondpass
* [7924] : Fix useloglambdasampling
* [7891] : Fix abs mtm tplratio recording
* [7744] : Use InternalException to catch internally
* [7996] : Fix catch exceptions in AbstractOutput.get_attributes
* [7990] : Fix polynome under lines overlapping

## 0.40.0 (2022-12-14)
### Changes
* [7648] : Clean IGM extinction curves and rebuild hires curves
* [7500] : Changed continuum positivity constraint ( amplitude_snr > null_threshold instead of amplitude > 0)
### Fixed
* [7658] : Refactor compressed pdf format and zgrid computation + bug fix
* [7678] : Extend python exceptions handling
* [7717] : Fix error handling at context initialization

## 0.38.0 (2022-11-14)
### Added
* [7460] : Wrap warning codes enum
### Changed
* [6181] : Spectrum processing continues if one solver processing among (galaxy/qso/star) fails
* [7511] : Prevent direct call to `Log.LogError` (raise an error instead) or `Log.LogWarning` (issue a warning flag instead)
* [7532] : PDF compact storage: store linemodel PDF with coarse grid outside 2nd pass windows + parameters to be able to rebuild regular zgrid
### Fixed
* [7481] : Fix incoherence in `continuumComponent` members between `lineRatioManager` and `ContinuumManager`/`linemodelfitting`
* [7510] : Fix side effect bug: chi2 recomputed with same model was different
* [7460] : Check z grid compatibility with neural network learned for reliability
* [7345] : Fix Lya profile using high resolution Meiksin extinction convolved by LSF

## 0.36.0 (2022-09-02)
### Changed
* [7470] : Correct bug on lines selected for amplitude fitting (has impact on qso and linemeas results)
### Removed
* [7434] : Remove `lmfit` fitting method
### Fixed
* [7397] : Refit the continuum.count best templates in precomputeContinuumFit
* [7396] : Fix linemodel exception when continuumfit.count >1 and secondpass.continuumfit="fromfirstpass"
* [5806] : Select best continuum for linemodel free/rules
* [7437] : Remove major quality issues from fully unit-tested classes
* [7445] : Reduce cognitive complexity
* [7459] : Refactor `linemodel` method

## 0.34.0 (2022-06-23)
### Added
* [7106] : Add sub-classification for template ratio method
* [7304] : Extend unit test coverage
* [6564] : Nullify non-significant continuum, based on `nullthreshold` param
### Changed
* [7225] : Refactor API output classes
### Fixed
* [7370] : Fix useless member in reliability

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
