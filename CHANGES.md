# 0.6.0 (2019-04-11)

## New Features

* Added the new option `tplfitauto` to `continuumcomponent` parameter. With this option, the code selects the best method (`linemodel` or `fullmodel`) for redshift determination depending on the signal to noise ratio (SNR) of all template fitting. If the best SNR is greater than 50 then the `fullmodel` is selected. Otherwise the `linemodel` is selected
* Added the new option `svdlcp2` to `secondpasslcfittingmethod`. Enable second degree polynomial on line support on second pass (correct continuum removal residue)
* Implemented new fit options for LyAlpha line: `lyafit.asymfitmin`, `lyafit.asymfitmax`, `lyafit.asymfitstep`, `lyafit.widthfitmin`, `lyafit.widthfitmax`, `lyafit.widthfitstep`, `lyafit.deltafitmin`, `lyafit.deltafitmax`, `lyafit.deltafitstep`

## API Changes
None

## Bug Fixes
None

## Other Changes and Additions
None
