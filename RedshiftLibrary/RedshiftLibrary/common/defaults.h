#ifndef _DEFAULT_H
#define _DEFAULT_H

#include <gsl/gsl_const_mksa.h>

#include "RedshiftLibrary/common/datatypes.h"

namespace NSEpic {
static const Float64 N_SIGMA_SUPPORT = 8.;
static const Float64 N_SIGMA_SUPPORT_DI = 6.;

static const Int32 NOT_OVERLAP_VALUE = 20;
// for CExtremum
static const Int32 PEAKS_MIN_THRESHOLD = 3;
static const Int32 PEAKS_SMOOTH_LIMIT = 20;
static const Int32 undefIdx = -1;
static const std::string undefStr = "undefined";

static const Int32 MEDIAN_FAST_OR_BEERS_THRESHOLD = 1000;
static const Float64 SPEED_OF_LIGHT_IN_VACCUM =
    GSL_CONST_MKSA_SPEED_OF_LIGHT / 1000.0; // km.s^-1

// static const Float64 INSTRUMENT_RESOLUTION_FACTOR =
//     230.0 / 325.0 / 2.35; // empirical factor set by A Schmitt
static const Float64 INSTRUMENT_RESOLUTION_FACTOR = 1.0 / 2.355;

static const Float64 RESTLAMBDA_LYA = 1216.;
static const Int32 IGM_OVERSAMPLING = 1;
static const Float64 IGM_RAW_STEP =
    0.05; //  wavelength step of input extinction curves (in Angstrom)

static const Float64 OVERLAP_THRES_HYBRID_FIT =
    0.15; // 15% seemed necessary for Ha/SII complex when lines are very
// wide (either because of PSF or source size)
// mainly for hybrid fitting
static const Int32 MIN_GRID_COUNT = 10;

static const Float64 LSF_MIN_LAMBDA = 200.0;
static const Float64 LSF_MAX_LAMBDA = 30000.0;

static const Float64 OVERLAP_THRESHOLD_PDF_INTEGRATION = 0.3;

static const Int32 MIN_SAMPLE_NUMBER_CONTINUUMM_UNCERTAINTY = 10;

static const Float64 SNR_THRESHOLD_FOR_NLINESOVER = 3.0;

static const Float64 MAX_LAMBDA_OFFSET = 400.0; // km/s
static const Float64 LAMBDA_OFFSET_STEP = 25.0; // km/s;

// For QSO power law calculation, lambda at which the power law coefs changes
static const Float64 POWER_LOW_WAVELENGTH_CUT = 5400;
static const Float64 POWER_LAW_N_SAMPLES_MIN_FOR_CONTINUUM_FIT = 10;

} // namespace NSEpic
#endif
