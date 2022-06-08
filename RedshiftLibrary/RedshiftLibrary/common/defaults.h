#ifndef _DEFAULT_H
#define _DEFAULT_H

#include "RedshiftLibrary/common/datatypes.h"
#include <gsl/gsl_const_mksa.h>
namespace NSEpic {
static const Float64 N_SIGMA_SUPPORT = 8.;
static const Float64 N_SIGMA_SUPPORT_DI = 6.;

static const Int32 LOG_WORKING_BUFFER_SIZE = 4096;
static const Int32 NOT_OVERLAP_VALUE = 20;
// for CExtremum
static const Int32 PEAKS_MIN_THRESHOLD = 3;
static const Int32 PEAKS_SMOOTH_LIMIT = 20;
static const Int32 undefIdx = -1;

static const Int32 MEDIAN_FAST_OR_BEERS_THRESHOLD = 1000;
static const Float64 SPEED_OF_LIGHT_IN_VACCUM =
    GSL_CONST_MKSA_SPEED_OF_LIGHT / 1000.0; // km.s^-1

static constexpr Float64 INSTRUMENT_RESOLUTION_EMPERICALFACTOR =
    230.0 / 325.0 / 2.35;
static const Float64 RESTLAMBDA_LYA = 1216.;
static const Int32 IGM_OVERSAMPLING = 1;
} // namespace NSEpic
#endif
