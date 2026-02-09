// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
//
// https://www.lam.fr/
//
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
//
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use,
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info".
//
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability.
//
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or
// data to be ensured and,  more generally, to use and operate it in the
// same conditions as regards security.
//
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================

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
static const Int32 allIdx = -9;
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

static const Int32 MIN_GRID_COUNT = 10;

static const Float64 LSF_MIN_LAMBDA = 200.0;
static const Float64 LSF_MAX_LAMBDA = 30000.0;

static const Float64 OVERLAP_THRESHOLD_PDF_INTEGRATION = 0.3;

static const Int32 RMS_MIN_SAMPLE_NUMBER = 10;

static const Float64 SNR_THRESHOLD_FOR_NLINESOVER = 3.0;

static const Float64 MAX_LAMBDA_OFFSET = 400.0; // km/s
static const Float64 LAMBDA_OFFSET_STEP = 25.0; // km/s;

// For QSO power law calculation, lambda at which the power law coefs changes
static const Float64 POWER_LOW_WAVELENGTH_CUT = 5400;

// For CContinuumIrregularSamplingMedian, number of smoothing cycles
static const Float64 N_SAMPLING_SMOOTH_CYCLES = 5;

} // namespace NSEpic
#endif
