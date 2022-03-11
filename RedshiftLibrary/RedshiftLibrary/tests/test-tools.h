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
#ifndef TEST_TOOLS_H
#define TEST_TOOLS_H

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include <boost/filesystem.hpp>

namespace CPFTest
{

enum TGenerationKind
{ FULL, HA, DUMB, // Spectra
  GOOD_SIMPLE, GOOD_FULL, ELRATIO, SUPERSTRONG, // LineCatalog
  BAD_1, BAD_2, BAD_3,
  RANDOM, FLAT // Spectra
};

void generate_spectrum(NSEpic::CSpectrum &spectrum, NSEpic::Int32 size,
                       NSEpic::Float64 lambda_min, NSEpic::Float64 lambda_max,
                       NSEpic::Float64 mean=0.34,
                       TGenerationKind kind=RANDOM);
boost::filesystem::path generate_spectrum_fits(NSEpic::Int32 size,
                                               NSEpic::Float64 lambda_min,
                                               NSEpic::Float64 lambda_max,
                                               NSEpic::Float64 mean=0,
                                               TGenerationKind kind=RANDOM);
boost::filesystem::path generate_noise_fits(NSEpic::Int32 size,
                                            NSEpic::Float64 lambda_min,
                                            NSEpic::Float64 lambda_max,
                                            NSEpic::Float64 mean=0,
                                            TGenerationKind kind=RANDOM);
boost::filesystem::path generate_bogus_fits_file();
boost::filesystem::path generate_linecatalog_file(TGenerationKind kind);
void generate_template(NSEpic::CTemplate &_template, NSEpic::Int32 size,
                       NSEpic::Float64 lambda_min, NSEpic::Float64 lambda_max);
void generate_template_catalog(NSEpic::CTemplateCatalog &catalog,
                               NSEpic::Int32 template_size,
                               NSEpic::Float64 lambda_min,
                               NSEpic::Float64 lambda_max);
boost::filesystem::path generate_calibration_dir();
boost::filesystem::path generate_line_catalog(TGenerationKind kind);

} // namespace CPFTest

#endif
