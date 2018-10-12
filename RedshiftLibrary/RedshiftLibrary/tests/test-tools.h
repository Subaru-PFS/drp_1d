#ifndef TEST_TOOLS_H
#define TEST_TOOLS_H

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <boost/filesystem.hpp>

namespace CPFTest
{

enum TGenerationKind
{ FULL, HA, DUMB, // Spectra
  GOOD_SIMPLE, GOOD_FULL, ELRATIO, SUPERSTRONG, // RayCatalog
  BAD_1, BAD_2, BAD_3,
  RANDOM, FLAT // Spectra
};

void generate_spectrum(NSEpic::CSpectrum &spectrum, NSEpic::UInt32 size,
                       NSEpic::Float64 lambda_min, NSEpic::Float64 lambda_max,
                       NSEpic::Float64 mean=0.34,
                       TGenerationKind kind=RANDOM);
boost::filesystem::path generate_spectrum_fits(NSEpic::UInt32 size,
                                               NSEpic::Float64 lambda_min,
                                               NSEpic::Float64 lambda_max,
                                               NSEpic::Float64 mean=0,
                                               TGenerationKind kind=RANDOM);
boost::filesystem::path generate_noise_fits(NSEpic::UInt32 size,
                                            NSEpic::Float64 lambda_min,
                                            NSEpic::Float64 lambda_max,
                                            NSEpic::Float64 mean=0,
                                            TGenerationKind kind=RANDOM);
boost::filesystem::path generate_bogus_fits_file();
boost::filesystem::path generate_linecatalog_file(TGenerationKind kind);
void generate_template(NSEpic::CTemplate &_template, NSEpic::UInt32 size,
                       NSEpic::Float64 lambda_min, NSEpic::Float64 lambda_max);
void generate_template_catalog(NSEpic::CTemplateCatalog &catalog,
                               NSEpic::UInt32 template_size,
                               NSEpic::Float64 lambda_min,
                               NSEpic::Float64 lambda_max);
boost::filesystem::path generate_calibration_dir();
boost::filesystem::path generate_ray_catalog(TGenerationKind kind);

} // namespace CPFTest

#endif
