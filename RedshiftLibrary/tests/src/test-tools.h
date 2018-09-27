#ifndef TEST_TOOLS_H
#define TEST_TOOLS_H

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <boost/filesystem.hpp>

namespace CPFTest
{

enum TLineCatalogKind { FULL, HA, DUMB };

void generate_spectrum(NSEpic::CSpectrum &spectrum, NSEpic::UInt32 size,
                       NSEpic::Float64 lambda_min, NSEpic::Float64 lambda_max);
boost::filesystem::path generate_spectrum_fits(NSEpic::UInt32 size,
                                               NSEpic::Float64 lambda_min,
                                               NSEpic::Float64 lambda_max);
boost::filesystem::path generate_noise_fits(NSEpic::UInt32 size,
                                            NSEpic::Float64 lambda_min,
                                            NSEpic::Float64 lambda_max);
boost::filesystem::path generate_linecatalog_file(TLineCatalogKind kind);
void generate_template(NSEpic::CTemplate &_template, NSEpic::UInt32 size,
                       NSEpic::Float64 lambda_min, NSEpic::Float64 lambda_max);
void generate_template_catalog(NSEpic::CTemplateCatalog &catalog,
                               NSEpic::UInt32 template_size,
                               NSEpic::Float64 lambda_min,
                               NSEpic::Float64 lambda_max,
                               const std::string &category);

} // namespace CPFTest

#endif
