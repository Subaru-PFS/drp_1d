#include <boost/random/normal_distribution.hpp>
#include <boost/random/taus88.hpp>
#include <ctime>
#include <fitsio.h>
#include <fstream>

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/io/fitswriter.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/spectrum/template/template.h>

#include "test-tools.h"

namespace bfs = boost::filesystem;
using namespace std;
using namespace NSEpic;
using namespace CPFTest;

/**
 * Generate a spectrum
 *
 * Based on spectrum1_z_1.2299.fits
 *
 */
void CPFTest::generate_spectrum(CSpectrum &spectrum, UInt32 size,
                                Float64 lambda_min, Float64 lambda_max)
{
    CSpectrumSpectralAxis spectralAxis(size, false);
    CSpectrumFluxAxis fluxAxis(size);

    // boost::random::random_device rng;
    boost::random::taus88 rng;
    boost::random::normal_distribution<double> distribution(0.34, 1.35);

    rng.seed(time(0));
    for (UInt32 i = 0; i < size; i++)
    {
        spectralAxis[i] = lambda_min + double(i) * (lambda_max - lambda_min) / double(size);
        fluxAxis[i] = distribution(rng) * 10;
        fluxAxis.GetError()[i] = 1;
    }
    spectrum.GetSpectralAxis() = spectralAxis;
    spectrum.GetFluxAxis() = fluxAxis;
}

/**
 * Generate a noise fits file
 */
bfs::path CPFTest::generate_noise_fits(UInt32 size,
                                       NSEpic::Float64 lambda_min,
                                       NSEpic::Float64 lambda_max)
{
    CSpectrum noise;
    int status = 0;
    fitsfile *fptr = NULL;
    const char *ttype[2] = {"WAVE", "ERR"};
    const char *tform[2] = {"E", "E"};

    bfs::path tempfile = bfs::unique_path("tst_%%%%%%%%%%-Err.fits");

    generate_spectrum(noise, size, 3800, 12600);

    if (fits_create_file(&fptr, tempfile.c_str(), &status))
    {
        cerr << "Can't fits_create_file : " << status;
        throw runtime_error("Can't fits_create_file");
    }

    if (fits_create_tbl(fptr, BINARY_TBL, noise.GetSampleCount(), 2,
                        (char **)ttype, (char **)tform, NULL, NULL, &status))
    {
        throw runtime_error("Can't fits_create_tbl");
    }

    if (fits_write_col(fptr, TDOUBLE, 1, 1, 1, noise.GetSampleCount(),
                       noise.GetSpectralAxis().GetSamples(), &status))
    {
        throw runtime_error("Can't fits_write_col");
    }

    if (fits_write_col(fptr, TDOUBLE, 2, 1, 1, noise.GetSampleCount(),
                       noise.GetFluxAxis().GetSamples(), &status))
    {
        throw runtime_error("Can't fits_write_col");
    }

    if (fits_close_file(fptr, &status))
    {
        throw runtime_error("Can't fits_close_file");
    }

    return tempfile;
}

/**
 * Generate a spectrum file
 */
bfs::path CPFTest::generate_spectrum_fits(UInt32 size, Float64 lambda_min,
                                          Float64 lambda_max)
{
    CSpectrum spectrum;
    CSpectrumIOFitsWriter writer;
    bfs::path temp = bfs::unique_path("tst_%%%%%%%%%%.fits");

    generate_spectrum(spectrum, size, lambda_min, lambda_max);

    if (!writer.Write(temp.c_str(), spectrum))
    {
        throw runtime_error("Can't create temporary fits file");
    }

    return temp;
}

/**
 * Generate a linecatalog file
 */
bfs::path CPFTest::generate_linecatalog_file(TLineCatalogKind kind)
{
    bfs::path tempfile = bfs::unique_path("tst_%%%%%%%%%%.txt");

    ofstream out(tempfile.c_str(), ofstream::out);

    string content = "#version:0.4.0\t(autoconverted)\n"
                     "#2017-03-06\n";

    switch (kind)
    {
    case FULL:
        content += "#Absorption lines list\n"
                   "#lambda\tName\ttype\tforce\tprofile\tgroup\tnominal_ampl\n"
                   "12821.59\tP5A\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "10941.09\tP6A\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "10052.13\tP7A\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "9548.59\tP8A\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "9231.55\tP9A\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "9017.38\tP10A\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "8865.22\tP11A\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "8664.50\tCaII_t3A\tA\tS\tSYM\t-1\t-1\t-1\n"
                   "8544.42\tCaII_t2A\tA\tS\tSYM\t-1\t-1\t-1\n"
                   "8500.35\tCaII_t1A\tA\tS\tSYM\t-1\t-1\t-1\n"
                   "6564.61\tHalphaA\tA\tS\tSYM\t-1\t-1\t-1\n"
                   "5895.6\tNaD\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "5176.71\tMgI5175\tA\tS\tSYM\t-1\t-1\t-1\n"
                   "4862.72\tHbetaA\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "4341.58\tHgammaA\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "4304.57\tGBand\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "4102.81\tHdeltaA\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "3971.15\tHepsilonA\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "3969.55\tCaII_H\tA\tS\tSYM\tA_Ca\t22\t-1\n"
                   "3934.73\tCaII_K\tA\tS\tSYM\tA_Ca\t23\t-1\n"
                   "3890.11\tH8A\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "3836.43\tH9A\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "3798.93\tH10A\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "3771.65\tH11A\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "3582.19\tFeI\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "2853.73\tMgI2852\tA\tS\tSYM\t-1\t-1\t-1\n"
                   "2804.29\tMgII2803\tA\tS\tSYM\tA_MgII\t0.3054\t-1\n"
                   "2797.11\tMgII2796\tA\tS\tSYM\tA_MgII\t0.6123\t-1\n"
                   "2600.87\tFeII2600\tA\tS\tSYM\t-1\t-1\t-1\n"
                   "2587.35\tFeII2586\tA\tW\tSYM\t-1\t-1\t-1\n"
                   "#Emission lines list\n"
                   "6732.68\t[SII]6731\tE\tW\tSYM\t-1\t-1\t-1\n"
                   "6718.29\t[SII]6716\tE\tW\tSYM\t-1\t-1\t-1\n"
                   "6585.27\t[NII](doublet-1)\tE\tW\tSYM\tE_NII\t2.95\t-1\n"
                   "6564.61\tHalpha\tE\tS\tSYM\t-1\t-1\t-1\n"
                   "6549.86\t[NII](doublet-1/2.95)\tE\tW\tSYM\tE_NII\t1\t-1\n"
                   "5008.24\t[OIII](doublet-1)\tE\tS\tSYM\tE_OIII\t3\t-1\n"
                   "4960.29\t[OIII](doublet-1/3)\tE\tW\tSYM\tE_OIII\t1\t-1\n"
                   "4862.72\tHbeta\tE\tS\tSYM\t-1\t-1\t-1\n"
                   "4341.58\tHgamma\tE\tW\tSYM\t-1\t-1\t-1\n"
                   "4102.81\tHdelta\tE\tW\tSYM\t-1\t-1\t-1\n"
                   "3971.15\tHepsilon\tE\tW\tSYM\t-1\t-1\t-1\n"
                   "3890.11\tH8\tE\tW\tSYM\t-1\t-1\t-1\n"
                   "3869.05\tNeIII\tE\tW\tSYM\t-1\t-1\t-1\n"
                   "3836.43\tH9\tE\tW\tSYM\t-1\t-1\t-1\n"
                   "3798.93\tH10\tE\tW\tSYM\t-1\t-1\t-1\n"
                   "3771.65\tH11\tE\tW\tSYM\t-1\t-1\t-1\n"
                   "3729.88\t[OII]3729\tE\tS\tSYM\t-1\t-1\t-1\n"
                   "3727.09\t[OII]3726\tE\tS\tSYM\t-1\t-1\t-1\n"
                   "3426.73\t[NeVa]\tE\tW\tSYM\t-1\t-1\t-1\n"
                   "3346.81\t[NeVb]\tE\tW\tSYM\t-1\t-1\t-1\n"
                   "2799.12\tMgII\tE\tW\tSYM\t-1\t-1\t-1\n"
                   "2632.82\tFe2632\tE\tW\tSYM\t-1\t-1\t-1\n";
        break;
    case HA:
        content += "6564.61\tHalpha\tE\tS\tSYM\t-1\t-1\t-1\n";
        break;
    case DUMB:
        content += "50000.0\tDumb\tE\tS\tSYM\t-1\t-1\t-1\n";
        break;
    }

    out << content;
    out.close();

    return tempfile;
}

/**
 * Generate a spectrum template
 */
void CPFTest::generate_template(CTemplate &_template, UInt32 size,
                                Float64 lambda_min, Float64 lambda_max)
{
    generate_spectrum((CSpectrum &)_template, size, lambda_min, lambda_max);
}

/**
 * Generate a template catalog
 */
void CPFTest::generate_template_catalog(CTemplateCatalog &catalog,
                                        UInt32 template_size,
                                        Float64 lambda_min, Float64 lambda_max,
                                        const string &category)
{
    shared_ptr<CTemplate> _template;
    for (int i = 0; i < 3; i++)
    {
        string name = "test_template_" + to_string(i);
        _template = (shared_ptr<CTemplate>)(new CTemplate(name, category));
        generate_template(*_template, template_size, lambda_min, lambda_max);
        catalog.Add(_template);
    }
}
