#include <boost/random/normal_distribution.hpp>
#include <boost/random/taus88.hpp>
#include <ctime>
#include <fitsio.h>
#include <fstream>
#include <cmath>

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/io/fitswriter.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/spectrum/template/template.h>

#include "RedshiftLibrary/tests/test-tools.h"

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
                                Float64 lambda_min, Float64 lambda_max,
                                Float64 mean,
                                TGenerationKind kind)
{
    CSpectrumSpectralAxis spectralAxis(size, false);
    CSpectrumFluxAxis fluxAxis(size);

    boost::random::taus88 rng;
    boost::random::normal_distribution<double> distribution(mean, 1.35);


    switch (kind)
    {
    case RANDOM:
        rng.seed(time(0));
        for (UInt32 i = 0; i < size; i++)
        {
            spectralAxis[i] = lambda_min +
                double(i) * (lambda_max - lambda_min) / double(size);
            fluxAxis[i] = mean + distribution(rng) * 10;
            fluxAxis.GetError()[i] = 1;
        }
        break;
    case FLAT:
        for (UInt32 i = 0; i < size; i++)
        {
            spectralAxis[i] = lambda_min +
                double(i) * (lambda_max - lambda_min) / double(size);
            fluxAxis[i] = mean;
            fluxAxis.GetError()[i] = mean;
        }
        break;
    case SUPERSTRONG:
        for (UInt32 i = 0; i < size; i++)
        {
            spectralAxis[i] = lambda_min +
                double(i) * (lambda_max - lambda_min) / double(size);
            if (i<5042 && i>4975)
            {
                fluxAxis[i] = 30 * pow( (i-5000)/42, 2);
            } else {
                fluxAxis[i] = mean;
            }
            fluxAxis.GetError()[i] = mean;
        }
        break;
    }
    spectrum.GetSpectralAxis() = spectralAxis;
    spectrum.GetFluxAxis() = fluxAxis;

}

/**
 * Generate a noise fits file
 */
bfs::path CPFTest::generate_noise_fits(UInt32 size,
                                       NSEpic::Float64 lambda_min,
                                       NSEpic::Float64 lambda_max,
                                       NSEpic::Float64 mean,
                                       TGenerationKind kind)
{
    CSpectrum noise;
    int status = 0;
    fitsfile *fptr = NULL;
    const char *ttype[2] = {"WAVE", "ERR"};
    const char *tform[2] = {"E", "E"};

    bfs::path tempfile = bfs::unique_path("tst_%%%%%%%%%%-Err.fits");

    generate_spectrum(noise, size, lambda_min, lambda_max, mean, kind);

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
                                          Float64 lambda_max,
                                          Float64 mean,
                                          TGenerationKind kind)
{
    CSpectrum spectrum;
    CSpectrumIOFitsWriter writer;
    bfs::path temp = bfs::unique_path("tst_%%%%%%%%%%.fits");

    generate_spectrum(spectrum, size, lambda_min, lambda_max, mean, kind);

    if (!writer.Write(temp.c_str(), spectrum))
    {
        throw runtime_error("Can't create temporary fits file");
    }

    return temp;
}

/**
 * Generate a bogus fits file
 */
bfs::path CPFTest::generate_bogus_fits_file()
{
    bfs::path temp = bfs::unique_path("tst_%%%%%%%%%%.fits");

    ofstream out(temp.c_str(), ofstream::out);

    out << "Ceci n'est pas un fichier FITS";
    out.close();
    return temp;
}


/**
 * Generate a linecatalog file
 */
bfs::path CPFTest::generate_linecatalog_file(TGenerationKind kind)
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
 * Generate a template catalog category
 */
static void add_template_category(CTemplateCatalog &catalog,
                                  UInt32 count,
                                  UInt32 template_size,
                                  Float64 lambda_min, Float64 lambda_max,
                                  const string &category)
{
    shared_ptr<CTemplate> _template;
    for (int i = 0; i < count; i++)
    {
        string name = category + "_test_template_" + to_string(i) + ".txt";
        _template = (shared_ptr<CTemplate>)(new CTemplate(name, category));
        generate_template(*_template, template_size, lambda_min, lambda_max);
        catalog.Add(_template);
    }
}

/**
 * Generate a template full catalog
 */
void CPFTest::generate_template_catalog(CTemplateCatalog &catalog, UInt32 template_size,
                                        Float64 lambda_min, Float64 lambda_max)
{
    add_template_category(catalog, 3, template_size, lambda_min, lambda_max, "galaxy");
    add_template_category(catalog, 2, template_size, lambda_min, lambda_max, "star");
    add_template_category(catalog, 1, template_size, lambda_min, lambda_max, "qso");
    add_template_category(catalog, 2, template_size, lambda_min, lambda_max, "emission");
}


static void fake_calzetti(bfs::path path)
{
    ofstream out(path.c_str(), ofstream::out);
    out << "#\tlambda\tflux\n"
        "100.0\t50.0\n"
        "10000.0\t1.8776\n"
        "20000.0\t0.49492\n"
        "30000.0\t0.000136928880816\n"
        "40000.0\t0.000117367612128\n"
        "50000.0\t9.780634344e-05\n"
        "60000.0\t7.8245074752e-05\n"
        "70000.0\t5.8683806064e-05\n"
        "80000.0\t3.9122537376e-05\n"
        "90000.0\t1.9561268688e-05\n";
}

static void fake_meiksin(bfs::path path)
{
    vector<string> fileNamesList;
    fileNamesList.push_back("Meiksin_Var_curves_2.0.txt");
    fileNamesList.push_back("Meiksin_Var_curves_2.5.txt");
    fileNamesList.push_back("Meiksin_Var_curves_3.0.txt");
    fileNamesList.push_back("Meiksin_Var_curves_3.5.txt");
    fileNamesList.push_back("Meiksin_Var_curves_4.0.txt");
    fileNamesList.push_back("Meiksin_Var_curves_4.5.txt");
    fileNamesList.push_back("Meiksin_Var_curves_5.0.txt");
    fileNamesList.push_back("Meiksin_Var_curves_5.5.txt");
    fileNamesList.push_back("Meiksin_Var_curves_6.0.txt");
    fileNamesList.push_back("Meiksin_Var_curves_6.5.txt");
    fileNamesList.push_back("Meiksin_Var_curves_7.0.txt");

    for (UInt32 i=0; i<fileNamesList.size(); i++)
    {
        ofstream out((path / fileNamesList[i]).c_str(), ofstream::out);
        for (UInt32 lambda=200; lambda<=1300; lambda++)
        {
            out << lambda << " " <<
                min(1.0, 0.05 + exp(float(i*lambda)/1200)/exp(1)) << " " <<
                min(1.0, 0.15 + exp(float(i*lambda)/1200)/exp(1)) << " " <<
                min(1.0, 0.20 + exp(float(i*lambda)/1200)/exp(1)) << " " <<
                min(1.0, 0.25 + exp(float(i*lambda)/1200)/exp(1)) << " " <<
                min(1.0, 0.30 + exp(float(i*lambda)/1200)/exp(1)) << " " <<
                min(1.0, 0.35 + exp(float(i*lambda)/1200)/exp(1)) << "\n";
        }
    }
}


/**
 * Generate a calibration directory
 */
bfs::path CPFTest::generate_calibration_dir()
{
    bfs::path _path = bfs::unique_path("tst_%%%%%%%%%%");
    bfs::path meiksin_path = _path / "igm" /"IGM_variation_curves_meiksin";
    bfs::create_directories(_path / "ism");
    fake_calzetti(_path / "ism" / "SB_calzetti.dl1.txt");

    bfs::create_directories(meiksin_path);
    fake_meiksin(meiksin_path);

    return _path;
}

/**
 * Generate a ray catalog file
 */
bfs::path CPFTest::generate_ray_catalog(TGenerationKind kind)
{
    string good_simple =
        "#version:0.4.0\t(autoconverted)\n"
        "# Test comment\n"
        "3968.5\tCaH\tA\t-1\t-1\t-1\t-1\t-1\n"
        "5175    MgI\tA\t-1\t-1\t-1\t-1\t-1\n"
        "10320   [SII]\t\t\tE\tW\t-1\t-1\t-1\t-1\n"
        "6562.8  Halpha\t\t\tE\tS\t-1\t-1\t-1\t-1\n"
        "4101.7  Hdelta\t\t\tE\tW\t-1\t-1\t-1\t-1\n";

    string good_full =
        "#version:0.4.0\n"
        "# Absorption lines list\n"
        "3968.5\tCaH\tA\t-1\t-1\t-1\t-1\t-1\n"
        "3933.7\tCaK\tA\t-1\t-1\t-1\t-1\t-1\n"
        "4304.4\tGband\tA\t-1\t-1\t-1\t-1\t-1\n"
        "4861.3  Hbeta\tA\t-1\t-1\t-1\t-1\t-1\n"
        "4340.4  Hgamma\tA\t-1\t-1\t-1\t-1\t-1\n"
        "4101.7  Hdelta\tA\t-1\t-1\t-1\t-1\t-1\n"
        "5175    MgI\tA\t-1\t-1\t-1\t-1\t-1\n"
        "5269    CaFe\tA\t-1\t-1\t-1\t-1\t-1\n"
        "1549.0  CIV1548 A\t-1\t-1\t-1\t-1\t-1\n"
        "1303.0  OI1302  A\t-1\t-1\t-1\t-1\t-1\n"
        "1397.0  SiV     A\t-1\t-1\t-1\t-1\t-1\n"
        "2799.0  MgII\tA\t-1\t-1\t-1\t-1\t-1\n"
        "# Emission lines list\n"
        "10320   [SII]\t\t\tE\tW\t-1\t-1\t-1\t-1\n"
        "6725.0  [SII]doublet\t\tE\tW\t-1\t-1\t-1\t-1\n"
        "6562.8  Halpha\t\t\tE\tS\t-1\t-1\t-1\t-1\n"
        "5006.8  [OIIIa]\t\t\tE\tS\t-1\t-1\t-1\t-1\n"
        "4958.9  [OIIIb]\t\t\tE\tS\t-1\t-1\t-1\t-1\n"
        "4861.3  Hbeta\t\t\tE\tS\t-1\t-1\t-1\t-1\n"
        "3727.5  [OII]\t\t\tE\tS\t-1\t-1\t-1\t-1\n"
        "3425.8  [NeVa]\t\t\tE\tW\t-1\t-1\t-1\t-1\n"
        "3345.9  [NeVb]\t\t\tE\tW\t-1\t-1\t-1\t-1\n"
        "2799.0  MgII\t\t\tE\tW\t-1\t-1\t-1\t-1\n"
        "1909.0  [CIII]\t\t\tE\tW\t-1\t-1\t-1\t-1\n"
        "1549.0  CIV1548+1550\t\tE\tW\t-1\t-1\t-1\t-1\n"
        "1215.7  LyA\t\t\tE\tS\t-1\t-1\t-1\t-1\n";

    string elratio = "#version:0.4.0\t(autoconverted)\n"
        "#Emission\tlines\tlist\n"
        "#lambda\tName\ttype\tforce\tprofile\tgroup\tnominal_ampl\n"
        "7600.00\t[OII]3729\tE\tS\tSYM\t-1\t-1\t-1\n"
        "7400.00\t[OII]3726\tE\tS\tSYM\t-1\t-1\t-1\n";

    string superstrong = "#version:0.4.0\t(autoconverted)\n"
        "# Emission lines list\n"
        "5008.24\t[OIII](doublet-1)\tE\tS\tSYM\tE_OIII\t3\t-1\n"
        "4960.29\t[OIII](doublet-1/3)\tE\tS\tSYM\tE_OIII\t1\t-1\n"
        "3729.88\t[OII]3729\tE\tS\tSYM\t-1\t-1\t-1\n"
        "3727.09\t[OII]3726\tE\tS\tSYM\t-1\t-1\t-1\n";

    string bad_1 =
        "3sdqgdsfg968.5\tMgI\tA\n"
        "5175\tMgI\tA\n";

    string bad_2 =
        "#version:0.4.0"
        "3968.5"
        "5175\tMgI\tA";

    string bad_3 =
        "#invalid because no version found at the beginning of the file"
        "3968.5\tCaH\tA\t-1\t-1\t-1\t-1\t-1"
        "5175\tMgI\tA\t-1\t-1\t-1\t-1\t-1";
    string content;

    bfs::path _path = bfs::unique_path("tst_%%%%%%%%%%.txt");

    switch (kind)
    {
    case GOOD_SIMPLE:
        content = good_simple;
        break;
    case GOOD_FULL:
        content = good_full;
        break;
    case ELRATIO:
        content = elratio;
        break;
    case SUPERSTRONG:
        content = superstrong;
        break;
    case BAD_1:
        content = bad_1;
        break;
    case BAD_2:
        content = bad_2;
        break;
    case BAD_3:
        content = bad_3;
        break;
    }
    
    ofstream out(_path.c_str(), ofstream::out);
    out << content;
    out.close();

    return _path;

}
