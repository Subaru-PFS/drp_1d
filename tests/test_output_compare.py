"""
File: test_csv_comparator.py

Created on: 31/01/19
Author: CeSAM
"""

import os
import csv
import json
import shutil
import pytest
from tempfile import NamedTemporaryFile, mkdtemp
from astropy.io import fits
import numpy as np
from pyamazed.output_compare import (read_spectrumlist, float_cmp, int_cmp,
                                     str_cmp,
                                     CSVComparator, JSONComparator,
                                     FLOAT_PRECISION,
                                     SpeClassificationCatalogComparator,
                                     SpeZCatalogComparator,
                                     SpePdfCatalogComparator,
                                     SpeLineCatalogComparator,
                                     SpeAbsorptionCatalogComparator,
                                     SpeParametersCatalogComparator)

#
# Fixture functions
#


@pytest.fixture
def spectrum_list():
    """
    Create a spectrum list file for the needs of the tests.
    """
    spectrum_list_file = NamedTemporaryFile()
    with open(spectrum_list_file.name, mode='w') as f:
        csv_writer = csv.writer(f, delimiter='\t')
        csv_writer.writerow(['flux1', 'err1', '0001'])
        csv_writer.writerow(['flux2', 'err2', '0002'])
        csv_writer.writerow(['flux3', 'err3', '0003'])

    return spectrum_list_file


def dummy_json(param):
    """
    Create a json file for the needs of the tests.

    :args: param Chose an arg the create a different json file.
    """
    temp_dir = mkdtemp()
    data = {}
    with open(os.path.join(temp_dir, "config.json"), mode='w') as f:
        data["parameters_file"] = "tests/pfs6b3_127984_0/parameters.json",
        data["config"] = "tests/pfs6b3_127984_0/config.json",
        data["output_folder"] = "results/5e37479c/pfs6b3_127984_0",
        data["input_file"] = "tests/pfs6b3_127984_0/input.spectrumlist",
        data["error_file"] = None,
        data["spectrum_dir"] = "data/pfs6b3_127984_0",
        data["template_dir"] = "data/calibration/templates/BC03_sdss_tremonti21",
        data["linecatalog"] = "data/calibration/linecatalogs/data[linecatalogamazedvacuum_C1_noHepsilon.txt",
        data["calibration_dir"] = "data/calibration/",
        data["zclassifier_dir"] = "data/calibration/reliability/data[zclassifier_C6A3iS1nD9cS1nS1_20180404/",
        data["log_level"] = 70,
        data["linecatalog_convert"] = True,
        data["linemeascatalog"] = "",
        data["save_intermediate_results"] = param
        json.dump(data, f)
    return temp_dir


class DummyResultComparator(CSVComparator):
    filename = "dummy_result_test.csv"
    key = '#Key'
    header_mark = '#Key'
    tests = {
        "#Key": int_cmp,
        "float": float_cmp(FLOAT_PRECISION),
        "str": str_cmp,
        "int": int_cmp
    }


def dummy_csv_res(float_arg):
    temp_dir = mkdtemp()
    res_csv_file = open(os.path.join(temp_dir, "dummy_result_test.csv"), "w")
    csv_writer = csv.writer(res_csv_file, delimiter='\t')
    csv_writer.writerow(['#Comment'])
    csv_writer.writerow(['#Key',
                         'float',
                         'str',
                         'int'])
    csv_writer.writerow(['0',
                         float_arg,
                         'string', #must be a string
                         '2'])
    csv_writer.writerow(['1',
                         '2.0',
                         'tpl_NEW-Im-g.txt',
                         '3'])
    res_csv_file.close()

    return temp_dir


def dummy_fits_file(data):
    """Create a dummy fits file"""
    fits_file = NamedTemporaryFile()
    hdul = fits.HDUList()
    hdul.append(fits.PrimaryHDU())
    for name in data.keys():
        columns = []
        for col in data[name]:
            columns.append(fits.Column(**col))
        hdu = fits.BinTableHDU.from_columns(columns)
        hdu.name = name
        hdul.append(hdu)
    hdul.writeto(fits_file)
    return fits_file

#
# Tests functions
#


def test_read_spectrumlist(spectrum_list):
    """
    Tests the spectrum list reader function.
    """
    reader = read_spectrumlist(spectrum_list.name)
    assert type(reader) == list
    for element in reader:
        assert type(element) == str


def test_compare():
    """
    Tests JSONComparator.compare() method.
    """
    s_val1 = 'all'
    s_val2 = 'none'
    f1 = dummy_json(s_val1)
    f2 = dummy_json(s_val2)
    jsonObj = JSONComparator("config.json")

    # Compare two identical json
    equal = jsonObj.compare(f1, f1)
    assert equal == []

    # Compare two identical json
    diff = jsonObj.compare(f1, f2)
    assert diff == ['config.json/save_intermediate_results : value mismatch ['\
                    + str(s_val1) + '] / [' + str(s_val2) + ']']

    shutil.rmtree(f1)
    shutil.rmtree(f2)


def test_CSVComparator():
    """
    Tests CSVComparator() class.
    """
    f_val1 = 0.001
    f_val2 = 0.002
    f1 = dummy_csv_res(f_val1)
    f2 = dummy_csv_res(f_val2)
    cmp = DummyResultComparator()

    # Compare two identical csv
    r_equal = cmp.compare(f1, f1)
    assert r_equal == []

    # Compare two different csv
    r_diff = cmp.compare(f1, f2)
    assert r_diff == ['0/float : [' + str(f_val1) + '] != [' \
                                    + str(f_val2) + ']']

    shutil.rmtree(f1)
    shutil.rmtree(f2)


def test_true_cmp():
    """
    # NOTE: true_cmp() does not need to be tested
    """
    pass


def test_float_cmp():
    """
    # NOTE: float_cmp() does not need to be tested
    """
    pass


def test_int_cmp():
    """
    # NOTE: int_cmp() does not need to be tested
    """
    pass


def test_speclassificationcatalogcomparator():
    comparator = SpeClassificationCatalogComparator()

    data = {'SPE_CLASSIFICATION_CAT': [
        {'name': 'OBJECT_ID',
         'array': np.array([0, 0, 1, 1]), 'format': 'K'},
        {'name': 'SPE_STAR_PROB',
         'array': np.array([0.5, 0.5, 0.5, 0.5]), 'format': 'E'},
        {'name': 'SPE_GAL_PROB',
         'array': np.array([0.5, 0.5, 0.5, 0.5]), 'format': 'E'},
        {'name': 'SPE_QSO_PROB',
         'array': np.array([0.5, 0.5, 0.5, 0.5]), 'format': 'E'}]}
    f1 = dummy_fits_file(data)
    r = comparator.compare(f1.name, f1.name)
    assert r == []

    data2 = {'SPE_CLASSIFICATION_CAT': [
        {'name': 'OBJECT_ID',
         'array': np.array([0, 0, 1, 1]), 'format': 'K'},
        {'name': 'SPE_STAR_PROB',
         'array': np.array([0.5, 0.5, 0.5, 0.5]), 'format': 'E'},
        {'name': 'SPE_GAL_PROB',
         'array': np.array([0.6, 0.7, 0.8, 0.9]), 'format': 'E'},
        {'name': 'SPE_QSO_PROB',
         'array': np.array([0.5, 0.5, 0.5, 0.5]), 'format': 'E'}]}
    f2 = dummy_fits_file(data2)
    r = comparator.compare(f1.name, f2.name)
    assert r[0].startswith('HDU SPE_CLASSIFICATION_CAT, column SPE_GAL_PROB, '
                           'line 0 : 0.5 != 0.')
    assert r[1].startswith('HDU SPE_CLASSIFICATION_CAT, column SPE_GAL_PROB, '
                           'line 1 : 0.5 != 0.')
    assert r[2].startswith('HDU SPE_CLASSIFICATION_CAT, column SPE_GAL_PROB, '
                           'line 2 : 0.5 != 0.')
    assert r[3] == '...'


def test_spezcatalogcomparator():
    comparator = SpeZCatalogComparator()

    data = {'SPE_REDSHIFT_CAT': [
        {'name': 'OBJECT_ID',
         'array': np.array([1234, 1234, 1234, 1234]), 'format': 'K'},
        {'name': 'SPE_RANK',
         'array': np.array([0, 1, 2, 3]), 'format': 'K'},
        {'name': 'SPE_Z',
         'array': np.array([2.2, 0.47, 0.43, 0.52]), 'format': 'E'},
        {'name': 'SPE_Z_REL',
         'array': np.array([5.1e-2, 3.2e-2, 2.2e-2, 1.2e-2]), 'format': 'E'}]}
    f1 = dummy_fits_file(data)
    r = comparator.compare(f1.name, f1.name)
    assert r == []

    data2 = {'SPE_REDSHIFT_CAT': [
        {'name': 'OBJECT_ID',
         'array': np.array([1234, 1234, 1234, 1234]), 'format': 'K'},
        {'name': 'SPE_RANK',
         'array': np.array([0, 1, 2, 3]), 'format': 'K'},
        {'name': 'SPE_Z',
         'array': np.array([2.2, 0.47, 0.43, 0.52]), 'format': 'E'},
        {'name': 'SPE_Z_REL',
         'array': np.array([5.1e-2, 3.2e-2, 2.3e-2, 1.2e-2]), 'format': 'E'}]}
    f2 = dummy_fits_file(data2)
    r = comparator.compare(f1.name, f2.name)
    assert r[0].startswith('HDU SPE_REDSHIFT_CAT, column SPE_Z_REL, '
                           'line 2 : 0.02')


def test_spepdfcatalogcomparator():
    comparator = SpePdfCatalogComparator()

    data = {'SPE_PDF_CAT': [
        {'name': 'OBJECT_ID',
         'array': np.array([1234, ]), 'format': 'K'},
        {'name': 'SPE_PDF_CLASS',
         'array': np.array(['6']), 'format': 'A'},
        {'name': 'SPE_PDF',
         'array': np.array([[0.1, 1.1, 2.2, 3.3, 4.4]]), 'format': '5E'}]}
    f1 = dummy_fits_file(data)
    r = comparator.compare(f1.name, f1.name)
    assert r == []

    data2 = {'SPE_PDF_CAT': [
        {'name': 'OBJECT_ID',
         'array': np.array([1234, ]), 'format': 'K'},
        {'name': 'SPE_PDF_CLASS',
         'array': np.array(['6']), 'format': 'A'},
        {'name': 'SPE_PDF',
         'array': np.array([[0.1, 1.1, 2.2, 3.3, 4.5]]), 'format': '5E'}]}
    f2 = dummy_fits_file(data2)
    r = comparator.compare(f1.name, f2.name)
    assert len(r) == 1
    assert r[0].startswith('HDU SPE_PDF_CAT, column SPE_PDF, line 0 : [0.1')

    data3 = {'SPE_PDF_CAT': [
        {'name': 'OBJECT_ID',
         'array': np.array([1234, ]), 'format': 'K'},
        {'name': 'SPE_PDF_CLASS',
         'array': np.array(['6']), 'format': 'A'},
        {'name': 'SPE_PDF',
         'array': np.array([[0.1, 1.1 + FLOAT_PRECISION / 1.1,
                             2.2, 3.3,
                             4.4 + FLOAT_PRECISION / 2]]), 'format': '5E'}]}
    f3 = dummy_fits_file(data3)
    r = comparator.compare(f1.name, f3.name)
    assert r == []


def test_spelinecatalogcomparator():
    comparator = SpeLineCatalogComparator()

    data = {'SPE_LINE_FEATURES_CAT': [
        {'name': 'OBJECT_ID',
         'array': np.array([0, 0, 0, 0]), 'format': 'K'},
        {'name': 'SPE_RANK',
         'array': np.array([1337, 1337, 1337, 1337]), 'format': 'I'},
        {'name': 'SPE_LINE_ID',
         'array': np.array([1337, 1337, 1337, 1337]), 'format': 'I'},
        {'name': 'SPE_LINE_NAME',
         'array': np.array([1, 1, 1, 1]), 'format': 'A'},
        {'name': 'SPE_LINE_CENTRAL_WL',
         'array': np.array([1.0, 1.1, 1.2, 1.3]), 'format': 'E'},
        {'name': 'SPE_LINE_CENTRAL_WL_ERR',
         'array': np.array([1.0e-3, 1.1e-3, 1.2e-3, 1.3e-3]), 'format': 'E'},
        {'name': 'SPE_LINE_FLUX',
         'array': np.array([1.0, 1.1, 1.2, 1.3]), 'format': 'E'},
        {'name': 'SPE_LINE_FLUX_ERR',
         'array': np.array([1.0e-3, 1.1e-3, 1.2e-3, 1.3e-3]), 'format': 'E'},
        {'name': 'SPE_LINE_EW',
         'array': np.array([1.0, 1.1, 1.2, 1.3]), 'format': 'E'},
        {'name': 'SPE_LINE_EW_ERR',
         'array': np.array([1.0e-3, 1.1e-3, 1.2e-3, 1.3e-3]), 'format': 'E'},
        {'name': 'SPE_LINE_FWHM',
         'array': np.array([1.0, 1.1, 1.2, 1.3]), 'format': 'E'},
        {'name': 'SPE_LINE_FWHM_ERR',
         'array': np.array([1.0e-3, 1.1e-3, 1.2e-3, 1.3e-3]), 'format': 'E'}]}
    f1 = dummy_fits_file(data)
    r = comparator.compare(f1.name, f1.name)
    assert r == []

    data2 = {'SPE_LINE_FEATURES_CAT': [
        {'name': 'OBJECT_ID',
         'array': np.array([0, 0, 0, 0]), 'format': 'K'},
        {'name': 'SPE_RANK',
         'array': np.array([1337, 1337, 1337, 1337]), 'format': 'I'},
        {'name': 'SPE_LINE_ID',
         'array': np.array([1337, 1337, 1337, 1337]), 'format': 'I'},
        {'name': 'SPE_LINE_NAME',
         'array': np.array([1, 1, 1, 1]), 'format': 'A'},
        {'name': 'SPE_LINE_CENTRAL_WL',
         'array': np.array([1.0, 1.1, 1.2, 1.3]), 'format': 'E'},
        {'name': 'SPE_LINE_CENTRAL_WL_ERR',
         'array': np.array([1.0e-3, 1.1e-3, 1.2e-3, 1.3e-3]), 'format': 'E'},
        {'name': 'SPE_LINE_FLUX',
         'array': np.array([1.0, 1.1, 1.2, 1.3]), 'format': 'E'},
        {'name': 'SPE_LINE_FLUX_ERR',
         'array': np.array([1.0e-3, 1.1e-4, 1.2e-3, 1.3e-3]), 'format': 'E'},
        {'name': 'SPE_LINE_EW',
         'array': np.array([1.0, 1.1, 1.2, 1.3]), 'format': 'E'},
        {'name': 'SPE_LINE_EW_ERR',
         'array': np.array([1.0e-3, 1.1e-3, 1.2e-3, 1.3e-3]), 'format': 'E'},
        {'name': 'SPE_LINE_FWHM',
         'array': np.array([1.0, 1.1, 1.2, 1.3]), 'format': 'E'},
        {'name': 'SPE_LINE_FWHM_ERR',
         'array': np.array([1.0e-3, 1.1e-3, 1.2e-3, 1.3e-3]), 'format': 'E'}]}

    f2 = dummy_fits_file(data2)
    r = comparator.compare(f1.name, f2.name)
    assert len(r) == 1
    assert r[0].startswith('HDU SPE_LINE_FEATURES_CAT, column '
                           'SPE_LINE_FLUX_ERR, line 1 :')


def test_speabsorptioncatalogcomparator():
    comparator = SpeAbsorptionCatalogComparator()

    data = {'SPE_SPECTRAL_FEATURES_CAT': [
        {'name': 'OBJECT_ID',
         'array': np.array([0, 0, 1, 1]), 'format': 'K'},
        {'name': 'SPE_RANK',
         'array': np.array([1337, 1337, 1337, 1337]), 'format': 'I'},
        {'name': 'SPE_ABSORPTION_NAME',
         'array': np.array(['foo', 'bar', 'baz', 'quux']), 'format': '8A'},
        {'name': 'SPE_ABSORPTION',
         'array': np.array([1.0, 1.1, 1.2, 1.3]), 'format': 'E'},
        {'name': 'SPE_ABSORPTION_ERR',
         'array': np.array([1.0e-3, 1.1e-3, 1.2e-3, 1.3e-3]), 'format': 'E'}]}
    f1 = dummy_fits_file(data)
    r = comparator.compare(f1.name, f1.name)
    assert r == []

    data2 = {'SPE_SPECTRAL_FEATURES_CAT': [
        {'name': 'OBJECT_ID',
         'array': np.array([0, 0, 1, 1]), 'format': 'K'},
        {'name': 'SPE_RANK',
         'array': np.array([1337, 1337, 1337, 1337]), 'format': 'I'},
        {'name': 'SPE_ABSORPTION_NAME',
         'array': np.array(['foo', 'bar', 'baz', 'quux']), 'format': '8A'},
        {'name': 'SPE_ABSORPTION',
         'array': np.array([1.0, 1.1, 1.2, 1.4]), 'format': 'E'},
        {'name': 'SPE_ABSORPTION_ERR',
         'array': np.array([1.0e-3, 1.1e-3, 1.2e-3, 1.3e-3]), 'format': 'E'}]}
    f2 = dummy_fits_file(data2)
    r = comparator.compare(f1.name, f2.name)
    assert len(r) == 2
    assert r[0].startswith('HDU SPE_SPECTRAL_FEATURES_CAT, '
                           'column SPE_LINE_ID, line 2 : 1337 != 1338')
    assert r[1].startswith('HDU SPE_SPECTRAL_FEATURES_CAT, column '
                           'SPE_ABSORPTION, line 3 : 1.')


def test_speparamernscatalogcomparator():
    comparator = SpeParametersCatalogComparator()

    data = {'SPE_REST_FRAME_PARAMETERS_CAT': [
        {'name': 'OBJECT_ID',
         'array': np.array([0, 0, 1, 1]), 'format': 'K'},
        {'name': 'SPE_RANK',
         'array': np.array([1337, 1337, 1337, 1337]), 'format': 'I'},
        {'name': 'SPE_SFR',
         'array': np.array([1.0, 1.1, 1.2, 1.3]), 'format': 'E'}]}
    f1 = dummy_fits_file(data)
    r = comparator.compare(f1.name, f1.name)
    assert r == []

    data2 = {'SPE_REST_FRAME_PARAMETERS_CAT': [
        {'name': 'OBJECT_ID',
         'array': np.array([0, 1, 1, 1]), 'format': 'K'},
        {'name': 'SPE_RANK',
         'array': np.array([1337, 1337, 1337, 1337]), 'format': 'I'},
        {'name': 'SPE_SFR',
         'array': np.array([1.0, 1.1, 1.2, 1.3]), 'format': 'E'}]}
    f2 = dummy_fits_file(data2)
    r = comparator.compare(f1.name, f2.name)
    assert len(r) == 1
    assert r[0].startswith('HDU SPE_REST_FRAME_PARAMETERS_CAT, '
                           'column OBJECT_ID, line 1 : 0 != 1')
