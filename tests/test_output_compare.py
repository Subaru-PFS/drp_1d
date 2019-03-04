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
from tempfile import NamedTemporaryFile, TemporaryDirectory, mkdtemp
from pyamazed.output_compare import read_spectrumlist, float_cmp, int_cmp, \
                                    str_cmp, main, \
                                    CSVComparator, JSONComparator,\
                                    FLOAT_PRECISION, PDF_PRECISION

###############################################################################
#                               fixture functions                             #
###############################################################################

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

###############################################################################
#                              tests functions                                #
###############################################################################


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
