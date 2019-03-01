"""
File: test_csv_comparator.py

Created on: 31/01/19
Author: CeSAM
"""

import os
import csv
import json
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
        data['cpf-redshift-version'] = param
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
    jsonObj = JSONComparator("config.json")
    r_diff = jsonObj.compare(dummy_json('a1b2'), dummy_json('c3d4'))
    assert r_diff != []
    assert type(r_diff) == list

    r_equal = jsonObj.compare(dummy_json('a1b2'), dummy_json('a1b2'))
    assert r_equal == []
    assert type(r_equal) == list

def test_CSVComparator():
    """
    Tests CSVComparator() class.
    """
    f_val1 = 0.001
    f_val2 = 0.002
    cmp = DummyResultComparator()

    # Compare two identical csv
    r_equal = cmp.compare(dummy_csv_res(f_val1), dummy_csv_res(f_val1))
    assert r_equal == []

    # Compare two different csv
    r_diff = cmp.compare(dummy_csv_res(f_val1), dummy_csv_res(f_val2))
    assert r_diff != []
    assert r_diff == ['0/float : [' + str(f_val1) + '] != [' \
                                    + str(f_val2) + ']']

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
