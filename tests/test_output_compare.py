"""
File: test_csv_comparator.py

Created on: 31/01/19
Author: CeSAM
"""

import os
import csv
import json
import pytest
from tempfile import NamedTemporaryFile, TemporaryDirectory
from pyamazed.output_compare import read_spectrumlist, float_cmp, int_cmp, \
                                    str_cmp, \
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


@pytest.fixture
def dummy_version_json():
    json_file1 = NamedTemporaryFile()
    data = {}
    with open(json_file1.name, 'w') as v1:
        data['cpf-redshift-version'] = '3f9ca6b'
        json.dump(data, v1)
    return json_file1


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


@pytest.fixture
def dummy_csv_res():
    d = TemporaryDirectory()
    redshiftresult_csv_file = open(os.path.join(d.name, "dummy_result_test.csv"), "w")
    csv_writer = csv.writer(redshiftresult_csv_file, delimiter='\t')
    csv_writer.writerow(['#Comment'])
    csv_writer.writerow(['#Key',
                         'float',
                         'str',
                         'int'])
    csv_writer.writerow(['0',
                         '1.0',
                         'tpl_NEW-Im-extended-blue_TF_catalog.txt',
                         '2'])
    csv_writer.writerow(['1',
                         '2.0',
                         'tpl_NEW-Im-g.txt',
                         '3'])
    redshiftresult_csv_file.close()

    return d

###############################################################################
#                              tests functions                                #
###############################################################################

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


def test_read_spectrumlist(spectrum_list):
    """
    Tests the spectrum list reader.
    """
    reader = read_spectrumlist(spectrum_list.name)
    assert type(reader) == list
    for element in reader:
        assert type(element) == str


def test_compare(dummy_version_json):
    jsonObj = JSONComparator(dummy_version_json.name)
    r_equal = jsonObj.compare(dummy_version_json.name, dummy_version_json.name)
    assert r_equal == []
    assert type(r_equal) == list


def test_CSVComparator(dummy_csv_res):
    cmp = DummyResultComparator()
    r = cmp.compare(dummy_csv_res.name, dummy_csv_res.name)

    assert r == []


def test_main():
    pass
