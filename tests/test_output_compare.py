"""
File: test_csv_comparator.py

Created on: 31/01/19
Author: CeSAM
"""

import csv
import pytest
import tempfile
import csv_comparator
from csv_comparator import read_spectrumlist

def test_int_cmp():
    pass

def test_float_cmp():
    pass

def test_reader():
    """
    This function tests the spectrum list reader.
    """
    spectrumlist = NamedTemporaryFile()
    table = [["flux1", "err1", "0001"],
             ["flux2", "err2", "0002"],
             ["flux3", "err3", "0003"]]
    with open(spectrumlist, 'w') as f:
        fake_writer = csv.writer(f)
        fake_writer.writerows(table)
    spectrumlist.close()
    #TODO: suite

class TestOutputDirComparator(object):

    def test_compare(self):
        pass


class TestCSVComparator(object):

    def test_read_csv(self):
        pass

    def test_compare(self):
        pass
