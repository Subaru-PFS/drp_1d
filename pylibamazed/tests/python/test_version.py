"""
File: test_version.py

Created on: 22/08/19
Author: CeSAM
"""

from pylibamazed.redshift import get_version


def test_version():
    assert get_version()
    print(get_version())
