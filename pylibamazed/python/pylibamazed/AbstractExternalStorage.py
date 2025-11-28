# ============================================================================
#
# This file is part of: AMAZED
#
# Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
#
# https://www.lam.fr/
#
# This software is a computer program whose purpose is to estimate the
# spectrocopic redshift of astronomical sources (galaxy/quasar/star)
# from there 1D spectrum.
#
# This software is governed by the CeCILL-C license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL-C
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL-C license and that you accept its terms.
# ============================================================================
from abc import ABCMeta, abstractmethod

import h5py
import pandas as pd
from astropy.io import fits
from astropy.table import Table

from pylibamazed.DocDecorator import doc_method

READER_CLASSES = dict()


@doc_method
def register_storage(storage_name, storage):
    READER_CLASSES[storage_name] = storage


@doc_method
def get_storage_from_name(storage_name):
    return READER_CLASSES[storage_name]


class AbstractExternalStorage(metaclass=ABCMeta):
    """
    Class dedicated to opening spectrum files and return their data for the readers
    to load it into themselves.

    Can only work with existing pyamazed readers, which names are listed in the
    constant `READER_CLASSES` in the same module as this class.
    """

    def __init__(self, config):
        if config.reader not in READER_CLASSES:
            raise Exception(f"Reader class must be one of the following: {READER_CLASSES}")
        self.config = config
        self.spectrum_infos = dict()
        self.global_infos = dict()

    #  to be used as context manager
    def __enter__(self):
        spectrum_id, path, obs_id = self._call_params
        del self._call_params
        self.resource = self.read(spectrum_id, path, obs_id)
        return self.resource

    def __call__(self, spectrum_id, path: str = "", obs_id: str = ""):
        # store read parameters
        self._call_params = (spectrum_id, path, obs_id)

        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.close(self.resource)
        self.resource = None
        return False

    @abstractmethod
    @doc_method
    def read(self, spectrum_id: str, path: str, obs_id: str = ""):
        """
        Read a spectrum file and return its data.

        :param spectrun_id: id of the source
        :type spectrum_id: str
        :param path: path or anything else neded to acquire the resource
        :type path: str
        :param obs_id: id of the observation, for multiple observations of the same source
        :type obs_id: str

        :return: resource
        """
        raise NotImplementedError("Implement in derived class")

    @abstractmethod
    @doc_method
    def close(self, resource):
        """
        Close a resource file (if opened).

        :param resource: path of the file to close
        :type resource: any
        """
        raise NotImplementedError("Implement in derived class")

    def _read_fits(self, filepath: str) -> fits.HDUList:
        """
        Read a FITS file and return its HDUList.

        :param filepath: path of the file to open
        :type filepath: str

        :return resource: hdulist contained in the FITS
        :type resource: HDUList
        """
        return fits.open(filepath)

    def _read_fits_table(self, filepath: str, hdu: int) -> Table:
        """
        Read a FITS file and return

        :param filepath: path of the file to open
        :type filepath: str

        :return resource: astropy Table containing fits table at hdu <hdu>
        :type resource: Table
        """
        return Table.read(filepath, hdu=hdu, unit_parse_strict="silent")

    def _read_ascii(self, filepath: str, **kwargs) -> pd.DataFrame:
        """
        Read an ASCII file and return its content as a DataFrame.

        :param filepath: path of the file to open
        :type filepath: str
        :param kwargs: specific arguments for ascii files
        :type kwargs: dict

        :return spectrum: data contained in the ASCII file
        :type spectrum: DataFrame
        """
        spectrum = pd.read_table(filepath, delimiter="\t", **kwargs)
        return spectrum

    def _read_hdf5(self, filepath: str) -> h5py.File:
        """
        Read an HDF5 file and return its content as a h5py File.

        :param filepath: path of the file to open
        :type filepath: str

        :return spectrum: data contained in the HDF5 file
        :type spectrum: h5py.File
        """
        spectrum = h5py.File(filepath, "r")
        return spectrum
