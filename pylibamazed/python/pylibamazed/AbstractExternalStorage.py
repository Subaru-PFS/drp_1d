import os

import h5py
import pandas as pd
from astropy.io import fits
from astropy.table import Table

READER_CLASSES = dict()


def register_storage(storage_name, storage):
    READER_CLASSES[storage_name] = storage


def get_storage_from_name(storage_name):
    return READER_CLASSES[storage_name]


class AbstractExternalStorage:
    """
    Class dedicated to opening spectrum files and return their data for the readers
    to load it into themselves.

    Can only work with existing pyamazed readers, which names are listed in the
    constant `READER_CLASSES` in the same module as this class.
    """

    def __init__(self, config, spectrum_id):
        if config.reader not in READER_CLASSES:
            raise Exception(f"Reader class must be one of the following: {READER_CLASSES}")
        self.config = config
        self.spectrum_id = spectrum_id
        self.spectrum_infos = dict()
        self.global_infos = dict()

    #  to be used as context manager
    def __enter__(self):
        obs_id, kwargs = self._read_param
        self.resource = self.read(obs_id, **kwargs)
        return self.resource

    def __call__(self, obs_id="", **kwargs):
        # store read parameters
        self._read_param = (obs_id, kwargs)
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self._read_param = None
        self.close(self.resource)
        self.resource = None
        return False

    def set_spectrum_id(self, spectrum_id):
        self.spectrum_id = spectrum_id

    def read(
        self,
        obs_id: str = "",
        **kwargs,
    ):
        """
        Read a spectrum file and return its data.

        :param obs_id: id of the observation
        :type obs_id: str
        :param kwargs: additional keyword arguments
        :type kwargs: dict

        :return: HDUList or DataFrame
        """
        raise NotImplementedError("Implement in derived class")

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
        return Table.read(filepath, hdu=hdu)

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

    def _get_spectrum_path(self, obs_id=""):
        if self.config.spectrum_path_col:
            if obs_id:
                s_filename = obs_id
            else:
                s_filename = self.spectrum_id.Path
        elif obs_id:
            s_filename = s_filename = (
                self.config.spectrum_prefix
                + self.spectrum_id.ProcessingID
                + "_"
                + obs_id
                + self.config.spectrum_suffix
            )
        else:
            s_filename = (
                self.config.spectrum_prefix + self.spectrum_id.ProcessingID + self.config.spectrum_suffix
            )
        return os.path.join(self.config.spectrum_dir, s_filename)
