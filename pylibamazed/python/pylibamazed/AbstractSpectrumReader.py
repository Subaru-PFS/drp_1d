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

from typing import Dict

import numpy as np
import pandas as pd
from pylibamazed.Spectrum import Spectrum
from pylibamazed.Container import Container
from pylibamazed.Exception import APIException
from pylibamazed.Filter import FilterList
from pylibamazed.FilterLoader import AbstractFilterLoader, ParamJsonFilterLoader
from pylibamazed.Parameters import Parameters
from pylibamazed.lsf import LSFParameters
from pylibamazed.redshift import (
    CFlagWarning,
    CLog,
    ErrorCode,
    WarningCode,
)

zlog = CLog.GetInstance()
zflag = CFlagWarning.GetInstance()

reader_classes = dict()


def register_reader(reader_name: str, reader):
    reader_classes[reader_name] = reader


def get_reader_from_name(reader_name):
    return reader_classes[reader_name]


class AbstractSpectrumReader:
    """
    Base class for spectrum reader, it handles at least wavelengths, flux and error (variance). It also
    handles Light Spread Function (LSF) whether its spectrum dependent or general. It can also handle
    photometric data if present.

    :param observation_id: Observation ID of the spectrum, there can be multiple observation id for the same
        source_id
    :type observation_id: str
    :param parameters: Parameters of the amazed run. It should at least have information about LSF and lambda
        range
    :type parameters: dict
    :param calibration_library: Calibration library object, containing all calibration data necessary,
        according to parameters
    :type calibration_library: class:`CalibrationLibrary`
    :param source_id: Astronomical source ID
    :type source_id: str

    """

    def __init__(
        self,
        parameters: Parameters,
        calibration_library,
        source_id: str,
        filter_loader_class: AbstractFilterLoader = ParamJsonFilterLoader,
    ):
        """Constructor method"""
        self.parameters = parameters
        self.calibration_library = calibration_library
        self.filter_loader_class = filter_loader_class
        self.source_id = str(source_id)

        # Initial loaded data
        self.waves = Container[np.ndarray]()
        self.fluxes = Container[np.ndarray]()
        self.errors = Container[np.ndarray]()
        self.others: Dict[str, Container[np.ndarray]] = dict()
        self.lsf_type = None
        self.lsf_data = Container[np.ndarray]()
        self.photometric_data = []
        self.w_frame = "vacuum"

        # pandas dataframe to instanciate Spectrum
        self.spectra_dataframe = pd.DataFrame()

    # setup context manger to automatize cleaning after getting sepectrum
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.clean()
        return False

    def load_wave(self, resource, obs_id=""):
        """Append the spectral axis in self.wave , units are in Angstrom by default

        :param resource: resource where the wave can be found, no restriction for type (can be a path,
            a file handler, an hdf5 node,...)

        """
        raise NotImplementedError("Implement in derived class")

    def load_flux(self, resource, obs_id=""):
        """Append the flux in self.flux , units are in erg.cm-2 by default

        :param resource: resource where the flux can be found, no restriction for type (can be a path, a file
            handler, an hdf5 node,...)

        """
        raise NotImplementedError("Implement in derived class")

    def load_error(self, resource, obs_id=""):
        """Append the variance in self.error , units are in erg.cm-2 by default

        :param resource: resource where the error can be found, no restriction for type (can be a path, a
            file handler, an hdf5 node,...)

        """
        raise NotImplementedError("Implement in derived class")

    def load_lsf(self, resource, obs_id=""):
        """Append the spectral axis in self.flux , units are in erg.cm-2 by default

        :param resource: resource where the error can be found, no restriction for type (can be a path, a
           file handler, an hdf5 node,...)

        """
        raise NotImplementedError("Implement in derived class")

    def load_photometry(self, resource):
        """Append the photometric data in self.photometric_data , units are in erg.cm-2 by default

        :param resource: resource where the error can be found, no restriction for type (can be a path, a
            file handler, an hdf5 node,...)

        """
        pass  # implementation not mandatory

    def load_others(self, resource, obs_id=""):
        """Appends other data in self.others

        :param resource: resource where other data vector indexed by the wavelengths
        can be found, no restriction for type (can be a path, a file handler, an hdf5 node,...)
        """
        pass  # implementation not mandatory

    def set_air_or_vaccum(self, resource):
        """
        Set w_frame to "air" or "vaccum" (default is vacuum).
        frame should be deduced from resource.
        """
        pass  # implemenation not mandatory

    def load_all(self, resource, obs_id_list=[""]) -> Spectrum:
        """
        Load all components of the spectrum. Reimplement this if resources are different

        :param resource: resource where wave, flux, error, lsf and photometry can be found
               obs_id_list: list of obs id, will loop on them to load all observations,
                            note: only usefull when resource does not need to be updated between observations
        """
        #  on first observation only:
        if self.waves.size() == 0:
            self.set_air_or_vaccum(resource)
            self.load_photometry(resource)
        for obs_id in obs_id_list:
            self.load_wave(resource, obs_id)
            self.load_flux(resource, obs_id)
            self.load_error(resource, obs_id)
            self.load_others(resource, obs_id)
            self.load_lsf(resource, obs_id)

    def load_and_get_spectrum(self, resource, obs_id_list=[""]) -> Spectrum:
        """
        Load all components of the spectrum, build and return Spectrum, then clean memory (re load necessary)
        """
        self.load_all(resource, obs_id_list)
        spectrum = self.get_spectrum()
        self.clean()
        return spectrum

    def get_spectrum(self) -> Spectrum:
        self._check_spectrum_is_loaded()
        if not self._sizes_are_consistent():
            sizes: str = ", ".join([str(container.size()) for container in self._all_containers()])
            raise APIException(
                ErrorCode.INVALID_SPECTRUM,
                f"Number of error, wavelength, flux and other arrays should be the same: {sizes}",
            )

        if not self._obs_ids_are_consistent():
            raise APIException(
                ErrorCode.INVALID_SPECTRUM,
                "Names of error, wavelength, flux and other arrays should be the same",
            )

        if not self.parameters.get_multiobs_method():
            # Add check names if multiobs type is null
            if list(self._get_observation_ids()) != [""]:
                raise APIException(ErrorCode.INVALID_NAME, "Non multi obs observations cannot be named")

        self._merge_spectrum_in_dataframe()

        self._apply_filters(self._get_filters())
        self._check_wavelengths()

        spectrum = Spectrum(
            self.source_id,
            self.parameters,
            self.spectra_dataframe,
            self._select_lsf(),
            self.photometric_data,
            self.w_frame,
        )
        return spectrum

    def clean(self):
        """
        clean the Containers
        """
        self.__init__(
            self.parameters,
            self.calibration_library,
            self.source_id,
            self.filter_loader_class,
        )

    def _get_filters(self):
        return self.filter_loader_class().get_filters(self.parameters)

    def _apply_filters(self, filters: FilterList) -> None:
        if filters is None:
            return
        filters.apply(self.spectra_dataframe)

    def _check_wavelengths(self):
        """Looks if lambda range specified in parameters is contained in spectrum range."""
        for obs_id in self.waves.keys():
            [params_lambda_min, params_lambda_max] = self.parameters.get_lambda_range(obs_id)
            obs_waves = np.array(self.waves.get(obs_id))
            if obs_waves.size == 0:
                raise APIException(ErrorCode.INVALID_SPECTRUM, "Filtered spectrum is empty")
            spectrum_lambda_min = obs_waves[0]
            spectrum_lambda_max = obs_waves[-1]

            params_lambda_range_in_spectrum_lambdas = (
                spectrum_lambda_min <= params_lambda_min and spectrum_lambda_max >= params_lambda_max
            )
            if not params_lambda_range_in_spectrum_lambdas:
                zflag.warning(
                    WarningCode.SPECTRUM_WAVELENGTH_TIGHTER_THAN_PARAM,
                    f"Parameters lambda range ([{params_lambda_min}, {params_lambda_max}])is not "
                    f"contained in spectrum wavelength ([{spectrum_lambda_min}, {spectrum_lambda_max}])",
                )

    def _all_containers(self):
        return [container for container in self.others.values()] + [self.waves, self.fluxes, self.errors]

    def _sizes_are_consistent(self) -> bool:
        expected_size = self.waves.size()
        return all(container.size() == expected_size for container in self._all_containers())

    def _obs_ids_are_consistent(self) -> bool:
        """Checks that all keys from waves container are present in all other containers.
        Must be combined with size test to be sure that all keys are equal
        """
        expected_keys = self._get_observation_ids()
        for obs_id in expected_keys:
            for container in self._all_containers():
                if obs_id not in container.keys():
                    return False
            # check lsf keys
            if (self.lsf_data.size()) and (obs_id not in self.lsf_data.keys()):
                return False
        return True

    def _select_lsf(self):
        parameter_lsf_type = self.parameters.get_lsf_type()

        if parameter_lsf_type == "fromSpectrumData":
            selected_lsf = {"type": self.lsf_type}
            lsf_obs_ids = self.lsf_data.keys()
            if not lsf_obs_ids:
                raise APIException(
                    ErrorCode.LSF_NOT_LOADED,
                    "No LSF loaded in reader, " "lsftype=fromSpectrumData " "parameter cannot be applied",
                )
            obs_id = next(iter(lsf_obs_ids))
            if len(lsf_obs_ids) > 1:
                zflag.warning(
                    WarningCode.MULTI_OBS_ARBITRARY_LSF,
                    f"lsf of observation {obs_id} chosen, other lsf ignored",
                )
            param_name = LSFParameters[selected_lsf["type"]]
            if selected_lsf["type"] == "gaussianVariableWidth":
                selected_lsf["data"] = {
                    "wave": self.lsf_data.get(obs_id)["wave"],
                    "width": self.lsf_data.get(obs_id)["width"],
                }
            else:
                selected_lsf["data"] = {param_name: self.lsf_data.get(obs_id)["width"][0]}

        else:
            selected_lsf = {"type": parameter_lsf_type}
            if parameter_lsf_type == "gaussianVariableWidth":
                selected_lsf["data"] = {
                    "wave": self.calibration_library.lsf["wave"],
                    "width": self.calibration_library.lsf["width"],
                }
            else:
                selected_lsf["data"] = self.parameters.get_lsf()
        return selected_lsf

    def _check_spectrum_is_loaded(self):
        if not self.waves.size():
            raise APIException(ErrorCode.SPECTRUM_NOT_LOADED, "Spectrum not loaded, load_all first")

    def _merge_spectrum_in_dataframe(self):
        frames = dict()
        for obs_id in self._get_observation_ids():
            # Creates dataframe with mandatory columns
            spectrum = pd.DataFrame(
                {
                    "wave": self.waves.get(obs_id),
                    "flux": self.fluxes.get(obs_id),
                    "error": self.errors.get(obs_id),
                }
            )
            for col_key in self.others:
                obs_others = self.others[col_key].get(obs_id)
                if obs_others is not None:
                    spectrum[col_key] = obs_others

            # sort by wavelength
            spectrum.sort_values(["wave"], inplace=True)

            # check unicity of wavelength
            if len(spectrum["wave"].unique()) != spectrum.index.size:
                raise APIException(
                    ErrorCode.UNALLOWED_DUPLICATES, f"Duplicated wavelength values in obs_id {obs_id}"
                )
            frames[obs_id] = spectrum

        self.spectra_dataframe = pd.concat(frames, names=["obs_id", "index"])

    def _get_observation_ids(self):
        return self.waves.keys()
