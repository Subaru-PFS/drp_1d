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
from pylibamazed.Container import Container
from pylibamazed.Exception import APIException
from pylibamazed.Filter import FilterList
from pylibamazed.FilterLoader import (AbstractFilterLoader,
                                      ParamJsonFilterLoader)
from pylibamazed.lsf import LSFParameters, TLSFArgumentsCtor
from pylibamazed.Parameters import Parameters
from pylibamazed.redshift import (CFlagWarning, CLog, CLSFFactory,
                                  CPhotometricData, CProcessFlowContext,
                                  CSpectrum, CSpectrumFluxAxis_withError,
                                  CSpectrumSpectralAxis, ErrorCode,
                                  TLSFGaussianVarWidthArgs, WarningCode)

zlog = CLog.GetInstance()
zflag = CFlagWarning.GetInstance()


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
        observation_id: str,
        parameters: Parameters,
        calibration_library,
        source_id: str,
        filter_loader_class: AbstractFilterLoader = ParamJsonFilterLoader
    ):
        """Constructor method
        """
        self.observation_id = observation_id

        # Initial loaded data
        self.waves = Container[np.ndarray]()
        self.fluxes = Container[np.ndarray]()
        self.errors = Container[np.ndarray]()
        self.others: Dict[str, Container[np.ndarray]] = dict()

        self.lsf_type = None
        self.lsf_data = Container[np.ndarray]()
        self.photometric_data = []
        self._spectra = dict()
        self.w_frame = 'vacuum'
        self.parameters = parameters
        self.calibration_library = calibration_library
        self.source_id = str(source_id)
        self.editable_spectra = Container[pd.DataFrame]()
        self.filter_loader_class = filter_loader_class

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

    def load_photometry(self, resource, obs_id=""):
        """Append the photometric data in self.photometric_data , units are in erg.cm-2 by default

        :param resource: resource where the error can be found, no restriction for type (can be a path, a
            file handler, an hdf5 node,...)

        """
        raise NotImplementedError("Implement in derived class")

    def load_others(self, resource, obs_id=""):
        """Appends other data in self.others

        :param resource: resource where the error can be found, no restriction for type (can be a path, a
            file handler, an hdf5 node,...)
        """

    def load_all(self, resource):
        """
        Load all components of the spectrum. Reimplement this if resources are different

        :param resource: resource where wave, flux, error, lsf and photometry can be found
        """
        self.load_wave(resource)
        self.load_flux(resource)
        self.load_error(resource)
        self.load_lsf(resource)
        self.load_photometry(resource)
        self.load_others(resource)

    def get_spectrum(self, obs_id=""):
        """
        Get spectrum c++ object

        :return: CSpectrum object, fully loaded
        """
        self._check_spectrum_is_loaded()
        return self._spectra[obs_id]

    def get_wave(self, obs_id=""):
        """
        :raises: Spectrum not loaded
        :return: wavelength
        :rtype: np.array
        """
        self._check_spectrum_is_loaded()
        return self._spectra[obs_id].GetSpectralAxis().GetSamplesVector().to_numpy()

    def get_flux(self, obs_id=""):
        """
        :raises: Spectrum not loaded
        :return: wavelength
        :rtype: np.array
        """
        self._check_spectrum_is_loaded()
        return self._spectra[obs_id].GetFluxAxis().GetSamplesVector().to_numpy()

    def get_error(self, obs_id=""):
        """
        :raises: Spectrum not loaded
        :return: error
        :rtype: np.array
        """
        self._check_spectrum_is_loaded()
        return self._spectra[obs_id].GetErrorAxis().GetSamplesVector().to_numpy()

    def get_lsf(self, obs_id=""):
        """
        :raises: Spectrum not loaded
        :return: lsf
        :rtype: np.array
        """
        self._check_spectrum_is_loaded()
        return self.lsf_data.get(obs_id)

    def set_air(self):
        """
        Set w_frame to "air" (default is vacuum). Client should use this method if input spectra are in air
        frame
        """
        self.w_frame = "air"

    def get_loaded_lsf_observation_ids(self):
        return list(self.lsf_data.keys())

    def init(self):
        if not self._sizes_are_consistent():
            sizes: str = ', '.join([str(container.size()) for container in self._all_containers()])
            raise APIException(
                ErrorCode.INVALID_SPECTRUM,
                f"Number of error, wavelength, flux and other arrays should be the same: {sizes}"
            )

        if not self._obs_ids_are_consistent():
            raise APIException(
                ErrorCode.INVALID_SPECTRUM,
                "Names of error, wavelength, flux and other arrays should be the same"
            )

        self._merge_spectrum_in_dataframe()
        self._apply_filters(self._get_filters())
        self._check_wavelengths()
        self._add_cspectra()
        lsf_factory = CLSFFactory.GetInstance()
        lsf_args = self._lsf_args()
        lsf_type = self.parameters.get_lsf_type()
        if self.parameters.get_lsf_type() == "fromSpectrumData":
            lsf_type = self.lsf_type
        lsf = lsf_factory.Create(lsf_type, lsf_args)
        for obs_id in self.parameters.get_observation_ids():
            self._spectra[obs_id].SetLSF(lsf)
        self._add_photometric_data()

    def _get_filters(self):
        return self.filter_loader_class().get_filters(self.parameters)

    def _apply_filters(self, filters: FilterList) -> None:
        if filters is None:
            return
        for spectrum_key in self.editable_spectra.keys():
            self.editable_spectra.data[spectrum_key] = filters.apply(self.editable_spectra.data[spectrum_key])

    def _check_wavelengths(self):
        """Looks if lambda range specified in parameters is contained in spectrum range."""
        for obs_id in self.waves.keys():
            [params_lambda_min, params_lambda_max] = self.parameters.get_lambda_range(obs_id)
            obs_waves = np.array(self.waves.get(obs_id))
            if obs_waves.size == 0:
                raise APIException(
                    ErrorCode.INVALID_SPECTRUM,
                    "Filtered spectrum is empty"
                )
            spectrum_lambda_min = obs_waves[0]
            spectrum_lambda_max = obs_waves[-1]

            params_lambda_range_in_spectrum_lambdas = spectrum_lambda_min <= params_lambda_min \
                and spectrum_lambda_max >= params_lambda_max
            if not params_lambda_range_in_spectrum_lambdas:
                zflag.warning(
                    WarningCode.SPECTRUM_WAVELENGTH_TIGHTER_THAN_PARAM.value,
                    f"Parameters lambda range ([{params_lambda_min}, {params_lambda_max}])is not "
                    f"contained in spectrum wavelength ([{spectrum_lambda_min}, {spectrum_lambda_max}])"
                )

    def _add_photometric_data(self):
        if len(self.photometric_data) > 0 and len(self.photometric_data[0]) > 0:
            names = tuple(self.photometric_data[0].Name)
            flux = tuple([float(f) for f in self.photometric_data[0].Flux])
            fluxerr = tuple([float(f) for f in self.photometric_data[0].Error])
            for obs_id in self.parameters.get_observation_ids():
                self._spectra[obs_id].SetPhotData(CPhotometricData(names, flux, fluxerr))

    def _all_containers(self):
        return [container for container in self.others.values()] \
            + [self.waves, self.fluxes, self.errors]

    def _sizes_are_consistent(self) -> bool:
        expected_size = self.waves.size()
        return all(container.size() == expected_size for container in self._all_containers())

    def _obs_ids_are_consistent(self) -> bool:
        """Checks that all keys from waves conatiner are present in all other containers.
        Must be combined with size test to be sure that all keys are equal
        """
        expected_keys = self.waves.keys()
        for obs_id in expected_keys:
            for container in self._all_containers():
                if obs_id not in container.keys():
                    return False
        return True

    def _editable_waves(self, obs_id=""):
        return self.editable_spectra.get(obs_id)["waves"]

    def _editable_fluxes(self, obs_id=""):
        return self.editable_spectra.get(obs_id)["fluxes"]

    def _editable_errors(self, obs_id=""):
        return self.editable_spectra.get(obs_id)["errors"]

    def get_filtered_others(self, obs_id: str = "") -> pd.DataFrame:
        """
        Return a dataframe with the filtered non-mandatory columns of the spectrum.

        :param obs_id: name of the observation
        :return: dataframe with the data of the other columns of the spectrum
        """
        others_data = pd.DataFrame()
        for col in self.others:
            others_data[col] = self.editable_spectra.get(obs_id)[col]
        return others_data

    def _add_cspectra(self):
        airvacuum_method = self._corrected_airvacuum_method()
        multiobs_type = self.parameters.get_multiobs_method()
        if not multiobs_type:

            # Add check names if multiobs type is null
            if list(self.waves.keys()) != ['']:
                raise APIException(
                    ErrorCode.INVALID_NAME,
                    "Non multi obs observations cannot be named"
                )

            spectralaxis = CSpectrumSpectralAxis(self._editable_waves(), airvacuum_method)
            signal = CSpectrumFluxAxis_withError(self._editable_fluxes(), self._editable_errors())
            self._add_cspectrum(spectralaxis, signal)

        elif multiobs_type == "merge":
            wses = []
            for i in self.waves.keys():
                wse_ = pd.DataFrame()
                wse_["wave"] = self._editable_waves(i)
                wse_["flux"] = self._editable_fluxes(i)
                wse_["error"] = self._editable_errors(i)
                wses.append(wse_)
            wse = pd.concat(wses)
            wse.sort_values(["wave"], inplace=True)
            codes, uniques = pd.factorize(wse["wave"])
            epsilon = np.concatenate([i for i in map(np.arange, np.bincount(codes))]) * 1e-10
            wse["wave"] = wse["wave"] + epsilon

            if len(wse["wave"].unique()) != wse.index.size:
                raise APIException(ErrorCode.UNALLOWED_DUPLICATES, "Duplicates in wavelengths")

            if not (np.diff(wse["wave"]) > 0).all():
                raise APIException(ErrorCode.UNSORTED_ARRAY, "Wavelenghts are not sorted")

            spectralaxis = CSpectrumSpectralAxis(np.array(wse["wave"]), airvacuum_method)
            signal = CSpectrumFluxAxis_withError(np.array(wse["flux"]),
                                                 np.array(wse["error"]))
            self._add_cspectrum(spectralaxis, signal)

        elif multiobs_type == "full":
            for obs_id in self.waves.keys():
                spectralaxis = CSpectrumSpectralAxis(self._editable_waves(obs_id),
                                                     airvacuum_method)
                signal = CSpectrumFluxAxis_withError(self._editable_fluxes(obs_id),
                                                     self._editable_errors(obs_id))
                self._add_cspectrum(spectralaxis, signal, obs_id)

    def _add_cspectrum(self, spectralaxis, signal, obs_id=""):
        self._spectra[obs_id] = CSpectrum(spectralaxis, signal)
        self._spectra[obs_id].SetName(self.source_id)
        self._spectra[obs_id].setObsID(obs_id)

        ctx = CProcessFlowContext.GetInstance()
        ctx.addSpectrum(self._spectra[obs_id])

    def _corrected_airvacuum_method(self):
        airvacuum_method = self.parameters.get_airvacuum_method()

        if airvacuum_method == "default":
            airvacuum_method = "morton2000"

        if airvacuum_method == "" and self.w_frame == "air":
            airvacuum_method = "morton2000"

        elif airvacuum_method != "" and self.w_frame == "vacuum":
            zflag.warning(
                WarningCode.AIR_VACUUM_CONVERSION_IGNORED.value,
                f"Air vacuum method {airvacuum_method} ignored, spectrum already in vacuum"
            )
            airvacuum_method = ""

        return airvacuum_method

    def _lsf_args(self):
        ctx = CProcessFlowContext.GetInstance()
        parameter_lsf_type = self.parameters.get_lsf_type()

        if parameter_lsf_type == "fromSpectrumData":
            lsf_obs_ids = self.get_loaded_lsf_observation_ids()
            if not lsf_obs_ids:
                raise APIException(ErrorCode.LSF_NOT_LOADED,
                                   "No LSF loaded in reader, "
                                   "lsftype=fromSpectrumData "
                                   "parameter cannot be applied")
            obs_id = lsf_obs_ids[0]
            if len(lsf_obs_ids) > 1:
                zflag.warning(WarningCode.MULTI_OBS_ARBITRARY_LSF.value,
                              f"lsf of observation {obs_id} chosen, other lsf ignored")
            if self.lsf_type != "gaussianVariableWidth":
                self.parameters.set_lsf_param(LSFParameters[self.lsf_type],
                                              self.lsf_data.get(obs_id)["width"][0])
                parameter_store = ctx.LoadParameterStore(self.parameters.to_json())
                lsf_args = TLSFArgumentsCtor[self.lsf_type](parameter_store)
            else:
                lsf_args = TLSFGaussianVarWidthArgs(self.lsf_data.get(obs_id)["wave"],
                                                    self.lsf_data.get(obs_id)["width"])
        else:
            if parameter_lsf_type != "gaussianVariableWidth":
                parameter_store = ctx.LoadParameterStore(self.parameters.to_json())
                lsf_args = TLSFArgumentsCtor[parameter_lsf_type](parameter_store)
            else:
                lsf_args = TLSFGaussianVarWidthArgs(self.calibration_library.lsf["wave"],
                                                    self.calibration_library.lsf["width"])
        return lsf_args

    def _check_spectrum_is_loaded(self):
        if not self._spectra:
            raise APIException(
                ErrorCode.SPECTRUM_NOT_LOADED,
                "Spectrum not loaded, launch init first"
            )

    def _merge_spectrum_in_dataframe(self):
        for obs_id in self.waves.keys():
            # Creates dataframe with mandatory columns
            full_spectrum = pd.DataFrame({
                "waves": self.waves.get(obs_id),
                "fluxes": self.fluxes.get(obs_id),
                "errors": self.errors.get(obs_id),
            })
            for col_key in self.others:
                obs_others = self.others[col_key].get(obs_id)
                if obs_others is not None:
                    full_spectrum[col_key] = obs_others

            self.editable_spectra.append(full_spectrum, obs_id)

    def get_observation_ids(self):
        return self.waves.keys()
