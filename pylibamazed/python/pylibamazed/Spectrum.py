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

import numpy as np
import pandas as pd
from pylibamazed.Exception import APIException
from pylibamazed.lsf import get_lsf_args_from_parameters
from pylibamazed.Parameters import Parameters
from pylibamazed.redshift import (
    CFlagWarning,
    CLog,
    CLSFFactory,
    CPhotometricData,
    CProcessFlowContext,
    CSpectrum,
    CSpectrumFluxAxis_withError,
    CSpectrumSpectralAxis,
    ErrorCode,
    WarningCode,
)
from pylibamazed.Filter import FilterList
from pylibamazed.FilterLoader import AbstractFilterLoader, ParamJsonFilterLoader

zlog = CLog.GetInstance()
zflag = CFlagWarning.GetInstance()


class Spectrum:
    """
    class for spectrum interface
    """

    def __init__(
        self,
        source_id: str,
        parameters: Parameters,
        spectra_dataframe,
        lsf,
        photometric_data,
        wave_frame="vacuum",
        filter_loader_class: AbstractFilterLoader = ParamJsonFilterLoader,
    ):
        self.source_id = str(source_id)
        self.parameters = parameters

        self.filter_loader_class = filter_loader_class

        # Data members
        self._dataframe = spectra_dataframe
        self._lsf = lsf
        self._photometric_data = photometric_data
        self.w_frame = wave_frame
        self.masks = dict()

    @property
    def _is_multiobs_merge(self) -> bool:
        return self.observation_number > 1 and self.parameters.get_multiobs_method() == "merge"

    def _merge_if_not_available(self):
        if "wave_merged" in self._dataframe.columns:
            return
        self._merge_spectra()

    def _is_obs_id_merge(self, obs_id) -> bool:
        if obs_id != "":
            return False
        if self._is_multiobs_merge:
            return True
        return False

    @property
    def observation_ids(self):
        return self._dataframe.index.levels[0]

    @property
    def observation_number(self):
        return len(self.observation_ids)

    def get_dataframe(self, obs_id=None, filtered_only=True) -> pd.DataFrame:
        """
        Get spectra pandas dataframe
            * all spectra if obs_id=None or merge
            * if filtered_only return only filtered samples,
              otherwise return all input samples

        :return: pandas dataframe
        """
        if obs_id is None:
            df = self._dataframe
        elif self._is_obs_id_merge(obs_id):
            self._merge_if_not_available()
            df = self._dataframe
        else:
            df = self._dataframe.loc[obs_id]
            if filtered_only and "amazed_mask" in df:
                return df[df["amazed_mask"]]
            else:
                return self._dataframe.loc[obs_id]
        if filtered_only and "amazed_mask" in df:
            return df[df["amazed_mask"]]
        else:
            return df

    def get_index(self, obs_id="", filtered_only=True) -> pd.Index:
        df = self.get_dataframe(obs_id, filtered_only)
        return df.index

    def get_samples_number(self, obs_id="", filtered_only=True) -> int:
        return len(self.get_index(obs_id, filtered_only))

    def get_wave(self, obs_id="", filtered_only=True, vacuum=True) -> pd.Series:
        """
        :return: wavelength
        : if obs_id is None return the unmerged wavelength of all observations
        : if obs_id is "" return the merged wavelength of all observations/
        :type: pandas.series
        """
        wave_column = "wave"

        if self.w_frame == "air":
            if vacuum:
                self._convert_air_to_vacuum()
            else:
                wave_column = "wave_air"
        if self._is_obs_id_merge(obs_id):
            wave_column = "wave_merged"
        spectrum = self.get_dataframe(obs_id, filtered_only)
        return spectrum[wave_column]

    def get_flux(self, obs_id="", filtered_only=True) -> pd.Series:
        """
        :return: wavelength
        :rtype: pandas.series
        """
        spectrum = self.get_dataframe(obs_id, filtered_only)

        return spectrum["flux"]

    def get_error(self, obs_id="", filtered_only=True) -> pd.Series:
        """
        :return: error
        :rtype: pandas.series
        """
        spectrum = self.get_dataframe(obs_id, filtered_only)
        return spectrum["error"]

    def get_others(self, obs_id: str = "", filtered_only=True) -> pd.DataFrame:
        """
        Return a dataframe with the filtered non-mandatory columns of the spectrum.

        :param obs_id: name of the observation
        :return: dataframe with the data of the other columns of the spectrum
        """
        spectrum = self.get_dataframe(obs_id, filtered_only)
        col = spectrum.columns != "amazed_mask"
        col &= spectrum.columns != "wave_air"
        col &= spectrum.columns != "wave_merged"
        col &= spectrum.columns != "wave"
        col &= spectrum.columns != "flux"
        col &= spectrum.columns != "error"

        return spectrum.loc[:, col]

    def get_lsf(self, obs_id=""):
        """
        :return: lsf
        :rtype: np.array
        """
        return self._lsf

    def get_photometric_data(self):
        return self._photometric_data

    def push_in_context(self):
        ctx = CProcessFlowContext.GetInstance()
        cpp_spectra = self._make_cspectra()
        for cpp_spectrum in cpp_spectra.values():
            ctx.addSpectrum(cpp_spectrum)

    def _make_clsf(self):
        lsf_factory = CLSFFactory.GetInstance()
        lsf_args = self._lsf_args()
        lsf_type = self._lsf["lsfType"]
        return lsf_factory.Create(lsf_type, lsf_args)

    def _make_photometric_data(self):
        if len(self._photometric_data) > 0 and len(self._photometric_data[0]) > 0:
            names = tuple(self._photometric_data[0].Name)
            flux = tuple([float(f) for f in self._photometric_data[0].Flux])
            fluxerr = tuple([float(f) for f in self._photometric_data[0].Error])
            cpp_phot = CPhotometricData(names, flux, fluxerr)
            return cpp_phot

    def _get_obs_ids(self) -> list:
        multiobs_type = self.parameters.get_multiobs_method()
        if not multiobs_type:
            obs_ids = [""]

        elif self._is_multiobs_merge:
            self._merge_if_not_available()
            obs_ids = [""]

        elif multiobs_type == "full":
            obs_ids = self.observation_ids

        return obs_ids

    def _make_cspectra(self) -> dict:
        cpp_spectra = dict()

        cpp_lsf = self._make_clsf()
        cpp_phot = self._make_photometric_data()

        obs_ids = self._get_obs_ids()

        for obs_id in obs_ids:
            spectralaxis = CSpectrumSpectralAxis(self.get_wave(obs_id))
            signal = CSpectrumFluxAxis_withError(self.get_flux(obs_id), self.get_error(obs_id))
            cpp_spectra[obs_id] = self._make_cspectrum(spectralaxis, signal, cpp_lsf, cpp_phot, obs_id)

        return cpp_spectra

    def _make_cspectrum(self, spectralaxis, signal, cpp_lsf, cpp_phot, obs_id="") -> CSpectrum:
        cpp_spectrum = CSpectrum(spectralaxis, signal)
        cpp_spectrum.SetName(self.source_id)
        cpp_spectrum.setObsID(obs_id)
        cpp_spectrum.SetLSF(cpp_lsf)
        cpp_spectrum.SetPhotData(cpp_phot)
        return cpp_spectrum

    @property
    def airvacuum_method(self):
        try:
            return self._airvacuum_method
        except AttributeError:
            return self._corrected_airvacuum_method()

    def _corrected_airvacuum_method(self):
        self._airvacuum_method = self.parameters.get_airvacuum_method()

        if self._airvacuum_method == "default":
            self._airvacuum_method = "morton2000"

        if self._airvacuum_method == "" and self.w_frame == "air":
            self._airvacuum_method = "morton2000"

        elif self._airvacuum_method != "" and self.w_frame == "vacuum":
            zflag.warning(
                WarningCode.AIR_VACUUM_CONVERSION_IGNORED,
                f"Air vacuum method {self._airvacuum_method} ignored, spectrum already in vacuum",
            )
            self._airvacuum_method = ""

        return self._airvacuum_method

    def _convert_air_to_vacuum(self):
        if (self.w_frame == "vacuum") or ("wave_air" in self._dataframe.columns):
            return
        try:
            self.w_frame = "vacuum"  # temporary set to vacuum to get wave column
            wave_air = self.get_wave(obs_id=None)
        finally:
            self.w_frame = "air"
        spectralAxis = CSpectrumSpectralAxis(wave_air, self.airvacuum_method)
        wave_vacuum = spectralAxis.GetSamplesVector().to_numpy()
        self._dataframe.rename(columns={"wave": "wave_air"}, inplace=True)
        self._dataframe["wave"] = wave_vacuum

    def _lsf_args(self):
        return get_lsf_args_from_parameters(self._lsf)

    def _merge_spectra(self):
        self._dataframe.sort_values(["wave"], inplace=True)
        #  add k*epsilon at wavelength of each n-duplicates, k in [0,n[
        codes, uniques = pd.factorize(self._dataframe["wave"])
        epsilon = np.concatenate([i for i in map(np.arange, np.bincount(codes))]) * 1e-10
        self._dataframe["wave_merged"] = self._dataframe["wave"] + epsilon

        if len(self._dataframe["wave_merged"].unique()) != len(self._dataframe.index):
            raise APIException(ErrorCode.UNALLOWED_DUPLICATES, "Duplicates in multi-obs merged wavelengths")

        if not (np.diff(self._dataframe["wave_merged"]) > 0).all():
            raise APIException(ErrorCode.UNSORTED_ARRAY, "Wavelenghts are not sorted")

    def _check_wavelengths(self):
        """Looks if lambda range specified in parameters is contained in spectrum range."""
        for obs_id in self.parameters.get_observation_ids():
            [params_lambda_min, params_lambda_max] = self.parameters.get_lambda_range(obs_id)
            obs_waves = np.array(self.get_wave(obs_id))
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

    def _get_filters(self, obs_id: str):
        return self.filter_loader_class().get_filters(self.parameters, obs_id=obs_id)

    def _apply_filters(self, filters: FilterList, obs_id: str) -> None:
        if filters is None:
            return
        self.masks[obs_id] = filters.apply(self._dataframe.loc[obs_id])

    def init(self):
        """
        Does three things :
         - Check if airvaccum conversion specified in parameters is legitimate
         - Apply filtering specified in parameters
         - Check if spectrum(s) wavelength(s) is/are correct:
            not empty, contains lambda range(s) specified in parameters
        """

        self._corrected_airvacuum_method()

        for obs_id in self._get_obs_ids():
            self._apply_filters(self._get_filters(obs_id), obs_id)
        if self.masks:
            self._dataframe["amazed_mask"] = [
                self.masks[obs_id][i]
                for obs_id, i in zip(
                    self._dataframe.index.get_level_values("obs_id"),
                    self._dataframe.index.get_level_values("index"),
                )
            ]
        self._check_wavelengths()
