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
import os,json

import pandas as pd
from pylibamazed.redshift import (CSpectrumSpectralAxis,
                                  CSpectrumFluxAxis_withError,
                                  CSpectrum,
                                  PC_Get_AxisSampleList,
                                  CProcessFlowContext,
                                  TLSFGaussianVarWidthArgs,
                                  CLSFFactory,
                                  CPhotometricData,
                                  CLog,
                                  CFlagWarning,
                                  ErrorCode
                                 )
from pylibamazed.lsf import LSFParameters, TLSFArgumentsCtor

zlog = CLog.GetInstance()
zflag = CFlagWarning.GetInstance()
from pylibamazed.Exception import APIException
class AbstractSpectrumReader:
    """
    Base class for spectrum reader, it handles at least wavelengths, flux and error (variance). It also handles
    Light Spread Function (LSF) whether its spectrum dependent or general. It can also handle photometric data
    if present.

    :param observation_id: Observation ID of the spectrum, there can be multiple observation id for the same source_id
    :type observation_id: str
    :param parameters: Parameters of the amazed run. It should at least have information about LSF and lambda range
    :type parameters: dict
    :param calibration_library: Calibration library object, containing all calibration data necessary,
        according to parameters
    :type calibration_library: class:`CalibrationLibrary`
    :param source_id: Astronomical source ID
    :type source_id: str

    """

    def __init__(self, observation_id, parameters, calibration_library, source_id):
        """Constructor method
        """
        self.observation_id = observation_id
        self.waves = []
        self.fluxes = []
        self.errors = []
        self.lsf_type = None
        self.lsf_data = []
        self.photometric_data = []
        self._spectra = []
        self.w_frame = 'vacuum'
        self.parameters = parameters.copy()
        self.calibration_library = calibration_library
        self.source_id = source_id

    def load_wave(self, resource):
        """Append the spectral axis in self.wave , units are in Angstrom by default

        :param resource: resource where the wave can be found, no restriction for type (can be a path,
            a file handler, an hdf5 node,...)

        """
        raise NotImplementedError("Implement in derived class")

    def load_flux(self, resource):
        """Append the flux in self.flux , units are in erg.cm-2 by default

        :param resource: resource where the flux can be found, no restriction for type (can be a path, a file handler,
            an hdf5 node,...)

        """
        raise NotImplementedError("Implement in derived class")

    def load_error(self, resource):
        """Append the variance in self.error , units are in erg.cm-2 by default

        :param resource: resource where the error can be found, no restriction for type (can be a path, a file handler,
            an hdf5 node,...)

        """
        raise NotImplementedError("Implement in derived class")

    def load_lsf(self, resource):
        """Append the spectral axis in self.flux , units are in erg.cm-2 by default

         :param resource: resource where the error can be found, no restriction for type (can be a path, a file handler,
            an hdf5 node,...)

        """
        raise NotImplementedError("Implement in derived class")

    def load_photometry(self, resource):
        """Append the photometric data in self.photometric_data , units are in erg.cm-2 by default

        :param resource: resource where the error can be found, no restriction for type (can be a path, a file handler,
            an hdf5 node,...)

        """
        raise NotImplementedError("Implement in derived class")

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


    def get_spectrum(self):
        """
        Get spectrum c++ object

        :return: CSpectrum object, fully loaded
        """
        if self._spectra[0] is not None:
            return self._spectra[0]
        else:
            raise  APIException(ErrorCode.SPECTRUM_NOT_LOADED,"Spectrum not loaded")
    def get_wave(self):
        """

        :raises: Spectrum not loaded
        :return: wavelength
        :rtype: np.array
        """
        if self._spectra[0] is not None:
            return PC_Get_AxisSampleList(self._spectra[0].GetSpectralAxis().GetSamplesVector())
        else:
            raise  APIException(ErrorCode.SPECTRUM_NOT_LOADED,"Spectrum not loaded")

    def get_flux(self):
        """

        :raises: Spectrum not loaded
        :return: wavelength
        :rtype: np.array
        """
        if self._spectra[0] is not None:
            return PC_Get_AxisSampleList( self._spectra[0].GetFluxAxis().GetSamplesVector())
        else:
            raise  APIException(ErrorCode.SPECTRUM_NOT_LOADED,"Spectrum not loaded")

    def get_error(self):
        """

        :raises: Spectrum not loaded
        :return: error
        :rtype: np.array
        """
        if self._spectra[0] is not None:
            return PC_Get_AxisSampleList(self._spectra[0].GetErrorAxis().GetSamplesVector())
        else:
            raise  APIException(ErrorCode.SPECTRUM_NOT_LOADED,"Spectrum not loaded")

    def get_lsf(self):
        """
        :raises: Spectrum not loaded
        :return: lsf
        :rtype: np.array
        """
        return self.lsf_data[0]

    def set_air(self):
        """
        Set w_frame to "air" (default is vacuum). Client should use this method if input spectra are in air frame
        """
        self.w_frame = "air"

    def init(self):
        if len(self.waves) != len(self.fluxes) or len(self.waves) != len(self.errors):
            raise  APIException(ErrorCode.INVALID_SPECTRUM,"Number of error, wavelength and flux arrays should be the same:{0} {1} {2}".format(str(len(self.waves)), str(len(self.errors)), str(len(self.fluxes))))
        if len(self.lsf_data) > 1:
            raise  APIException(ErrorCode.MULTILSF_NOT_HANDELED,"Multiple LSF not handled")
        airvacuum_method = self.parameters.get("airvacuum_method", "")
        if airvacuum_method == "default":
            airvacuum_method = "Morton2000"
        if airvacuum_method == "" and self.w_frame == "air":
            airvacuum_method = "Morton2000"
        elif airvacuum_method != "" and self.w_frame == "vacuum":
            zflag.warning(zflag.AIR_VACCUM_CONVERSION_IGNORED, "Air vaccum method " + airvacuum_method + " ignored, spectrum already in vacuum")
            airvacuum_method = ""

        if len(self.waves) == 1:
            spectralaxis = CSpectrumSpectralAxis(self.waves[0], airvacuum_method)
            signal = CSpectrumFluxAxis_withError(self.fluxes[0], self.errors[0])
        else:
            wse = pd.DataFrame()
            wse["wave"] = self.waves[0]
            wse["flux"] = self.fluxes[0]
            wse["error"] = self.errors[0]

            for i in range(1, len(self.waves)):
                wse_ = pd.DataFrame()
                wse_["wave"] = self.waves[i]
                wse_["flux"] = self.fluxes[i]
                wse_["error"] = self.errors[i]
                wse = wse.append(wse_)
            wse.sort_values(["wave"], inplace=True)
            codes, uniques = pd.factorize(wse["wave"])
            epsilon = np.concatenate([i for i in map(np.arange, np.bincount(codes))]) * 1e-10
            wse["wave"] = wse["wave"] + epsilon
            if len(wse["wave"].unique()) != wse.index.size:
                raise  APIException(ErrorCode.UNALLOWED_DUPLICATES,"Duplicates in wavelengths")
            if not (np.diff(wse["wave"]) > 0).all():
                raise  APIException(ErrorCode.UNSORTED_ARRAY,"Wavelenghts are not sorted")
            spectralaxis = CSpectrumSpectralAxis(np.array(wse["wave"]), airvacuum_method)
            signal = CSpectrumFluxAxis_withError(np.array(wse["flux"]),
                                                 np.array(wse["error"]))

        self._spectra.append(CSpectrum(spectralaxis, signal))
        self._spectra[0].SetName(self.source_id)

        ctx = CProcessFlowContext.GetInstance()
        ctx.setSpectrum(self._spectra[0])
        parameter_lsf_type = self.parameters["LSF"]["LSFType"]
        if parameter_lsf_type == "FROMSPECTRUMDATA":
            self.parameters["LSF"]["LSFType"] = self.lsf_type
            if self.lsf_type != "GaussianVariableWidth":
                self.parameters["LSF"][LSFParameters[self.lsf_type]] = self.lsf_data[0]
                parameter_store = ctx.LoadParameterStore(json.dumps(self.parameters))
                lsf_args = TLSFArgumentsCtor[self.lsf_type](parameter_store)
            else:
                lsf_args = TLSFGaussianVarWidthArgs(self.lsf_data[0]["wave"], self.lsf_data[0]["width"])
        else:
            if parameter_lsf_type != "GaussianVariableWidth":
                parameter_store = ctx.LoadParameterStore(json.dumps(self.parameters))
                lsf_args = TLSFArgumentsCtor[parameter_lsf_type](parameter_store)
            else:
                lsf_args = TLSFGaussianVarWidthArgs(self.calibration_library.lsf["wave"],
                                                    self.calibration_library.lsf["width"])

        lsf_factory = CLSFFactory.GetInstance()
        lsf = lsf_factory.Create(self.parameters["LSF"]["LSFType"], lsf_args)
        self._spectra[0].SetLSF(lsf)
        if len(self.photometric_data) > 0 and len(self.photometric_data[0]) > 0:
            names = tuple(self.photometric_data[0].Name)
            flux = tuple([float(f) for f in self.photometric_data[0].Flux])
            fluxerr = tuple([float(f) for f in self.photometric_data[0].Error])
            self._spectra[0].SetPhotData(CPhotometricData(names, flux, fluxerr))

