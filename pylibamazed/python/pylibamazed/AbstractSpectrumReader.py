import numpy as np
import os,json
from pylibamazed.redshift import (CSpectrumSpectralAxis,
                                  CSpectrumFluxAxis_withError,
                                  CSpectrum,
                                  PC_Get_AxisSampleList,
                                  CProcessFlowContext,
                                  TLSFGaussianVarWidthArgs,
                                  CLSFFactory,
                                  CPhotometricData)
from .lsf import LSFParameters, TLSFArgumentsCtor


class AbstractSpectrumReader:

    def __init__(self, observation_id, parameters, calibration_library, source_id):
        '''
        :param observation_id:
        :param parameters: dict containing the parameters
        :param calibration_dir:
        :param source_id: Astronomical source ID
        '''
        self.observation_id = observation_id
        self.waves = []
        self.fluxes = []
        self.errors = []
        self.lsf_type = None
        self.lsf_data = []
        self.photometric_data = []
        self._spectra = []
        self.w_frame = 'vacuum'
        self.parameters = parameters
        self.calibration_library = calibration_library
        self.source_id = source_id

    def load_wave(self, arg):
        """
        Load the spectral axis in self.wave , units are in Angstrom by default
        :param location: location of the resource where the wave can be found
        """
        raise NotImplementedError("Implement in derived class")

    def load_flux(self, arg):
        """
        Load the spectral axis in self.flux , units are in erg.cm-2 by default
        :param location: location of the resource where the wave can be found
        """
        raise NotImplementedError("Implement in derived class")

    def load_error(self, arg):
        raise NotImplementedError("Implement in derived class")

    def load_lsf(self, arg):
        raise NotImplementedError("Implement in derived class")

    def load_photometry(self, location):
        raise NotImplementedError("Implement in derived class")

    def load_all(self, location):
        '''
        Load all components of the spectrum. Reimplement this if locations are different
        :param location:
        @return:
        '''
        self.load_wave(location)
        self.load_flux(location)
        self.load_error(location)
        self.load_lsf(location)
        self.load_photometry(location)
        self.init()

    def get_spectrum(self):
        if self._spectra[0] is not None:
            return self._spectra[0]
        else:
            raise Exception("Spectrum not loaded")

    def get_wave(self):
        if self._spectra[0] is not None:

            return PC_Get_AxisSampleList(self._spectra[0].GetSpectralAxis().GetSamplesVector())
        else:
            raise Exception("Spectrum not loaded")

    def get_flux(self):
        if self._spectra[0] is not None:
            return PC_Get_AxisSampleList( self._spectra[0].GetFluxAxis().GetSamplesVector())
        else:
            raise Exception("Spectrum not loaded")

    def get_error(self):
        if self._spectra[0] is not None:
            return PC_Get_AxisSampleList(self._spectra[0].GetFluxAxis().GetSamplesVector())
        else:
            raise Exception("Spectrum not loaded")

    def get_lsf(self):
        return self.lsf_data[0]

    def set_air(self):
        self.w_frame = "air"

    def init(self):
        if len(self.lsf_data) > 1:
            raise Exception("Multiple LSF not handled")
        airvacuum_method = self.parameters.get("airvacuum_method", "")
        if airvacuum_method == "" and self.w_frame == "air":
            airvacuum_method = "Morton2000"
        elif airvacuum_method != "" and self.w_frame == "vacuum":
            print("Air vaccum method " + airvacuum_method + " ignored, spectrum already in vacuum")
            airvacuum_method = ""

        spectralaxis = CSpectrumSpectralAxis(self.waves[0], airvacuum_method)
        signal = CSpectrumFluxAxis_withError(self.fluxes[0], self.errors[0])
        self._spectra.append(CSpectrum(spectralaxis, signal))
        self._spectra[0].SetName(self.source_id)

        ctx = CProcessFlowContext()
        parameter_lsf_type = self.parameters["LSF"]["LSFType"]
        if parameter_lsf_type == "FROMSPECTRUMDATA":
            self.parameters["LSF"]["LSFType"] = self.lsf_type
            if self.lsf_type != "GaussianVariableWidth":
                self.parameters["LSF"][LSFParameters[self.lsf_type]] = self.lsf_data
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
        if len(self.photometric_data[0]) > 0:
            names = tuple(self.photometry.Name)
            flux = tuple([float(f) for f in self.photometry.Flux])
            fluxerr = tuple([float(f) for f in self.photometry.Error])
            self._spectra[0].SetPhotData(CPhotometricData(names, flux, fluxerr))

        del ctx


