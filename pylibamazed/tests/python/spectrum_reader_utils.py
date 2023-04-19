import random
import numpy as np
from pylibamazed.CalibrationLibrary import CalibrationLibrary
import tempfile

from pylibamazed.AbstractSpectrumReader import AbstractSpectrumReader, Container

class TestSpectrumReaderUtils:

    def make_parameters(self, **kwargs):
        parameters = dict()
        parameters["LSF"] = dict()
        parameters["LSF"]["LSFType"] = kwargs.get("lsf_type", "FROMSPECTRUMDATA")
        parameters["airvacuum_method"] = kwargs.get("airvacuum_method", "")
        parameters["objects"]=[]
        parameters["multiobsmethod"]=kwargs.get("multiobsmethod", "")
        return parameters

    def full_load(self, fsr, **kwargs):
        fsr.load_wave([0,10], kwargs.get('obs_id', ''))
        fsr.load_flux([0,10], kwargs.get('obs_id', ''))
        fsr.load_error([0,10], kwargs.get('obs_id', ''))

    def initialize_fsr_with_data(self, **kwargs):
        parameters = self.make_parameters(**kwargs)
        cl = CalibrationLibrary(parameters, tempfile.mkdtemp())
        fsr = FakeSpectrumReader("000", parameters, cl, "000","range")
        self.full_load(fsr, **kwargs)
        fsr.load_lsf(None)
        return fsr

class FakeSpectrumReader(AbstractSpectrumReader):

    def __init__(self, observation_id, parameters, calibration_library, source_id, lambda_type):
        AbstractSpectrumReader.__init__(self, observation_id,parameters,calibration_library,source_id)
        self.lambda_type = lambda_type

    def load_wave(self, l_range, obs_id=""):
        if self.lambda_type == "random":
            self.waves.append(
                np.array([float(i) + random.uniform(0, 0.0001) for i in range(10)]),
                obs_id=obs_id,
            )
        elif self.lambda_type == "range":
            self.waves.append(
                np.array([float(i) for i in range(l_range[0], l_range[1])]),
                obs_id=obs_id,
            )
        else:
            raise Exception("unknown lambda type")

    def load_flux(self, l_range, obs_id=""):
        if l_range is not None:
            self.fluxes.append(
                np.array([random.random() for i in range(l_range[0],l_range[1])]),
                obs_id=obs_id,
            )
        else:
            self.fluxes.append(
                np.array([random.random() for i in range(10)]),
                obs_id=obs_id,
            )

    def load_error(self, l_range, obs_id=""):
        if l_range is not None:
            self.errors.append(
                np.array([random.random() for i in range(l_range[0],l_range[1])]),
                obs_id=obs_id,
            )
        else:
            self.errors.append(
                np.array([random.random() for i in range(10)]),
                obs_id=obs_id,
            )
    
    def load_others(self, data_name: str, obs_id=""):
        if not self.others.get(data_name):
            self.others[data_name] = Container()
        self.others[data_name].append(
            np.array([random.random() for i in range(10)]),
            obs_id=obs_id,
        )

    def load_lsf(self, location):
        self.lsf_type = "GaussianConstantWidth"
        self.lsf_data.append(3.)

    def load_photometry(self, location):
        pass