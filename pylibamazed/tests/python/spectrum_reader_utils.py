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
import random
import tempfile

import numpy as np
from pylibamazed.AbstractSpectrumReader import (AbstractSpectrumReader,
                                                Container)
from pylibamazed.CalibrationLibrary import CalibrationLibrary
from pylibamazed.Parameters import Parameters
from tests.python.fake_parameters_checker import FakeParametersChecker


class TestSpectrumReaderUtils:

    def make_parameters_dict(self, **kwargs):
        params_dict = dict()
        params_dict["LSF"] = dict()
        params_dict["LSF"]["LSFType"] = kwargs.get("lsf_type", "FROMSPECTRUMDATA")
        params_dict["airvacuum_method"] = kwargs.get("airvacuum_method", "")
        params_dict["objects"] = []
        params_dict["multiobsmethod"] = kwargs.get("multiobsmethod", "")
        if params_dict["multiobsmethod"] == "full":
            params_dict["lambdarange"] = kwargs.get("parameters_lambdarange", {"": [0, 10]})
        else:
            params_dict["lambdarange"] = kwargs.get("parameters_lambdarange", [0, 10])
        if kwargs.get("filters"):
            params_dict["filters"] = kwargs.get("filters")
        if kwargs.get("width"):
            params_dict["LSF"]["width"] = kwargs.get("width")
        return params_dict

    def full_load(self, fsr, **kwargs):
        obs_id = kwargs.get('obs_id', '')
        if kwargs.get("multiobsmethod") == "full":
            wave_ranges = kwargs.get("spectrum_wave_range", {obs_id: [0, 10]})
            wave_range = wave_ranges[obs_id]
        else:
            wave_range = kwargs.get("spectrum_wave_range", [0, 10])
        fsr.load_wave(wave_range, kwargs.get('obs_id', ''))
        fsr.load_flux(wave_range, kwargs.get('obs_id', ''))
        fsr.load_error(wave_range, kwargs.get('obs_id', ''))
        fsr.load_lsf(None, kwargs.get('obs_id', ''))

    def initialize_fsr_with_data(self, **kwargs):
        params_dict = self.make_parameters_dict(**kwargs)
        params = Parameters(params_dict, FakeParametersChecker)
        cl = CalibrationLibrary(params, tempfile.mkdtemp())
        fsr = FakeSpectrumReader("000", params, cl, "000", "range")
        self.full_load(fsr, **kwargs)
        return fsr


class FakeSpectrumReader(AbstractSpectrumReader):

    def __init__(self, observation_id, parameters: Parameters, calibration_library, source_id, lambda_type):
        AbstractSpectrumReader.__init__(self, observation_id, parameters, calibration_library, source_id)
        self.lambda_type = lambda_type
        super().__init__(observation_id, parameters, calibration_library, source_id)

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
                np.array([random.random() for i in range(l_range[0], l_range[1])]),
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
                np.array([random.random() for i in range(l_range[0], l_range[1])]),
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

    def load_lsf(self, location, obs_id=""):
        self.lsf_type = "GaussianConstantWidth"
        lsf = np.zeros(np.dtype([("width", '<f8')]))
        lsf["width"][0] = 3.0
        self.lsf_data.append(lsf, obs_id)

    def load_photometry(self, location):
        pass
