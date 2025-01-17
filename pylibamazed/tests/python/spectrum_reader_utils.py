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
import tempfile

import numpy as np
from pylibamazed.AbstractSpectrumReader import AbstractSpectrumReader, Container
from pylibamazed.CalibrationLibrary import CalibrationLibrary
from pylibamazed.Parameters import Parameters


class TestSpectrumReaderUtils:
    def make_parameters_dict(self, **kwargs):
        params_dict = dict({"version": 2})
        params_dict["lsf"] = dict()
        params_dict["lsf"]["lsfType"] = kwargs.get("lsfType", "fromSpectrumData")
        params_dict["airVacuumMethod"] = kwargs.get("airVacuumMethod", "")
        params_dict["spectrumModels"] = []
        params_dict["multiObsMethod"] = kwargs.get("multiObsMethod", "")
        if params_dict["multiObsMethod"] == "full":
            params_dict["lambdaRange"] = kwargs.get("parameters_lambdaRange", {"": [0, 10]})
        else:
            params_dict["lambdaRange"] = kwargs.get("parameters_lambdaRange", [0, 10])
        if kwargs.get("filters"):
            params_dict["filters"] = kwargs.get("filters")
        if kwargs.get("width"):
            params_dict["lsf"]["width"] = kwargs.get("width")
        params_dict["nbSamplesMin"] = kwargs.get("nbSamplesMin", 0)
        return params_dict

    def full_load(self, fsr, **kwargs):
        multiObsMethod = kwargs.get("multiObsMethod")
        for obs_id in kwargs.get("obs_ids", [""]):
            if multiObsMethod == "full" or multiObsMethod == "merge":
                wave_ranges = kwargs.get("spectrum_wave_range", {obs_id: [0, 10]})
                wave_range = wave_ranges[obs_id]
            else:
                wave_range = kwargs.get("spectrum_wave_range", [0, 10])
            fsr.load_wave(wave_range, obs_id)
            fsr.load_flux(wave_range, obs_id)
            fsr.load_error(wave_range, obs_id)
            fsr.load_lsf(None, obs_id)

    def initialize_empty_fsr(self, **kwargs):
        params_dict = self.make_parameters_dict(**kwargs)
        params = Parameters(params_dict, make_checks=False)
        cl = CalibrationLibrary(params, tempfile.mkdtemp())
        fsr = FakeSpectrumReader(params, cl, "000", "range")
        return fsr

    def initialize_fsr_with_data(self, **kwargs):
        fsr = self.initialize_empty_fsr(**kwargs)
        self.full_load(fsr, **kwargs)
        return fsr


class FakeSpectrumReader(AbstractSpectrumReader):
    def __init__(self, parameters: Parameters, calibration_library, source_id, lambda_type):
        self.lambda_type = lambda_type
        super().__init__(parameters, calibration_library, source_id)

    def load_wave(self, l_range, obs_id=""):
        if self.lambda_type == "random":
            self.waves.append(
                np.arange(10, dtype=float) + np.random.uniform(0, 0.0001, 10),
                obs_id=obs_id,
            )
        elif self.lambda_type == "range":
            self.waves.append(
                np.arange(l_range[0], l_range[1], dtype=float),
                obs_id=obs_id,
            )
        else:
            raise Exception("unknown lambda type")

    def load_flux(self, l_range, obs_id=""):
        if l_range is not None:
            self.fluxes.append(
                np.random.uniform(0, 1, l_range[1] - l_range[0]),
                obs_id=obs_id,
            )
        else:
            self.fluxes.append(
                np.random.random(10),
                obs_id=obs_id,
            )

    def load_error(self, l_range, obs_id=""):
        if l_range is not None:
            self.errors.append(
                np.random.random(l_range[1] - l_range[0]),
                obs_id=obs_id,
            )
        else:
            self.errors.append(
                np.random.random(10),
                obs_id=obs_id,
            )

    def load_others(self, data_name: str, obs_id=""):
        if not self.others.get(data_name):
            self.others[data_name] = Container()
        self.others[data_name].append(
            np.random.random(10),
            obs_id=obs_id,
        )

    def load_lsf(self, location, obs_id=""):
        self.lsf_type = "gaussianConstantWidth"
        lsf = np.ndarray((1,), dtype=np.dtype([("width", "<f8")]))
        lsf["width"][0] = 3.0
        self.lsf_data.append(lsf, obs_id)
