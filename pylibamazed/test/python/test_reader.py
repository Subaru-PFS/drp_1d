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

from pylibamazed.AbstractSpectrumReader import AbstractSpectrumReader
from pylibamazed.CalibrationLibrary import CalibrationLibrary

import random
import numpy as np
import tempfile


class FakeSpectrumReader(AbstractSpectrumReader):

    def __init__(self, observation_id, parameters, calibration_library, source_id, lambda_type):
        AbstractSpectrumReader.__init__(self, observation_id,parameters,calibration_library,source_id)
        self.lambda_type = lambda_type

    def load_wave(self, l_range):
        if self.lambda_type == "random":
            self.waves.append(np.array([float(i) + random.uniform(0, 0.0001) for i in range(10)]))
        elif self.lambda_type == "range":
            self.waves.append(np.array([float(i) for i in range(l_range[0], l_range[1])]))
        else:
            raise Exception("unknown lambda type")

    def load_flux(self, l_range):
        if l_range is not None:
            self.fluxes.append(np.array([random.random() for i in range(l_range[0],l_range[1])]))
        else:
            self.fluxes.append(np.array([random.random() for i in range(10)]))

    def load_error(self, l_range):
        if l_range is not None:
            self.errors.append(np.array([random.random() for i in range(l_range[0],l_range[1])]))
        else:
            self.errors.append(np.array([random.random() for i in range(10)]))

    def load_lsf(self, location):
        self.lsf_type = "GaussianConstantWidth"
        self.lsf_data.append(3.)

    def load_photometry(self, location):
        pass
        #self.photometric_data.append


def test_reader():
    parameters = dict()
    parameters["LSF"] = dict()
    parameters["LSF"]["LSFType"] = "FROMSPECTRUMDATA"
    parameters["airvacuum_method"] = ""
    parameters["objects"]=[]
    cl = CalibrationLibrary(parameters, "/tmp")
    fsr = FakeSpectrumReader("000", parameters, cl, "000","random")
    fsr.load_all(None)
    fsr.init()
    s = fsr.get_spectrum()


def test_multi_obs():
    parameters = dict()
    parameters["LSF"] = dict()
    parameters["LSF"]["LSFType"] = "FROMSPECTRUMDATA"
    parameters["airvacuum_method"] = ""
    parameters["objects"]=[]
    cl = CalibrationLibrary(parameters, tempfile.mkdtemp())
    fsr = FakeSpectrumReader("000", parameters, cl, "000","range")
    fsr.load_wave([0,10])
    fsr.load_flux([0,10])
    fsr.load_error([0,10])
    fsr.load_wave([8,20])
    fsr.load_flux([8,20])
    fsr.load_error([8,20])
    fsr.load_wave([7,12])
    fsr.load_flux([7,12])
    fsr.load_error([7,12])
    fsr.load_lsf(None)
    fsr.init()
    s = fsr.get_spectrum()
