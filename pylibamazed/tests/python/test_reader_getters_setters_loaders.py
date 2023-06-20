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

from pylibamazed.CalibrationLibrary import CalibrationLibrary
from pylibamazed.Parameters import Parameters
from tests.python.spectrum_reader_utils import TestSpectrumReaderUtils
from pylibamazed.Exception import APIException
import pytest
from pylibamazed.AbstractSpectrumReader import AbstractSpectrumReader
import tempfile

class TestReaderGettersSettersLoaders(TestSpectrumReaderUtils):

    def test_get_wave(self):
        fsr = self.initialize_fsr_with_data()
        fsr.init()
        fsr.get_wave()

    def test_get_flux(self):
        fsr = self.initialize_fsr_with_data()
        fsr.init()
        fsr.get_flux()

    def test_get_error(self):
        fsr = self.initialize_fsr_with_data()
        fsr.init()
        fsr.get_error()
    
    def test_get_lsf(self):
        fsr = self.initialize_fsr_with_data()
        fsr.init()
        fsr.get_lsf()
    
    def test_get_spectrum(self):
        fsr = self.initialize_fsr_with_data()
        fsr.init()
        fsr.get_spectrum()
    
    def test_set_air(self):
        fsr = self.initialize_fsr_with_data()
        fsr.set_air()
    
    def test_load_all(self):
        fsr = self.initialize_fsr_with_data()
        fsr.load_others = lambda *args : None
        fsr.load_all([1,10])
    
    def test_load_others(self):
        # Make test from the abstract spectrum reader
        parameters = Parameters(self.make_parameters_dict())
        cl = CalibrationLibrary(parameters, tempfile.mkdtemp())
        fsr = AbstractSpectrumReader("000", parameters, cl, "000")
        fsr.load_others('data name')
    
class SpectrumReader(AbstractSpectrumReader):
    def __init__(self):
        self._spectra = []
        pass

class TestSpectrumReaderGetErrors:

    spectrumReader = SpectrumReader()

    def test_get_spectrum_error(self):
        with pytest.raises(APIException, match=r"SPECTRUM_NOT_LOADED"):
            self.spectrumReader.get_spectrum()
    
    def test_get_wave_error(self):
        with pytest.raises(APIException, match=r"SPECTRUM_NOT_LOADED"):
            self.spectrumReader.get_wave()

    def test_get_flux_error(self):
        with pytest.raises(APIException, match=r"SPECTRUM_NOT_LOADED"):
            self.spectrumReader.get_flux()

    def test_get_error_error(self):
        with pytest.raises(APIException, match=r"SPECTRUM_NOT_LOADED"):
            self.spectrumReader.get_error()
    
    def test_get_lsf_error(self):
        with pytest.raises(APIException, match=r"SPECTRUM_NOT_LOADED"):
            self.spectrumReader.get_lsf()