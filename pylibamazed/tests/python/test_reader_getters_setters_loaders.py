from .spectrum_reader_utils import TestSpectrumReaderUtils
from pylibamazed.Exception import APIException
import pytest
from pylibamazed.AbstractSpectrumReader import AbstractSpectrumReader

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
        fsr.load_all([1,10])

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