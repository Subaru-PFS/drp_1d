from pylibamazed.AbstractSpectrumReader import AbstractSpectrumReader
import pytest

class SpectrumReader(AbstractSpectrumReader):
    def __init__(self):
        pass

class TestSpectrumReaderImplementationErrors:

    spectrumReader = SpectrumReader()

    def test_load_wave_not_implemented_error(self):
        with pytest.raises(NotImplementedError):
            self.spectrumReader.load_wave('', '')

    def test_load_flux_not_implemented_error(self):
        with pytest.raises(NotImplementedError):
            self.spectrumReader.load_flux('', '')

    def test_load_error_not_implemented_error(self):
        with pytest.raises(NotImplementedError):
            self.spectrumReader.load_error('', '')

    def test_load_lsf_not_implemented_error(self):
        with pytest.raises(NotImplementedError):
            self.spectrumReader.load_lsf('', '')
    
    def test_load_photometry_not_implemented_error(self):
        with pytest.raises(NotImplementedError):
            self.spectrumReader.load_photometry('', '')