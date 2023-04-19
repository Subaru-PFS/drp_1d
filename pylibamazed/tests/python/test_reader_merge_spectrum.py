from pylibamazed.AbstractSpectrumReader import Container
from pylibamazed.CalibrationLibrary import CalibrationLibrary

from .spectrum_reader_utils import FakeSpectrumReader, TestSpectrumReaderUtils
import numpy as np
import pandas as pd
import tempfile


class TestMergeSpectrumInDataframe(TestSpectrumReaderUtils):
    # Some default data
    waves = [1,2,3]
    fluxes = [11,12,13]
    errors = [111,112,113]
    lsf_data = [1111, 1112, 1113]
    photometric_data = [11111, 11112, 11113]
    lambdas = [4, 5, 6]

    # Base of the expected spectrum
    expected_full_scpectrum = pd.DataFrame({
        "waves": waves,
        "fluxes": fluxes,
        "errors": errors,
        "lsf_data": lsf_data,
        "photometric_data": photometric_data,
    })

    def _initialize_fsr(self):
        # Initializes empty fake spectrum reader
        parameters = self.make_parameters()
        cl = CalibrationLibrary(parameters, tempfile.mkdtemp())
        fsr = FakeSpectrumReader("000", parameters, cl, "000","range")

        # Adds some data
        fsr.waves.append(np.array(self.waves))
        fsr.fluxes.append(np.array(self.fluxes))
        fsr.errors.append(np.array(self.errors))
        fsr.lsf_data = [self.lsf_data]
        fsr.photometric_data = [self.photometric_data]

        return fsr

    def test_without_others(self):
        # TODO in other commit change list for container
        fsr = self._initialize_fsr()
        fsr._merge_spectrum_in_dataframe()

        assert fsr.full_spectra.size() == 1
        assert self.expected_full_scpectrum.equals(fsr.full_spectra.get())

    def test_with_others(self):
        fsr = self._initialize_fsr()

        # Adds a "other" colum
        fsr.others["lambdas"] = Container(**{"": np.array(self.lambdas) })

        # Updates expected spectrum
        self.expected_full_scpectrum["lambdas"] = self.lambdas
        
        fsr._merge_spectrum_in_dataframe()
        assert fsr.full_spectra.size() == 1
        assert self.expected_full_scpectrum.equals(fsr.full_spectra.get())

    def test_with_multi_obs(self):

        waves2 = [21,22,23]
        fluxes2 = [211,212,213]
        errors2  = [2111,2112,2113]
        lambdas2 = [24, 25, 26]

        fsr = self._initialize_fsr()

        # Adds multi obs + other column
        fsr.waves.append(np.array(waves2), "2")
        fsr.fluxes.append(np.array(fluxes2), "2")
        fsr.errors.append(np.array(errors2), "2")

        fsr.others["lambdas"] = Container(**{
            "": np.array(self.lambdas),
            "2": np.array(lambdas2)
        })
        
        # Expected spectrum
        fsr._merge_spectrum_in_dataframe()
        self.expected_full_scpectrum["lambdas"] = self.lambdas
        expected_full_scpectrum2 = pd.DataFrame({
            "waves": waves2,
            "fluxes": fluxes2,
            "errors": errors2,
            "lsf_data": self.lsf_data,
            "photometric_data": self.photometric_data,
            "lambdas": lambdas2,
        })

        assert fsr.full_spectra.size() == 2
        assert self.expected_full_scpectrum.equals(fsr.full_spectra.get())
        assert expected_full_scpectrum2.equals(fsr.full_spectra.get("2"))
