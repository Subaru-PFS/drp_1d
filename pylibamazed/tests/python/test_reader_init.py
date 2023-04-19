from pylibamazed.Exception import APIException
import pytest

from pylibamazed.redshift import GlobalException
from .spectrum_reader_utils import TestSpectrumReaderUtils
import numpy as np
import pandas as pd


class TestReaderInit(TestSpectrumReaderUtils):
    def test_input_fluxes_size_error(self):
        fsr = self.initialize_fsr_with_data()
        # Adds a wave -> incoherent size
        fsr.load_wave([20,30], "2")
        with pytest.raises(APIException, match=r"INVALID_SPECTRUM"):
            fsr.init()
    
    def test_input_fluxes_key_error(self):
        fsr = self.initialize_fsr_with_data()
        # Adds an error with different obs_id -> incoherent keys
        fsr.load_wave([20,30], "2")
        fsr.load_flux([20,30], "2")
        fsr.load_error([20,30], "3")
        with pytest.raises(APIException, match=r"INVALID_SPECTRUM"):
            fsr.init()

    def test_multilsf_not_handled(self):
        fsr = self.initialize_fsr_with_data()
        fsr.lsf_data.append('exceeding lsf')
        with pytest.raises(APIException, match=r"MULTILSF_NOT_HANDLED"):
            fsr.init()

    def test_non_multi_obs_naming_restrictions(self):
        fsr = self.initialize_fsr_with_data(**{"obs_id": "name that shouldn't be here"})
        with pytest.raises(APIException, match=r"INVALID_NAME"):
            fsr.init()

    def test_wavelength_duplicates_error(self):
        fsr = self.initialize_fsr_with_data(**{"multiobsmethod": "merge"})
        self.full_load(fsr, **{"obs_id": "2"})
        fsr.waves.append(
            np.array([float(i) + 1e-10 for i in range(0, 10)]),
            obs_id="3",
        )
        fsr.fluxes.append(
            np.array([float(i) + 1e-10 for i in range(0, 10)]),
            obs_id="3",
        )
        fsr.errors.append(
            np.array([float(i) + 1e-10 for i in range(0, 10)]),
            obs_id="3",
        )
        with pytest.raises(APIException, match=r"UNALLOWED_DUPLICATES"):
            fsr.init()

    def test_wavelength_not_sorted_error(self):
        fsr = self.initialize_fsr_with_data(**{"multiobsmethod": "merge"})
        self.full_load(fsr, **{"obs_id": "2"})
        fsr.waves.append(
            np.array([float(i) + 1e-11 for i in range(0, 10)]),
            obs_id="3",
        )
        fsr.fluxes.append(
            np.array([float(i) + 1e-11 for i in range(0, 10)]),
            obs_id="3",
        )
        fsr.errors.append(
            np.array([float(i) + 1e-11 for i in range(0, 10)]),
            obs_id="3",
        )
        with pytest.raises(APIException, match=r"UNSORTED_ARRAY"):
            fsr.init()

    def test_no_multi_obs(self):
        fsr = self.initialize_fsr_with_data()
        fsr.init()

    def test_multi_obs_merge(self):
        fsr = self.initialize_fsr_with_data(**{"multiobsmethod": "merge"})
        fsr.load_wave([8,20], "2")
        fsr.load_flux([8,20], "2")
        fsr.load_error([8,20], "2")
        fsr.load_wave([7,12], "3")
        fsr.load_flux([7,12], "3")
        fsr.load_error([7,12], "3")
        fsr.init()

    def test_multi_obs_full(self):
        fsr = self.initialize_fsr_with_data(**{"multiobsmethod": "full"})
        fsr.init()

    def test_add_photometric_data(self):
        fsr = self.initialize_fsr_with_data()
        fsr.photometric_data = [pd.DataFrame([["a",2,3]], columns=['Name', 'Flux', 'Error'])]
        fsr.init()
    
    def test_airvacuum_method_settings(self):
        # No meaning - for test coverage only
        fsr = self.initialize_fsr_with_data(**{"airvacuum_method": "default"})
        fsr.init()

        fsr = self.initialize_fsr_with_data()
        fsr.set_air()
        with pytest.raises(GlobalException):
            fsr.init()

        fsr = self.initialize_fsr_with_data(**{"airvacuum_method": "default"})
        fsr.set_air()
        with pytest.raises(GlobalException):
            fsr.init()
    
    def test_lsf_args(self):
        # No meaning - for test coverage only
        fsr = self.initialize_fsr_with_data(**{"lsf_type": "GaussianNISPSIM2016"})
        with pytest.raises(BaseException):
            fsr.init()

        fsr = self.initialize_fsr_with_data(**{"lsf_type": "GaussianVariableWidth"})
        with pytest.raises(KeyError):
            fsr.init()
        
        fsr = self.initialize_fsr_with_data()
        fsr.lsf_type = "GaussianVariableWidth"
        with pytest.raises(TypeError):
            fsr.init()

    def test_sizes_are_consistent(self):
        # Without "others"
        # All sizes consistent
        fsr = self.initialize_fsr_with_data()
        assert fsr._sizes_are_consistent() == True

        # One different size
        fsr.load_wave([20,30], "2")
        assert fsr._sizes_are_consistent() == False
        
        # With others
        fsr = self.initialize_fsr_with_data()
        fsr.load_others("lambda")
        assert fsr._sizes_are_consistent() == True
        
        # One different size
        fsr.load_others("lambda", "2")
        assert fsr._sizes_are_consistent() == False

    def test_obs_ids_are_consistent(self):
        # Without "others"
        # All keys consistent, no multi obs
        fsr = self.initialize_fsr_with_data()
        assert fsr._obs_ids_are_consistent() == True

        # All key consistent, with multi obs
        TestSpectrumReaderUtils().full_load(fsr, **{"obs_id": "2"})
        assert fsr._obs_ids_are_consistent() == True
        
        # One key inconsistent, with multi obs
        fsr.load_wave([20,30], "3")
        fsr.load_flux([20,30], "3")
        fsr.load_error([20,30], "other")
        assert fsr._obs_ids_are_consistent() == False

        # With others
        # All keys consistent, no multi obs
        fsr = self.initialize_fsr_with_data()
        fsr.load_others("lambdas")
        assert fsr._obs_ids_are_consistent() == True
        
        # All key consistent, with multi obs
        TestSpectrumReaderUtils().full_load(fsr, **{"obs_id": "2"})
        fsr.load_others("lambdas", "2")
        assert fsr._obs_ids_are_consistent() == True
        
        # One key inconsistent, with multi obs
        TestSpectrumReaderUtils().full_load(fsr, **{"obs_id": "3"})
        fsr.load_others("lambdas", "other")
        assert fsr._obs_ids_are_consistent() == False
