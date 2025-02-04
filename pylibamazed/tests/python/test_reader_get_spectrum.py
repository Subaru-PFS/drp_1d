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

import numpy as np
import pandas as pd
import pytest
from pylibamazed.Exception import APIException
from pylibamazed.redshift import WarningCode
from tests.python.spectrum_reader_utils import TestSpectrumReaderUtils
from tests.python.utils import WarningUtils


class TestReaderGetSpectrum(TestSpectrumReaderUtils):
    def test_input_fluxes_size_error(self):
        fsr = self.initialize_fsr_with_data()
        # Adds a wave -> incoherent size
        fsr.load_wave([20, 30], "2")
        with pytest.raises(APIException, match=r"INVALID_SPECTRUM"):
            fsr.get_spectrum()

    def test_input_fluxes_key_error(self):
        fsr = self.initialize_fsr_with_data()
        # Adds an error with different obs_id -> incoherent keys
        fsr.load_wave([20, 30], "2")
        fsr.load_flux([20, 30], "2")
        fsr.load_error([20, 30], "3")
        with pytest.raises(APIException, match=r"INVALID_SPECTRUM"):
            fsr.get_spectrum()

    def test_non_multi_obs_naming_restrictions(self):
        fsr = self.initialize_fsr_with_data(**{"obs_ids": ["name that shouldn't be here"]})
        with pytest.raises(APIException, match=r"INVALID_NAME"):
            fsr.get_spectrum()

    def test_wavelength_duplicates_error(self):
        fsr = self.initialize_fsr_with_data()
        self.full_load(fsr)
        fsr.waves.append(np.array([0, 1, 1, 3, 4, 5, 6, 7, 8, 9], dtype=float))
        with pytest.raises(APIException, match=r"UNALLOWED_DUPLICATES"):
            fsr.get_spectrum()

    def test_no_multi_obs(self):
        fsr = self.initialize_fsr_with_data()
        fsr.get_spectrum()

    def test_multi_obs_merge(self):
        fsr = self.initialize_fsr_with_data(
            **{
                "multiObsMethod": "merge",
                "obs_ids": ["1", "2", "3"],
                "wave_range": {"1": [0, 10], "2": [8, 20], "3": [7, 12]},
            }
        )
        fsr.get_spectrum()

    def test_multi_obs_full(self):
        fsr = self.initialize_fsr_with_data(**{"multiObsMethod": "full"})
        fsr.get_spectrum()

    def test_add_photometric_data(self):
        fsr = self.initialize_fsr_with_data()
        fsr.photometric_data = [pd.DataFrame([["a", 2, 3]], columns=["Name", "Flux", "Error"])]
        fsr.get_spectrum()

    def test_lsf_args(self):
        # No meaning - for test coverage only
        fsr = self.initialize_fsr_with_data(**{"lsfType": "gaussianNISPSIM2016"})
        fsr.get_spectrum()

        fsr = self.initialize_fsr_with_data(**{"lsfType": "gaussianVariableWidth"})
        with pytest.raises(APIException):
            fsr.get_spectrum()

    def test_sizes_are_consistent(self):
        # Without "others"
        # All sizes consistent
        fsr = self.initialize_fsr_with_data()
        assert fsr._sizes_are_consistent() is True

        # One different size
        fsr.load_wave([20, 30], "2")
        assert fsr._sizes_are_consistent() is False

        # With others
        fsr = self.initialize_fsr_with_data()
        fsr.load_others("lambda")
        assert fsr._sizes_are_consistent() is True

        # One different size
        fsr.load_others("lambda", "2")
        assert fsr._sizes_are_consistent() is False

    def test_obs_ids_are_consistent(self):
        # Without "others"
        # All keys consistent, no multi obs
        fsr = self.initialize_fsr_with_data(**{"obs_ids": ["1"]})
        assert fsr._obs_ids_are_consistent() is True

        # All key consistent, with multi obs
        TestSpectrumReaderUtils().full_load(fsr, **{"obs_ids": ["2"]})
        assert fsr._obs_ids_are_consistent() is True

        # One key inconsistent, with multi obs
        fsr.load_wave([20, 30], "3")
        fsr.load_flux([20, 30], "3")
        fsr.load_error([20, 30], "other")
        assert fsr._obs_ids_are_consistent() is False

        # With others
        # All keys consistent, no multi obs
        fsr = self.initialize_fsr_with_data()
        fsr.load_others("lambdas")
        assert fsr._obs_ids_are_consistent() is True

        # All key consistent, with multi obs
        TestSpectrumReaderUtils().full_load(fsr, **{"obs_ids": ["2"]})
        fsr.load_others("lambdas", "2")
        assert fsr._obs_ids_are_consistent() is True

        # One key inconsistent, with multi obs
        TestSpectrumReaderUtils().full_load(fsr, **{"obs_ids": ["3"]})
        fsr.load_others("lambdas", "other")
        assert fsr._obs_ids_are_consistent() is False

    def test_input_output_size_coherence_with_filtering(self):
        fsr = self.initialize_fsr_with_data(
            **{
                "obs_ids": ["1", "2"],
                "parameters_lambdaRange": {"1": [0, 1], "2": [3, 4]},
                "multiObsMethod": "full",
                "filters": {
                    "1": [{"key": "wave", "instruction": "<=", "value": 3}],
                    "2": [{"key": "wave", "instruction": "<=", "value": 3}],
                },
            }
        )

        # Default data: 0 -> 9 => 4 items correspond to the filter
        spectra = fsr.get_spectrum()
        spectra.init()

        # Check size is correct for spectrum first obs
        spectrum_id1 = spectra.get_dataframe("1")

        assert len(spectrum_id1.index) == 4

        # Check size is correct for second obs
        spectrum_id2 = spectra.get_dataframe("2")
        assert len(spectrum_id2.index) == 4

    def test_nb_samples_min(self):
        fsr = self.initialize_fsr_with_data(
            **{"filters": [{"key": "wave", "instruction": "<=", "value": 3}], "nbSamplesMin": 5}
        )

        spectra = fsr.get_spectrum()
        with pytest.raises(APIException, match=r"INVALID_SPECTRUM"):
            spectra.init()

    def test_nb_samples_min_full_multiobs(self):
        fsr = self.initialize_fsr_with_data(
            **{
                "obs_ids": ["1", "2"],
                "parameters_lambdaRange": {"1": [0, 1], "2": [3, 4]},
                "multiObsMethod": "full",
                "filters": {
                    "1": [{"key": "wave", "instruction": "<=", "value": 3}],
                    "2": [{"key": "wave", "instruction": "<=", "value": 3}],
                },
                "nbSamplesMin": 5,
            }
        )

        spectra = fsr.get_spectrum()
        with pytest.raises(APIException, match=r"INVALID_SPECTRUM"):
            spectra.init()

    class TestSpectrumLambdaRange(TestSpectrumReaderUtils):
        def _init_fsr(self, spectrum_wave_range, parameters_lambda_range):
            fsr = self.initialize_fsr_with_data(
                **{
                    "spectrum_wave_range": spectrum_wave_range,
                    "parameters_lambdaRange": parameters_lambda_range,
                }
            )
            fsr.get_spectrum().init()

        def test_warning_if_parameters_lambdarange_is_outside_spectrum(self, zflag):
            fsr = self.initialize_fsr_with_data(
                **{"spectrum_wave_range": [1, 3], "parameters_lambdaRange": [5, 40]}
            )
            fsr.get_spectrum().init()

            assert WarningUtils.has_any_warning()

        def test_warning_if_parameters_lambdarange_is_outside_spectrum_multiobs(self, zflag):
            fsr = self.initialize_fsr_with_data(
                **{
                    "obs_ids": ["1", "2"],
                    "multiObsMethod": "full",
                    "spectrum_wave_range": {"1": [10, 40], "2": [1, 4]},
                    "parameters_lambdaRange": {"1": [20, 30], "2": [0, 3]},
                }
            )
            fsr.get_spectrum().init()

            assert WarningUtils.has_warning(WarningCode.SPECTRUM_WAVELENGTH_TIGHTER_THAN_PARAM)

        def test_warning_if_parameters_lambdarange_begins_lower_than_spectrum(self, zflag):
            self._init_fsr(spectrum_wave_range=[10, 30], parameters_lambda_range=[5, 25])
            assert WarningUtils.has_any_warning()

        def test_warning_if_parameters_lambdarange_ends_higher_than_spectrum(self, zflag):
            self._init_fsr(spectrum_wave_range=[10, 30], parameters_lambda_range=[11, 31])
            assert WarningUtils.has_any_warning()

        def test_OK_if_parameters_lambdarange_boundaries_are_contained_in_spectrum(self, zflag):
            self._init_fsr(spectrum_wave_range=[10, 30], parameters_lambda_range=[15, 20])
            assert not WarningUtils.has_any_warning()

        def test_OK_if_parameters_lambdarange_boundaries_are_contained_in_spectrum_multiobs(self, zflag):
            fsr = self.initialize_fsr_with_data(
                **{
                    "obs_ids": ["1", "2"],
                    "multiObsMethod": "full",
                    "lsf_type": "gaussianConstantWidth",
                    "width": 0.01,
                    "spectrum_wave_range": {"1": [10, 40], "2": [1, 4]},
                    "parameters_lambdaRange": {"1": [20, 30], "2": [2, 3]},
                }
            )
            fsr.get_spectrum()

            assert not WarningUtils.has_warning(WarningCode.SPECTRUM_WAVELENGTH_TIGHTER_THAN_PARAM)

        def test_OK_if_parameters_lambdarange_boundaries_are_equal_to_spectrum(self, zflag):
            self._init_fsr(
                spectrum_wave_range=[10, 30], parameters_lambda_range=[10, 29]  # range => last item is 29
            )
            assert not WarningUtils.has_any_warning()
