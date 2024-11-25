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
from pylibamazed.redshift import AmzException, WarningCode
from tests.python.spectrum_reader_utils import TestSpectrumReaderUtils
from tests.python.utils import WarningUtils


class TestMergeSpectrum(TestSpectrumReaderUtils):
    def test_wavelength_duplicates_error(self):
        fsr = self.initialize_fsr_with_data(**{"multiObsMethod": "merge", "obs_id": "1"})
        self.full_load(fsr, **{"obs_id": "2"})
        fsr.waves.append(
            np.arange(10, dtype=float) + 1e-10,
            obs_id="3",
        )
        fsr.fluxes.append(
            np.arange(10, dtype=float) + 1e-10,
            obs_id="3",
        )
        fsr.errors.append(
            np.arange(10, dtype=float) + 1e-10,
            obs_id="3",
        )
        fsr.load_lsf(None, "3")
        spectrum = fsr.get_spectrum()
        with pytest.raises(APIException, match=r"UNALLOWED_DUPLICATES"):
            spectrum._merge_spectra()

    def test_no_multi_obs(self):
        fsr = self.initialize_fsr_with_data()
        spectrum = fsr.get_spectrum()
        spectrum.get_wave()

    def test_multi_obs_merge(self):
        fsr = self.initialize_fsr_with_data(**{"multiObsMethod": "merge", "obs_id": "1"})
        fsr.load_wave([8, 20], "2")
        fsr.load_flux([8, 20], "2")
        fsr.load_error([8, 20], "2")
        fsr.load_lsf(None, "2")
        fsr.load_wave([7, 12], "3")
        fsr.load_flux([7, 12], "3")
        fsr.load_error([7, 12], "3")
        fsr.load_lsf(None, "3")
        fsr.get_spectrum()._merge_spectra()

    def test_multi_obs_full(self):
        fsr = self.initialize_fsr_with_data(**{"multiObsMethod": "full"})
        fsr.get_spectrum()._merge_spectra()

    def test_airvacuum_method_settings(self):
        # No meaning - for test coverage only
        fsr = self.initialize_fsr_with_data(**{"airVacuumMethod": "default"})
        fsr.get_spectrum().get_wave()

        fsr = self.initialize_fsr_with_data()
        fsr.w_frame = "air"
        with pytest.raises(AmzException):
            fsr.get_spectrum().get_wave()

        fsr = self.initialize_fsr_with_data(**{"airVacuumMethod": "default"})
        fsr.w_frame = "air"
        with pytest.raises(AmzException):
            fsr.get_spectrum().get_wave()

    def test_lsf_args(self):
        # No meaning - for test coverage only
        fsr = self.initialize_fsr_with_data(**{"lsfType": "gaussianNISPSIM2016"})
        spectrum = fsr.get_spectrum()
        spectrum.get_lsf()
        spectrum._lsf_args()

        fsr = self.initialize_fsr_with_data(**{"lsfType": "gaussianVariableWidth"})
        with pytest.raises(
            KeyError
        ):  # error raised since missing parameter lsf.gaussianVariableWidthFileName
            fsr.get_spectrum().get_lsf()

    # TODO to make a real test we should add an lsf compatible with spectrum
    # fsr = self.initialize_fsr_with_data()
    # fsr.lsf_type = "gaussianVariableWidth"
    # with pytest.raises(TypeError):
    #     fsr.init()


##  TODO add make_cspectra testing
