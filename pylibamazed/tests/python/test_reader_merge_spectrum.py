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
import pandas as pd
from pylibamazed.AbstractSpectrumReader import Container
from pylibamazed.CalibrationLibrary import CalibrationLibrary
from pylibamazed.Parameters import Parameters
from tests.python.fake_parameters_checker import FakeParametersChecker
from tests.python.spectrum_reader_utils import FakeSpectrumReader, TestSpectrumReaderUtils


class TestMergeSpectrumInDataframe(TestSpectrumReaderUtils):
    # Some default data
    waves = [1, 2, 3]
    fluxes = [11, 12, 13]
    errors = [111, 112, 113]
    lambdas = [4, 5, 6]

    # Base of the expected spectrum
    expected_spectrum_dataframe = pd.DataFrame(
        {
            "wave": waves,
            "flux": fluxes,
            "error": errors,
        }
    )

    def _initialize_fsr(self):
        # Initializes empty fake spectrum reader
        parameters = Parameters(self.make_parameters_dict(), Checker=FakeParametersChecker)
        cl = CalibrationLibrary(parameters, tempfile.mkdtemp())
        fsr = FakeSpectrumReader("000", parameters, cl, "000", "range")

        # Adds some data
        fsr.waves.append(np.array(self.waves))
        fsr.fluxes.append(np.array(self.fluxes))
        fsr.errors.append(np.array(self.errors))

        return fsr

    def test_without_others(self):
        fsr = self._initialize_fsr()
        fsr._merge_spectrum_in_dataframe()

        assert len(fsr.spectra_dataframe.index.levels[0]) == 1
        assert self.expected_spectrum_dataframe.equals(fsr.spectra_dataframe.loc[""])

    def test_with_others(self):
        fsr = self._initialize_fsr()

        # Adds a "other" colum
        fsr.others["lambda"] = Container(**{"": np.array(self.lambdas)})

        # Updates expected spectrum
        self.expected_spectrum_dataframe["lambda"] = self.lambdas

        fsr._merge_spectrum_in_dataframe()
        assert len(fsr.spectra_dataframe.index.levels[0]) == 1
        assert self.expected_spectrum_dataframe.equals(fsr.spectra_dataframe.loc[""])

    def test_with_multi_obs(self):
        waves2 = [21, 22, 23]
        fluxes2 = [211, 212, 213]
        errors2 = [2111, 2112, 2113]
        lambdas2 = [24, 25, 26]

        fsr = self._initialize_fsr()

        # Adds multi obs + other column
        fsr.waves.append(np.array(waves2), "2")
        fsr.fluxes.append(np.array(fluxes2), "2")
        fsr.errors.append(np.array(errors2), "2")

        fsr.others["lambda"] = Container(**{"": np.array(self.lambdas), "2": np.array(lambdas2)})

        # Expected spectrum
        fsr._merge_spectrum_in_dataframe()
        self.expected_spectrum_dataframe["lambda"] = self.lambdas
        expected_spectrum_dataframe2 = pd.DataFrame(
            {
                "wave": waves2,
                "flux": fluxes2,
                "error": errors2,
                "lambda": lambdas2,
            }
        )

        assert len(fsr.spectra_dataframe.index.levels[0]) == 2
        assert self.expected_spectrum_dataframe.equals(fsr.spectra_dataframe.loc[""])
        assert expected_spectrum_dataframe2.equals(fsr.spectra_dataframe.loc["2"])
