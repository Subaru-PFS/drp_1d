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

import pytest
from pylibamazed.AbstractSpectrumReader import AbstractSpectrumReader


class SpectrumReader(AbstractSpectrumReader):
    def __init__(self):
        pass

    def load_wave(self, resource, obs_id=""):
        return super().load_wave(resource, obs_id)

    def load_flux(self, resource, obs_id=""):
        return super().load_flux(resource, obs_id)

    def load_error(self, resource, obs_id=""):
        return super().load_error(resource, obs_id)

    def load_lsf(self, resource, obs_id=""):
        return super().load_lsf(resource, obs_id)


class TestSpectrumReaderImplementationErrors:
    spectrumReader = SpectrumReader()

    def test_load_wave_not_implemented_error(self):
        with pytest.raises(NotImplementedError):
            self.spectrumReader.load_wave(None)

    def test_load_flux_not_implemented_error(self):
        with pytest.raises(NotImplementedError):
            self.spectrumReader.load_flux(None)

    def test_load_error_not_implemented_error(self):
        with pytest.raises(NotImplementedError):
            self.spectrumReader.load_error(None)

    def test_load_others_not_implemented_error(self):
        self.spectrumReader.load_others(None)

    def test_load_lsf_not_implemented_error(self):
        with pytest.raises(NotImplementedError):
            self.spectrumReader.load_lsf(None)

    def test_load_photometry_not_implemented_error(self):
        self.spectrumReader.load_photometry(None)

    def test_set_air_or_vaccum_not_implemented_error(self):
        self.spectrumReader.set_air_or_vaccum(None)
