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
from pylibamazed.Container import Container
from pylibamazed.redshift import CFlagWarning, CLog

zlog = CLog.GetInstance()
zflag = CFlagWarning.GetInstance()


class ASCIISpectrumReader(AbstractSpectrumReader):
    """
    Child class for spectrum reader, it handles at least wavelengths, flux and error (variance). It also
    handles Light Spread Function (LSF) whether its spectrum dependent or general. It can also handle
    photometric data if present.

    Reads a CSpectrum
    :param observation_id: Observation ID of the spectrum, there can be multiple observation id for the same
        source_id
    :type observation_id: str
    :param parameters: Parameters of the amazed run. It should at least have information about LSF and lambda
        range
    :type parameters: dict
    :param calibration_library: Calibration library object, containing all calibration data necessary,
        according to parameters
    :type calibration_library: class:`CalibrationLibrary`
    :param source_id: Astronomical source ID
    :type source_id: str

    """

    def load_wave(self, spectrum, obs_id=""):
        self.waves.append(spectrum.wave)

    def load_flux(self, spectrum, obs_id=""):
        self.fluxes.append(spectrum.flux)

    def load_error(self, spectrum, obs_id=""):
        self.errors.append(spectrum.err)

    def load_lsf(self, spectrum, obs_id=""):
        pass
        # self.lsf_type = "no_lsf"

    def add_photometry(self, phot):
        self.photometric_data.append(phot)

    def load_others(self, spectrum, obs_id: str = ""):
        additional_cols = self.parameters.get_additional_cols()
        if additional_cols is None:
            return

        for col_name in additional_cols:
            col_data = spectrum[col_name]
            if self.others.get(col_name) is None:
                self.others[col_name] = Container(**{obs_id: col_data})
            else:
                self.others[col_name].append(col_data, obs_id)
