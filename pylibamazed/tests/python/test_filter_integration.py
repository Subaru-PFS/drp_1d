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

from pylibamazed.ASCIISpectrumReader import ASCIISpectrumReader
from pylibamazed.Parameters import Parameters
from pylibamazed.ProcessFlow import ProcessFlow
from tests.python.test_ITlike import get_observation, get_parameters, get_spectra, make_config


class TestFilterIntegration:
    def test_filter_params(self):
        # Checks that parameters containing columns names & filters are correctly accessed and used

        # Creates a "real" configuration
        config = make_config(**{"config_filename": "config_filters.json"})
        param = Parameters(get_parameters(config["parameters_file"]), make_checks=False)
        process_flow = ProcessFlow(config, param)  # vars returns the dict version of config
        observation = get_observation(config["input_file"])

        # Read and load spectra using spectra reader
        spectra = get_spectra(config, observation)
        reader = ASCIISpectrumReader(
            parameters=param,
            calibration_library=process_flow.calibration_library,
            source_id=observation.ProcessingID[0],
        )

        spectrum = reader.load_and_get_spectrum(spectra)
        spectrum.init()

        # Checks that the number of waves kept has decreased (6 to 3) with filtering
        assert len(spectrum.get_wave()) == 3
