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

import json
import os
import tempfile

import h5py
import pandas as pd
from pylibamazed.ASCIISpectrumReader import ASCIISpectrumReader
from pylibamazed.Context import Context
from pylibamazed.H5Writer import H5Writer
from pylibamazed.Parameters import Parameters
from tests.python.config import test_dir
from tests.python.fake_parameters_checker import FakeParametersChecker


def read_photometry_fromfile(fname):
    df = pd.DataFrame(
        [
            ["riz", 2e-29, 3e-31],
            ["Y", 3e-29, 4e-31],
            ["J", 4e-29, 5e-31],
            ["H", 5e-29, 3e-31],
        ],
        columns=["Name", "Flux", "Error"],
    )
    return df


def accessOutputData(output):
    return output.get_attributes(
        [
            "classification",
            "error",
            "ContextWarningFlags",
            "WarningFlags",
            "test.WarningFlags",
        ],
        None,
    )


def make_config(**kwargs):
    config_path = os.path.join(
        test_dir,
        kwargs.get("config_filename", "config.json"),
    )

    with open(os.path.expanduser(config_path), "r") as f:
        config = json.load(f)

    config["calibration_dir"] = os.path.join(test_dir, config["calibration_dir"])

    return config


def get_parameters(parameters_file_path):
    parameters_json_path = os.path.join(
        test_dir,
        parameters_file_path,
    )
    with open(parameters_json_path) as f:
        param = json.load(f)

    return param


def get_observation(input_file_path):
    input_spectra_path = os.path.join(test_dir, input_file_path)
    observation = pd.read_table(
        input_spectra_path, delimiter=" ", names=["ProcessingID"]
    )
    return observation


def get_spectra(config, observation):
    s_filename = (
        config["spectrum_prefix"]
        + str(observation.ProcessingID[0])  # only one spectra
        + config["spectrum_suffix"]
    )
    s_filename = os.path.join(test_dir, config["spectrum_dir"], s_filename)

    # read and load spectra using spectra reader
    spectra = pd.read_table(s_filename, delimiter="\t")
    return spectra


def add_photometry_to_reader(config, observation, reader):
    phot_fname = (
        os.path.join(test_dir, config["spectrum_dir"])
        + "/"
        + str(observation.ProcessingID[0])
        + "_phot.txt"
    )
    if os.path.exists(phot_fname):
        phot = read_photometry_fromfile(phot_fname)
        reader.load_photometry(phot)


def save_output(output, config, observation):
    # save in temporary file
    tf = tempfile.TemporaryFile()
    output_file = h5py.File(tf, "w")
    writer = H5Writer(output)
    writer.excluded_datasets = config["excluded_datasets"]
    writer.write_hdf5(output_file, str(observation.ProcessingID[0]))
    output_file.close()


def test_ITLikeTest():
    config = make_config()
    param = Parameters(get_parameters(config["parameters_file"]), Checker=FakeParametersChecker)
    context = Context(config, param)  # vars returns the dict version of config
    observation = get_observation(config["input_file"])

    # read and load spectra using spectra reader
    spectra = get_spectra(config, observation)

    reader = ASCIISpectrumReader(
        observation_id=observation.ProcessingID[0],
        parameters=param,
        calibration_library=context,
        source_id=observation.ProcessingID[0],
    )

    reader.load_all(spectra)
    add_photometry_to_reader(config, observation, reader)

    output = context.run(reader)  # passing spectra reader to launch amazed

    # check results (no errors)
    for spectrum_model, stage in (("", "init"),
                                  ("galaxy", "redshiftSolver"),
                                  ("galaxy", "linemeas_catalog_load"),
                                  ("galaxy", "lineMeasSolver"),
                                  ("galaxy", "subClassifSolver"),
                                  ("galaxy", "reliabilitySolver"),
                                  ("", "classification"),
                                  ("", "load_result_store")):
        if output.has_error(spectrum_model, stage):
            print("object_type", spectrum_model, "stage", stage, output.get_error(spectrum_model, stage))
        assert output.has_error(spectrum_model, stage) is False

    # add calls to output
    # accessOutputData(output)

    save_output(output, config, observation)
