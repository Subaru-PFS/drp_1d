
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

import os

from pylibamazed.CalibrationLibrary import CalibrationLibrary
from pylibamazed.Parameters import Parameters
from tests.python.fake_parameters_checker import FakeParametersChecker

module_root_dir = os.path.split(__file__)[0]
calibration_dir = os.path.join(module_root_dir, "..", "..", "auxdir", "pylibamazed", "test", "calibration")


def make_parameters() -> Parameters:
    parameters_dict = dict()
    parameters_dict["objects"] = ["galaxy"]
    parameters_dict["ebmv"] = dict()
    parameters_dict["ebmv"]["count"] = 3
    parameters_dict["ebmv"]["start"] = 0.0
    parameters_dict["ebmv"]["step"] = 0.5
    parameters_dict["LSF"] = {"LSFType": "GaussianVariableWidth",
                                         "GaussianVariablewidthFileName": "LSF/EuclidNISPVSSPSF201707.fits"}
    parameters_dict["galaxy"] = dict()
    parameters_dict["galaxy"]["template_dir"] = "templates/BC03_sdss_tremonti21"
    parameters_dict["galaxy"]["LineModelSolve"] = dict()
    parameters_dict["galaxy"]["LineModelSolve"]["linemodel"] = dict()
    parameters_dict["galaxy"]["LineModelSolve"]["linemodel"]["linecatalog"] = \
        "linecatalogs/linecatalogamazedvacuum_H0.tsv"
    parameters_dict["galaxy"]["LineModelSolve"]["linemodel"]["tplratio_catalog"] = \
        "lineratiocataloglists/lineratiocatalogs_v16/"
    parameters_dict["galaxy"]["LineModelSolve"]["linemodel"]["tplratio_ismfit"] = True
    parameters_dict["galaxy"]["LineModelSolve"]["linemodel"]["nsigmasupport"] = 8
    parameters_dict["galaxy"]["LineModelSolve"]["linemodel"]["igmfit"] = True
    parameters_dict["photometryTransmissionDir"] = "photometric_transmission/EL-COSMOSv2/"
    parameters_dict["photometryBand"] = ["H", "J", "Y", "riz"]
    return Parameters(parameters_dict, Checker=FakeParametersChecker)


def test_calibration_linecatalog():
    parameters = make_parameters()
    cl = CalibrationLibrary(parameters, calibration_dir)
    cl.load_linecatalog("galaxy", "LineModelSolve")


def test_calibration_lineratiocatalog():
    parameters = make_parameters()
    cl = CalibrationLibrary(parameters, calibration_dir)
    cl.load_linecatalog("galaxy", "LineModelSolve")
    cl.load_line_ratio_catalog_list("galaxy")


def test_calibraton_meiksin():
    parameters = make_parameters()
    cl = CalibrationLibrary(parameters, calibration_dir)
    cl.load_Meiksin()


def test_calibraton_calzetti():
    parameters = make_parameters()
    cl = CalibrationLibrary(parameters, calibration_dir)
    cl.load_calzetti()


def test_calibration_load_templates_catalog():
    parameters = make_parameters()
    cl = CalibrationLibrary(parameters, calibration_dir)
    cl.load_templates_catalog("galaxy")


def test_calibration_load_lsf():
    parameters = make_parameters()
    cl = CalibrationLibrary(parameters, calibration_dir)
    cl.load_lsf()


def test_calibration_load_photometric_bands():
    parameters = make_parameters()
    cl = CalibrationLibrary(parameters, calibration_dir)
    cl.load_photometric_bands()


def test_calibration_load_all():
    parameters = make_parameters()
    cl = CalibrationLibrary(parameters, calibration_dir)
    cl.load_all()
