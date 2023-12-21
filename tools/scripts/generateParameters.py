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
from pylibamazed.get_json_parameter import get_parameter
import collections.abc
from pylibamazed.Paths import module_root_dir


def deep_update(original, overrider):
    for key, overrideValue in overrider.items():
        if isinstance(overrideValue, collections.abc.Mapping):
            original[key] = deep_update(original.get(key, {}), overrideValue)
        else:
            original[key] = overrideValue
    return original


def generate_parameters(inputString: str, specificJson: dict = None, locations=[], toCpp=False):
    root = inputString.split('|')[0]
    bricks = inputString.split('|')[1:]

    finalParam = get_parameter(root, locations)

    for brick in bricks:
        spectrumModel = brick.split(".")[0]
        method = brick.split(".")[1]

        print("\nobject type\n", spectrumModel)
        print("\nmethod\n", method)

        paramsToUpdate = {f"spectrumModel_{spectrumModel}": {"stages": []}}
        finalParam["spectrumModels"].append(spectrumModel)
        if method == "lineMeasSolve":
            paramsToUpdate[f"spectrumModel_{spectrumModel}"]["stages"].append("lineMeasSolve")
        deep_update(paramsToUpdate[f"spectrumModel_{spectrumModel}"],
                    get_parameter(f"{spectrumModel}.{method}", locations))
        deep_update(finalParam, paramsToUpdate)
    if specificJson is not None:
        deep_update(finalParam, specificJson)

    if toCpp:
        return json.dumps(finalParam, indent=0).replace('"', "\\\"").replace("\n", "\"\n\"")
    return json.dumps(finalParam, indent=4)


if __name__ == "__main__":
    params = generate_parameters(
        "linemeas_solve|galaxy.lineMeasSolve",
        locations=[os.path.join(module_root_dir, "resources/parameters/test-cpp/")],
        toCpp=True,
        specificJson={
            "spectrumModel_galaxy": {
                "lineMeasSolve": {
                    "lineMeasSolve": {
                        "lineModel": {
                            "nSigmaSupport": 14,
                            "fittingMethod": "lbfgsb",
                            "velocityFit": True
                        }
                    }
                }
            }
        }
    )
    print(params)
