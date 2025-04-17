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
from pylibamazed.Exception import APIException
from pylibamazed.JsonParametersChecker import JsonSchemaFileAccessor, JsonParametersChecker
from pylibamazed.Paths import module_root_dir


class TestParametersChecker:
    def test_json_schema_error(self):
        param_dict = {}
        with pytest.raises(APIException, match=r"INVALID_PARAMETER_FILE"):
            JsonParametersChecker(param_dict, 1).check()


class TestJsonSchemaFileAccessor:
    version = 2
    json_file_accessor = JsonSchemaFileAccessor(version)

    def test_json_schema_filename(self):
        assert (
            self.json_file_accessor.json_schema_filename()
            == module_root_dir + "/resources/jsonschema-v2/general.json"
        )
