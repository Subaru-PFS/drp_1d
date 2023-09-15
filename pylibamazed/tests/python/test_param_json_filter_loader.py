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
from pylibamazed.Filter import FilterList, SpectrumFilterItem
from pylibamazed.FilterLoader import ParamJsonFilterLoader
from pylibamazed.Parameters import Parameters


class TestParamJsonFilterLoader:

    def test_load_check_on_json_format(self, mocker):
        jsonFilterLoader = ParamJsonFilterLoader()

        mocker.patch(
            'pylibamazed.Parameters.Parameters.check_params',
            return_value=True
        )
        # Json is not a list
        params = Parameters({"filters": {"i should be a list": ""}})
        with pytest.raises(ValueError, match=r"must be a list"):
            jsonFilterLoader.get_filters(params)

        # Wrong key
        params = Parameters({"filters": [{"wrong key": ""}]})
        with pytest.raises(ValueError, match=r"must have exactly the following keys"):
            jsonFilterLoader.get_filters(params)

        # Missing key
        params = Parameters({"filters": [{"key": "", "instruction": ""}]})
        with pytest.raises(ValueError, match=r"must have exactly the following keys"):
            jsonFilterLoader.get_filters(params)

    def test_load_returns_expected_filters(self, mocker):
        mocker.patch(
            'pylibamazed.Parameters.Parameters.check_params',
            return_value=True
        )
        jsonFilterLoader = ParamJsonFilterLoader()
        params = Parameters({"filters": [
            {"key": "col1", "instruction": "<", "value": 2},
            {"key": "col2", "instruction": ">=", "value": 2}
        ]})
        assert jsonFilterLoader.get_filters(params) == FilterList([
            SpectrumFilterItem("col1", "<", 2),
            SpectrumFilterItem("col2", ">=", 2)
        ])
