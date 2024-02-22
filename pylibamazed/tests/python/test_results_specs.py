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


import pandas as pd
from pylibamazed.OutputSpecifications import ResultsSpecifications

rs = ResultsSpecifications()


def test_get_dataframe_by_criteria():
    REDSHIFT_INDEX = [0]
    MODEL_INDEXES = [50, 51]
    METHOD_INDEX = [129]
    name_data = {
        "name": "Redshift",
        "dataset": "model_parameters",
        "extended_results": False,
        "level": "candidate",
        "OperatorResult_name": "Redshift",
        "ResultStore_key": "extrema_results",
        "hdf5_type": "f8"
    }
    name_df = pd.DataFrame(name_data, index=REDSHIFT_INDEX)
    test_name_df = rs.get_df_by_name("Redshift")
    assert test_name_df.equals(name_df)

    dataset_data = {
        "name": ["ModelLambda", "ModelFlux"],
        "dataset": ["model<ObsID>", "model<ObsID>"],
        "extended_results": [False, False],
        "level": ["candidate", "candidate"],
        "OperatorResult_name": ["ModelLambda[obs_id]", "ModelFlux[obs_id]"],
        "ResultStore_key": ["extrema_results", "extrema_results"],
        "hdf5_type": ["f8", "f8"]
    }
    dataset_df = pd.DataFrame(dataset_data, index=MODEL_INDEXES)
    test_dataset_df = rs.get_df_by_dataset("model<ObsID>")
    assert test_dataset_df.equals(dataset_df)

    level_data = {
        "name": "<MethodType>WarningFlags",
        "dataset": "warningFlag",
        "extended_results": False,
        "level": "method",
        "OperatorResult_name": "flagValue",
        "ResultStore_key": "warningFlag",
        "hdf5_type": "i"
    }
    level_df = pd.DataFrame(level_data, index=METHOD_INDEX)
    test_level_df = rs.get_df_by_level("method")
    assert test_level_df.equals(level_df)
