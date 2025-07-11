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

import numpy as np
import sklearn
from collections import OrderedDict
from pylibamazed.Exception import APIException
from pylibamazed.redshift import ErrorCode
from pylibamazed.redshift import CLog
from pylibamazed.ResultStoreOutput import ResultStoreOutput
from pylibamazed.AbstractReliabilitySolver import AbstractReliabilitySolver, register_reliability_solver

zlog = CLog.GetInstance()

class SkLearnSolve(AbstractReliabilitySolver):
    def Compute(self, source):
        output = ResultStoreOutput(
            source.GetResultStore(),
            self.parameters,
            auto_load=False,
            extended_results=False,
        )
        zlog.LogInfo(f"SkLearnSolver: Compute reliability for {self.object_type}")
        classifier = self.calibration_library.reliability["sklearn"][self.object_type]["classifier"]
        classes = self.calibration_library.reliability["sklearn"][self.object_type]["classes"]
        success = classes[-1]
        return self.get_probas(output, classifier, classes)[success]

    def get_probas(self, output, classifier, classes):
        output.load_object_level(self.object_type)
        attributes = OrderedDict( {
            'A_IMAGE':"",
            'ELLIPTICITY':"",
            'POINT_LIKE_PROB':"",
            'MAG_VIS':"",
            'MER_Y_MAG':"",
            'MER_J_MAG':"",
            'MER_H_MAG':"",
            'NDITH':"",
            'LSF_SIG':"",
            'Z':"galaxy.Redshift",
            'Z_ERR':"galaxy.RedshiftUncertainty",
            'Z_PROB':"",
            'HA_FLUX':"galaxy.lfHaNII",
            'HA_SNR':"galaxy.snrHaNII",
            'OII_FLUX':"galaxy.lfOII",
            'OII_SNR':"galaxy.snrOII",
            'VEL_EMI':"galaxy.VelocityEmission",
            'RELIABILITY':"",
            'SPEC_COLOR':"",
            'SNR_MEAN':"",
            'SNR_STD':"",
            'NDITH_MEAN':"",
            'NDITH_STD':"" 
            }
        )
        col_used = ['LSF_SIG', 'Z', 'Z_ERR', 'HA_FLUX', 'HA_SNR', 'OII_FLUX', 'OII_SNR', 'VEL_EMI']

        v = np.ndarray([len(col_used)])
        idx = 0 
        for k,att in attributes.items():
            if k in col_used and att!="":
                v[idx] =  output.get_attribute_short(self.object_type, att)
                idx += 1
        ret = dict()
        probas = classifier.predict_proba(v.reshape(1,-1))
        for i,c in enumerate(classes):
            ret[c] = float(probas[0, i])
        zlog.LogInfo(f"SkLearnSolver: probas are {ret}")
        return ret

register_reliability_solver("skLearnSolver", SkLearnSolve, "sk")
