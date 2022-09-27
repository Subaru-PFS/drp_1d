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

from pylibamazed.AbstractOutput import AbstractOutput,RootStages,ObjectStages
from pylibamazed.Exception import AmazedError
from pylibamazed.redshift import ErrorCode

import pandas as pd

def _create_dataset_from_dict(h5_node, name, source, compress=False):
    df = pd.DataFrame(source)
    records = df.to_records(index=False)
    h5_node.create_dataset(name,
                           len(records),
                           records.dtype,
                           records,
                           compression="lzf")


class H5Writer():

    def __init__(self, output):
        self.output = output

    def write_hdf5_root(self, hdf5_spectrum_node):
        for ds in self.output.get_available_datasets("root"):
            dsg = hdf5_spectrum_node.create_group(ds)
            for attr_name,attr in self.output.get_dataset(ds).items():
                dsg.attrs[attr_name] = attr

    def write_hdf5_object_level(self, object_type, object_results_node):
        level = "object"
        for ds in self.output.get_available_datasets("object",
                                                     object_type=object_type):
            ds_size = self.output.get_dataset_size(object_type, ds)
            dataset = self.output.get_dataset(ds,
                                              object_type=object_type)
            if ds_size > 1:
                _create_dataset_from_dict(object_results_node,
                                          ds,
                                          dataset)
            else:
                object_results_node.create_group(ds)
                for attr_name,attr in dataset.items():
                    object_results_node.get(ds).attrs[attr_name] = attr

    def write_hdf5_candidate_level(self, object_type, object_results_node):
        if self.output.parameters.get_solve_method(object_type):
            candidates = object_results_node.create_group("candidates")
            level = "candidate"
            nb_candidates = self.output.get_nb_candidates(object_type)
            for rank in range(nb_candidates):
                candidate = candidates.create_group(self.output.get_candidate_group_name(rank))
                for ds in self.output.get_available_datasets("candidate",object_type = object_type):
                    dataset = self.output.get_dataset(ds,
                                               object_type=object_type,
                                               rank=rank)
                    ds_dim = self.output.get_dataset_size(object_type, ds, rank)
                    if ds_dim == 1:
                        candidate.create_group(ds)
                        for attr_name, attr in dataset.items():
                            candidate.get(ds).attrs[attr_name] = attr
                    else:
                        _create_dataset_from_dict(candidate,
                                                  ds,
                                                  dataset)
        
    def write_hdf5(self,hdf5_root,spectrum_id):
        try:
            obs = hdf5_root.create_group(spectrum_id)
            self.write_hdf5_root(obs)

            for stage in RootStages:
                if self.output.has_error(None,stage):
                    self.write_error(obs, None, stage)
            if self.output.has_error(None,"init"):
                return
            for object_type in self.output.object_types:
                object_results = obs.create_group(object_type) #h5
                self.write_hdf5_object_level(object_type, object_results)
#                self.write_hdf5_method_level(object_type, object_results)
                if self.output.get_nb_candidates(object_type) > 0 :
                    self.write_hdf5_candidate_level(object_type, object_results)
                for stage in ObjectStages:
                    if self.output.has_error(object_type, stage):
                        self.write_error(object_results,object_type, stage)
                    
        except Exception as e:
            raise AmazedError(ErrorCode.EXTERNAL_LIB_ERROR,"Failed writing h5:".format(e))

    def write_error(self, hdf5_object_node, object_type, stage):
        error = hdf5_object_node.create_group(stage + "_error")
        for k,v in self.output.get_error(object_type,stage).items():
            error.attrs[k]= v
