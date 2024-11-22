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
from pathlib import Path

from jsonschema.exceptions import ValidationError
from pylibamazed.Exception import APIException
from pylibamazed.FilterLoader import ParamJsonFilterLoader
from pylibamazed.Paths import module_root_dir
from pylibamazed.redshift import CFlagWarning, ErrorCode
from pylibamazed.ParametersChecker import ParametersChecker
from referencing import Registry, Resource
from jsonschema import Draft202012Validator


zflag = CFlagWarning.GetInstance()


class JsonSchemaFileAccessor:
    def __init__(self, version: int = 2):
        self.version = version

    def get_json_schema(self) -> dict:
        with open(self.json_schema_filename()) as schemaFile:
            jsonSchema = json.load(schemaFile)
        return jsonSchema

    def json_schema_filename(self) -> str:
        return os.path.join(self._json_schema_path(), "general.json")

    def _json_schema_path(self) -> str:
        return os.path.join(module_root_dir, "resources", f"jsonschema-v{self.version}")

    def retrieve_from_filesystem(self, uri: str):
        schema_path = Path(self._json_schema_path())
        path = schema_path / Path(uri)
        contents = json.loads(path.read_text())
        return Resource.from_contents(contents)


class JsonParametersChecker(ParametersChecker):
    def __init__(
        self,
        parameters: dict,
        version: int = 2,
        FilterLoader=ParamJsonFilterLoader,
        FileAccessor=JsonSchemaFileAccessor,
    ):
        self.filter_loader = FilterLoader()
        self.file_accessor = FileAccessor(version)
        self.parameters = parameters

    def check(self) -> None:
        # Parameters dict must be at "raw parameters dict"
        # How to link the different json schema files
        json_schema = self.file_accessor.get_json_schema()

        # Add the base schema to the registry
        registry = Registry(retrieve=self.file_accessor.retrieve_from_filesystem)
        validator = Draft202012Validator(json_schema, registry=registry)
        try:
            validator.validate(self.parameters)
        except ValidationError as e:
            raise APIException(ErrorCode.INVALID_PARAMETER_FILE, e.message)
