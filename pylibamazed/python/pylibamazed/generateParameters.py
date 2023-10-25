import json
import os
from pylibamazed.get_json_parameter import get_parameter_base, get_parameter
import collections.abc
from pylibamazed.Paths import module_root_dir

# TODO add license


def deep_update(original, overrider):
    # print("\nkey", original, "\nvalue", overrider)
    for key, overrideValue in overrider.items():
        if isinstance(overrideValue, collections.abc.Mapping):
            # print('\nin if')
            # print("\nkey", key, "\nold value", original.get(key, "not found"), "\nnew value", overrideValue)
            original[key] = deep_update(original.get(key, {}), overrideValue)
        else:
            # print("\nin else", "\nkey", key, "\nold value", original.get(
            # key, "not found"), "\nnew value", overrideValue)
            original[key] = overrideValue
    # print("AAAAAAAAAAAA overrider", overrider)
    return original

# definition example : se8|galaxy.LineModelSolve


def generate_parameters(inputString: str, specificJson: dict = None, locations=[], toCpp=False):
    root = inputString.split('|')[0]
    bricks = inputString.split('|')[1:]

    # baseRoot = get_parameter_base("root")
    finalParam = get_parameter(root, locations)

    # finalParam = deep_update(baseRoot, datasetParam)

    # TODO here limited to one instance per method (e.g. only once redshift solver etc)
    for brick in bricks:
        # elts = o.split(".")[2:]
        spectrumModel = brick.split(".")[0]
        method = brick.split(".")[1]

        print("\nobject type\n", spectrumModel)
        print("\nmethod\n", method)

        paramsToUpdate = {f"spectrumModel_{spectrumModel}": {"stages": []}}
        # paramsToUpdate[spectrumModel] = get_parameter(method, locations)
        finalParam["spectrumModels"].append(spectrumModel)
        if method == "lineMeasSolve":
            paramsToUpdate[f"spectrumModel_{spectrumModel}"]["stages"].append("lineMeasSolve")
        deep_update(paramsToUpdate[f"spectrumModel_{spectrumModel}"],
                    get_parameter(f"{spectrumModel}.{method}", locations))
        # update(op[object_type][method], get_parameter(f"{object_type}.{method}.{root}", locations))
        # for v in elts:
        #     o_v = get_parameter(v, locations)
        #     update(op[object_type][method], o_v)
        deep_update(finalParam, paramsToUpdate)
    if specificJson is not None:
        deep_update(finalParam, specificJson)

    # return json.dumps(finalParam, indent=4)
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
