import json
import os

# TODO add license

module_root_dir = os.path.split(__file__)[0]

base_param_dir = os.path.join(module_root_dir, "resources", "parameters")


def get_json(directory, name):
    path = os.path.join(directory, f"{name}.json")
    print(f"try open {path}")
    with open(path) as f:
        try:
            json_ = json.load(f)
        except Exception as e:
            raise Exception(f"Failed to open {path} : {e}")
    return json_


def get_parameter_base(name):
    return get_json(base_param_dir, name)


def get_parameter(name, locations=[]):
    directories = locations.copy()
    directories.append(base_param_dir)
    for dir_ in directories:
        try:
            return get_json(dir_, name)
        except FileNotFoundError:
            print(f"{name} not found in {dir_}")
            continue
        except Exception as e:
            raise Exception(f"{dir_}/{name}.json malformed : {e}")
        else:
            print(f"OK file {name} found in {dir_}")
        finally:
            print("yyoyoyo")
