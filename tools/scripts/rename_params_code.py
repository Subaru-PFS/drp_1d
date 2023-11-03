#!/usr/bin/env python

import glob
import os
import pandas as pd

project_path = "/home/fdufresne/Projects/refacto_parameters/cpf-redshift"

renaming_filename = os.path.join(project_path, "pylibamazed/auxdir/pylibamazed/rename_params_v1_to_treed.csv")

folders_to_process = [
    "pylibamazed/auxdir/pylibamazed/jsonschema-v2",
    "pylibamazed/python",
    "pylibamazed/tests",
    "RedshiftLibrary"
]

list_to_exclude = [
    "pylibamazed/python/pylibamazed/resources",
    "redshift.py",
    "egg-info",
    ".bak",
    ".pyc"
]


def list_files_to_rename():
    to_rename = []
    for folder in folders_to_process:
        filenames = glob.glob(
            os.path.join(project_path,folder, "**/*"),  recursive=True)
        for filename in filenames:
            keep_file = True
            for to_exclude in list_to_exclude:
                if to_exclude in filename:
                    keep_file = False
                    break
            if keep_file:
                to_rename.append(filename)

    return to_rename


def make_renaming_df():
    df = pd.read_csv(renaming_filename)
    return df


def massive_rename():
    files_to_rename = list_files_to_rename()
    renaming_df = make_renaming_df()

    for filename in files_to_rename:
        full_filename = os.path.join(project_path, filename)

        if os.path.isfile(full_filename):
            print("\nfile ", filename)
            with open(full_filename, 'r') as file:
                file_data = file.read()
            for _, row in renaming_df.iterrows():
                file_data = file_data.replace(f"{row.loc['old_name']}\"", f"{row.loc['new_name']}\"")
                file_data = file_data.replace(f"\"{row.loc['old_name']}", f"\"{row.loc['new_name']}")
                file_data = file_data.replace(f"\'{row.loc['old_name']}\'", f"\'{row.loc['new_name']}\'")
                file_data = file_data.replace(
                    f"\\\"{row.loc['old_name']}\\\"", f"\\\"{row.loc['new_name']}\\\"")
                file_data = file_data.replace(
                    f"\\\'{row.loc['old_name']}\\\'", f"\\\'{row.loc['new_name']}\\\'")
                file_data = file_data.replace(f".{row.loc['old_name']}.", f".{row.loc['new_name']}.")

            with open(full_filename, "w") as file:
                file.write(file_data)


if __name__ == "__main__":
    massive_rename()
