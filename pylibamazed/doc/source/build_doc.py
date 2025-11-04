#!/usr/bin/env python3
import os
import shutil
import subprocess
import sys

OUTPUT_DIR = "pylibamazed/doc/build"
MAIN_VERSION = "develop"


def deduplicate_static(output_dir: str, main_version: str):
    # To avoid taking to much space from duplicated static files

    print(f"[deduplicate_static] Running in {output_dir}")
    versions = [
        d
        for d in os.listdir(output_dir)
        if os.path.isdir(os.path.join(output_dir, d)) and not d.startswith("_")
    ]
    if not versions:
        print("No versions found")
        return

    main_static = os.path.join(output_dir, main_version, "_static")
    shared_static = os.path.join(output_dir, "_static")

    if not os.path.exists(main_static):
        main_static = os.path.join(output_dir, versions[0], "_static")

    if not os.path.exists(main_static):
        print("No _static found")
        return

    if not os.path.exists(shared_static):
        print(f"Copying {main_static} → {shared_static}")
        shutil.copytree(main_static, shared_static)

    for version in versions:
        static_path = os.path.join(output_dir, version, "_static")
        if os.path.exists(static_path) and not os.path.islink(static_path):
            shutil.rmtree(static_path)
            os.symlink("../_static", static_path)
            print(f"Linked {static_path} → ../_static")


def main():
    output_dir = sys.argv[1] if len(sys.argv) > 1 else OUTPUT_DIR
    subprocess.run(
        ["sphinx-multiversion", "pylibamazed/doc/source", output_dir, "--dump-metadata"], check=True
    )
    subprocess.run(["sphinx-multiversion", "pylibamazed/doc/source", output_dir], check=True)
    deduplicate_static(output_dir, MAIN_VERSION)


if __name__ == "__main__":
    main()
