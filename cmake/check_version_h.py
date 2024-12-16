#!/usr/bin/env python
import os,subprocess

version_h_ref="RedshiftLibrary/RedshiftLibrary/version.h.ref"
version_h    ="RedshiftLibrary/RedshiftLibrary/version.h"

print (f"{version_h_ref=}")
print (f"{version_h=}")

if (os.path.exists(version_h_ref)):
    print("File H_REF exist, run diff H_REF H")
    status=subprocess.run(["diff",version_h,version_h_ref], stdout=subprocess.PIPE)
    print(f"{status.stdout.decode()=}")
    if (status.stdout.decode()=="") :
        print(f"same, Copy version H_REF to H to reduce compilation time")
        status=subprocess.run(["cp","-p",version_h_ref,version_h], stdout=subprocess.PIPE)
        print(f"status cp : {status.stdout.decode()}")
    else :
        print(f"New H, Copy version H to H_REF")
        status=subprocess.run(["cp",version_h,version_h_ref], stdout=subprocess.PIPE)
else:
    print("File H_REF does not exist")
    if (os.path.exists(version_h)):
        print(f">>> {version_h} exist")
        print(f"First time Copy version H to H_REF")
        status=subprocess.run(["cp","-p",version_h,version_h_ref], stdout=subprocess.PIPE)
        print(f"status cp : {status.stdout.decode()}")
    else:
        print(">>> RedshiftLibrary/RedshiftLibrary/version.h does not exist")


    