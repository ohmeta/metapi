#!/usr/bin/env python
import shutil
import sys


program_list = ["sickle", "bwa", "samtools", "megahit", "metabat2", "pigz"]


def checking_dependencies(program_list):
    install = []
    exit = False
    for program in program_list:
        where = shutil.which(program)
        if where is None:
            exit = True
            install.append(program)
            print(program + ":\tno")
        else:
            print(program + ":\tyes")

    if exit:
        if "metabat2" not in install:
            print("\npelase use conda to install these program:")
            print("conda install %s" % (" ".join(install)))
        else:
            install_info = " ".join(install).replace("metabat2", "")
            print("\npelase use conda to install these program:")
            if install_info != "":
                print("conda install %s" % install_info)
            print("conda install -c ursky metabat2")
        sys.exit()


checking_dependencies(program_list)
