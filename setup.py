#!/usr/bin/env python
import os
import sys
from setuptools import setup

exec(open("metapi/__about__.py").read())

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

with open("README.md") as f:
    long_description = f.read()

packages = ["metapi"]
package_data = {"metapi": ["metapi/*.yaml", "metapi/Snakefile", "metapi/rules/*.smk"]}
data_files = [(".", ["LICENSE", "README.md"])]

entry_points = {"console_scripts": ["metapi=metapi.corer:main"]}

requires = [
    "numpy",
    "pandas",
    "openpyxl",
    "snakemake",
    "ruamel.yaml",
    "biopython>=1.73",
    "InSilicoSeq",
    "multiqc",
    "quast",
    "checkm-genome",
    "gtdbtk",
    "drep",
]

classifiers = [
    "Development Status :: 3 - Alpha",
    "Environment :: Console",
    "Intended Audience :: Developers",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Natural Language :: English",
    "Operating System :: OS Independent",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name="metapi",
    version=__version__,
    description="a pipeline to construct a genome catalogue from metagenomics data",
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=packages,
    package_data=package_data,
    data_files=data_files,
    entry_points=entry_points,
    install_requires=requires,
    author=__author__,
    author_email="alienchuj@gmail.com",
    url="https://github.com/ohmeta/metapi",
    license="GPLv3",
    classifiers=classifiers,
)
