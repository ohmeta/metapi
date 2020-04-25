#!/usr/bin/env python
import setuptools

with open("README.md") as f:
    long_description = f.read()

version = "0.7.1"
print(
    """------------------------
Installing metapi version {}
------------------------
""".format(
        version
    )
)

setuptools.setup(
    name="metapi",
    version=version,
    url="https://github.com/ohmeta/metapi",
    license="GPLv3",
    author="Jie Zhu",
    author_email="zhujie@genomics.cn",
    description="a pipeline to construct a genome catalogue from metagenomics data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=["metapi"],
    package_data={
        "metapi": ["metapi/*.yaml", "metapi/Snakefile", "metapi/rules/*.smk",]
    },
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
        "openpyxl",
        "snakemake",
        "ruamel.yaml",
        "biopython",
        "InSilicoSeq",
        "multiqc",
        "quast",
        "checkm-genome",
        "gtdbtk",
        "drep",
    ],
    zip_safe=False,
    entry_points={"console_scripts": ["metapi = metapi.corer:main"]},
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.6",
)

print(
    """
---------------------------------
metapi installation complete!
---------------------------------
For help in running metapi, please see the documentation available
at https://github.com/ohmeta/metapi or run metapi --help
"""
)
