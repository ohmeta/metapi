import io
import os

import setuptools

long_description = "a metagenomics data processing pipeline to help research"


def get_version(relpath):
    """read version info from a file without importing it"""
    for line in io.open(os.path.join(os.path.dirname(__file__), relpath)):
        if "__version__" in line:
            if '"' in line:
                return line.split('"')[1]
            elif "'" in line:
                return line.split("'")[1]


setuptools.setup(
    name='metapipe',
    version=get_version("metapipe/__init__.py"),
    url='https://github.com/ohmeta/metapipe',
    license='MIT',
    author='Jie Zhu',
    author_email='zhujie@genomics.cn',
    description='a metagenomics data processing pipeline to help research',
    long_description=long_description,
    packages=['metapipe'],
    package_data={'': ['metapipe/Snakefile',
                       'metapipe/rules/step.smk',
                       'metapipe/rules/simulation.smk',
                       'metapipe/rules/fastqc.smk',
                       'metapipe/rules/trimming.smk',
                       'metapipe/rules/rmhost.smk',
                       'metapipe/rules/assembly.smk',
                       'metapipe/rules/alignment.smk',
                       'metapipe/rules/binning.smk',
                       'metapipe/rules/checkm.smk',
                       'metapipe/rules/dereplication.smk',
                       'metapipe/rules/classification.smk',
                       'metapipe/rules/annotation.smk',
                       'metapipe/rules/profilling.smk',
                       'metapipe/__init__.py',
                       'metapipe/metapipe.py',
                       'metapipe/metacheck.py',
                       'metapipe/metaconfig.py',
                       'metapipe/metasample.py',
                       'metapipe/metareport.py',
                       'metapipe/metacluster.yaml',
                       'metapipe/metaconfig.yaml']},
    include_package_data=True,
    install_requires=["pandas", "ruamel.yaml"],
    entry_points={},
    classifiers=["Topic :: Scientific/Engineering:: Bio-Informatics"],
)
