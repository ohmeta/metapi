import io
from os.path import dirname, join
from setuptools import setup

long_description = open('README.md').read()

def get_version(relpath):
    """read version info from a file without importing it"""
    for line in io.open(join(dirname(__file__), relpath), encoding="cp437"):
        if "__version__" in line:
            if '"' in line:
                return line.split('"')[1]
            elif "'" in line:
                return line.split("'")[1]

setup(
    name='metapipe',
    version=get_version("metapipe/__init__.py"),
    url='https://github.com/ohmeta/metapipe',
    license='MIT',
    author='Jie Zhu and Jiahui Zhu',
    author_email='zhujie@genomics.cn, zhujiahui@genomics.cn',
    description='a metagenomics data processing pipeline to help research',
    long_description=long_description,
    packages=['metapipe'],
    package_data={'': ['metapipe/Snakefile',
                       'metapipe/rules/fastqc.smk',
                       'metapipe/rules/trim.smk',
                       'metapipe/rules/rmhost.smk',
                       'metapipe/rules/qcreport.smk',
                       'metapipe/rules/assembly.smk',
                       'metapipe/rules/alignment.smk',
                       'metapipe/rules/binning.smk',
                       'metapipe/rules/checkm.smk',
                       'metapipe/rules/dereplicate.smk',
                       'metapipe/rules/classification.smk',
                       'metapipe/rules/annotation.smk',
                       'metapipe/__init__.py',
                       'metapipe/metapipe.py',
                       'metapipe/project.py',
                       'metapipe/report.py']},
    include_package_data=True,
    install_requires=[],
    entry_points={},
    classifiers=["Topic :: Scientific/Engineering:: Bio-Informatics"],
)