#!/usr/bin/env python
import setuptools

with open("README.md") as f:
    long_description = f.read()

version = "0.5.0"
print("""------------------------
Installing metapi version {}
------------------------
""".format(version))

setuptools.setup(
    name='metapi',
    version=version,
    url='https://github.com/ohmeta/metapi',
    license='GPLv3',
    author='Jie Zhu',
    author_email='zhujie@genomics.cn',
    description='a pipeline to construct a genome catalogue from metagenomics data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=['metapi'],
    package_data={'metapi': ['metapi/Snakefile',
                             'metapi/metaconfig.yaml'
                             'metapi/metacluster.yaml',
                             'metapi/__init__py',
                             'metapi/metacheck.py',
                             'metapi/metaconfig.py',
                             'metapi/metapi.py',
                             'metapi/metareport.py',
                             'metapi/metasample.py',
                             'metapi/rules/alignment.smk',
                             'metapi/rules/annotation.smk',
                             'metapi/rules/assembly.smk',
                             'metapi/rules/binning.smk',
                             'metapi/rules/burst.smk',
                             'metapi/rules/checkm.smk',
                             'metapi/rules/classification.smk',
                             'metapi/rules/coassembly.smk',
                             'metapi/rules/cobinning.smk',
                             'metapi/rules/dereplication.smk',
                             'metapi/rules/fastqc.smk',
                             'metapi/rules/metaquast.smk',
                             'metapi/rules/prediction.smk',
                             'metapi/rules/profilling.smk',
                             'metapi/rules/rmhost.smk',
                             'metapi/rules/simulation.smk',
                             'metapi/rules/step.smk',
                             'metapi/rules/trimming.smk',
                             ]},
    include_package_data=True,
    install_requires=['pandas', 'ruamel.yaml'],
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'metapi = metapi.metapi:main']
    },
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
)

print("""
---------------------------------
metapi installation complete!
---------------------------------
For help in running metapi, please see the documentation available
at https://github.com/ohmeta/metapi or run metapi --help
""")
