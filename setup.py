#!/usr/bin/env python
import setuptools

with open("README.md") as f:
    long_description = f.read()

version = "0.1.3"
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
    description='a metagenomics data processing pipeline to help research',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=['metapi'],
    package_data={'metapi': ['metapi/Snakefile',
                             'metapi/metaconfig.yaml'
                             'metapi/metacluster.yaml',
                             'metapi/rules/step.smk',
                             'metapi/rules/simulation.smk',
                             'metapi/rules/fastqc.smk',
                             'metapi/rules/trimming.smk',
                             'metapi/rules/rmhost.smk',
                             'metapi/rules/assembly.smk',
                             'metapi/rules/alignment.smk',
                             'metapi/rules/binning.smk',
                             'metapi/rules/checkm.smk',
                             'metapi/rules/dereplication.smk',
                             'metapi/rules/classification.smk',
                             'metapi/rules/annotation.smk',
                             'metapi/rules/profilling.smk',
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
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Visualization'
    ],
)

print("""
---------------------------------
metapi installation complete!
---------------------------------
For help in running metapi, please see the documentation available
at https://github.com/ohmeta/metapi or run metapi --help
""")
