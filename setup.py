#!/usr/bin/env python
"""
This is a prototype MultiQC module for Atropos.
"""

from setuptools import setup, find_packages

version = '0.1'

setup(
    name = 'multiqc_atropos',
    version = version,
    author = 'John Didion',
    author_email = 'john.didion@nih.gov',
    description = "MultiQC plugin for Atropos",
    long_description = __doc__,
    keywords = 'bioinformatics',
    url = 'https://github.com/jdidion/atropos-multiqc',
    download_url = 'https://github.com/jdidion/atropos-multiqc',
    license = 'CC0',
    packages = find_packages(),
    include_package_data = True,
    entry_points = {
        'multiqc.modules.v1': [
            'atropos = multiqc_atropos.modules.atropos.atropos:MultiqcModule'
        ]
    },
    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)
