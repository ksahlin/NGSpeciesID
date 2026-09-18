"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# using example setup file from https://github.com/pypa/sampleproject/blob/master/setup.py

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(

    name='NGSpeciesID',  # Required
    version='0.4.0',  # Required
    description='Reconstructs viral consensus sequences from a set of ONT reads.',  # Required
    long_description=long_description,  # Optional
    url='https://github.com/ksahlin/NGSpeciesID',  # Optional
    author='Kristoffer Sahlin',  # Optional
    author_email='ksahlin@math.su.se',  # Optional

    # Classifiers help users find your project by categorizing it.
    #
    # For a list of valid classifiers, see
    # https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        #'Intended Audience :: Developers',
        #'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish
        #'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
    ],

    keywords='viral sequeces ONT Oxford Nanopore Technologies long reads',  # Optional

    # You can just specify package directories manually here if your project is
    # simple. Or you can use find_packages().
    #
    # Alternatively, if you just want to distribute a single Python file, use
    # the `py_modules` argument instead as follows, which will expect a file
    # called `my_module.py` to exist:
    #
    #   py_modules=["my_module"],
    #
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),  # Required

    # The code uses f-strings throughout, so 3.6 is the floor on syntax alone.
    # 3.10 is the floor in practice: it is the oldest interpreter for which
    # bioconda still ships parasail-python and python-edlib builds. Reproducible
    # results additionally want >=3.12 -- see the README.
    python_requires='>=3.10, <4',
    # This field lists other packages that your project depends on to run.
    # Any package you put here will be installed by pip when your project is
    # installed, so they must be valid existing projects.
    #
    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    # parasail 1.2.4 is from 2018 and publishes no wheel for any ARM platform,
    # so pinning it made `pip install NGSpeciesID` compile it from source on
    # Apple Silicon and ARM Linux -- a build that fails with
    # "RuntimeError: autoreconf -fi failed". 1.3.4 is no better on ARM (also no
    # aarch64 wheel), so the README installs both libraries from bioconda
    # (medaka pulls in parasail-python and python-edlib) and then runs
    # `pip install --no-deps NGSpeciesID`. These bounds are what a bare
    # `pip install` should ask for when it has to resolve them itself.
    install_requires=['parasail>=1.3.4',
                      'edlib>=1.1.2'],  # Optional
    # dependency_links=[], # Optional
    # List additional groups of dependencies here (e.g. development
    # dependencies). Users will be able to install these using the "extras"
    # syntax, for example:
    #
    #   $ pip install sampleproject[dev]
    #
    # Similar to `install_requires` above, these must be valid existing
    # projects.
    # extras_require={  # Optional
    #     'dev': ['check-manifest'],
    #     'test': ['coverage'],
    # },

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # `pip` to create the appropriate form of executable for the target
    # platform.
    # entry_points={  # Optional
    #     'console_scripts': [
    #         'NGSpeciesID=NGSpeciesID:main()',
    #     ],
    # },
    scripts=['NGSpeciesID'],
)
