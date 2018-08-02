from __future__ import absolute_import
# Setup script for the pyampr package
# $Id: setup.py,v 1.0 2014/07/28 tjlang Exp $
#
# Usage: python setup.py install
#
from distutils.core import setup

try:
    # add download_url syntax to distutils
    from distutils.dist import DistributionMetadata
    DistributionMetadata.classifiers = None
    DistributionMetadata.download_url = None
except:
    pass

VERSION = '1.6'

DESCRIPTION = "The Python Advanced Microwave Precipitation Radiometer " + \
    "Data Toolkit (PyAMPR) - a package to read, analyze, and display AMPR data"

LONG_DESCRIPTION = """The Advanced Microwave Precipitation Radiometer (AMPR)
is an airborne passive microwave radiometer managed by NASA Marshall Space
Flight Center. Download AMPR data from http://ghrc.nsstc.nasa.gov.
AMPR brightness temperature data from NASA field projects
are in ASCII or netCDF format. This python script defines a class that will
read in single file from an individual aircraft flight and pull out
timing, brightness temperatures from each channel, geolocation, and
other information and store them as attributes using numpy
arrays of the appropriate type. The file is read and the data are populated
when the class is instantiated with the full path and name of an AMPR file.
Numerous visualization methods are provided, including track plots,
strip charts, and Google Earth KMZs. In addition, polarization
deconvolution is available."""

setup(
    name="pyampr",
    version=VERSION,
    author="Timothy J. Lang",
    author_email="timothy.j.lang@nasa.gov",
    url="http://github.com/nasa/PyAMPR/",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    download_url="http://github.com/nasa/PyAMPR/",
    license="LICENSE.md",
    packages=["pyampr"],
    platforms="Python 2.7, 3.4",
    classifiers=["""
        Development Status :: Beta,
        Programming Language :: Python :: 2.7, 3.4
        Topic :: Scientific/Engineering
        Topic :: Scientific/Engineering :: Atmospheric Science
        Operating System :: Unix
        Operating System :: POSIX :: Linux
        Operating System :: MacOS
        """]
    )
