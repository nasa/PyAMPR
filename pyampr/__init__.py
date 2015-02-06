#!/usr/bin/python
#
#  The Python Advanced Microwave Precipitation Radiometer Data Toolkit (PyAMPR)
#  A package to read, analyze, and display AMPR data
#  
#  Download AMPR data from http://ghrc.nsstc.nasa.gov.
#  AMPR brightness temperature data from NASA field projects
#  are in ASCII format. This python script defines a class that will
#  read in single file from an individual aircraft flight and pull out
#  timing, brightness temperatures from each channel, geolocation, and
#  other information and store them as attributes using numpy
#  arrays of the appropriate type. The approach is to ingest the entire
#  file as a text string and then parse and type convert as necessary.
#  The file is read and the data are populated when the class is
#  instantiated with the full path and name of an AMPR file.
#  Numerous visualization methods are provided, including track plots,
#  strip charts, and Google Earth KMZs. In addition, polarization
#  deconvolution is available.
# 
__author__ = "Timothy J. Lang"

__email__ = "timothy.j.lang@nasa.gov"

__version__ = "1.3.2"

__doc__ = """pyampr v%s by %s

PyAMPR is a package to read, analyze, and display AMPR data

Please e-mail bug reports to: %s""" % (__version__, __author__,__email__)

from pyampr import (AmprTb, _get_timestring_and_sod, _get_sod,
               _method_footer_printout, _method_header_printout,
               _print_times_not_valid)


