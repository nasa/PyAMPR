from __future__ import absolute_import
# !/usr/bin/python
#
# The Python Advanced Microwave Precipitation Radiometer Data Toolkit (PyAMPR)
# A package to read, analyze, and display AMPR data
#
#
from .pyampr import AmprTb
# from .pyampr import (AmprTb, _get_timestring_and_sod, _get_sod,
#                      _method_footer_printout, _method_header_printout,
#                      _print_times_not_valid)
from .misc_tools import read_aircraft_nav_into_awot

__author__ = "Timothy J. Lang"

__email__ = "timothy.j.lang@nasa.gov"

__version__ = "1.5.2"

__doc__ = """pyampr v%s by %s

PyAMPR is a package to read, analyze, and display AMPR data

Please e-mail bug reports to: %s""" % (__version__, __author__, __email__)
