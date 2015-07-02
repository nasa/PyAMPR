Title/Version
-------------
Python AMPR Data Toolkit (PyAMPR) v1.4.2
Last changed 07/02/2015


Lead Author
-----------
Timothy Lang
NASA MSFC
timothy.j.lang@nasa.gov
(256) 961-7861


Contributing Authors
--------------------
Brent Roberts
NASA MSFC
jason.b.roberts@nasa.gov
(256) 961-7477


Overview
--------
The Advanced Microwave Precipitation Radiometer (AMPR) is an airborne 
passive microwave radiometer managed by NASA Marshall Space Flight Center.
Download AMPR data from http://ghrc.nsstc.nasa.gov.
AMPR brightness temperature data from NASA field projects
are in ASCII or netCDF format. This python script defines a class that will 
read in single file from an individual aircraft flight and pull out
timing, brightness temperatures from each channel, geolocation, and
other information and store them as attributes using numpy 
arrays of the appropriate type. The file is read and the data are populated when
the class is instantiated with the full path and name of an AMPR file.
Numerous visualization methods are provided, including track plots,
strip charts, and Google Earth KMZs. In addition, polarization
deconvolution is available.


Installation and Use
--------------------
Dependencies: Python 2.7, NumPy, matplotlib, Basemap,
              PyLab, os, time, simplekml, datetime, calendar, 
              codecs, gzip, netCDF4
Most of these are provided with standard Python distributions.
You may need to install Basemap via your Python distribution's
package directory. The simplekml package can be found at 
https://code.google.com/p/simplekml/
Since Google Code is being shut down in the near future, I have also
migrated this package to https://github.com/tjlang/simplekml.

In the same directory as this README is setup.py, to install this
package enter the following command at the prompt:
python setup.py install

Then to import, in your python program include:
import pyampr

To read an AMPR TB file type:
ampr_data = pyampr.AmprTb('FILE_NAME_HERE', project='PROJECT_NAME_HERE')

Then the ampr_data object will have access to all the plotting and analysis 
methods. Use 'help(pyampr.AmprTb)' to find out more.

In particular, 'help(pyampr.AmprTb.read_ampr_tb_level2b)' will give a full 
rundown on the data structure.

A demonstration IPython notebook can be found in the notebooks directory.

A simple interactive testing notebook is available in the test directory.

The following conference presentation describes PyAMPR (among other modules):
https://ams.confex.com/ams/95Annual/webprogram/Paper262779.html
