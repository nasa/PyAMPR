"""
Program: udf_cmap.py


Purpose: Define colormaps to be used for plotting AMPR data.


Description: A colormap is defined with relatively sharp transitions to
             help highlight transitions in observed data.


Author: Brent Roberts, jason.b.roberts@nasa.gov
Contributor: Timothy Lang, timothy.j.lang@nasa.gov

Rev 1.1, 07.08.2015 - Added amprQC_cmap
Rev 1.0, 05.04.2014 - Created amprTB_cmap

"""
from __future__ import absolute_import
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

# Specify colors in list:
rgb = []
rgb.append([50, 50, 50])
rgb.append([81, 81, 81])
rgb.append([107, 107, 107])
rgb.append([139, 139, 139])
rgb.append([160, 160, 160])
rgb.append([178, 178, 178])
rgb.append([200, 200, 200])

rgb.append([226, 219, 254])
rgb.append([190, 178, 252])
rgb.append([125, 106, 224])
rgb.append([115, 100, 215])
rgb.append([74, 54, 212])
rgb.append([37, 34, 177])

rgb.append([24, 104, 233])
rgb.append([43, 132, 244])
rgb.append([55, 144, 250])
rgb.append([76, 160, 253])
rgb.append([144, 214, 255])

rgb.append([183, 252, 163])
rgb.append([142, 242, 133])
rgb.append([80, 238, 80])
rgb.append([33, 183, 33])
rgb.append([13, 161, 13])

rgb.append([254, 249, 163])
rgb.append([254, 230, 121])
rgb.append([244, 194, 56])
rgb.append([253, 158, 14])
rgb.append([254, 94, 0])
rgb.append([254, 48, 0])
rgb.append([235, 15, 0])
rgb.append([188, 2, 0])
rgb.append([166, 0, 0])

rgb.append([96, 63, 46])
rgb.append([114, 78, 61])
rgb.append([137, 101, 86])
rgb.append([155, 119, 104])
rgb.append([177, 141, 126])
rgb.append([198, 161, 153])
rgb.append([215, 196, 175])
rgb.append([223, 151, 167])
rgb.append([231, 128, 138])
rgb.append([226, 96, 107])
rgb.append([218, 75, 84])
rgb.append([209, 55, 46])

rgbarray = np.asarray(rgb)
# normalize rgb array.
rgbarray = rgbarray / 255.0

# Use "from_list"
# amprTB_cmap=LinearSegmentedColormap.from_list('amprTB',rgbarray);

nrows = rgbarray.shape[0]
# create y values.
xvals = np.linspace(0, 1, nrows)
# Now create dictionary of list of tuples (x,y0,y1), x=R, G, or B
cdict = {}
cdict['red'] = [tuple([xvals[x], rgbarray[x, 0], rgbarray[x, 0]])
                for x in np.arange(nrows)]
cdict['green'] = [tuple([xvals[x], rgbarray[x, 1], rgbarray[x, 1]])
                  for x in np.arange(nrows)]
cdict['blue'] = [tuple([xvals[x], rgbarray[x, 2], rgbarray[x, 2]])
                 for x in np.arange(nrows)]

amprTB_cmap = LinearSegmentedColormap('testcmap', cdict)

# amprQC_cmap definition - TJL
qc_colors = ['#000000', '#33CC33', '#66FF33', '#FFFF00',
             '#FFCC00', '#FF0000']
amprQC_cmap = ListedColormap(qc_colors)
