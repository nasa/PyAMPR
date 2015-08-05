"""
google_earth_tools v1.3
Collator and Editor - Timothy J. Lang of NASA MSFC (timothy.j.lang@nasa.gov)
Library of tools for creating KMZ files. Amalgamated from several sources and
edited. Requires simplekml from https://code.google.com/p/simplekml/.
Also requires numpy, and matplotlib. Tested with Python 2.7 and 3.4.
Last edited - 07/08/2015

Change Log
----------
v1.3 Major Changes (07/08/2015)
1. Made code pep8 and Python 3.4 compatible.

v1.2 Major Changes (09/24/2014):
1. Deleted kml_contour, kml_begin, kml_end, line begin, line_end,
   and place_label functions as the write_ampr_kml method was removed
   from pyampr.AmprTb, given the superior write_ampr_kmz method.

v1.1 Major Changes:
1. Added timespan stamping capability via passing the times string list as
   an argument to make_kml().

"""
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
try:
    from simplekml import (Kml, OverlayXY, ScreenXY, Units, RotationXY,
                           AltitudeMode, Camera)
    SIMPLEKML_FLAG = True
except ImportError:
    SIMPLEKML_FLAG = False


def gearth_fig(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, pixels=1024):
    """
    Return a Matplotlib `fig` and `ax` handles for a Google-Earth Image.
    TJL - Obtained from
    http://ocefpaf.github.io/python4oceanographers/blog/2014/03/10/gearth/

    """
    aspect = np.cos(np.mean([llcrnrlat, urcrnrlat]) * np.pi/180.0)
    xsize = np.ptp([urcrnrlon, llcrnrlon]) * aspect
    ysize = np.ptp([urcrnrlat, llcrnrlat])
    aspect = ysize / xsize

    if aspect > 1.0:
        figsize = (10.0 / aspect, 10.0)
    else:
        figsize = (10.0, 10.0 * aspect)

    if False:
        plt.ioff()  # Make `True` to prevent KML components from popping up
    fig = plt.figure(figsize=figsize, frameon=False, dpi=pixels//10)
    # KML friendly image. If using basemap try: `fix_aspect=False`.
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(llcrnrlon, urcrnrlon)
    ax.set_ylim(llcrnrlat, urcrnrlat)
    return fig, ax


def make_kml(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat,
             figs, colorbar=None, times=None, **kw):
    """
    TODO: LatLon bbox, list of figs, optional colorbar figure,
    and several simplekml kw ...
    TJL - Obtained from
    http://ocefpaf.github.io/python4oceanographers/blog/2014/03/10/gearth/

    """
    if not SIMPLEKML_FLAG:
        print('***ERROR!***')
        print('simplekml not installed, download from',
              'https://pypi.python.org/pypi/simplekml/')
        return
    kml = Kml()
    altitude = kw.pop('altitude', 2e7)
    roll = kw.pop('roll', 0)
    tilt = kw.pop('tilt', 0)
    altitudemode = kw.pop('altitudemode', AltitudeMode.relativetoground)
    camera = Camera(latitude=np.mean([urcrnrlat, llcrnrlat]),
                    longitude=np.mean([urcrnrlon, llcrnrlon]),
                    altitude=altitude, roll=roll, tilt=tilt,
                    altitudemode=altitudemode)

    kml.document.camera = camera
    draworder = 0

    for fig in figs:  # NOTE: Overlays are limited to the same bbox.
        draworder += 1
        ground = kml.newgroundoverlay(name='GroundOverlay')
        ground.draworder = draworder
        ground.visibility = kw.pop('visibility', 1)
        ground.name = kw.pop('name', 'overlay')
        ground.color = kw.pop('color', '9effffff')
        ground.atomauthor = kw.pop('author', 'ocefpaf')
        ground.latlonbox.rotation = kw.pop('rotation', 0)
        ground.description = kw.pop('description', 'Matplotlib figure')
        ground.gxaltitudemode = kw.pop('gxaltitudemode',
                                       'clampToSeaFloor')
        if times:
            ground.timespan.begin = times[0]
            ground.timespan.end = times[1]
        ground.icon.href = fig
        # TJL - swapping west/east to match LL vs. UR properly
        ground.latlonbox.west = llcrnrlon
        ground.latlonbox.south = llcrnrlat
        ground.latlonbox.north = urcrnrlat
        ground.latlonbox.east = urcrnrlon

    if colorbar:  # Options for colorbar are hard-coded (to avoid a big mess).
        screen = kml.newscreenoverlay(name='ScreenOverlay')
        screen.icon.href = colorbar
        screen.overlayxy = OverlayXY(x=0, y=0,
                                     xunits=Units.fraction,
                                     yunits=Units.fraction)
        screen.screenxy = ScreenXY(x=0.015, y=0.075,
                                   xunits=Units.fraction,
                                   yunits=Units.fraction)
        screen.rotationXY = RotationXY(x=0.5, y=0.5,
                                       xunits=Units.fraction,
                                       yunits=Units.fraction)
        screen.size.x = 0
        screen.size.y = 0
        screen.size.xunits = Units.fraction
        screen.size.yunits = Units.fraction
        screen.visibility = 1

    kmzfile = kw.pop('kmzfile', 'overlay.kmz')
    kml.savekmz(kmzfile)
