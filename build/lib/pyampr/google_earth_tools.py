"""
google_earth_tools v1.1
Collator and Editor - Timothy J. Lang of NASA MSFC (timothy.j.lang@nasa.gov)
Library of tools for creating KML/KMZ files. Amalgamated from several sources.
Requires simplekml from https://code.google.com/p/simplekml/
Also requires os, pylab, numpy, and matplotlib. Tested with Python 2.7.
Last edited - 7/28/2014

Change Log
----------
v1.1 major changes:
1. Added timespan stamping capability via passing the times string list as
   an argument to make_kml().

"""

SIMPLEKML_FLAG = True
try:
    from simplekml import (Kml, OverlayXY, ScreenXY, Units, RotationXY,
                           AltitudeMode, Camera)
except ImportError:
    SIMPLEKML_FLAG = False

import numpy as np
import matplotlib.pyplot as plt

############################################################
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
        plt.ioff()  # Make `True` to prevent the KML components from popping-up.
    fig = plt.figure(figsize=figsize, frameon=False, dpi=pixels//10)
    # KML friendly image.  If using basemap try: `fix_aspect=False`.
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(llcrnrlon, urcrnrlon)
    ax.set_ylim(llcrnrlat, urcrnrlat)
    return fig, ax

############################################################
def make_kml(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat,
             figs, colorbar=None, times=None, **kw):
    """
    TODO: LatLon bbox, list of figs, optional colorbar figure,
    and several simplekml kw...
    TJL - Obtained from 
    http://ocefpaf.github.io/python4oceanographers/blog/2014/03/10/gearth/
    
    """
    if not SIMPLEKML_FLAG:
        print '***ERROR!***'
        print 'simplekml not installed, download from',\
              'https://pypi.python.org/pypi/simplekml/'
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

    for fig in figs: #NOTE: Overlays are limited to the same bbox.
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
            ground.timespan.end   = times[1]
        ground.icon.href = fig
        #TJL - swapping west/east to match LL vs. UR properly
        ground.latlonbox.west =  llcrnrlon
        ground.latlonbox.south = llcrnrlat
        ground.latlonbox.north = urcrnrlat
        ground.latlonbox.east =  urcrnrlon
                    
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

############################################################
#TJL - Changed arguments to allow user-defined filename
def kml_contour(lon, lat, z, full_path_and_name, *args, **kwargs):
    """
    Overlay Python contour lines onto Google Earth
    
    Syntax:
    KML_CONTOUR(LON,LAT,Z) writes contour lines in the same format as
    pylab's CONTOUR(LON,LAT,Z).
    KML_CONTOUR(LON,LAT,Z,N) draws N contour lines, overriding the
    automatic value
    KML_CONTOUR(LON,LAT,Z,V) draws LEN(V) contour lines at the values
    specified in the vector V
    KML_CONTOUR(LON,LAT,Z,[v v]) computes a single contour at the level v
    
    Input:
    LON: This can be either a matrix the same size as Z or a vector with
    length the same as the number of columns in Z.
    LAT: This can be either a matrix the same size as Z or a vector with
    length the same as the number of rows in Z.
    Z: Matrix of elevations
    
    Output:
    This function creates a kml file with the name and path
    it is passed
    
    Cameron Sparr - Nov. 14, 2011
    cameronsparr@gmail.com
    Minor edits by Timothy Lang (4/19/2014) - timothy.j.lang@nasa.gov
    """
    import pylab as pl
    import os
        
    # STYLE:
    # Edit the color, width, and labelsize as you see fit.
    # Colors: (NOTICE order is blue, green, red)
    #     white: FFFFFF      red: 0000FF         green: 00FF00
    #     blue: FF0000       cyan: FFFF00        magenta: FF00FF
    #     black: 000000
    # Width: int or float (default is 1)
    # Labelsize: int or float text line label size (default is 0.9).
    color = 'FFFFFF'
    width = 1
    labelsize = 0.9
        
    #TJL trying to reduce the number of labels
    labellimit = round(pl.sqrt(pl.sqrt(pl.size(z))))*10
    #TJL - Allowing user-defined filename
    kmlfile = open(full_path_and_name, 'w')
    fname = os.path.basename(full_path_and_name)
    c = pl.contour(lon, lat, z, *args, **kwargs)
    pl.close(pl.gcf())
    zs = c.levels
        
    kml_begin(kmlfile, color, width, labelsize, fname)
        
    for i in pl.arange(0,len(zs)):
        zz = zs[i]
        coll = c.collections[i]
        for p in coll.get_paths():
            s = len(p.vertices)
            if(s > 1):
                plon = p.vertices[0][0]
                plat = p.vertices[0][1]
                    
                if(s > labellimit):
                    place_label(kmlfile, plon, plat, zz)
                    
                # place contour line labels first:
                clon = p.vertices[:,0]
                clat = p.vertices[:,1]
                for j in pl.arange(0, pl.size(clon)):
                    if((j % labellimit*5) == 0):
                        place_label(kmlfile, clon[j], clat[j], zz)
                    
                #begin drawing a contour line:
                line_begin(kmlfile, zz)
                for k in pl.arange(0, pl.size(clon)):
                    lon = str(clon[k])
                    lat = str(clat[k])
                    kmlfile.write(lon+ ','+ lat+ ','+str(zz)+' ')
                line_end(kmlfile)
        
    kml_end(kmlfile)

############################################################
def kml_begin(kmlfile, color, width, labelsize, fname):
    kmlfile.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    kmlfile.write('<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n')
    kmlfile.write('<Document>\n')
    kmlfile.write('	<name>'+fname+'</name>\n')
    kmlfile.write('	<Style id="sn_noicon">\n')
    kmlfile.write('      <IconStyle>\n')
    kmlfile.write('          <Icon>\n')
    kmlfile.write('          </Icon>\n')
    kmlfile.write('      </IconStyle>\n')
    kmlfile.write('      <LabelStyle>\n')
    kmlfile.write('          <scale>')
    kmlfile.write(str(labelsize))
    kmlfile.write('</scale>\n')
    kmlfile.write('      </LabelStyle>\n')
    kmlfile.write('	</Style>\n')
    kmlfile.write(' <Style id="linestyle">\n')
    kmlfile.write('      <LineStyle>\n')
    kmlfile.write('          <color>#FF')
    kmlfile.write(color)
    kmlfile.write('</color>\n')
    kmlfile.write('          <width>')
    kmlfile.write(str(width))
    kmlfile.write('</width>\n')
    kmlfile.write('      </LineStyle>\n')
    kmlfile.write(' </Style>\n')

############################################################
def kml_end(kmlfile):
    kmlfile.write('</Document>\n')
    kmlfile.write('</kml>\n')
    kmlfile.close()

############################################################
def line_begin(kmlfile, zz):
    kmlfile.write('	<Placemark>\n')
    kmlfile.write('		<name>'+str(zz)+'</name>\n')
    kmlfile.write('		<styleUrl>#linestyle</styleUrl>\n')
    kmlfile.write('		<LineString>\n')
    kmlfile.write('			<tessellate>1</tessellate>\n')
    kmlfile.write('			<altitudeMode>clampToSeaFloor</altitudeMode>\n')
    kmlfile.write('			<gx:altitudeMode>clampToSeaFloor</gx:altitudeMode>\n')
    kmlfile.write('			<coordinates>\n')
    kmlfile.write('				')

############################################################
def line_end(kmlfile):
    kmlfile.write('\n')
    kmlfile.write('			</coordinates>\n')
    kmlfile.write('		</LineString>\n')
    kmlfile.write('	</Placemark>\n')

############################################################
def place_label(kmlfile, plon, plat, z):
    z = int(round(z))
    kmlfile.write('	<Placemark>\n')
    kmlfile.write('		<name>' + str(z) + '</name>\n')
    kmlfile.write('		<styleUrl>#sn_noicon</styleUrl>\n')
    kmlfile.write('		<Point>\n')
    kmlfile.write('			<altitudeMode>clampToSeaFloor</altitudeMode>\n')
    kmlfile.write('			<gx:altitudeMode>clampToSeaFloor</gx:altitudeMode>\n')
    kmlfile.write('			<coordinates>'+str(plon)+','+str(plat)+',0</coordinates>\n')
    kmlfile.write('		</Point>\n')
    kmlfile.write('	</Placemark>\n')
