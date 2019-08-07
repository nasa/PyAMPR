from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import time
import datetime
import calendar
import gzip
from netCDF4 import Dataset, num2date, date2num
import codecs
import cartopy
from .google_earth_tools import gearth_fig, make_kml
from .misc_tools import _FourPanelTrack
from .defaults import (
    DEFAULT_CLEVS, DEFAULT_GRID_DEL, DEFAULT_VAR, DEFAULT_CHAN_LIST,
    DEFAULT_SWATH_SIZE, DEFAULT_NAV_SIZE, DEFAULT_BAD_DATA,
    DEFAULT_SWATH_LEFT, DEFAULT_PROJECT_NAME, VERSION)
from .udf_cmap import amprTB_cmap, amprQC_cmap

CMAP_FLAG = True
FREQS = ['10', '19', '37', '85']
POLS = ['H', 'V']
CHANS = ['A', 'B']

#######################
# Main class definition
#######################


class AmprTb(object):

    def __init__(self, full_path_and_filename=None,
                 project=DEFAULT_PROJECT_NAME):
        """
If passed a filename, call the read_ampr_tb_level2b() method,
otherwise just instance the class with nothing
        """
        if full_path_and_filename is None:
            print('Class instantiated,',
                  'call read_ampr_tb_level2b() to populate')
        else:
            self.read_ampr_tb_level2b(full_path_and_filename, project=project)

    #########################################

    def read_ampr_tb_level2b(self, full_path_and_filename,
                             project=DEFAULT_PROJECT_NAME):
        """
Reads Level 2B AMPR data in text or netCDF files provided from
http://ghrc.msfc.nasa.gov. Tested and working with all AMPR
data from this site, as well as IPHEX. Depending on project,
some variables are unused or are duplicates. Most notably,
pre-MC3E there are no B channels.

Currently available field projects: CAMP2EX, ORACLES, OLYMPEX, IPHEX, MC3E,
TC4, TCSP, JAX90, COARE, CAMEX1, CAMEX2, CAMEX3, CAMEX4, TRMMLBA, KWAJEX,
TEFLUNA, FIRE3ACE, CAPE
If you read one project's data while mistakenly telling PyAMPR the data
are from a different project, then errors are likely.

Notable attributes in output data class
(Note - Order in documentation does not necessarily match order in data files)
---------------------------------------
nscans = Number of scans (depends on file)
swath_size = 50 (hard coded)
nav_size =   18 (hard coded)

shape = (nscans)
*****
Scan - Individual scan record number
       (usually thousands of scans per flight)
Year, Month, Day, Hour, Minute, Second, Day_of_Year, Second _of_Day -
    Scan timing info (UTC)
Icon - QC flag (not currently used by PyAMPR)

shape = (nscans, swath_size)
*****
TB10A, TB10B - 10 GHz brightness temperatures
    (A: Left V -> Right H, B: Left H -> Right V, units: K)
TB19A, TB19B - 19 GHz brightness temperatures (V->H, H->V, K)
TB37A, TB37B - 37 GHz brightness temperatures (V->H, H->V, K)
TB85A, TB85B - 85 GHz brightness temperatures (V->H, H->V, K)
Latitude, Longitude -  Geolocation for the AMPR beam (degrees)
Land_Fraction10 - Estimated fraction of land in 10/19 GHz pixel
Land_Fraction19 - Estimated fraction of land in 19 GHz pixel (IPHEX only)
Land_Fraction37 - Estimated fraction of land in 37 GHz pixel
Land_Fraction85 - Estimated fraction of land in 85 GHz pixel
    (0 = All Ocean, 1 = All Land)
Elevation - Topographic elevation (m MSL)

shape = (nscans, nav_size)
*****
Aircraft_Nav - Python dict of Aircraft navigation info: key (units)
GPS Latitude       (deg)
GPS Longitude      (deg)
GPS Altitude       (m MSL)
Pitch              (deg, + is nose up)
Roll               (deg, + is right wing down)
Yaw                (deg from N)
Heading            (deg from N)
Ground Speed       (m/s)
Air Speed          (m/s)
Static Pressure    (hPa)
Total Pressure     (hPa)
Total Temperature  (C)
Static Temperature (C)
Wind Speed         (m/s)
Wind Direction     (deg from N)
INS Latitude       (deg)
INS Longitude      (deg)
INS Altitude       (m MSL)

Version 1.4.0: Added support for Level 2B netCDF files. Starting with IPHEX,
processed AMPR instrument files will be provided in a netCDF-4 format.
        """
        _method_header_printout()
        print('read_ampr_tb_level2b(): Reading', full_path_and_filename)

        try:
            self._read_level2b_netcdf(full_path_and_filename, project=project)
        except IOError:
            print('Not netCDF file, trying ASCII read ...')
            self._read_level2b_ascii(full_path_and_filename, project=project)
        _method_footer_printout()

    #########################################

    def help(self):
        """AmprTb.help() = help(AmprTb)"""
        help(self)

    #########################################

    def plot_ampr_track(
            self, var=DEFAULT_VAR, latrange=None, lonrange=None,
            parallels=DEFAULT_GRID_DEL, meridians=DEFAULT_GRID_DEL,
            title=None, clevs=DEFAULT_CLEVS, cmap=None,
            save=None, show_track=False, maneuver=True,
            scanrange=None, show_grid=True, equator=False,
            timerange=None, return_flag=False, show_qc=False,
            resolution='50m', projection=cartopy.crs.PlateCarree(),
            ax=None, fig=None, title_flag=True,
            colorbar_label=True, verbose=False,
            show_borders=True, colorbar_flag=True,
            wmts_layer=None):

        """
This method plots geolocated AMPR data, along with the Aircraft track if
requested. matplotlib.pyplot.pcolormesh() on a Cartopy basemap is the workhorse
plotting routine.

var = String with channel number and letter (e.g., 10A for 10 GHz (A) channel
latrange = List with lat range defined. Order and size (>= 2) is irrelevant
           as max and min are retrieved.
lonrange = List with lon range defined. Order and size (>= 2) is irrelevant
           as max and min are retrieved.
parallels = Scalar spacing (deg) for parallels (i.e., constant latitude)
meridians = Scalar spacing (deg) for meridians (i.e., constant longitude)
ptitle = Plot title as string
clevs = List with contour levels. Only max and min values are used.
cmap = Colormap desired. See, e.g.,
       https://matplotlib.org/3.1.0/gallery/color/colormap_reference.html
save = Filename+path as string to save plot to. Type determined from suffix.
       Careful - .ps/.eps/.pdf files can get huge!
show_track = Boolean to plot Aircraft track along with AMPR data.
             Plotted in black with white highlights for significant maneuvers
             (abs(Aircraft_Nav['Roll']) > 5).
scanrange = List of scan numbers (from AmprTb.Scan) to plot.
            Only max/min are used.
show_grid = Set to False to turn off gridlines
equator = Boolean to consider 0s/-1s in Latitude/Longitude as good geolocations
          (e.g., flight crosses Equator or Prime Meridian).
          Default is bad geolocations.
maneuver = Set to False to suppress the plotting of swaths during significant
           aircraft maneuvers.
timerange = Time range to plot. Overrides scanrange if both are set.
            Format: timerange = ['hh:mm:ss', 'HH:MM:SS']
return_flag = Set to True to return figure, plot axes, etc. Order of items
              returned is fig (figure instance), ax (main axis instance),
              cbar (colorbar instance).
show_qc = Set to True to show QC flags instead of TB variables.
resolution = Resolution of Cartopy map ('10m', '50m', '110m')
projection = Cartopy map projection to use.
ax, fig = Matplotlib Axes and Figure objects. Either both must be None
          or both must be valid objects for the plot to work right.
title_flag = Set to False to suppress a title
colorbar_label = Set to False to suppress the colorbar label
verbose = Set to True for text output
show_borders = False removes coastlines and state/national borders
wmts_layer = Use this to add a WMTS map layer from
             https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi, like
             'ASTER_GDEM_Color_Shaded_Relief'
        """

        # plt.close()  # mpl seems buggy if you don't clean up old windows
        if verbose:
            _method_header_printout()
            print('plot_ampr_track():')

        # 10 GHz (A) channel plotted by default if mistake made
        if not isinstance(var, str):
            var = DEFAULT_VAR

        # Check to make sure data exist!
        if not hasattr(self, 'TB'+var.upper()):
            self._missing_channel_printout()
            if verbose:
                _method_footer_printout()
            return

        # Check that QC data exist
        show_qc, clevs = self._check_qc(show_qc, clevs)

        if self.Year[0] < 2011 and verbose is True:
            print('Warning: Older projects commonly had bad or missing',
                  'geolocation data.')
            print('If there are plotting problems, try a strip chart with',
                  'plot_ampr_channels(),')
            print('or try adjusting scanrange, lonrange, or latrange.')

        # Adjustable scan range limits
        # Fairly robust - will go down to a width of 10 scans or so before
        # plotting artifacts begin to occur.
        ind1, ind2 = self._get_scan_indices(scanrange, timerange, verbose)
        plon, plat, zdata = self._get_data_subsection(
            var, ind1, ind2, maneuver, show_qc, verbose)
        plon, plat, zdata = self._filter_bad_geolocations(
            plon, plat, zdata, equator, verbose)
        enough_data = self._check_for_enough_data_to_plot(plon, plat)
        if not enough_data:
            return

        latrange, lonrange = self._get_latrange_lonrange(
            plat, plon, latrange, lonrange)
        self._check_aspect_ratio(latrange, lonrange, verbose)
        ax, fig = parse_ax_fig(ax, fig, projection=projection)

        llcrnrlon = np.min(lonrange)
        urcrnrlon = np.max(lonrange)
        llcrnrlat = np.min(latrange)
        urcrnrlat = np.max(latrange)
        ax.set_extent([llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat])

        if wmts_layer is not None:
            url = 'https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi'
            ax.add_wmts(url, wmts_layer)

        if show_borders:
            ax.coastlines(resolution=resolution)
            countries = cartopy.feature.NaturalEarthFeature(
                category='cultural', name='admin_0_boundary_lines_land',
                scale=resolution, facecolor='none')
            ax.add_feature(countries, edgecolor='black')
            states_provinces = cartopy.feature.NaturalEarthFeature(
                category='cultural', name='admin_1_states_provinces_lines',
                scale=resolution, facecolor='none')
            ax.add_feature(states_provinces, edgecolor='black')

        if show_grid:
            gl = ax.gridlines(draw_labels=True, linestyle='--')
            gl.xlabels_top = False
            gl.ylabels_right = False
            vmeridians = np.arange(-180, 180, meridians)
            gl.xlocator = mticker.FixedLocator(vmeridians)
            vparallels = np.arange(-90, 90, parallels)
            gl.ylocator = mticker.FixedLocator(vparallels)
            del vparallels, vmeridians

        # Draw filled contours
        cmap = self._get_colormap(cmap, CMAP_FLAG, show_qc)
        cs = ax.pcolormesh(plon, plat, zdata, vmin=np.min(clevs),
                           vmax=np.max(clevs), cmap=cmap, zorder=2,
                           transform=projection)

        # Add Aircraft track
        if show_track:
            aplon = self.Aircraft_Nav['GPS Longitude'][ind1:ind2]
            aplat = self.Aircraft_Nav['GPS Latitude'][ind1:ind2]
            # Black dots during normal flight
            ax.plot(aplon, aplat, 'k.', transform=projection)
            indices = np.where(
                np.abs(self.Aircraft_Nav['Roll'][ind1:ind2]) >= 5)
            # White dots during maneuvers
            ax.plot(aplon[indices[0]], aplat[indices[0]], 'w.',
                    transform=projection)

        # Plot title & display
        if title_flag:
            if title is None:
                title = str(self._get_ampr_title(var)) + '\n' + \
                    str(self._get_date_string(ind1)) + \
                    str(', ') + str(self.Time_String[ind1]) + str('-') + \
                    str(self.Time_String[ind2-1]) + str(' UTC')
            ax.set_title(title)

        # Add colorbar
        # Colorbar independent of Basemap
        # cax = fig.add_axes([0.20, 0.07, 0.60, 0.02])
        # cbar = plt.colorbar(cs, cax=cax, orientation='horizontal',
        #                     extend='both')
        if colorbar_flag:
            cbar = plt.colorbar(cs, ax=ax, orientation='vertical',
                                extend='both', shrink=0.75)
            cbar = self._adjust_colorbar(cbar, show_qc, colorbar_label)

        # Save image to file
        if save is not None:
            plt.savefig(save)

        if verbose:
            _method_footer_printout()
        if return_flag:
            return fig, ax, cbar  # cax, cbar

    #########################################

    def plot_ampr_track_4panel(self, **kwargs):

        """
This method plots 4 panels of geolocated AMPR data, along with the aircraft
track if requested. matplotlib.pyplot.pcolormesh() on a Cartopy basemap is the
workhorse plotting routine.

chan = String with channel letter. Default is 'A'
latrange = List with lat range defined. Order and size (>= 2) is irrelevant
           as max and min are retrieved.
lonrange = List with lon range defined. Order and size (>= 2) is irrelevant
           as max and min are retrieved.
parallels = Scalar spacing (deg) for parallels (i.e., constant latitude)
meridians = Scalar spacing (deg) for meridians (i.e., constant longitude)
ptitle = Plot title as string
clevs = List with contour levels. Only max and min values are used.
cmap = Colormap desired. See, e.g.,
       https://matplotlib.org/3.1.0/gallery/color/colormap_reference.html
save = Filename+path as string to save plot to. Type determined from suffix.
       Careful - .ps/.eps/.pdf files can get huge!
show_track = Boolean to plot Aircraft track along with AMPR data.
             Plotted in black with white highlights for significant maneuvers
             (abs(Aircraft_Nav['Roll']) > 5).
scanrange = List of scan numbers (from AmprTb.Scan) to plot.
            Only max/min are used.
show_grid = Set to False to turn off gridlines
equator = Boolean to consider 0s/-1s in Latitude/Longitude as good geolocations
          (e.g., flight crosses Equator or Prime Meridian).
          Default is bad geolocations.
maneuver = Set to False to suppress the plotting of swaths during significant
           aircraft maneuvers.
timerange = Time range to plot. Overrides scanrange if both are set.
            Format: timerange = ['hh:mm:ss', 'HH:MM:SS']
return_flag = Set to True to return the _FourPanelTrack instance to make
              additional figure modifications
show_qc = Set to True to show QC flags instead of TB variables.
resolution = Resolution of Cartopy map ('10m', '50m', '110m')
projection = Cartopy map projection to use.
title_flag = Set to False to suppress a title
colorbar_label = Set to False to suppress the colorbar label
verbose = Set to True for text output
wmts_layer = Use this to add a WMTS map layer from
             https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi, like
             'ASTER_GDEM_Color_Shaded_Relief'
        """
        fourpan = _FourPanelTrack(self, **kwargs)
        if 'return_flag' in kwargs.keys():
            if kwargs['return_flag']:
                return fourpan

    #########################################

    def plot_ampr_channels(self, scanrange=None, cmap=None,
                           clevs=DEFAULT_CLEVS, show_qc=False,
                           save=None, show_pol=False, timerange=None):
        """
This method plots a strip chart akin to those seen here:
ftp://gpm.nsstc.nasa.gov/gpm_validation/mc3e/ampr/browse/
matplotlib.pyplot.pcolormesh() is the workhorse plotting routine

clevs = List with contour levels. Only max and min values are used.
cmap = Colormap desired.
       See http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
       and dir(cm) for more
save = Filename+path as string to save plot to. Type determined from suffix.
       Careful - .ps/.eps/.pdf files can get huge!
scanrange = List of scan numbers (from AmprTb.Scan) to plot.
            Only max/min are used.
show_pol = Set to True to show deconvolved H & V polarizations. Will call
           calc_polarization() beforehand if these channels are missing.
show_qc = Set to True to show QC flags instead of TB variables.
          Overrides show_pol.
timerange = Time range to plot. Overrides scanrange if both are set.
            Format: timerange = ['hh:mm:ss', 'HH:MM:SS']
        """

        plt.close()  # mpl seems buggy if you don't clean up old windows
        _method_header_printout()
        print('plot_ampr_channels():')

        # Check that QC data exist
        show_qc, clevs = self._check_qc(show_qc, clevs)

        tb_list = self._get_list_of_channels_to_plot(show_pol, show_qc)
        if show_pol is True and show_qc is False:
            pol_flag = self._check_for_pol_data()
            if not pol_flag:
                _method_footer_printout()
                return

        # Set up plot - 10 rows, 1 column, 8 channels for 8 rows,
        # one row for colorbar, and one row for Aircraft_Nav info
        fig, axes = plt.subplots(nrows=10, ncols=1, sharex=False)
        fig.set_size_inches(11, 8.5)

        # Adjustable scan range limits
        # Fairly robust - will go down to a width of 10 scans before plotting
        # artifacts begin to occur.
        ind1, ind2 = self._get_scan_indices(scanrange, timerange)
        xran = [self.Scan[ind1], self.Scan[ind2-1]]
        plt.xlim(xran)

        # The following will put the title above the x-axis tick labels
        ptitle = 'AMPR ' + self._get_date_string(ind1)
        plt.text(0.5, 1.60, ptitle, horizontalalignment='center', fontsize=14,
                 transform=axes[0].transAxes)

        colormap = self._get_colormap(cmap, CMAP_FLAG, show_qc)

        # Loop and plot available channels
        itest = 0
        for index, chan in enumerate(tb_list):
            if hasattr(self, 'TB'+chan):
                itest = itest + 1
                if show_qc:
                    var = np.transpose(getattr(self, 'qctb'+chan.lower()))
                else:
                    var = np.transpose(getattr(self, 'TB'+chan))
                # AMPR data arranged L-to-R in PyAMPR (index 0 to index 49
                # in a scan),
                # so need to reverse order to have L on top in strip charts.
                im = axes[index].pcolormesh(
                    self.Scan, (self.swath_size-1)-np.arange(self.swath_size),
                    var, vmin=np.min(clevs), vmax=np.max(clevs), cmap=colormap)
                axes[index].set_xlim(xran)
                axes[index].set_ylim([0, self.swath_size-1])
                axes[index].set_ylabel(chan)
                axes[index].yaxis.set_ticks([2, self.swath_size-5])
                axes[index].tick_params(axis='y', labelsize=7)
                axes[index].tick_params(axis='x', labelbottom=False, top=True,
                                        length=2.5)
                axes[index].tick_params(axis='y', which='both', left=False,
                                        right=False)
                # Shows how V and H vary with beam position
                if show_pol is False or show_qc is True:
                    if index % 2 == 0:
                        axes[index].set_yticklabels(['H', 'V'])
                    else:
                        axes[index].set_yticklabels(['V', 'H'])
                else:
                    axes[index].set_yticklabels(['R', 'L'])

                # Get scan timing and plot it as tick labels
                if index == 0:
                    axes[index].tick_params(
                        axis='x', labelsize=7, labeltop=True,
                        direction='out', labelbottom=False, pad=0)
                    locs, labels = plt.xticks()
                    new_labels = []
                    for t in locs:
                        indt = np.where(self.Scan == int(t))
                        tmpst = np.array(self.Time_String)[indt[0]]
                        if not isinstance(tmpst, str):
                            if len(tmpst) > 0:
                                tmpst = tmpst[0]
                            else:
                                tmpst = ''
                        new_labels.append(''.join(str(tmpst)))
                    axes[index].set_xticklabels(new_labels)

        # Check if missing too much data!
        if itest == 0:
            print('No data to plot, try reading in a file')
            _method_footer_printout()
            return

        # Add colorbar
        axes[8].axis('off')
        cax = fig.add_axes([0.126, 0.245, 0.775, 0.01])
        cbar = plt.colorbar(im, cax=cax, orientation='horizontal',
                            extend='both')
        cbar = self._adjust_colorbar(cbar, show_qc)
        cbar.ax.tick_params(labelsize=7)

        # Add Aircraft roll angle time series
        if hasattr(self, 'Aircraft_Nav'):
            axes[9].plot(self.Scan, self.Aircraft_Nav['Roll'], 'b-')
            axes[9].plot(self.Scan, 0.0*self.Aircraft_Nav['Roll'], 'k:')
            axes[9].plot(self.Scan, 5.0+0.0*self.Aircraft_Nav['Roll'], 'k:')
            axes[9].plot(self.Scan, -5.0+0.0*self.Aircraft_Nav['Roll'], 'k:')
            axes[9].set_xlabel('Scan Number or Time (UTC)', fontsize=10)
            axes[9].tick_params(axis='x', labelsize=9, top=False,
                                direction='out')
            axes[9].set_ylabel('Aircraft Roll (deg)', fontsize=10)
            axes[9].tick_params(axis='y', labelsize=7)
            axes[9].set_xlim(xran)
            axes[9].set_ylim([-10, 10])
            axes[9].yaxis.set_ticks([-10, 0, 10])

        # Save the plot and clean up
        if save is not None:
            plt.savefig(save)
        _method_footer_printout()

    #########################################

    def calc_polarization(self, simple=False, force_match=True,
                          chan_list=DEFAULT_CHAN_LIST):

        """
*** THIS METHOD IS EXPERIMENTAL ***

This method calculates H & V given the mixed-pol A & B channels.
Solves Equation 1 in Vivekanandan et al. (1993) for Tb in H and V.
Where A or B are not good will be populated with AmprTb.bad_data.
If successful, TB10H, TB10V, TB19H, TB19V, TB37H, TB37V, TB85H, TB85V
will now be attributes of the AmprTb instance. Missing channels will
not be processed. Calculation performed via methodology of Brent Roberts.

Creates attributes called TB##_offset, where ## is channel
frequency (in GHz). This is Channel A TB - Channel B TB (in K)
at nadir scan angle.

Author: Brent Roberts w/ adjustments by Timothy Lang

simple = Boolean flag to swap between simple linear substitution or
         constrained linear inversion.
         (Default = constrained linear inversion)
force_match = Boolean flag to force A & B channels to match at nadir
              (Default = match forced)
chan_list = List of strings to enable individual freqs to be deconvolved
            using different methodologies (Default = All 4 freqs).

        """
        begin_time = time.time()
        _method_header_printout()
        print('calc_polarization():')

        if self.Year[0] < 2011:
            print('Pre-2011, AMPR only had one channel per frequency.')
            print('Thus, PyAMPR cannot deconvolve polarization for')
            print('this project\'s data. Sorry!')
            _method_footer_printout()
            return

        for chan in chan_list:
            if hasattr(self, 'TB'+chan+'A') and hasattr(self, 'TB'+chan+'B'):

                print('Calculating for', chan, 'GHz channel')
                # Use dummy variables and setattr to get the right-sized arrays
                T1 = 1.0 * getattr(self, 'TB'+chan+'A')
                T2 = 1.0 * getattr(self, 'TB'+chan+'B')

                # Get angular argument and convert to radians.
                # angle = np.deg2rad(data.Incidence_Angle - 45.0);
                angle = np.deg2rad(
                    np.linspace(self.swath_left,
                                -1.0*self.swath_left, self.swath_size) - 45.0)

                # There appears to be an offset between the mixed-pol
                # brightness temperatures at 0-deg incidence angle.
                # Theoretically these should be the same.
                T2 = self._compute_nadir_offset_and_apply_if_desired(
                    angle, T1, T2, force_match, chan)

                if simple:
                    tbv, tbh = \
                        self._solve_using_simple_linear_substitution(
                            angle, T1, T2)

                # Solve using Tikhonov regularization
                else:
                    tbv, tbh = \
                        self._solve_using_constrained_linear_inversion(
                            angle, T1, T2)

                # Finalize V & H attributes and clean up
                setattr(self, 'TB'+chan+'V', tbv)
                setattr(self, 'TB'+chan+'H', tbh)

            else:
                print('TB' + chan, 'does not have both A and B channels,',
                      'read in a file to obtain.')
                print('Scene H & V not produced for this channel')

        print(time.time() - begin_time, 'seconds to calculate H & V')
        print('If successful, following attributes are now available:')
        for chan in chan_list:
            print('TB'+chan+'H', 'TB'+chan+'V', end=' ')
        print('')
        _method_footer_printout()

    #########################################

    def write_ampr_kmz(self, var=DEFAULT_VAR, latrange=None, lonrange=None,
                       clevs=DEFAULT_CLEVS, cmap=None, timerange=None,
                       file_path=None, file_name=None, scanrange=None,
                       show_legend=True, equator=False, maneuver=True,
                       show_qc=False):
        """
        This method plots geolocated AMPR data as a filled color Google Earth
        kmz. Qualitatively similar plot to plot_ampr_track() but for Google
        Earth.
        Will produce overlay.png and, if a legend is created,
        legend.png as temporary image files in the current working
        directory.

        var = AMPR channel to plot (Default = DEFAULT_VAR)
        scanrange = List of scan numbers (from AmprTb.Scan) to plot.
                    Only max/min are used.
        file_path = Desired path to kmz file (Default = '')
        file_name = Desired name for kmz file (Default = YYYYMMDD_TB###.kmz,
                    ### = channel)
        latrange = List with lat range defined. Order and size (>= 2)
                   is irrelevant as max and min are retrieved
        lonrange = List with lon range defined. Order and size (>= 2)
                   is irrelevant as max and min are retrieved
        clevs = List with contour levels. Only max and min values are used.
                (Default = DEFAULT_CLEVS)
        cmap = Colormap desired.
               See http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
               and dir(cm) for more
        equator = Boolean to consider 0s/-1s in Latitude/Longitude as
                  good geolocations
                  (e.g., flight crosses Equator or Prime Meridian).
                  Default is bad geolocations.
        show_legend = Set to False to suppress the color bar
        maneuver = Set to False to suppress plotting of swaths during
                   significant aircraft maneuvers.
        timerange = Time range to plot. Overrides scanrange if both are set.
                    Format: timerange = ['hh:mm:ss', 'HH:MM:SS']
        show_qc = Set to True to show QC flags instead of TB variables.
        """
        plt.close()  # mpl seems buggy if multiple windows are left open
        _method_header_printout()
        print('write_ampr_kmz():')

        # 10 GHz (A) channel plotted by default
        if var is None or isinstance(var, str) is False:
            var = DEFAULT_VAR

        # Check to make sure data exist!
        if not hasattr(self, 'TB'+var.upper()):
            self._missing_channel_printout()
            _method_footer_printout()
            return

        # Check that QC data exist
        show_qc, clevs = self._check_qc(show_qc, clevs)

        # Adjustable scan range limits
        # Fairly robust - will go down to a width of 25 scans before plotting
        # artifacts begin to occur.
        ind1, ind2 = self._get_scan_indices(scanrange, timerange)
        plon, plat, zdata = self._get_data_subsection(
            var, ind1, ind2, maneuver, show_qc)
        plon, plat, zdata = self._filter_bad_geolocations(plon, plat,
                                                          zdata, equator)
        enough_data = self._check_for_enough_data_to_plot(plon, plat)
        if not enough_data:
            return

        latrange, lonrange = self._get_latrange_lonrange(
            plat, plon, latrange, lonrange)
        times = self._get_timestamps_for_gearth(ind1, ind2)

        # Set file info
        if file_path is None:
            file_path = ''
        if file_name is None:
            file_name = self._get_gearth_file_name(var, ind1, '.kmz')

        # Google Earth image production
        fig, ax = gearth_fig(np.min(lonrange), np.min(latrange),
                             np.max(lonrange), np.max(latrange))
        cmap = self._get_colormap(cmap, CMAP_FLAG, show_qc)
        cs = ax.pcolormesh(plon, plat, zdata,
                           vmin=np.min(clevs), vmax=np.max(clevs), cmap=cmap)
        ax.set_axis_off()
        fig.savefig('overlay.png', transparent=True, format='png')

        # Now we convert to KMZ
        if show_legend is True:
            fig = plt.figure(figsize=(1.0, 4.0), facecolor=None, frameon=False)
            ax = fig.add_axes([0.0, 0.05, 0.2, 0.9])
            cb = fig.colorbar(cs, cax=ax)
            cbytick_obj = plt.getp(cb.ax.axes, 'yticklabels')
            plt.setp(cbytick_obj, color='w', weight='bold')
            ptitle = self._get_ampr_title(var)
            if show_qc:
                cb.set_ticks([1, 2, 3, 4, 5])
                clabel = 'QC Flag'
            else:
                clabel = 'TB [K]'
            cb.set_label(ptitle+clabel, rotation=-90, color='w',
                         labelpad=20, weight='bold')
            fig.savefig('legend.png', transparent=True, format='png')
            make_kml(np.min(lonrange), np.min(latrange), np.max(lonrange),
                     np.max(latrange), figs=['overlay.png'],
                     kmzfile=str(file_path+file_name), colorbar='legend.png',
                     times=times)
        else:
            make_kml(np.min(lonrange), np.min(latrange), np.max(lonrange),
                     np.max(latrange), figs=['overlay.png'],
                     kmzfile=str(file_path+file_name), times=times)

        print('Google Earth image hopefully written to:',
              file_path + file_name)
        _method_footer_printout()

    #################################################################
    # Internal methods below here. Average user can stop reading now.
    #################################################################

    def _read_level2b_netcdf(self, inputFile, project=DEFAULT_PROJECT_NAME):
        """
        Internal method to read L2 netCDF-format AMPR data files.
        Accounts for whether the data are from older netCDF formats
        (2014-2017) or the new CF-compliant format (2017-).

        Parameters
        ----------
        inputFile : str
            Name of input AMPR data file.
        """
        # Open the data
        level2b = Dataset(inputFile, mode="r")

        # Set bad data and nav_size
        self.bad_data = DEFAULT_BAD_DATA
        self.nav_size = DEFAULT_NAV_SIZE

        # Assigning project name (self.Project) based on user input
        self._assign_project_name(project)
        self.keylist = list(level2b.variables.keys())

        if 'TB' in self.keylist:
            self.CF_flag = True
        else:
            self.CF_flag = False

        # Check for navigation
        if 'lat' in self.keylist or 'Lat' in self.keylist:
            print('Found Navigation Data!')
            self.hasNav = True
        else:
            self.hasNav = False
            print('No navigation data, track plots unavailable ...')

        self._initialize_vars_netcdf(level2b)
        self._populate_time_vars_netcdf()
        self._fill_epoch_time()  # defines and populates self.Epoch_Time

        if self.hasNav:  # Add the geospatial information
            if self.CF_flag:
                self.Latitude = level2b.variables['Lat'][:, :]
                self.Longitude = level2b.variables['Lon'][:, :]
                if 'IncidenceAngle' in self.keylist:
                    self.Incidence_Angle = level2b.variables[
                        'IncidenceAngle'][:, :]
                if 'RelativeAzimuth' in self.keylist:
                    self.Relative_Azimuth = level2b.variables[
                        'RelativeAzimuth'][:, :]
            else:
                self.Latitude = level2b.variables['lat'][:, :]
                self.Longitude = level2b.variables['lon'][:, :]
                self.Incidence_Angle = level2b.variables[
                    'incidence_angle'][:, :]
                self.Relative_Azimuth = level2b.variables[
                    'relative_azimuth'][:, :]

        # Add the Calibrated TBs
        self._assign_tbs_netcdf(level2b)

        # Other Variables -- Set to bad data for now.
        self._set_old_vars_to_bad_netcdf(level2b)

        # Aircraft navigation info
        self._consider_aircraft_nav_netcdf(level2b)

        # QC flags and FOV (IPHEx V2 Release)
        self._consider_qc_flags_netcdf(level2b)
        self._consider_land_frac_netcdf(level2b)
        self._remove_nan_inf()

    #########################################

    def _initialize_vars_netcdf(self, level2b):
        """
        Helper method to _read_level2b_netcdf()
        Initializes miscellaneous variables based on level2b netCDF input.
        Accounts for whether the data are from older netCDF formats
        (2014-2017) or the new CF-compliant format (2017-).

        Parameters
        ----------
        level2b : netCDF4.Dataset object
            The Dataset object read in from the AMPR data file.
        """
        # Handle the times, convert to a datetime object
        # Add the scan and swath information
        if self.CF_flag:
            self.netcdfTimes = level2b.variables['Time']
            self.dateTimes = num2date(self.netcdfTimes[:],
                                      self.netcdfTimes.units)
            self.nscans = len(level2b.dimensions['AlongTrackDim'])
            self.ncross = len(level2b.dimensions['CrossTrackDim'])
            self.swath_size = self.ncross
            self.Scan = np.arange(self.nscans, dtype='int') + 1
            self.Scan_Position = np.arange(self.ncross, dtype='int') + 1
            if 'ScanAngle' in self.keylist:
                self.swath_angle = level2b.variables['ScanAngle'][:]
                self.swath_left = self.swath_angle[0]
            else:
                self.swath_left = DEFAULT_SWATH_LEFT
                self.swath_angle = self.swath_left - \
                    (2.0 * self.swath_left / float(self.swath_size - 1.0)) * \
                    np.arange(self.swath_size)
        else:
            self.netcdfTimes = level2b.variables['time']
            self.dateTimes = num2date(self.netcdfTimes[:],
                                      self.netcdfTimes.units)
            self.nscans = len(level2b.dimensions['scan_number'])
            self.ncross = len(level2b.dimensions['scan_position'])
            self.swath_size = self.ncross
            self.Scan = level2b.variables['scan_number'][:]
            self.Scan_Position = level2b.variables['scan_position'][:]
            self.swath_angle = level2b.variables['scan_angle'][:]
            self.swath_left = self.swath_angle[0]

        # Common to both formats
        self._initialize_time_fields()

    #########################################

    def _assign_tbs_netcdf(self, level2b):
        """
        Helper method to _read_level2b_netcdf()
        Assigns all brightness temperature variables to their appropriate
        attributes. Accounts for whether the data are from older netCDF
        formats (2014-2017) or the new CF-compliant format (2017-).

        Parameters
        ----------
        level2b : netCDF4.Dataset object
            The Dataset object read in from the AMPR data file.
        """
        if self.CF_flag:
            tbdata = np.array(level2b.variables['TB'][:, :, :, :])
            nchan = np.shape(tbdata)[0]
            self.hasDeconvolvedHV = False
            if nchan == 1:  # CF files may be from pre-dual-pol era
                chanlist = [CHANS[0]]
            elif nchan == 2:
                chanlist = CHANS
            else:
                chanlist = np.concatenate([CHANS, POLS])
                self.hasDeconvolvedHV = True
            for j, pol in enumerate(chanlist):
                for i, freq in enumerate(FREQS):
                    setattr(self, 'TB' + freq + pol, tbdata[j, i, :, :])
        else:
            if 'tbs_10h' in self.keylist:
                self.hasDeconvolvedHV = True
                chanlist = np.concatenate([CHANS, POLS])
            else:
                self.hasDeconvolvedHV = False
                chanlist = CHANS  # Older netCDF always dual-pol
            for freq in FREQS:
                for pol in chanlist:
                    setattr(
                        self, 'TB' + freq + pol,
                        level2b.variables['tbs_'+freq+pol.lower()][:, :])
        if self.hasDeconvolvedHV:
            print('Found deconvolved H & V data!')

    #########################################

    def _set_old_vars_to_bad_netcdf(self, level2b):
        """
        Helper method to _read_level2b_netcdf()
        Checks for presence of older variables in the netCDF file, and
        if they are not present sets the corresponding attributes to bad.
        Accounts for whether the data are from older netCDF
        formats (2014-2017) or the new CF-compliant format (2017-).

        Parameters
        ----------
        level2b : netCDF4.Dataset object
            The Dataset object read in from the AMPR data file.
        """
        # Icon (whatever that is ...)
        if 'Icon' in self.keylist:
            self.Icon = level2b.variables['Icon'][:]
        else:
            self.Icon = self.bad_data * np.ones(self.nscans, dtype=np.int32)
        # Noise
        if 'Noise' in self.keylist:
            for j, freq in enumerate(FREQS):
                setattr(self, 'Noise' + freq, level2b.variables['Noise'][j, :])
        else:
            for freq in FREQS:
                setattr(self, 'Noise' + freq,
                        self.bad_data * np.ones(self.nscans, dtype='float'))
        # Land Fraction
        for freq in FREQS:
            setattr(self, 'Land_Fraction' + freq,
                    self.bad_data * np.ones((self.nscans, self.ncross),
                                            dtype='float'))

    #########################################

    def _consider_aircraft_nav_netcdf(self, level2b):
        """
        Helper method to _read_level2b_netcdf()
        Assigns all aircraft navigation variables to their appropriate
        attributes. Accounts for whether the data are from older netCDF
        formats (2014-2017) or the new CF-compliant format (2017-).

        Parameters
        ----------
        level2b : netCDF4.Dataset object
            The Dataset object read in from the AMPR data file.
        """
        self._initialize_aircraft_dict()
        for var in self.Aircraft_varlist:
            self.Aircraft_Nav[var][:] = self.bad_data
        # Now, return a structured array of Aircraft_Nav parameters.
        if self.hasNav:
            if self.CF_flag:
                self.Aircraft_Nav['GPS Latitude'] = \
                    level2b.variables['GPSLatitude'][:]
                self.Aircraft_Nav['GPS Longitude'] = \
                    level2b.variables['GPSLongitude'][:]
                self.Aircraft_Nav['GPS Altitude'] = \
                    level2b.variables['GPSAltitude'][:]
                self.Aircraft_Nav['Pitch'] = level2b.variables['Pitch'][:]
                self.Aircraft_Nav['Roll'] = level2b.variables['Roll'][:]
                self.Aircraft_Nav['Yaw'] = level2b.variables['Yaw'][:]
                self.Aircraft_Nav['Heading'] = level2b.variables['Head'][:]
                self.Aircraft_Nav['Ground Speed'] = \
                    level2b.variables['GroundSpeed'][:]
                self.Aircraft_Nav['Air Speed'] = \
                    level2b.variables['AirSpeed'][:]
                if 'WindSpeed' in self.keylist:
                    self.Aircraft_Nav['Wind Speed'] = \
                        level2b.variables['WindSpeed'][:]
                if 'WindDirection' in self.keylist:
                    self.Aircraft_Nav['Wind Direction'] = \
                        level2b.variables['WindDirection'][:]
                if 'Pressure' in self.keylist:
                    self.Aircraft_Nav['Static Pressure'] = \
                        level2b.variables['Pressure'][:]
                if 'Temperature' in self.keylist:
                    self.Aircraft_Nav['Total Temperature'] = \
                        level2b.variables['Temperature'][:]
                # INS Lat/Lon ignored in CF datasets
            else:
                self.Aircraft_Nav['GPS Latitude'] = \
                    level2b.variables['gLat'][:]
                self.Aircraft_Nav['GPS Longitude'] = \
                    level2b.variables['gLon'][:]
                self.Aircraft_Nav['GPS Altitude'] = \
                    level2b.variables['gAlt'][:]
                self.Aircraft_Nav['Pitch'] = level2b.variables['pitch'][:]
                self.Aircraft_Nav['Roll'] = level2b.variables['roll'][:]
                self.Aircraft_Nav['Yaw'] = level2b.variables['track_angle'][:]
                self.Aircraft_Nav['Heading'] = level2b.variables['heading'][:]
                self.Aircraft_Nav['Ground Speed'] = \
                    level2b.variables['groundSpeed'][:]
                self.Aircraft_Nav['Air Speed'] = \
                    level2b.variables['airSpeed'][:]
                self.Aircraft_Nav['Static Pressure'] = \
                    level2b.variables['staticPressure'][:]
                self.Aircraft_Nav['Total Temperature'] = \
                    level2b.variables['totalTemp'][:]
                self.Aircraft_Nav['Wind Speed'] = \
                    level2b.variables['iWindSpeed'][:]
                self.Aircraft_Nav['Wind Direction'] = \
                    level2b.variables['iWindDir'][:]
                self.Aircraft_Nav['INS Latitude'] = \
                    level2b.variables['iLat'][:]
                self.Aircraft_Nav['INS Longitude'] = \
                    level2b.variables['iLon'][:]

    #########################################

    def _consider_qc_flags_netcdf(self, level2b):
        """
        Helper method to _read_level2b_netcdf()
        Assigns QC variables to their appropriate attributes.
        Accounts for whether the data are from older netCDF
        formats (2014-2017) or the new CF-compliant format (2017-).

        Parameters
        ----------
        level2b : netCDF4.Dataset object
            The Dataset object read in from the AMPR data file.
        """
        if self.CF_flag:
            if 'QC' in self.keylist:
                qcdata = np.array(level2b.variables['QC'])
                if np.ndim(qcdata) != 1:  # Filter out pre-dual-pol QC
                    self.qcIncidence = \
                        level2b.variables['IncidenceAngleQC'][:, :]
                    for j, chan in enumerate(CHANS):
                        for i, freq in enumerate(FREQS):
                            setattr(self, 'qctb' + freq + chan.lower(),
                                    qcdata[j, i, :, :])
        else:
            if 'qctb10a' in self.keylist:
                # If one there, assumes all are in file
                self.qcIncidence = level2b.variables['qcIncidence'][:, :]
                for freq in FREQS:
                    for chan in CHANS:
                        setattr(
                            self, 'qctb' + freq + chan,
                            level2b.variables['qctb'+freq+chan.lower()][:, :])

    #########################################

    def _consider_land_frac_netcdf(self, level2b):
        """
        Helper method to _read_level2b_netcdf()
        Assigns land fraction variables to their appropriate attributes.
        Accounts for whether the data are from older netCDF
        formats (2014-2017) or the new CF-compliant format (2017-).

        Parameters
        ----------
        level2b : netCDF4.Dataset object
            The Dataset object read in from the AMPR data file.
        """
        if self.CF_flag:
            if 'LandFraction' in self.keylist:
                lf = np.array(level2b.variables['LandFraction'])
                if np.ndim(lf) == 2:
                    self.Land_Fraction10 = lf
                else:
                    for j, freq in enumerate(FREQS):
                        setattr(self, 'Land_Fraction' + freq, lf[j, :, :])
        else:
            if 'FovWaterFrac10' in level2b.variables:
                self.Land_Fraction10 = 1.0 - \
                    level2b.variables['FovWaterFrac10'][:, :]
                self.Land_Fraction19 = 1.0 - \
                    level2b.variables['FovWaterFrac19'][:, :]
                self.Land_Fraction37 = 1.0 - \
                    level2b.variables['FovWaterFrac37'][:, :]
                self.Land_Fraction85 = 1.0 - \
                    level2b.variables['FovWaterFrac85'][:, :]

    #########################################

    def _populate_time_vars_netcdf(self):
        """Extracts needed info from dateTime attribute"""
        for icount in np.arange(self.nscans):
            self.Year[icount] = np.int32(self.dateTimes[icount].year)
            self.Month[icount] = np.int32(self.dateTimes[icount].month)
            self.Day[icount] = np.int32(self.dateTimes[icount].day)
            self.Day_of_Year[icount] = \
                np.int32(self.dateTimes[icount].timetuple().tm_yday)
            self.Hour[icount] = np.int32(self.dateTimes[icount].hour)
            self.Minute[icount] = np.int32(self.dateTimes[icount].minute)
            self.Second[icount] = np.int32(self.dateTimes[icount].second)
            self._set_timestring_and_sod(icount)

    #########################################

    def _set_timestring_and_sod(self, index):
        """Create a Time String and get Second of Day for each scan"""
        ts, self.Second_of_Day[index] = \
            _get_timestring_and_sod(self.Hour[index], self.Minute[index],
                                    self.Second[index])
        self.Time_String.append(ts)

    #########################################

    def _read_level2b_ascii(self, full_path_and_filename,
                            project=DEFAULT_PROJECT_NAME):

        # Following puts entire file contents into self.ampr_string
        read_successful = self._read_ampr_ascii_file(full_path_and_filename)
        if not read_successful:
            return

        # Assigning project name (self.Project) based on user input
        self._assign_project_name(project)

        # Determine number of scans
        self.nscans = np.size(self.ampr_string)
        print('Number of scans =', self.nscans)

        # Hard-coding sizes of AMPR arrays and dictionaries
        self._hard_code_ampr_array_sizes_and_other_metadata()
        self._declare_ampr_variables()

        # Populate the variables with the file's data.
        # Begin master loop
        for index, line in enumerate(self.ampr_string):
            line_split = line.split()

            # Header info
            if self.Project in ['CAMP2EX', 'ORACLES', 'OLYMPEX',
                                'IPHEX', 'MC3E']:
                # 2011+, header info placed before TBs
                self._fill_2011on_header_info(line_split, index)
            else:
                # All pre-MC3E projects, aircraft data placed with header
                # info before TBs
                self._fill_pre2011_header_and_aircraft_info(line_split,
                                                            index)
            self._set_timestring_and_sod(index)

            # Get TBs, Latitudes, Longitudes, etc.
            for i in np.arange(self.swath_size):

                if self.Project == 'IPHEX':
                    self._fill_2011on_ampr_variables(line_split, index, i)
                    # Below are IPHEX-specific fixes
                    # 37 GHZ A&B accidentally swapped during IPHEX
                    # Remains unfixed in ASCII data (fixed in netCDF)
                    self.TB37B[index, i] = float(line_split[i + 9 +
                                                 4 * self.swath_size])
                    self.TB37A[index, i] = float(line_split[i + 9 +
                                                 5 * self.swath_size])
                    # Note: Currently (May 2014) Land_Fraction## is set to
                    # bad data in IPHEX ASCII files
                    # Terrain elevation data not recorded, instead incidence
                    # angle is in its position
                    self.Elevation[index, i] = self.bad_data
                    # IPHEX incidence angle currently ignored by PyAMPR
                    # self.Incidence_Angle[index, i] = float(line_split[i
                    # +27+13*self.swath_size])

                elif self.Project in ['CAMP2EX', 'ORACLES', 'OLYMPEX', 'MC3E']:
                    self._fill_2011on_ampr_variables(line_split, index, i)

                else:  # Pre-MC3E projects
                    self._fill_pre2011_ampr_variables(line_split, index, i)

            if self.Project in ['CAMP2EX', 'ORACLES', 'OLYMPEX',
                                'IPHEX', 'MC3E']:
                # Aircraft_Nav out of order to improve efficiency
                # for 2011+ projects
                self._fill_2011on_aircraft_info(line_split, index)

        self._fill_epoch_time()  # defines and populates self.Epoch_Time

        # Replace unfilled data.
        self._set_unfilled_data_to_bad()
        if self.Year[0] < 2011:
            print('Only A channels available in this project\'s data')

    #########################################

    def _filter_bad_geolocations(self, plon=None, plat=None, zdata=None,
                                 equator=False, verbose=True):
        """
        Internal method to filter bad geolocation data.
        Called by following:
        plot_ampr_track()
        write_ampr_kmz()
        """
        # Attempt to deal with bad geolocation data
        # (e.g., Latitude/Longitude=bad_data)
        cond1 = np.logical_or(plon < -180, plon > 180)
        cond2 = np.logical_or(plat < -90, plat > 90)
        condition = np.logical_or(cond1, cond2)
        indices = np.where(condition)
        if np.shape(indices)[1] > 0:
            if verbose:
                print(np.shape(indices)[1],
                      'bad geolocation(s) (e.g., ' + str(self.bad_data) +
                      's) exist, attempting correction.')
            zdata = np.delete(zdata, indices[0], axis=0)
            plon = np.delete(plon, indices[0], axis=0)
            plat = np.delete(plat, indices[0], axis=0)

        # Attempt to deal with bad 0s/-1s in Latitude/Longitude
        if not equator:
            cond1 = np.logical_or(plon == 0, plat == 0)
            cond2 = np.logical_or(plon == -1, plat == -1)
            condition = np.logical_or(cond1, cond2)
            indices = np.where(condition)
            if np.shape(indices)[1] > 0:
                if verbose:
                    print(np.shape(indices)[1],
                          'bad geolocation(s) (0s/-1s) exist,',
                          'attempting correction.')
                    print('If aircraft crossed Equator or Prime Meridian,',
                          'try keyword equator=True.')
                zdata = np.delete(zdata, indices[0], axis=0)
                plon = np.delete(plon, indices[0], axis=0)
                plat = np.delete(plat, indices[0], axis=0)

        # Attempt to deal with remaining wildly varying
        # Latitude/Longitude in single scan
        plat_max = np.amax(plat, axis=1)
        plat_min = np.amin(plat, axis=1)
        plon_max = np.amax(plon, axis=1)
        plon_min = np.amin(plon, axis=1)
        cond1 = plat_max - plat_min > 5  # Usually 2 deg is sufficient,
        cond2 = plon_max - plon_min > 5  # being extra safe here.
        condition = np.logical_or(cond1, cond2)
        indices = np.where(condition)
        if np.shape(indices)[1] > 0:
            if verbose:
                print(np.shape(indices)[1],
                      'scan(s) with high intra-scan geolocation variance',
                      'exist, attempting correction.')
            zdata = np.delete(zdata, indices[0], axis=0)
            plon = np.delete(plon, indices[0], axis=0)
            plat = np.delete(plat, indices[0], axis=0)
        return plon, plat, zdata

    #########################################

    def _get_scan_indices(self, scanrange=None, timerange=None, verbose=True):

        """
        Internal method to get scan indices. Used by:
        plot_ampr_track()
        write_ampr_kml()
        write_ampr_kmz()
        plot_ampr_channels()
        Prioritizes timerange over scanrange. Currently cannot handle date
        changes between times. User must manually break down timerange by
        dates and call for each date separately.

        """
        if verbose:
            print('Available scans =',
                  np.min(self.Scan), 'to', np.max(self.Scan))
            print('Available times =', str(self.Time_String[0]), str('-'),
                  str(self.Time_String[self.nscans-1]))
            bt = time.time()

        if not scanrange and not timerange:
            ind1, ind2 = self._get_min_max_indices()

        elif scanrange is not None and timerange is None:
            indices = np.where(self.Scan == np.min(scanrange))
            if np.shape(indices[0])[0] == 0:
                if verbose:
                    print('Scan number too small,',
                          'using first scan for beginning')
                ind1 = 0
            else:
                ind1 = indices[0][0]
            indices = np.where(self.Scan == np.max(scanrange))
            if np.shape(indices[0])[0] == 0:
                if verbose:
                    print('Scan number too high, using last scan for end')
                ind2 = self.nscans
            else:
                ind2 = indices[0][0]

        else:
            try:
                t1 = _get_sod(
                    float(timerange[0][0:2]), float(timerange[0][3:5]),
                    float(timerange[0][6:8]))
                t2 = _get_sod(
                    float(timerange[1][0:2]), float(timerange[1][3:5]),
                    float(timerange[1][6:8]))
                cond1 = self.Second_of_Day >= np.min([t1, t2])
                cond2 = self.Second_of_Day <= np.max([t1, t2])
                condition = np.logical_and(cond1, cond2)
                indices = np.where(condition)
                if np.shape(indices[0])[0] > 0:
                    ind1 = np.min(indices[0])
                    ind2 = np.max(indices[0])
                else:
                    if verbose:
                        _print_times_not_valid()
                    ind1, ind2 = self._get_min_max_indices()
            except:
                if verbose:
                    _print_times_not_valid()
                ind1, ind2 = self._get_min_max_indices()

        if verbose:
            print(time.time()-bt, 'seconds to process scan indices')
        return ind1, ind2

    #########################################

    def _adjust_colorbar(self, cbar, show_qc, colorbar_label=True):
        if show_qc:
            cbar.set_ticks([1, 2, 3, 4, 5])
            if colorbar_label:
                cbar.set_label('QC Flag')
        else:
            if colorbar_label:
                cbar.set_label('Brightness Temperature (K)')
        return cbar

    #########################################

    def _check_qc(self, show_qc, clevs):
        """Used by strip charts, track plots, and Google Earth maps"""
        if show_qc:
            if not hasattr(self, 'qctb10a'):
                print('No QC flags available, plotting TB data instead')
                show_qc = False
            else:
                clevs = [0, 6]
        return show_qc, clevs

    #########################################

    def _initialize_time_fields(self):
        """
        Common to ASCII and netCDF reads.
        Before looping through the datetime object, initialize numpy arrays.
        """
        self.Year = np.zeros(self.nscans, dtype=np.int32)
        self.Month = np.zeros(self.nscans, dtype=np.int32)
        self.Day = np.zeros(self.nscans, dtype=np.int32)
        self.Day_of_Year = np.zeros(self.nscans, dtype=np.int32)
        self.Hour = np.zeros(self.nscans, dtype=np.int32)
        self.Minute = np.zeros(self.nscans, dtype=np.int32)
        self.Second = np.zeros(self.nscans, dtype=np.int32)
        self.Second_of_Day = np.zeros(self.nscans, dtype=np.int32)
        self.Time_String = []

    #########################################

    def _declare_ampr_variables(self):
        """Define variables to be populated"""
        # Timing, Icon, and Noise
        self._initialize_time_fields()
        self.Scan = np.zeros(self.nscans, dtype=np.int32)
        self.Icon = np.zeros(self.nscans, dtype=np.int32)
        self.Noise10 = np.zeros(self.nscans, dtype='float')
        self.Noise19 = np.zeros(self.nscans, dtype='float')
        self.Noise37 = np.zeros(self.nscans, dtype='float')
        self.Noise85 = np.zeros(self.nscans, dtype='float')

        # Geolocation and land information
        self.Latitude = np.zeros((self.nscans, self.swath_size),
                                 dtype='float')
        self.Longitude = np.zeros((self.nscans, self.swath_size),
                                  dtype='float')
        self.Land_Fraction10 = np.zeros((self.nscans, self.swath_size),
                                        dtype='float')
        self.Land_Fraction37 = np.zeros((self.nscans, self.swath_size),
                                        dtype='float')
        self.Land_Fraction85 = np.zeros((self.nscans, self.swath_size),
                                        dtype='float')
        self.Elevation = np.zeros((self.nscans, self.swath_size),
                                  dtype='float')

        # Brightness Temperatures
        self.TB10A = np.zeros((self.nscans, self.swath_size), dtype='float')
        self.TB10B = np.zeros((self.nscans, self.swath_size), dtype='float')
        self.TB19A = np.zeros((self.nscans, self.swath_size), dtype='float')
        self.TB19B = np.zeros((self.nscans, self.swath_size), dtype='float')
        self.TB37A = np.zeros((self.nscans, self.swath_size), dtype='float')
        self.TB37B = np.zeros((self.nscans, self.swath_size), dtype='float')
        self.TB85A = np.zeros((self.nscans, self.swath_size), dtype='float')
        self.TB85B = np.zeros((self.nscans, self.swath_size), dtype='float')
        # Aircraft Nav
        self._initialize_aircraft_dict()

    #########################################

    def _initialize_aircraft_dict(self):
        """Aircraft Navgation info"""
        self.Aircraft_varlist = [u'GPS Latitude', u'GPS Longitude',
                                 u'GPS Altitude', u'Pitch', u'Roll',
                                 u'Yaw', u'Heading', u'Ground Speed',
                                 u'Air Speed', u'Static Pressure',
                                 u'Total Pressure', u'Total Temperature',
                                 u'Static Temperature', u'Wind Speed',
                                 u'Wind Direction', u'INS Latitude',
                                 u'INS Longitude', u'INS Altitude']
        self.Aircraft_Nav = {}
        for var in self.Aircraft_varlist:
            self.Aircraft_Nav[var] = np.zeros(self.nscans, dtype=float)

    #########################################

    def _fill_epoch_time(self):
        """Calculate Epoch_Time attribute"""
        self.Epoch_Time = 0 * self.Second
        for i in np.arange(self.nscans):
            self.Epoch_Time[i] = calendar.timegm(
                (int(self.Year[i]), int(self.Month[i]), int(self.Day[i]),
                 int(self.Hour[i]), int(self.Minute[i]), int(self.Second[i])))

    #########################################

    def _set_unfilled_data_to_bad(self):
        """
        Vectorizing remaining data assignments for unused/duplicate variables
        Project dependent which variables are or are not used
        """
        if self.Project in ['CAMP2EX', 'ORACLES', 'OLYMPEX', 'IPHEX', 'MC3E']:
            self.Noise10[:] = self.bad_data
            self.Noise19[:] = self.bad_data
            self.Noise37[:] = self.bad_data
            self.Noise85[:] = self.bad_data
        else:
            self.Aircraft_Nav['Static Pressure'][:] = self.bad_data
            self.Aircraft_Nav['Total Pressure'][:] = self.bad_data
            self.Aircraft_Nav['Static Temperature'][:] = self.bad_data
            self.Aircraft_Nav['Total Temperature'][:] = self.bad_data
            self.Aircraft_Nav['Wind Speed'][:] = self.bad_data
            self.Aircraft_Nav['Wind Direction'][:] = self.bad_data
            self.Aircraft_Nav['INS Latitude'][:] = self.bad_data
            self.Aircraft_Nav['INS Longitude'][:] = self.bad_data
            self.Aircraft_Nav['INS Altitude'][:] = self.bad_data
            self.TB10B[:, :] = self.bad_data
            self.TB19B[:, :] = self.bad_data
            self.TB37B[:, :] = self.bad_data
            self.TB85B[:, :] = self.bad_data
            # Only one Land_Fraction variable with old projects
            self.Land_Fraction37[:, :] = self.bad_data
            self.Land_Fraction85[:, :] = self.bad_data
        self._remove_nan_inf()

    #########################################

    def _remove_nan_inf(self):
        """
        Attempting to control for infinite and nan numbers.
        They can blow up plot_ampr_track() and write_ampr_kmz().
        """
        # infinite and nan check
        self.TB10A[np.isfinite(self.TB10A) is False] = self.bad_data
        self.TB19A[np.isfinite(self.TB19A) is False] = self.bad_data
        self.TB37A[np.isfinite(self.TB37A) is False] = self.bad_data
        self.TB85A[np.isfinite(self.TB85A) is False] = self.bad_data
        if hasattr(self, 'TB10B'):  # Assume if one, all are there
            self.TB10B[np.isfinite(self.TB10B) is False] = self.bad_data
            self.TB19B[np.isfinite(self.TB19B) is False] = self.bad_data
            self.TB37B[np.isfinite(self.TB37B) is False] = self.bad_data
            self.TB85B[np.isfinite(self.TB85B) is False] = self.bad_data
        else:  # Set all the B channels to bad, these are single-channel data
            for freq in FREQS:
                setattr(self, 'TB' + freq + 'B',
                        self.bad_data * np.ones((self.nscans, self.ncross),
                                                dtype='float'))
        # Address geolocations if available
        if hasattr(self, 'Latitude'):
            self.Latitude[np.isfinite(self.Latitude) is False] = self.bad_data
            self.Longitude[np.isfinite(self.Longitude) is False] = \
                self.bad_data

    #########################################

    def _get_gearth_file_name(self, var=None, index=0, suffix=None):
        """Obtains default file name for Google Earth KMZ"""
        mo, dy = self._get_month_and_day_string(index)
        timestamp = str(self.Time_String[index]).replace(':', '')
        return str(self.Year[index]) + str(mo) + str(dy) + str('_') + \
            str(timestamp) + str('z_TB') + str(var.upper()) + str(suffix)

    #########################################

    def _get_data_subsection(self, var=None, ind1=None, ind2=None,
                             maneuver=True, show_qc=False, verbose=True):
        """Subsections the data for later plotting"""
        if show_qc:
            zdata = 1.0 * getattr(self, 'qctb'+var.lower())
        else:
            zdata = 1.0 * getattr(self, 'TB'+var.upper())
        zdata = zdata[ind1:ind2]
        plon = 1.0 * getattr(self, 'Longitude')
        plon = plon[ind1:ind2]
        plat = 1.0 * getattr(self, 'Latitude')
        plat = plat[ind1:ind2]
        if not maneuver:
            if verbose:
                print('Filtering out significant aircraft maneuvers')
            roll = self.Aircraft_Nav['Roll'][ind1:ind2]
            indices = np.where(np.abs(roll) >= 5)
            if np.shape(indices)[1] > 0:
                zdata = np.delete(zdata, indices[0], axis=0)
                plon = np.delete(plon, indices[0], axis=0)
                plat = np.delete(plat, indices[0], axis=0)
        return plon, plat, zdata

    #########################################

    def _fill_2011on_header_info(self, line_split=None, index=None):
        """For ASCII data, fill 2011+ metadata variables"""
        self.Scan[index] = int(line_split[0])
        self.Year[index] = int(line_split[1])
        self.Month[index] = int(line_split[2])
        self.Day[index] = int(line_split[3])
        self.Day_of_Year[index] = int(line_split[4])
        self.Hour[index] = int(line_split[5])
        self.Minute[index] = int(line_split[6])
        self.Second[index] = int(line_split[7])
        self.Icon[index] = int(line_split[8])

    #########################################

    def _fill_pre2011_header_and_aircraft_info(self, line_split=None,
                                               index=None):
        """For ASCII data, fill < 2011 metadata variables"""
        # Begin theatrics to handle Year appearing only in first line,
        # but in different positions depending on project.
        # In JAX90, COARE, & CAPE files it wipes out Day_of_Year[0].
        # In other files it occupies the place of Scan[0].
        self.Scan[index] = index + 1
        if (self.Project == 'JAX90' or self.Project == 'CAPE' or
                self.Project == 'COARE') and index == 0:
            self.Year[:] = int(line_split[1])
        if (self.Project != 'JAX90' and self.Project != 'CAPE' and
                self.Project != 'COARE' and index == 0):
            self.Year[:] = int(line_split[0])
        self.Day_of_Year[index] = int(line_split[1])
        doy = int(line_split[1]) - 1
        sdate = str(datetime.date(self.Year[index], 1, 1) +
                    datetime.timedelta(doy))
        self.Month[index] = int(sdate[5:7])
        self.Day[index] = int(sdate[8:10])
        if (self.Project == 'JAX90' or self.Project == 'CAPE' or
                self.Project == 'COARE') and index == 1:
            # Assuming date change doesn't occur between lines 1 and 2
            self.Day_of_Year[0] = self.Day_of_Year[1]
            self.Month[0] = self.Month[1]
            self.Day[0] = self.Day[1]
        # End theatrics - WHEW!

        self.Hour[index] = int(line_split[2])
        self.Minute[index] = int(line_split[3])
        self.Second[index] = int(line_split[4])
        self.Icon[index] = int(line_split[5])
        self.Aircraft_Nav['GPS Latitude'][index] = float(line_split[6])
        self.Aircraft_Nav['GPS Longitude'][index] = float(line_split[7])
        self.Aircraft_Nav['GPS Altitude'][index] = float(line_split[8])
        self.Aircraft_Nav['Pitch'][index] = float(line_split[9])
        self.Aircraft_Nav['Roll'][index] = float(line_split[10])
        self.Aircraft_Nav['Yaw'][index] = float(line_split[11])
        self.Aircraft_Nav['Heading'][index] = float(line_split[12])
        self.Aircraft_Nav['Air Speed'][index] = float(line_split[13])
        self.Aircraft_Nav['Ground Speed'][index] = float(line_split[14])
        self.Noise10[index] = float(line_split[15])
        self.Noise19[index] = float(line_split[16])
        self.Noise37[index] = float(line_split[17])
        self.Noise85[index] = float(line_split[18])

    #########################################

    def _fill_2011on_ampr_variables(self, line_split=None, index=None,
                                    i=None):
        """
        Data written out left to right; i=0 is Left edge,
        i=49 is Right edge
        """
        self.TB10A[index, i] = float(line_split[i + 9])
        self.TB10B[index, i] = float(line_split[i + 9 + self.swath_size])
        self.TB19A[index, i] = float(line_split[i + 9 + 2 * self.swath_size])
        self.TB19B[index, i] = float(line_split[i + 9 + 3 * self.swath_size])
        self.TB37A[index, i] = float(line_split[i + 9 + 4 * self.swath_size])
        self.TB37B[index, i] = float(line_split[i + 9 + 5 * self.swath_size])
        self.TB85A[index, i] = float(line_split[i + 9 + 6 * self.swath_size])
        self.TB85B[index, i] = float(line_split[i + 9 + 7 * self.swath_size])
        self.Latitude[index, i] = float(
            line_split[i + 9 + 8 * self.swath_size])
        self.Longitude[index, i] = float(
            line_split[i + 9 + 9 * self.swath_size])
        self.Land_Fraction10[index, i] = float(
            line_split[i + 27 + 10 * self.swath_size])
        self.Land_Fraction37[index, i] = float(
            line_split[i + 27 + 11 * self.swath_size])
        self.Land_Fraction85[index, i] = float(
            line_split[i + 27 + 12 * self.swath_size])
        self.Elevation[index, i] = float(
            line_split[i + 27 + 13 * self.swath_size])

    #########################################

    def _fill_pre2011_ampr_variables(self, line_split=None, index=None,
                                     i=None):
        """
        Data written out left to right; i=0 is Left edge,
        i=49 is Right edge
        """
        self.TB10A[index, i] = float(line_split[i + 19])
        self.TB19A[index, i] = float(line_split[i + 19 + self.swath_size])
        self.TB37A[index, i] = float(line_split[i + 19 + 2 * self.swath_size])
        self.TB85A[index, i] = float(line_split[i + 19 + 3 * self.swath_size])
        self.Latitude[index, i] = float(
            line_split[i + 19 + 4 * self.swath_size])
        self.Longitude[index, i] = float(
            line_split[i + 19 + 5 * self.swath_size])
        self.Elevation[index, i] = float(
            line_split[i + 19 + 6 * self.swath_size])
        self.Land_Fraction10[index, i] = float(
            line_split[i + 19 + 7 * self.swath_size])

    #########################################

    def _check_for_enough_data_to_plot(self, plon=None, plat=None):
        """Test on amount of data available to plot, returns True/False"""
        cond1 = np.logical_and(plon >= -180, plon <= 180)
        cond2 = np.logical_and(plat >= -90, plat <= 90)
        condition = np.logical_and(cond1, cond2)
        indices = np.where(condition)
        if np.shape(indices)[1] < 100:
            print(np.shape(indices)[1], 'good gelocations,',
                  'need 100+ (i.e., 2+ scans).')
            print('Not enough good geolocation data to plot, returning.')
            return False
        else:
            return True

    #########################################

    def _check_aspect_ratio(self, latrange=None, lonrange=None,
                            verbose=True):
        """
        Provides a warning if the prospective plot's aspect ratio will
        cause colorbar to be far removed from plot window itself.
        """
        aspect_ratio = (float(np.max(latrange)) - float(np.min(latrange))) /\
                       (float(np.max(lonrange)) - float(np.min(lonrange)))
        if aspect_ratio < 0.5 or aspect_ratio > 2:
            if verbose:
                print('Warning: Your aspect ratio choice could lead to poor',
                      'plotting results.')
                print('Best results are obtained when latrange ~ lonrange.')

    #########################################

    def _hard_code_ampr_array_sizes_and_other_metadata(self):
        """Hard coding certain fixed AMPR characteristics"""
        self.swath_size = DEFAULT_SWATH_SIZE
        self.nav_size = DEFAULT_NAV_SIZE
        self.bad_data = DEFAULT_BAD_DATA
        self.swath_left = DEFAULT_SWATH_LEFT
        self.swath_angle = self.swath_left - \
            (2.0 * self.swath_left / float(self.swath_size - 1.0)) * \
            np.arange(self.swath_size)

    #########################################

    def _fill_2011on_aircraft_info(self, line_split, index):
        """For ASCII data, fill 2011+ aircraft navigation variables"""
        aircraft_i = 0
        for name in self.Aircraft_varlist:
            self.Aircraft_Nav[name][index] = \
                    float(line_split[aircraft_i + 9 + 10 * self.swath_size])
            aircraft_i += 1

    #########################################

    def _get_list_of_channels_to_plot(self, show_pol=False, show_qc=False):
        """Used by plot_ampr_channels() to figure out what variables to plot"""
        if show_qc:
            return ['10A', '10B', '19A', '19B', '37A', '37B', '85A', '85B']
        if show_pol:
            tb_list = ['10V', '10H', '19V', '19H', '37V', '37H', '85V', '85H']
            if (not hasattr(self, 'TB10V') or not hasattr(self, 'TB19V') or not
                    hasattr(self, 'TB37V') or not hasattr(self, 'TB85V')):
                print('Missing some pol channels,',
                      'trying calc_polarization() before plot')
                self.calc_polarization()
        else:
            tb_list = ['10A', '10B', '19A', '19B', '37A', '37B', '85A', '85B']
        return tb_list

    #########################################

    def _missing_channel_printout(self):
        """Simple warning message if requested channel is missing"""
        print('Channel doesn\'t exist, check typing or try reading in a file')
        print('Acceptable channels = 10A, 10B, 19A, 19B, 37A, 37B, 85A, 85B')
        print('If calculated, also = 10H, 10V, 19H, 19V, 37H, 37V, 85H, 85V')

    #########################################

    def _get_latrange_lonrange(self, plat=None, plon=None,
                               latrange=None, lonrange=None):
        """Determine domain of plot based on what user provided"""
        if latrange is None:
            latrange = [np.min(plat), np.max(plat)]
        if lonrange is None:
            lonrange = [np.min(plon), np.max(plon)]
        return latrange, lonrange

    #########################################

    def _get_colormap(self, cmap, flag, qc_flag=False):
        """Figure out colormap based on user input"""
        if cmap is None:
            if flag:
                if qc_flag:
                    cmap = amprQC_cmap
                else:
                    cmap = amprTB_cmap
            else:
                cmap = 'jet'
        return cmap

    #########################################

    def _get_date_string(self, index=0):
        """Get date string that is used in plot titles"""
        return str(self.Month[index]) + '/' + str(self.Day[index]) + '/' + \
            str(self.Year[index])

    #########################################

    def _get_ampr_title(self, var=None):
        """Get default"""
        return 'AMPR '+var[0:2]+' GHz ('+var[2].upper()+') '

    #########################################

    def _read_ampr_ascii_file(self, full_path_and_filename):
        """
        Ingest the AMPR ASCII file as a huge string array, which
        will be parsed for data later on. Support for gzipped files
        and error checking provided.
        """
        if full_path_and_filename[-3:] == '.gz':
            try:
                fileobj = gzip.open(full_path_and_filename)
                ascii = codecs.getreader('ASCII')
                fileobj = ascii(fileobj)
            except:
                print('Incorrect file or file doesn\'t exist, returning')
                return False
        else:
            try:
                fileobj = open(full_path_and_filename, 'r')
            except:
                print('Incorrect file or file doesn\'t exist, returning')
                return False
        try:
            contents = fileobj.readlines()
        except:
            print('File not ASCII format, returning')
            fileobj.close()
            return False
        fileobj.close()
        self.ampr_string = np.array(contents)
        return True

    #########################################

    def _assign_project_name(self, project=DEFAULT_PROJECT_NAME):
        """Check user-provided project name and keep track of it"""
        if not isinstance(project, str):
            print('Bad project name, provide actual string')
            print('Assuming', DEFAULT_PROJECT_NAME, 'data structure.')
            project = DEFAULT_PROJECT_NAME
        else:
            print('Assuming', project.upper(), 'data structure.')
        print('Change to proper project if incorrect, otherwise errors',
              'will occur.')
        print('Currently available field projects: CAMP2EX, ORACLES, OLYMPEX,'
              'IPHEX, MC3E, TC4, TCSP, JAX90, COARE,')
        print('CAMEX1, CAMEX2, CAMEX3, CAMEX4, TRMMLBA, KWAJEX, TEFLUNA,',
              'FIRE3ACE, CAPE')
        print('Default: project = \''+DEFAULT_PROJECT_NAME+'\'')
        self.Project = project.upper()

    #########################################

    def _get_min_max_indices(self):
        """Return all possible scan indices"""
        return 0, self.nscans

    #########################################

    def _solve_using_simple_linear_substitution(self, angle=None,
                                                T1=None, T2=None):
        """Lead author: Brent Robers"""
        tbv = (T1 - T1 * np.cos(angle)**2 - T2 * np.cos(angle)**2) /\
              (np.sin(angle)**2 - np.cos(angle)**2)
        tbh = T1 + T2 - tbv
        return tbv, tbh

    #########################################

    def _solve_using_constrained_linear_inversion(self, angle=None,
                                                  T1=None, T2=None):
        """Lead author: Brent Robers"""
        # Set regularization parameter
        gam = 10.0**(-1)
        # Solve
        # n.b. I solved this on paper so that a single set of equations
        # can be used and applied using matrix notation.
        # Define parameters needed.
        a = np.cos(angle)**2
        b = np.sin(angle)**2
        c = np.sin(angle)**2
        d = np.cos(angle)**2
        zeta = 1.0 / ((a**2 + c**2 + gam**2) *
                      (b**2 + d**2 + gam**2) - (a * b + c * d)**2)
        # Get tbv and tbh
        tbv = T1 * zeta * (b * (a**2 + c**2 + gam**2) - a * (a * b + c * d)) +\
            T2 * zeta * (d * (a**2 + c**2 + gam**2) - c * (a * b + c * d))
        tbh = T1 * zeta * (a * (b**2 + d**2 + gam**2) - b * (a * b + c * d)) +\
            T2 * zeta * (c * (b**2 + d**2 + gam**2) - d * (a * b + c * d))
        return tbv, tbh

    #########################################

    def _compute_nadir_offset_and_apply_if_desired(self, angle=None, T1=None,
                                                   T2=None, force_match=True,
                                                   chan=None):
        """Lead author: Brent Roberts"""
        scanang = np.rad2deg(angle)
        # xoff = np.empty(self.nscans) * np.nan
        xoff = np.zeros(self.nscans)
        for iscan in np.arange(self.nscans):
            # Get value.
            x1 = T1[iscan, :]
            x2 = T2[iscan, :]
            # Interpolate x1 and x2 to 0 degrees.
            x1_zero = np.interp(-45.0, scanang, x1)
            x2_zero = np.interp(-45.0, scanang, x2)
            # Compute the offset.
            xoffset = x1_zero - x2_zero
            # Adjust T2 to match T1.
            if force_match:
                T2[iscan, :] = T2[iscan, :] + xoffset
            xoff[iscan] = xoffset
        setattr(self, 'TB'+chan+'_offset', xoff)
        return T2

    #########################################

    def _check_for_pol_data(self):
        """Simple check on whether polarization deconvolved data exist"""
        if (not hasattr(self, 'TB10V') and not hasattr(self, 'TB19V') and not
                hasattr(self, 'TB37V') and not hasattr(self, 'TB85V')):
            print('No polarization data to plot! Returning ...')
            return False
        else:
            return True

    #########################################

    def _get_timestamps_for_gearth(self, ind1=None, ind2=None):
        """Format example: '1997-07-16T07:30:15Z'"""
        mo, dy = self._get_month_and_day_string(ind1)
        time1 = str(self.Year[ind1]) + str('-') + str(mo) + str('-') + \
            str(dy) + str('T') + str(self.Time_String[ind1]) + str('Z')
        mo, dy = self._get_month_and_day_string(ind2-1)
        time2 = str(self.Year[ind2-1]) + str('-') + str(mo) + str('-') + \
            str(dy) + str('T') + str(self.Time_String[ind2-1]) + str('Z')
        times = [time1, time2]
        return times

    #########################################

    def _get_month_and_day_string(self, index):
        """
        Return month and day as strings for use in creating file
        names. It places 0s in front single-digit numbers.
        """
        smo = str(self.Month[index])
        if self.Month[index] < 10:
            smo = '0' + smo
        sdy = str(self.Day[index])
        if self.Day[index] < 10:
            sdy = '0' + sdy
        return str(smo), str(sdy)

    #########################################

    #########################################

    #######################################
    # Add more attributes and methods here!
    #######################################

    #########################################

# End class AmprTb definition
###########################################################

# Stand-alone functions


def parse_ax_fig(ax, fig, projection=None):
    """ Parse and return ax and fig parameters. """
    if fig is None:
        fig = plt.figure(figsize=(8, 8))
    if ax is None:
        ax = fig.add_axes([0.1, 0.15, 0.8, 0.8], projection=projection)
#    if ax is None:
#        ax = plt.gca()
#    if fig is None:
#        fig = plt.gcf()
    return ax, fig


def _get_timestring_and_sod(hour=None, minute=None, second=None):
    """Time_String: Fill size gaps with 0s"""
    d = str(hour)
    if hour < 10:
        d = str('0' + d)
    e = str(minute)
    if minute < 10:
        e = str('0' + e)
    f = str(second)
    if second < 10:
        f = str('0' + f)
    return str(d + str(':') + e + str(':') + f), _get_sod(hour, minute, second)


def _get_sod(hour=None, minute=None, second=None):
    """Calculate second of day given hour, minute, second"""
    return 3600.0 * hour + 60.0 * minute + second


def _method_footer_printout():
    """Helps clarify text output"""
    print('********************')
    print('')


def _method_header_printout():
    """Helps clarify text output"""
    print('')
    print('********************')


def _print_times_not_valid():
    """Warning message if user provided bad timerange keyword"""
    print('Times not valid, just plotting everything')
    print('Next time try timerange=[\'hh:mm:ss\',\'HH:MM:SS\']')
