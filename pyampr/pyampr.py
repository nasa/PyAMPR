import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import time
import datetime
import calendar
import gzip
import codecs
from google_earth_tools import gearth_fig, make_kml
from netCDF4 import Dataset, num2date, date2num

try:
    from udf_cmap import amprTB_cmap
    CMAP_FLAG = True
except ImportError:
    CMAP_FLAG = False

VERSION = '1.4.0'

#Fixed constants used by PyAMPR set here
DEFAULT_CLEVS = [75, 325]
DEFAULT_GRID_DEL = 2.0
DEFAULT_VAR = '10A'
DEFAULT_CHAN_LIST = ['10', '19', '37', '85']
DEFAULT_SWATH_SIZE = 50
DEFAULT_NAV_SIZE = 18
DEFAULT_BAD_DATA = -999.0
DEFAULT_SWATH_LEFT = -44.1
DEFAULT_PROJECT_NAME = 'IPHEX'

######################
#Main class definition
######################

class AmprTb(object):

    def __init__(self, full_path_and_filename=None,
                 project=DEFAULT_PROJECT_NAME):
        """
If passed a filename, call the read_ampr_tb_level2b() method, 
otherwise just instance the class with nothing

        """
        if full_path_and_filename == None:
           print 'Class instantiated, call read_ampr_tb_level2b() to populate'
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

Currently available field projects: IPHEX, MC3E, TC4, TCSP, JAX90, COARE,
CAMEX1, CAMEX2, CAMEX3, CAMEX4, TRMMLBA, KWAJEX, TEFLUNA, FIRE3ACE, CAPE
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
        print 'read_ampr_tb_level2b(): Reading', full_path_and_filename
            
        try:
            self._read_level2b_netcdf(full_path_and_filename, project=project)
        except:
            print 'Not netCDF file, trying ASCII read ...'
            self._read_level2b_ascii(full_path_and_filename, project=project)
        _method_footer_printout()
    
    #########################################

    def help(self):
        help(self)

    #########################################
    
    def plot_ampr_track(self, var=DEFAULT_VAR, latrange=None, lonrange=None,
                        parallels=DEFAULT_GRID_DEL, meridians=DEFAULT_GRID_DEL,
                        title=None, clevs=DEFAULT_CLEVS, cmap=None,
                        save=None, show_track=False, maneuver=True,
                        scanrange=None, show_grid=True, equator=False,
                        timerange=None, return_flag=False):

        """
This method plots geolocated AMPR data, along with the Aircraft track if
requested. matplotlib.pyplot.pcolormesh() on a Basemap is the workhorse 
plotting routine.

var = String with channel number and letter (e.g., 10A for 10 GHz (A) channel
latrange = List with lat range defined. Order and size (>= 2) is irrelevant
           as max and min are retrieved
lonrange = List with lon range defined. Order and size (>= 2) is irrelevant
           as max and min are retrieved
parallels = Scalar spacing (deg) for parallels (i.e., constant latitude)
meridians = Scalar spacing (deg) for meridians (i.e., constant longitude)
ptitle = Plot title as string
clevs = List with contour levels. Only max and min values are used.
cmap = Colormap desired. See http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
       and dir(cm) for more
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
return_flag = Set to True to return Basemap, plot axes, etc. Order of items 
              returned is fig (figure instance), ax (main axis instance), 
              m (Basemap instance), cax (colorbar axis instance), 
              cbar (colorbar instance).
        """

        plt.close() #mpl seems buggy if you don't clean up old windows
        _method_header_printout()
        print 'plot_ampr_track():'

        #10 GHz (A) channel plotted by default if mistake made
        if not isinstance(var, str):
            var = DEFAULT_VAR
        
        #Check to make sure data exist!
        if not hasattr(self, 'TB'+var.upper()):
            self._missing_channel_printout()
            _method_footer_printout()
            return

        if self.Year[0] < 2011:
            print 'Warning: Older projects commonly had bad or missing',\
                  'geolocation data.'
            print 'If there are plotting problems, try a strip chart with',\
                  'plot_ampr_channels(),'
            print 'or try adjusting scanrange, lonrange, or latrange.'
        
        #Adjustable scan range limits
        #Fairly robust - will go down to a width of 10 scans or so before
        #plotting artifacts begin to occur.
        ind1, ind2 =        self._get_scan_indices(scanrange, timerange)
        plon, plat, zdata = self._get_data_subsection(var, ind1, ind2, maneuver)
        plon, plat, zdata = self._filter_bad_geolocations(plon, plat,
                                                         zdata, equator)
        enough_data =       self._check_for_enough_data_to_plot(plon, plat)
        if not enough_data:
            return
        
        latrange, lonrange = self._get_latrange_lonrange(plat, plon,
                                                     latrange, lonrange)
        self._check_aspect_ratio(latrange, lonrange)

        #Create Basemap instance
        lon_0 = np.median(plon)
        lat_0 = np.median(plat)
        llcrnrlat = np.min(latrange)
        urcrnrlat = np.max(latrange)
        llcrnrlon = np.min(lonrange)
        urcrnrlon = np.max(lonrange)
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_axes([0.1, 0.15, 0.8, 0.8])
        m = Basemap(projection='merc', lon_0=lon_0, lat_0=lat_0,\
                    llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,\
                    llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,\
                    resolution='l')
        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()

        if show_grid:
            #Draw parallels
            vparallels = np.arange(np.floor(np.min(latrange)),
                                   np.ceil(np.max(latrange)), parallels)
            m.drawparallels(vparallels, labels=[1,0,0,0], fontsize=10)
            
            #Draw meridians
            vmeridians = np.arange(np.floor(np.min(lonrange)),
                                   np.ceil(np.max(lonrange)), meridians)
            m.drawmeridians(vmeridians, labels=[0,0,0,1], fontsize=10)
            del vparallels, vmeridians

        #Compute map proj coordinates
        x, y = m(plon, plat)

        #Draw filled contours
        cmap = self._get_colormap(cmap, CMAP_FLAG)
        cs = m.pcolormesh(x, y, zdata, vmin=np.min(clevs),
                          vmax=np.max(clevs), cmap=cmap)

        #Add Aircraft track
        if show_track:
            x1,y1 = m(self.Aircraft_Nav['GPS Longitude'][ind1:ind2],
                      self.Aircraft_Nav['GPS Latitude' ][ind1:ind2])
            m.plot(x1, y1, 'k.') #Black dots during normal flight
            indices = np.where(np.abs(self.Aircraft_Nav['Roll'][ind1:ind2]) >= 5)
            #White dots during maneuvers
            m.plot(x1[indices[0]], y1[indices[0]], 'w.')

        #Plot title & display
        if title == None:
            title = self._get_ampr_title(var) + self._get_date_string(ind1) +\
                    ', '+self.Time_String[ind1]+'-'+self.Time_String[ind2-1]+\
                    ' UTC'
        plt.title(title)

        #Add colorbar
        #Colorbar independent of Basemap
        cax = fig.add_axes([0.20, 0.07, 0.60, 0.02])
        cbar = plt.colorbar(cs, cax=cax, orientation='horizontal', extend='both')
        cbar.set_label('Brightness Temperature (K)')

        #Save file
        if save != None:
            plt.savefig(save)

        _method_footer_printout()
        
        if return_flag:
            return fig, ax, m, cax, cbar

    #########################################

    def plot_ampr_channels(self, scanrange=None, cmap=None, clevs=DEFAULT_CLEVS,
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
timerange = Time range to plot. Overrides scanrange if both are set.
            Format: timerange = ['hh:mm:ss', 'HH:MM:SS']

        """

        plt.close() #mpl seems buggy if you don't clean up old windows
        _method_header_printout()
        print 'plot_ampr_channels():'
        
        tb_list = self._get_list_of_channels_to_plot(show_pol)
        if show_pol:
            pol_flag = self._check_for_pol_data()
            if not pol_flag:
                _method_footer_printout()
                return
    
        #Set up plot - 10 rows, 1 column, 8 channels for 8 rows,
        #one row for colorbar, and one row for Aircraft_Nav info
        fig, axes = plt.subplots(nrows=10, ncols=1, sharex=False)
        fig.set_size_inches(11, 8.5)

        #Adjustable scan range limits
        #Fairly robust - will go down to a width of 10 scans before plotting
        #artifacts begin to occur.
        ind1, ind2 = self._get_scan_indices(scanrange, timerange)
        xran = [self.Scan[ind1], self.Scan[ind2-1]]
        plt.xlim(xran)
        
        #The following will put the title above the x-axis tick labels
        ptitle = 'AMPR ' + self._get_date_string(ind1)
        plt.text(0.5, 1.60, ptitle, horizontalalignment='center', fontsize=14,
                 transform = axes[0].transAxes)
        
        colormap = self._get_colormap(cmap, CMAP_FLAG)

        #Loop and plot available channels
        itest = 0
        for index, chan in enumerate(tb_list):
            if hasattr(self, 'TB'+chan):
                itest = itest + 1
                var = np.transpose(getattr(self, 'TB'+chan))
                #AMPR data arranged L-to-R in PyAMPR (index 0 to index 49
                #in a scan),
                #so need to reverse order to have L on top in strip charts.
                im = axes[index].pcolormesh(self.Scan,
                     (self.swath_size-1)-np.arange(self.swath_size), var,
                     vmin=np.min(clevs), vmax=np.max(clevs), cmap=colormap)
                axes[index].set_xlim(xran)
                axes[index].set_ylim([0, self.swath_size-1])
                axes[index].set_ylabel(chan)
                axes[index].yaxis.set_ticks([2,self.swath_size-5])
                axes[index].tick_params(axis='y', labelsize=7)
                axes[index].tick_params(axis='x', labelbottom='off')
                axes[index].tick_params(axis='y', which='both', left='off', 
                                        right='off')
                #Shows how V and H vary with beam position
                if show_pol == False:
                    if index % 2 == 0:
                        axes[index].set_yticklabels(['H','V'])
                    else:
                        axes[index].set_yticklabels(['V','H'])
                else:
                    axes[index].set_yticklabels(['R','L'])

                #Get scan timing and plot it as tick labels
                if index == 0:
                    axes[index].tick_params(axis='x', labelsize=7, labeltop='on',
                                            direction='out', bottom='off', pad=0)
                    locs, labels = plt.xticks()
                    new_labels = []
                    for t in locs:
                        indt = np.where(self.Scan == t)
                        new_labels.append(''.join(self.Time_String[indt[0]]))
                    axes[index].set_xticklabels(new_labels)

        #Check if missing too much data!
        if itest == 0:
            print 'No data to plot, try reading in a file'
            _method_footer_printout()
            return

        #Add colorbar
        axes[8].axis('off')  
        cax = fig.add_axes([0.126, 0.245, 0.775, 0.01])
        cbar = plt.colorbar(im, cax=cax, orientation='horizontal', extend='both')
        cbar.ax.tick_params(labelsize=7)
        cbar.set_label('Brightness Temperature (K)', size=7)

        #Add Aircraft roll angle time series
        if hasattr(self, 'Aircraft_Nav'):
            axes[9].plot(self.Scan,          self.Aircraft_Nav['Roll'], 'b-')
            axes[9].plot(self.Scan,      0.0*self.Aircraft_Nav['Roll'], 'k:')
            axes[9].plot(self.Scan,  5.0+0.0*self.Aircraft_Nav['Roll'], 'k:')
            axes[9].plot(self.Scan, -5.0+0.0*self.Aircraft_Nav['Roll'], 'k:')
            axes[9].set_xlabel('Scan Number or Time (UTC)', fontsize=10)
            axes[9].tick_params(axis='x', labelsize=10, top='off',
                                direction='out')
            axes[9].set_ylabel('Aircraft Roll (deg)', fontsize=10)
            axes[9].tick_params(axis='y', labelsize=7)
            axes[9].set_xlim(xran)
            axes[9].set_ylim([-10,10])
            axes[9].yaxis.set_ticks([-10,0,10])

        #Save the plot and clean up
        if save != None:
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
        print 'calc_polarization():'

        if self.Year[0] < 2011:
            print 'Pre-2011, AMPR only had one channel per frequency.'
            print 'Thus, PyAMPR cannot deconvolve polarization for'
            print 'this project\'s data. Sorry!'
            _method_footer_printout()
            return
        
        for chan in chan_list:
            if hasattr(self, 'TB'+chan+'A') and hasattr(self, 'TB'+chan+'B'):

                print 'Calculating for', chan, 'GHz channel'
                #Use dummy variables and setattr to get the right-sized arrays
                T1 = 1.0 * getattr(self, 'TB'+chan+'A')
                T2 = 1.0 * getattr(self, 'TB'+chan+'B')
                
                #Get angular argument and convert to radians.
                #angle = np.deg2rad(data.Incidence_Angle - 45.0);
                angle = np.deg2rad(np.linspace(self.swath_left,
                                   -1.0*self.swath_left, self.swath_size) - 45.0)

                #There appears to be an offset between the mixed-pol
                #brightness temperatures at 0-deg incidence angle.
                #Theoretically these should be the same.
                T2 = self._compute_nadir_offset_and_apply_if_desired(angle,
                                                T1, T2, force_match, chan)

                if simple:
                    tbv, tbh = \
                    self._solve_using_simple_linear_substitution(angle, T1, T2)

                #Solve using Tikhonov regularization
                else:
                    tbv, tbh = \
                    self._solve_using_constrained_linear_inversion(angle, T1, T2)

                #Finalize V & H attributes and clean up
                setattr(self, 'TB'+chan+'V', tbv)
                setattr(self, 'TB'+chan+'H', tbh)

            else:
                print 'TB'+chan,'does not have both A and B channels,',\
                      'read in a file to obtain.'
                print 'Scene H & V not produced for this channel'

        print time.time() - begin_time, 'seconds to calculate H & V'
        print 'If successful, following attributes are now available:'
        for chan in chan_list:
            print 'TB'+chan+'H', 'TB'+chan+'V',
        print
        _method_footer_printout()

    #########################################

    def write_ampr_kmz(self, var=DEFAULT_VAR, latrange=None, lonrange=None,
                       clevs=DEFAULT_CLEVS, cmap=None, timerange=None,
                       file_path=None, file_name=None, scanrange=None,
                       show_legend=True, equator=False, maneuver=True):
    
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
        """
        plt.close() #mpl seems buggy if multiple windows are left open
        _method_header_printout()
        print 'write_ampr_kmz():'
    
        #10 GHz (A) channel plotted by default
        if var == None or isinstance(var, str) == False:
            var = DEFAULT_VAR
    
        #Check to make sure data exist!
        if not hasattr(self, 'TB'+var.upper()):
            self._missing_channel_printout()
            _method_footer_printout()
            return

        #Adjustable scan range limits
        #Fairly robust - will go down to a width of 25 scans before plotting
        #artifacts begin to occur.
        ind1, ind2 = self._get_scan_indices(scanrange, timerange)
        plon, plat, zdata = self._get_data_subsection(var, ind1, ind2, maneuver)
        plon, plat, zdata = self._filter_bad_geolocations(plon, plat,
                                                         zdata, equator)
        enough_data = self._check_for_enough_data_to_plot(plon, plat)
        if not enough_data:
            return

        latrange, lonrange = self._get_latrange_lonrange(plat, plon,
                                                     latrange, lonrange)
        times = self._get_timestamps_for_gearth(ind1, ind2)

        #Set file info
        if file_path == None:
            file_path = ''
        if file_name == None:
            file_name = self._get_gearth_file_name(var, ind1, '.kmz')

        #Google Earth image production
        fig, ax = gearth_fig(np.min(lonrange), np.min(latrange),
                             np.max(lonrange), np.max(latrange))
        cmap = self._get_colormap(cmap, CMAP_FLAG)
        cs = ax.pcolormesh(plon, plat, zdata,
                           vmin=np.min(clevs), vmax=np.max(clevs), cmap=cmap)
        ax.set_axis_off()
        fig.savefig('overlay.png', transparent=True, format='png')

        #Now we convert to KMZ
        if show_legend == True:
            fig = plt.figure(figsize=(1.0, 4.0), facecolor=None, frameon=False)
            ax = fig.add_axes([0.0, 0.05, 0.2, 0.9])
            cb = fig.colorbar(cs, cax=ax)
            cbytick_obj = plt.getp(cb.ax.axes, 'yticklabels')
            plt.setp(cbytick_obj, color='w', weight='bold')
            ptitle = self._get_ampr_title(var)
            cb.set_label(ptitle+'TB [K]', rotation=-90, color='w',
                         labelpad=20, weight='bold')
            fig.savefig('legend.png', transparent=True, format='png')
            make_kml(np.min(lonrange), np.min(latrange), np.max(lonrange),
                     np.max(latrange), figs = ['overlay.png'],
                     kmzfile=file_path+file_name, colorbar='legend.png',
                     times=times)
        else:
            make_kml(np.min(lonrange), np.min(latrange), np.max(lonrange),
                     np.max(latrange), figs = ['overlay.png'],
                     kmzfile=file_path+file_name, times=times)
    
        print 'Google Earth image hopefully written to:', file_path+file_name
        _method_footer_printout()
        
    ################################################################
    #Internal methods below here. Average user can stop reading now.
    ################################################################

    def _read_level2b_netcdf(self, inputFile, project=DEFAULT_PROJECT_NAME):
        """Internal method to read new L2 netCDF-format AMPR data files"""
        # Open the data
        level2b = Dataset(inputFile , mode="r")
    
        # Set bad data and nav_size
        self.bad_data = DEFAULT_BAD_DATA
        self.nav_size = DEFAULT_NAV_SIZE

        #Assigning project name (self.Project) based on user input
        self._assign_project_name(project)
    
        # Check for navigation
        if 'lat' in level2b.variables.keys():
            print "Found Navigation Data!"
            self.hasNav = True
        else:
            self.hasNav = False

        # Handle the times.
        # Convert to a datetime object
        self.netcdfTimes = level2b.variables['time']
        self.dateTimes = num2date(self.netcdfTimes[:], self.netcdfTimes.units)
        
        # Add the Scan and Swath information.
        self.nscans = len(level2b.dimensions['scan_number'])
        self.ncross = len(level2b.dimensions['scan_position'])
        self.Scan = level2b.variables['scan_number'][:]
        self.Scan_Position = level2b.variables['scan_position'][:]
        self.swath_angle = level2b.variables['scan_angle'][:]
        self.swath_left = self.swath_angle[0]
        self.swath_size = self.ncross
        
        # Before looping through the datetime object, initialize the numpy arrays
        self.Year = np.zeros(self.nscans, dtype=np.int32)
        self.Month = np.zeros(self.nscans, dtype=np.int32)
        self.Day = np.zeros(self.nscans, dtype=np.int32)
        self.Day_of_Year = np.zeros(self.nscans, dtype=np.int32)
        self.Hour = np.zeros(self.nscans, dtype=np.int32)
        self.Minute = np.zeros(self.nscans, dtype=np.int32)
        self.Second = np.zeros(self.nscans, dtype=np.int32)
        self.Second_of_Day = np.zeros(self.nscans, dtype=np.int32)
        self.Time_String = np.zeros(self.nscans, dtype='a8')
        
        for icount in np.arange(self.nscans):
            self.Year[icount] = np.int32(self.dateTimes[icount].year)
            self.Month[icount] = np.int32(self.dateTimes[icount].month)
            self.Day[icount] = np.int32(self.dateTimes[icount].day)
            self.Day_of_Year[icount] = \
                             np.int32(self.dateTimes[icount].timetuple().tm_yday)
            self.Hour[icount] = np.int32(self.dateTimes[icount].hour)
            self.Minute[icount] = np.int32(self.dateTimes[icount].minute)
            self.Second[icount] = np.int32(self.dateTimes[icount].second)
            # Create a Time String and get Second of Data
            self.Time_String[icount], self.Second_of_Day[icount] = \
                _get_timestring_and_sod(self.Hour[icount], self.Minute[icount],
                                        self.Second[icount])
            
        # Compute Epoch Time
        self._fill_epoch_time() #defines and populates self.Epoch_Time
            
        if self.hasNav:
            # Add the geomedateTimes information
            self.Latitude = level2b.variables['lat'][:,:]
            self.Longitude = level2b.variables['lon'][:,:]
            self.Incidence_Angle = level2b.variables['incidence_angle'][:,:]
            self.Relative_Azimuth = level2b.variables['relative_azimuth'][:,:]
        
        # Add the Calibrated TBs
        self.TB10A = level2b.variables['tbs_10a'][:,:]
        self.TB10B = level2b.variables['tbs_10b'][:,:]
        self.TB19A = level2b.variables['tbs_19a'][:,:]
        self.TB19B = level2b.variables['tbs_19b'][:,:]
        #37 GHz A & B accidentally swapped during IPHEX
        if self.Project == 'IPHEX':
            self.TB37A = level2b.variables['tbs_37b'][:,:]
            self.TB37B = level2b.variables['tbs_37a'][:,:]
        else:
            self.TB37A = level2b.variables['tbs_37a'][:,:]
            self.TB37B = level2b.variables['tbs_37b'][:,:]
        self.TB85A = level2b.variables['tbs_85a'][:,:]
        self.TB85B = level2b.variables['tbs_85b'][:,:]

        # Other Variables -- Set to bad data for now.
        self.Icon = self.bad_data * np.ones(self.nscans, dtype = np.int32)
        self.Noise10 = self.bad_data * np.ones(self.nscans, dtype = 'float')
        self.Noise19 = self.bad_data * np.ones(self.nscans, dtype = 'float')
        self.Noise37 = self.bad_data * np.ones(self.nscans, dtype = 'float')
        self.Noise85 = self.bad_data * np.ones(self.nscans, dtype = 'float')
        self.Land_Fraction10 = self.bad_data * np.ones((self.nscans,
                                      self.swath_size), dtype = 'float')
        self.Land_Fraction37 = self.bad_data * np.ones((self.nscans,
                                      self.swath_size), dtype = 'float')
        self.Land_Fraction85 = self.bad_data * np.ones((self.nscans,
                                      self.swath_size), dtype = 'float')
        self.Elevation = self.bad_data * np.ones((self.nscans,
                                      self.swath_size), dtype = 'float')
     
        self._initialize_aircraft_dict()
        for var in self.Aircraft_varlist:
            self.Aircraft_Nav[var][:] = self.bad_data
        if self.hasNav:
            # Now, return a structured array of Aircraft_Nav parameters.
            print np.shape(level2b.variables['gLat'])
            self.Aircraft_Nav['GPS Latitude'] = level2b.variables['gLat'][:]
            self.Aircraft_Nav['GPS Longitude'] = level2b.variables['gLon'][:]
            self.Aircraft_Nav['GPS Altitude'] = level2b.variables['gAlt'][:]
            self.Aircraft_Nav['Pitch'] = level2b.variables['pitch'][:]
            self.Aircraft_Nav['Roll'] = level2b.variables['roll'][:]
            self.Aircraft_Nav['Yaw'] = level2b.variables['track_angle'][:]
            self.Aircraft_Nav['Heading'] =  level2b.variables['heading'][:]
            self.Aircraft_Nav['Ground Speed'] = \
                                  level2b.variables['groundSpeed'][:]
            self.Aircraft_Nav['Air Speed'] = level2b.variables['airSpeed'][:]
            self.Aircraft_Nav['Static Pressure'] = \
                                  level2b.variables['staticPressure'][:]
            self.Aircraft_Nav['Total Temperature'] = \
                                  level2b.variables['totalTemp'][:]
            self.Aircraft_Nav['Wind Speed'] = \
                                  level2b.variables['iWindSpeed'][:]
            self.Aircraft_Nav['Wind Direction'] = \
                                  level2b.variables['iWindDir'][:]
            self.Aircraft_Nav['INS Latitude'] = level2b.variables['iLat'][:]
            self.Aircraft_Nav['INS Longitude'] = level2b.variables['iLon'][:]

    #########################################

    def _read_level2b_ascii(self, full_path_and_filename,
                             project=DEFAULT_PROJECT_NAME):
        
        #Following puts entire file contents into self.ampr_string
        read_successful = self._read_ampr_ascii_file(full_path_and_filename)
        if not read_successful:
            return
    
        #Assigning project name (self.Project) based on user input
        self._assign_project_name(project)
        
        #Determine number of scans
        self.nscans = np.size(self.ampr_string)
        print 'Number of scans =', self.nscans
        
        #Hard-coding sizes of AMPR arrays and dictionaries
        self._hard_code_ampr_array_sizes_and_other_metadata()
        self._declare_ampr_variables()
    
        #Populate the variables with the file's data.
        #Begin master loop
        for index, line in enumerate(self.ampr_string):
            line_split = line.split()
            
            #Header info
            if self.Project == 'IPHEX' or self.Project == 'MC3E':
                #2011+, header info placed before TBs
                self._fill_2011on_header_info(line_split, index)
            else:
                #All pre-MC3E projects, aircraft data placed with header
                #info before TBs
                self._fill_pre2011_header_and_aircraft_info(line_split,
                                                            index)
           
            self.Time_String[index], self.Second_of_Day[index] =\
                      _get_timestring_and_sod(self.Hour[index],
                                    self.Minute[index], self.Second[index])
            
            #Get TBs, Latitudes, Longitudes, etc.
            for i in xrange(self.swath_size):
                
                if self.Project == 'IPHEX':
                    self._fill_2011on_ampr_variables(line_split, index, i)
                    #Below are IPHEX-specific fixes
                    #37 GHZ A&B accidentally swapped during IPHEX
                    self.TB37B[index, i] = float(line_split[i+ 9+
                                                 4*self.swath_size])
                    self.TB37A[index, i] = float(line_split[i+ 9+
                                                 5*self.swath_size])
                    #Note: Currently (May 2014) Land_Fraction## is set to
                    #bad data in IPHEX files
                    #Terrain elevation data not recorded, instead incidence
                    #angle is in its position
                    self.Elevation[index, i] = self.bad_data
                    #IPHEX incidence angle currently ignored by PyAMPR
                    #self.Incidence_Angle[index, i] = float(line_split[i
                    #+27+13*self.swath_size])
                    
                elif self.Project == 'MC3E':
                    self._fill_2011on_ampr_variables(line_split, index, i)
                    
                else: #Pre-MC3E projects
                    self._fill_pre2011_ampr_variables(line_split, index, i)
           
            if self.Project == 'IPHEX' or self.Project == 'MC3E':
                #Aircraft_Nav out of order to improve efficiency for 2011+ projects
                self._fill_2011on_aircraft_info(line_split, index)
    
        self._fill_epoch_time() #defines and populates self.Epoch_Time

        # Replace unfilled data.        
        self._set_unfilled_data_to_bad()
        if self.Year[0] < 2011:
            print 'Only A channels available in this project\'s data'

    #########################################

    def _filter_bad_geolocations(self, plon=None, plat=None, zdata=None,
                                 equator=False):

        """
        Internal method to filter bad geolocation data. Called by the following:
        plot_ampr_track()
        write_ampr_kmz()

        """
        
        #Attempt to deal with bad geolocation data
        #(e.g., Latitude/Longitude=bad_data)
        cond1 = np.logical_or(plon < -180, plon > 180)
        cond2 = np.logical_or(plat <  -90, plat >  90)
        condition = np.logical_or(cond1, cond2)
        indices = np.where(condition)
        if np.shape(indices)[1] > 0:
            print np.shape(indices)[1],\
                 'bad geolocation(s) (e.g., '+str(self.bad_data)+'s) exist, attempting correction.'
            zdata = np.delete(zdata, indices[0], axis=0)
            plon  = np.delete( plon, indices[0], axis=0)
            plat  = np.delete( plat, indices[0], axis=0)
            
        #Attempt to deal with bad 0s/-1s in Latitude/Longitude
        if not equator:
            cond1 = np.logical_or(plon ==  0, plat ==  0)
            cond2 = np.logical_or(plon == -1, plat == -1)
            condition = np.logical_or(cond1, cond2)
            indices = np.where(condition)
            if np.shape(indices)[1] > 0:
                print np.shape(indices)[1],\
                      'bad geolocation(s) (0s/-1s) exist, attempting correction.'
                print 'If aircraft crossed Equator or Prime Meridian,'\
                      'try keyword equator=True.'
                zdata = np.delete(zdata, indices[0], axis=0)
                plon  = np.delete( plon, indices[0], axis=0)
                plat  = np.delete( plat, indices[0], axis=0)

        #Attempt to deal with remaining wildly varying
        #Latitude/Longitude in single scan
        plat_max = np.amax(plat, axis=1)
        plat_min = np.amin(plat, axis=1)
        plon_max = np.amax(plon, axis=1)
        plon_min = np.amin(plon, axis=1)
        cond1 = plat_max - plat_min > 5 #Usually 2 deg is more than sufficient,
        cond2 = plon_max - plon_min > 5 #being extra safe here.
        condition = np.logical_or(cond1, cond2)
        indices=np.where(condition)
        if np.shape(indices)[1] > 0:
            print np.shape(indices)[1],\
                  'scan(s) with high intra-scan geolocation variance exist,',\
                  'attempting correction.'
            zdata = np.delete(zdata, indices[0], axis=0)
            plon  = np.delete( plon, indices[0], axis=0)
            plat  = np.delete( plat, indices[0], axis=0)
        return plon, plat, zdata
    
    #########################################

    def _get_scan_indices(self, scanrange=None, timerange=None):
    
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
        print 'Available scans =', np.min(self.Scan), 'to', np.max(self.Scan)
        print 'Available times =', self.Time_String[0], '-',\
                                   self.Time_String[self.nscans-1]
        if not scanrange and not timerange:
            ind1, ind2 = self._get_min_max_indices()
   
        elif scanrange != None and timerange == None:
            indices = np.where(self.Scan == np.min(scanrange))
            if np.shape(indices[0])[0] == 0:
                print 'Scan number too small, using first scan for beginning'
                ind1 = 0
            else:
                ind1 = indices[0][0]
            indices=np.where(self.Scan == np.max(scanrange))
            if np.shape(indices[0])[0] == 0:
                print 'Scan number too high, using last scan for end'
                ind2 = self.nscans
            else:
                ind2 = indices[0][0]

        else:
            try:
                t1 = _get_sod(float(timerange[0][0:2]), float(timerange[0][3:5]),
                              float(timerange[0][6:8]))
                t2 = _get_sod(float(timerange[1][0:2]), float(timerange[1][3:5]),
                              float(timerange[1][6:8]))
                cond1 = self.Second_of_Day >= np.min([t1, t2])
                cond2 = self.Second_of_Day <= np.max([t1, t2])
                condition = np.logical_and(cond1, cond2)
                indices = np.where(condition)
                if np.shape(indices[0])[0] > 0:
                    ind1 = np.min(indices[0])
                    ind2 = np.max(indices[0])
                else:
                    _print_times_not_valid()
                    ind1, ind2 = self._get_min_max_indices()
            except:
                _print_times_not_valid()
                ind1, ind2 = self._get_min_max_indices()

        return ind1, ind2

    #########################################

    def _declare_ampr_variables(self):
    
        """Define variables to be populated"""
        
        #Timing, Icon, and Noise
        self.Scan =          np.zeros(self.nscans, dtype = np.int32)
        self.Year =          np.zeros(self.nscans, dtype = np.int32)
        self.Month =         np.zeros(self.nscans, dtype = np.int32)
        self.Day =           np.zeros(self.nscans, dtype = np.int32)
        self.Day_of_Year =   np.zeros(self.nscans, dtype = np.int32)
        self.Hour =          np.zeros(self.nscans, dtype = np.int32)
        self.Minute =        np.zeros(self.nscans, dtype = np.int32)
        self.Second =        np.zeros(self.nscans, dtype = np.int32)
        self.Icon =          np.zeros(self.nscans, dtype = np.int32)
        self.Second_of_Day = np.zeros(self.nscans, dtype = np.int32)
        self.Time_String =   np.zeros(self.nscans, dtype = 'a8')
        self.Noise10 =       np.zeros(self.nscans, dtype = 'float')
        self.Noise19 =       np.zeros(self.nscans, dtype = 'float')
        self.Noise37 =       np.zeros(self.nscans, dtype = 'float')
        self.Noise85 =       np.zeros(self.nscans, dtype = 'float')

        #Geolocation and land information
        self.Latitude  =       np.zeros((self.nscans, self.swath_size),
                                        dtype = 'float')
        self.Longitude =       np.zeros((self.nscans, self.swath_size),
                                        dtype = 'float')
        self.Land_Fraction10 = np.zeros((self.nscans, self.swath_size),
                                        dtype = 'float')
        self.Land_Fraction37 = np.zeros((self.nscans, self.swath_size),
                                        dtype = 'float')
        self.Land_Fraction85 = np.zeros((self.nscans, self.swath_size),
                                        dtype = 'float')
        self.Elevation =       np.zeros((self.nscans, self.swath_size),
                                        dtype = 'float')

        #Brightness Temperatures
        self.TB10A = np.zeros((self.nscans, self.swath_size), dtype = 'float')
        self.TB10B = np.zeros((self.nscans, self.swath_size), dtype = 'float')
        self.TB19A = np.zeros((self.nscans, self.swath_size), dtype = 'float')
        self.TB19B = np.zeros((self.nscans, self.swath_size), dtype = 'float')
        self.TB37A = np.zeros((self.nscans, self.swath_size), dtype = 'float')
        self.TB37B = np.zeros((self.nscans, self.swath_size), dtype = 'float')
        self.TB85A = np.zeros((self.nscans, self.swath_size), dtype = 'float')
        self.TB85B = np.zeros((self.nscans, self.swath_size), dtype = 'float')
        #Aircraft Nav
        self._initialize_aircraft_dict()

    #########################################

    def _initialize_aircraft_dict(self):
        """Aircraft Navgation info"""
        self.Aircraft_varlist = ['GPS Latitude', 'GPS Longitude',
                                 'GPS Altitude', 'Pitch', 'Roll',
                                 'Yaw', 'Heading', 'Ground Speed',
                                 'Air Speed', 'Static Pressure',
                                 'Total Pressure', 'Total Temperature',
                                 'Static Temperature', 'Wind Speed',
                                 'Wind Direction', 'INS Latitude',
                                 'INS Longitude', 'INS Altitude']
        dt = []
        for i in xrange(self.nav_size):
            dt.append((self.Aircraft_varlist[i], 'float'))
        self.Aircraft_Nav = np.zeros(self.nscans, dtype = dt)
    
    #########################################

    def _fill_epoch_time(self):
        
        self.Epoch_Time = 0 * self.Second
        for i in xrange(self.nscans):
            self.Epoch_Time[i] = calendar.timegm((int(self.Year[i]),
                                 int(self.Month[i]), int(self.Day[i]),
                                 int(self.Hour[i]), int(self.Minute[i]),
                                 int(self.Second[i])))

    #########################################

    def _set_unfilled_data_to_bad(self):
    
        """
        Vectorizing remaining data assignments for unused/duplicate variables
        Project dependent which variables are or are not used
        """
        if self.Project == 'IPHEX' or self.Project == 'MC3E':
            self.Noise10[:] = self.bad_data
            self.Noise19[:] = self.bad_data
            self.Noise37[:] = self.bad_data
            self.Noise85[:] = self.bad_data
        
        else:
            self.Aircraft_Nav['Static Pressure'   ][:] = self.bad_data
            self.Aircraft_Nav['Total Pressure'    ][:] = self.bad_data
            self.Aircraft_Nav['Static Temperature'][:] = self.bad_data
            self.Aircraft_Nav['Total Temperature' ][:] = self.bad_data
            self.Aircraft_Nav['Wind Speed'        ][:] = self.bad_data
            self.Aircraft_Nav['Wind Direction'    ][:] = self.bad_data
            self.Aircraft_Nav['INS Latitude'      ][:] = self.bad_data
            self.Aircraft_Nav['INS Longitude'     ][:] = self.bad_data
            self.Aircraft_Nav['INS Altitude'      ][:] = self.bad_data
            self.TB10B[:,:] = self.bad_data
            self.TB19B[:,:] = self.bad_data
            self.TB37B[:,:] = self.bad_data
            self.TB85B[:,:] = self.bad_data
            #Only one Land_Fraction variable with old projects
            self.Land_Fraction37[:,:] = self.bad_data
            self.Land_Fraction85[:,:] = self.bad_data

        #Attempting to control for infinite numbers.
        #They blow up plot_ampr_track() and write_ampr_kmz().
        self.TB10A[np.isinf(self.TB10A) == True] = self.bad_data
        self.TB10B[np.isinf(self.TB10A) == True] = self.bad_data
        self.TB19A[np.isinf(self.TB19A) == True] = self.bad_data
        self.TB19B[np.isinf(self.TB19B) == True] = self.bad_data
        self.TB37A[np.isinf(self.TB37A) == True] = self.bad_data
        self.TB37B[np.isinf(self.TB37B) == True] = self.bad_data
        self.TB85A[np.isinf(self.TB85A) == True] = self.bad_data
        self.TB85B[np.isinf(self.TB85B) == True] = self.bad_data
        self.Latitude [np.isinf(self.Latitude)  == True] = self.bad_data
        self.Longitude[np.isinf(self.Longitude) == True] = self.bad_data
        #Attempting to control for NaNs.
        #They blow up plot_ampr_track() and write_ampr_kmz().
        self.TB10A[np.isnan(self.TB10A) == True] = self.bad_data
        self.TB10B[np.isnan(self.TB10A) == True] = self.bad_data
        self.TB19A[np.isnan(self.TB19A) == True] = self.bad_data
        self.TB19B[np.isnan(self.TB19B) == True] = self.bad_data
        self.TB37A[np.isnan(self.TB37A) == True] = self.bad_data
        self.TB37B[np.isnan(self.TB37B) == True] = self.bad_data
        self.TB85A[np.isnan(self.TB85A) == True] = self.bad_data
        self.TB85B[np.isnan(self.TB85B) == True] = self.bad_data
        self.Latitude [np.isnan(self.Latitude)  == True] = self.bad_data
        self.Longitude[np.isnan(self.Longitude) == True] = self.bad_data

    #########################################

    def _get_gearth_file_name(self, var=None, index=0, suffix=None):
    
        mo, dy = self._get_month_and_day_string(index)
        timestamp = self.Time_String[index].replace(':', '')
        return str(self.Year[index]) + mo + dy + '_' +\
               timestamp + 'z_TB' + var.upper() + suffix

    #########################################

    def _get_data_subsection(self, var=None, ind1=None, ind2=None,
                             maneuver=True):
    
        zdata = 1.0 * getattr(self, 'TB'+var.upper())
        zdata = zdata[ind1:ind2]
        plon =  1.0 * getattr(self, 'Longitude')
        plon = plon[ind1:ind2]
        plat =  1.0 * getattr(self, 'Latitude')
        plat = plat[ind1:ind2]
        if not maneuver:
            print 'Filtering out significant aircraft maneuvers'
            roll = self.Aircraft_Nav['Roll'][ind1:ind2]
            indices = np.where(np.abs(roll) >= 5)
            if np.shape(indices)[1] > 0:
                zdata = np.delete(zdata, indices[0], axis=0)
                plon  = np.delete( plon, indices[0], axis=0)
                plat  = np.delete( plat, indices[0], axis=0)
        return plon, plat, zdata

    #########################################
    
    def _fill_2011on_header_info(self, line_split=None, index=None):
    
        self.Scan[index] =        long(line_split[0])
        self.Year[index] =        long(line_split[1])
        self.Month[index] =       long(line_split[2])
        self.Day[index] =         long(line_split[3])
        self.Day_of_Year[index] = long(line_split[4])
        self.Hour[index] =        long(line_split[5])
        self.Minute[index] =      long(line_split[6])
        self.Second[index] =      long(line_split[7])
        self.Icon[index] =        long(line_split[8])

    #########################################

    def _fill_pre2011_header_and_aircraft_info(self, line_split=None,
                                               index=None):
        
        #Begin theatrics to handle Year appearing but once in first line of file,
        #but in different positions depending on project.
        #In JAX90, COARE, & CAPE files it wipes out Day_of_Year[0].
        #In other files it occupies the place of Scan[0].
        self.Scan[index] = index + 1
        if (self.Project == 'JAX90' or self.Project == 'CAPE' or\
                                       self.Project == 'COARE') and index == 0:
            self.Year[:] = long(line_split[1])
        if self.Project != 'JAX90' and self.Project != 'CAPE' and\
                                       self.Project != 'COARE' and index == 0:
            self.Year[:] = long(line_split[0])
        self.Day_of_Year[index] = long(line_split[1])
        doy = int(line_split[1]) - 1
        sdate = str(datetime.date(self.Year[index], 1, 1) +
                    datetime.timedelta(doy))
        self.Month[index] = long(sdate[5:7])
        self.Day[index] =   long(sdate[8:10])
        if (self.Project == 'JAX90' or self.Project == 'CAPE' or\
                                       self.Project == 'COARE') and index == 1:
            #Assuming date change doesn't occur between first and second lines
            self.Day_of_Year[0] = self.Day_of_Year[1]
            self.Month      [0] = self.Month      [1]
            self.Day        [0] = self.Day        [1]
        #End theatrics - WHEW!
                
        self.Hour[index] =   long(line_split[2])
        self.Minute[index] = long(line_split[3])
        self.Second[index] = long(line_split[4])
        self.Icon[index] =   long(line_split[5])
        self.Aircraft_Nav['GPS Latitude' ][index] = float(line_split[6])
        self.Aircraft_Nav['GPS Longitude'][index] = float(line_split[7])
        self.Aircraft_Nav['GPS Altitude' ][index] = float(line_split[8])
        self.Aircraft_Nav['Pitch'        ][index] = float(line_split[9])
        self.Aircraft_Nav['Roll'         ][index] = float(line_split[10])
        self.Aircraft_Nav['Yaw'          ][index] = float(line_split[11])
        self.Aircraft_Nav['Heading'      ][index] = float(line_split[12])
        self.Aircraft_Nav['Air Speed'    ][index] = float(line_split[13])
        self.Aircraft_Nav['Ground Speed' ][index] = float(line_split[14])
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
        
        self.TB10A[index, i] = float(line_split[i+ 9                   ])
        self.TB10B[index, i] = float(line_split[i+ 9+   self.swath_size])
        self.TB19A[index, i] = float(line_split[i+ 9+ 2*self.swath_size])
        self.TB19B[index, i] = float(line_split[i+ 9+ 3*self.swath_size])
        self.TB37A[index, i] = float(line_split[i+ 9+ 4*self.swath_size])
        self.TB37B[index, i] = float(line_split[i+ 9+ 5*self.swath_size])
        self.TB85A[index, i] = float(line_split[i+ 9+ 6*self.swath_size])
        self.TB85B[index, i] = float(line_split[i+ 9+ 7*self.swath_size])
        self.Latitude [index, i] = float(line_split[i+ 9+ 8*self.swath_size])
        self.Longitude[index, i] = float(line_split[i+ 9+ 9*self.swath_size])
        self.Land_Fraction10[index, i] = float(line_split[i+27+10*
                                               self.swath_size])
        self.Land_Fraction37[index, i] = float(line_split[i+27+11*
                                               self.swath_size])
        self.Land_Fraction85[index, i] = float(line_split[i+27+12*
                                               self.swath_size])
        self.Elevation      [index, i] = float(line_split[i+27+13*
                                               self.swath_size])

    #########################################

    def _fill_pre2011_ampr_variables(self, line_split=None, index=None, i=None):
        """
        Data written out left to right; i=0 is Left edge, 
        i=49 is Right edge
        """
        self.TB10A    [index, i] = float(line_split[i+ 19                   ])
        self.TB19A    [index, i] = float(line_split[i+ 19+   self.swath_size])
        self.TB37A    [index, i] = float(line_split[i+ 19+ 2*self.swath_size])
        self.TB85A    [index, i] = float(line_split[i+ 19+ 3*self.swath_size])
        self.Latitude [index, i] = float(line_split[i+ 19+ 4*self.swath_size])
        self.Longitude[index, i] = float(line_split[i+ 19+ 5*self.swath_size])
        self.Elevation[index, i] = float(line_split[i+ 19+ 6*self.swath_size])
        self.Land_Fraction10[index, i] = float(line_split[i+19+ 7*
                                               self.swath_size])

    #########################################

    def _check_for_enough_data_to_plot(self, plon=None, plat=None):
        
        cond1 = np.logical_and(plon >= -180, plon <= 180)
        cond2 = np.logical_and(plat >=  -90, plat <=  90)
        condition = np.logical_and(cond1, cond2)
        indices = np.where(condition)
        if np.shape(indices)[1] < 100:
            print np.shape(indices)[1], 'good gelocations,',\
            'need 100+ (i.e., 2+ scans).'
            print 'Not enough good geolocation data to plot, returning.'
            return False
        else:
            return True

    #########################################

    def _check_aspect_ratio(self, latrange=None, lonrange=None):
    
        aspect_ratio = (float(np.max(latrange)) - float(np.min(latrange))) /\
                       (float(np.max(lonrange)) - float(np.min(lonrange)))
        if aspect_ratio < 0.5 or aspect_ratio > 2:
            print 'Warning: Your aspect ratio choice could lead to poor',\
            'plotting results.'
            print 'Best results are obtained when latrange ~ lonrange.'

    #########################################

    def _hard_code_ampr_array_sizes_and_other_metadata(self):
        
        self.swath_size =  DEFAULT_SWATH_SIZE
        self.nav_size =    DEFAULT_NAV_SIZE
        self.bad_data =    DEFAULT_BAD_DATA
        self.swath_left  = DEFAULT_SWATH_LEFT
        self.swath_angle = self.swath_left - (2.0 * self.swath_left /\
                           float(self.swath_size - 1.0)) *\
                           np.arange(self.swath_size)

    #########################################

    def _fill_2011on_aircraft_info(self, line_split, index):
    
        aircraft_i = 0L
        for name in self.Aircraft_varlist:
            self.Aircraft_Nav[name][index] = \
                    float(line_split[aircraft_i + 9 + 10 * self.swath_size])
            aircraft_i += 1

    #########################################

    def _get_list_of_channels_to_plot(self, show_pol=False):
    
        if show_pol:
            tb_list=['10V', '10H', '19V', '19H', '37V', '37H', '85V', '85H']
            if not hasattr(self, 'TB10V') or not hasattr(self, 'TB19V') or not\
                   hasattr(self, 'TB37V') or not hasattr(self, 'TB85V'):
                print 'Missing some pol channels, trying calc_polarization() before plot'
                self.calc_polarization()
        else:
            tb_list=['10A', '10B', '19A', '19B', '37A', '37B', '85A', '85B']
        return tb_list

    #########################################

    def _missing_channel_printout(self):
    
        print 'Channel doesn\'t exist, check typing or try reading in a file'
        print 'Acceptable channels = 10A, 10B, 19A, 19B, 37A, 37B, 85A, 85B'
        print 'If calculated, also = 10H, 10V, 19H, 19V, 37H, 37V, 85H, 85V'

    #########################################

    def _get_latrange_lonrange(self, plat=None, plon=None,
                               latrange=None, lonrange=None):
    
        if latrange == None:
            latrange = [np.min(plat), np.max(plat)]
        if lonrange == None:
            lonrange = [np.min(plon), np.max(plon)]
        return latrange, lonrange

    #########################################

    def _get_colormap(self, cmap, flag):
    
        if cmap == None:
            if flag:
                cmap = amprTB_cmap
            else:
                cmap = cm.GMT_wysiwyg
        return cmap

    #########################################

    def _get_date_string(self, index=0):

        return str(self.Month[index])+'/'+str(self.Day[index])+'/'+\
               str(self.Year[index])

    #########################################

    def _get_ampr_title(self, var=None):
    
        return 'AMPR '+var[0:2]+' GHz ('+var[2].upper()+') '

    #########################################
    
    def _read_ampr_ascii_file(self, full_path_and_filename):

        if full_path_and_filename[-3:] == '.gz':
            try:
                fileobj = gzip.open(full_path_and_filename)
                ascii = codecs.getreader('ASCII')
                fileobj = ascii(fileobj)
            except:
                print 'Incorrect file or file doesn\'t exist, returning'
                _method_footer_printout()
                return False
        else:
            try:
                fileobj = open(full_path_and_filename, 'r')
            except:
                print 'Incorrect file or file doesn\'t exist, returning'
                _method_footer_printout()
                return False
        try:
            contents = fileobj.readlines()
        except:
            print 'File not ASCII format, returning'
            _method_footer_printout()
            fileobj.close()
            return False
        fileobj.close()
        self.ampr_string = np.array(contents)
        return True

    #########################################
    
    def _assign_project_name(self, project=DEFAULT_PROJECT_NAME):
    
        if not isinstance(project, str):
            print 'Bad project name, provide actual string'
            print 'Assuming', DEFAULT_PROJECT_NAME, 'data structure.'
            project = DEFAULT_PROJECT_NAME
        else:
            print 'Assuming', project.upper(), 'data structure.'
        print 'Change to proper project if incorrect, otherwise errors',\
              'will occur.'
        print 'Currently available field projects: IPHEX, MC3E, TC4, TCSP,',\
              'JAX90, COARE,'
        print 'CAMEX1, CAMEX2, CAMEX3, CAMEX4, TRMMLBA, KWAJEX, TEFLUNA,',\
              'FIRE3ACE, CAPE'
        print 'Default: project = \''+DEFAULT_PROJECT_NAME+'\''
        self.Project = project.upper()
        
    #########################################
    
    def _get_min_max_indices(self):
    
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
        
        #Set regularization parameter
        gam = 10.0**(-1)

        #Solve
        #n.b. I solved this on paper so that a single set of equations
        #     can be used and applied using matrix notation.
        #Define parameters needed.
        a = np.cos(angle)**2
        b = np.sin(angle)**2
        c = np.sin(angle)**2
        d = np.cos(angle)**2
        zeta = 1.0 / ((a**2 + c**2 + gam**2) *\
                      (b**2 + d**2 + gam**2) - (a * b + c * d)**2)
                
        #Get tbv and tbh
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
        #xoff = np.empty(self.nscans) * np.nan
        xoff = np.zeros(self.nscans)
        for iscan in np.arange(self.nscans):
            #Get value.
            x1 = T1[iscan, :]
            x2 = T2[iscan, :]
            #Interpolate x1 and x2 to 0 degrees.
            x1_zero = np.interp(-45.0, scanang, x1)
            x2_zero = np.interp(-45.0, scanang, x2)
            #Compute the offset.
            xoffset = x1_zero - x2_zero
            #Adjust T2 to match T1.
            if force_match:
                T2[iscan, :] = T2[iscan, :] + xoffset
            xoff[iscan] = xoffset
        setattr(self, 'TB'+chan+'_offset', xoff)
        return T2
        
    #########################################
    
    def _check_for_pol_data(self):
    
        if not hasattr(self, 'TB10V') and not hasattr(self, 'TB19V') and not\
               hasattr(self, 'TB37V') and not hasattr(self, 'TB85V'):
            print 'No polarization data to plot! Returning ...'
            return False
        else:
            return True
    
    #########################################

    def _get_timestamps_for_gearth(self, ind1=None, ind2=None):

        #Format example: '1997-07-16T07:30:15Z'
        mo, dy = self._get_month_and_day_string(ind1)
        time1 = str(self.Year[ind1]) + '-' + mo + '-' + dy + 'T' +\
                self.Time_String[ind1] + 'Z'
        mo, dy = self._get_month_and_day_string(ind2-1)
        time2 = str(self.Year[ind2-1]) + '-' + mo + '-' + dy + 'T' +\
                self.Time_String[ind2-1] + 'Z'
        times = [time1, time2]
        return times

    #########################################
    
    def _get_month_and_day_string(self, index):
        
        smo = str(self.Month[index])
        if self.Month[index] < 10:
            smo = '0' + smo
        sdy = str(self.Day[index])
        if self.Day[index] < 10:
            sdy = '0' + sdy
        return smo, sdy
    
    #########################################
    
    #########################################
    
    ######################################
    #Add more attributes and methods here!
    ######################################

    #########################################

#End class AmprTb definition    
###########################################################

#Stand-alone functions

def _get_timestring_and_sod(hour=None, minute=None, second=None):
    """ Time_String: Fill size gaps with 0s """
    d = str(hour)
    if hour < 10:
        d = '0' + d
    e = str(minute)
    if minute < 10:
        e = '0' + e
    f = str(second)
    if second < 10:
        f = '0' + f
    return d + ':' + e + ':' + f, _get_sod(hour, minute, second)

def _get_sod(hour=None, minute=None, second=None):
    return 3600.0 * hour + 60.0 * minute + second

def _method_footer_printout():
    print '********************'
    print

def _method_header_printout():
    print
    print '********************'

def _print_times_not_valid():
    print 'Times not valid, just plotting everything'
    print 'Next time try timerange=[\'hh:mm:ss\',\'HH:MM:SS\']'


