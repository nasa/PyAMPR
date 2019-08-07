from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import cartopy
from .defaults import DEFAULT_CLEVS, DEFAULT_GRID_DEL


def read_aircraft_nav_into_awot(
        AmprTB, project='OLYMPEX', platform='NASA ER-2', flight_number=None):
    """
    Import the AmprTb aicraft nav into a format suitable for AWOT
    (https://github.com/nguy/AWOT). This simplifies plotting the track
    with time stamps, if the user has AWOT installed.

    Note: AWOT is not required for this function to work
    """

    if not hasattr(AmprTB, 'Aircraft_Nav'):
        print('No aircraft information in argument, failing ...')
        return

    flight = {}
    varlist = ['latitude', 'longitude', 'altitude', 'time']
    for var in varlist:
        flight[var] = {}
    flight['latitude']['data'] = AmprTB.Aircraft_Nav['GPS Latitude']
    flight['longitude']['data'] = AmprTB.Aircraft_Nav['GPS Longitude']
    flight['altitude']['data'] = AmprTB.Aircraft_Nav['GPS Altitude']

    ampr_datetime = []
    for et in AmprTB.Epoch_Time:
        ampr_datetime.append(dt.datetime(1970, 1, 1) +
                             dt.timedelta(seconds=np.float(et)))
    flight['time']['data'] = ampr_datetime

    for var in varlist:
        flight[var]['data'] = np.ma.masked_array(
            flight[var]['data'], mask=False)
    flight['flight_number'] = flight_number
    flight['project'] = project
    flight['platform'] = platform
    flight['Uwind'] = None
    flight['Vwind'] = None
    return flight


class _FourPanelTrack(object):

    def __init__(self, amprtb, **kwargs):
        self.kw_dict = {
            'latrange': None, 'lonrange': None, 'parallels': DEFAULT_GRID_DEL,
            'meridians': DEFAULT_GRID_DEL, 'clevs': DEFAULT_CLEVS,
            'cmap': None, 'save': None, 'show_track': False, 'scanrange': None,
            'show_grid': True, 'equator': False, 'maneuver': True,
            'timerange': None, 'show_qc': False, 'resolution': '50m',
            'projection': cartopy.crs.PlateCarree(), 'verbose': False,
            'wmts_layer': None}
        self.parse_kwargs(**kwargs)
        self.flag = self.construct_plot(amprtb)

    def parse_kwargs(self, **kwargs):
        if 'chan' in kwargs.keys():
            self.chan = kwargs['chan']
        else:
            self.chan = 'A'
        for key in self.kw_dict.keys():
            if key not in kwargs.keys():
                setattr(self, key, self.kw_dict[key])
            else:
                setattr(self, key, kwargs[key])

    def make_title(self, freq, amprtb, ind1, ind2):
        lab_dict = {'10': '(a) ', '19': '(b) ', '37': '(c) ', '85': '(d) '}
        title = lab_dict[freq] + \
            str(amprtb._get_ampr_title(freq+self.chan)) + \
            '\n' + str(amprtb._get_date_string(ind1)) + \
            str(', ') + str(amprtb.Time_String[ind1]) + str('-') + \
            str(amprtb.Time_String[ind2-1]) + str(' UTC')
        return title

    def construct_plot(self, amprtb):
        """
        This makes the actual 4-panel plot using repeated calls to the
        pyampr.AmprTb.plot_ampr_track() method.
        """
        self.fig, [[self.ax1, self.ax2], [self.ax3, self.ax4]] = \
            plt.subplots(2, 2, figsize=(10, 10),
                         subplot_kw={'projection': self.projection})
        ind1, ind2 = amprtb._get_scan_indices(
            self.scanrange, self.timerange, False)

        # 10 GHz plot
        stuff = amprtb.plot_ampr_track(
            var='10'+self.chan, latrange=self.latrange,
            lonrange=self.lonrange, parallels=self.parallels,
            meridians=self.meridians, title='', wmts_layer=self.wmts_layer,
            clevs=self.clevs, cmap=self.cmap, show_track=self.show_track,
            maneuver=self.maneuver, scanrange=self.scanrange,
            show_grid=self.show_grid, equator=self.equator,
            show_qc=self.show_qc, resolution=self.resolution,
            projection=self.projection, ax=self.ax1, fig=self.fig,
            verbose=self.verbose, timerange=self.timerange, return_flag=True)
        self.ax1.set_title(self.make_title('10', amprtb, ind1, ind2))

        # 19 GHz plot
        amprtb.plot_ampr_track(
            var='19'+self.chan, latrange=self.latrange,
            lonrange=self.lonrange, parallels=self.parallels,
            meridians=self.meridians, title='', wmts_layer=self.wmts_layer,
            clevs=self.clevs, cmap=self.cmap, show_track=self.show_track,
            maneuver=self.maneuver, scanrange=self.scanrange,
            show_grid=self.show_grid, equator=self.equator,
            show_qc=self.show_qc, resolution=self.resolution,
            projection=self.projection, ax=self.ax2, fig=self.fig,
            verbose=self.verbose, timerange=self.timerange)
        self.ax2.set_title(self.make_title('19', amprtb, ind1, ind2))

        # 37 GHz plot
        amprtb.plot_ampr_track(
            var='37'+self.chan, latrange=self.latrange,
            lonrange=self.lonrange, parallels=self.parallels,
            meridians=self.meridians, title='', wmts_layer=self.wmts_layer,
            clevs=self.clevs, cmap=self.cmap, show_track=self.show_track,
            maneuver=self.maneuver, scanrange=self.scanrange,
            show_grid=self.show_grid, equator=self.equator,
            show_qc=self.show_qc, resolution=self.resolution,
            projection=self.projection, ax=self.ax3, fig=self.fig,
            verbose=self.verbose, timerange=self.timerange)
        self.ax3.set_title(self.make_title('37', amprtb, ind1, ind2))

        # 85 GHz plot
        amprtb.plot_ampr_track(
            var='85'+self.chan, latrange=self.latrange,
            lonrange=self.lonrange, parallels=self.parallels,
            meridians=self.meridians, title='', wmts_layer=self.wmts_layer,
            clevs=self.clevs, cmap=self.cmap, show_track=self.show_track,
            maneuver=self.maneuver, scanrange=self.scanrange,
            show_grid=self.show_grid, equator=self.equator,
            show_qc=self.show_qc, resolution=self.resolution,
            projection=self.projection, ax=self.ax4, fig=self.fig,
            verbose=self.verbose, timerange=self.timerange)
        self.ax4.set_title(self.make_title('85', amprtb, ind1, ind2))

        # plt.tight_layout()
        return True
