from __future__ import print_function
import numpy as np
import datetime as dt


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
