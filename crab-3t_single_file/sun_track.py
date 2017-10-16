import astropy
from astropy.coordinates import AltAz
from astropy.coordinates import get_sun
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt

observing_time = Time('2017-09-25T9:52:00', format = 'isot', scale = 'utc')#use UTC time you want
observing_location = EarthLocation(lat='37d55.1m', lon='-122d9.4m', height=304*u.m) #leuschner location#
AA = AltAz(location=observing_location, obstime=observing_time)
sun	=	SkyCoord(get_sun(observing_time))
#Change to local Az and Alt#
LA_sun_at_obs = sun.transform_to(AA)
print LA_sun_at_obs
print '#######Postion for source in Observatory########'
print 'Azimuth:',LA_sun_at_obs.az,'Altitude:',LA_sun_at_obs.alt
