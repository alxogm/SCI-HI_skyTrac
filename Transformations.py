import numpy as np
import datetime
from sidereal import *

rad = np.pi/180.

def equa_sys(lon, lat, time):
    """
    lon, lat: in radians (pls)
    time: a datetime instance
    This definition calculates the right ascencion. For this particular case (a=90)
    the declination is given by the latitude
    """

    t_sidereal = valt = SiderealTime.fromDatetime(time)
    LSST = SiderealTime.lst(t_sidereal,lon)*rad*14.96
    #times = SiderealTime.lst(t_sidereal,lon)
    return LSST

def gal_sys(lon, lat, time):

    alpha_g=192.85*rad
    delta_g=27.128333*rad
    alpha_c=266.4*rad
    delta_c=-28.929656*rad
    dec = lat
    ra = equa_sys(lon, lat, time)

    B=np.arcsin(np.sin(dec)*np.sin(delta_g)+np.cos(dec)*np.cos(delta_g)*np.cos(ra-alpha_g))
    J=(np.sin(dec)*np.cos(delta_g)-np.cos(dec)*np.sin(delta_g)*np.cos(ra-alpha_g))/np.cos(B)
    K=np.arcsin(np.cos(dec)*np.sin(ra-alpha_g)/np.cos(B))
    Q=np.arccos(np.sin(delta_c)/np.cos(delta_g))

    if J<0.:
        L=Q/rad+K/rad-180.
    else:
        L=(Q/rad-K/rad)

    if L<0.:
        L=L+360.

    return L, B/rad #returns in degrees for some reason
