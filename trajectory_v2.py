#Este código hace transfomaciones de coordenasdas geograficas a coordenadas galacticas
#y las projecta sobre el plano galactico usando healpy

import numpy as np
import datetime
from sidereal import *
import matplotlib.pyplot as plt
import healpy as hp
import seaborn as sb
from Transformations import *
sb.set_style('white')
rad = np.pi/180.

#aqui se definen las horas de obsrvacion y en cuantos intervalos de tiempo (min)
#se esta integrando.
n_H = 24
delta = 10.

#estas son las longitu y latitud de Isla Guadalupe
glon = -118.3011 * rad
glat = 28.9733 * rad

h = datetime.timedelta(minutes=int(delta))
t = datetime.datetime(2013,6,1,23,0,0)

t_vec = int(n_H*60./delta) #number of intervals

H = np.arange(t_vec)

gal_coords=[]
for i in range(len(H)):
    time_int = t + H[i]*h
    g_l, g_b = gal_sys(glon,glat,time_int) #remember, return in degrees.
    g_b = 90. - g_b #colatitude
    gal_coords.append(np.array([g_l,g_b]))

gal_coords=np.transpose(gal_coords)
"""
RA = []
dec = np.zeros(len(H))
for i in range(len(H)):
    time_int1 = t + H[i]*h
    ra = equa_sys(glon, glat, time_int1)
    RA.append(np.array([ra]))
    dec[i] = glat/rad
RA = np.transpose(RA)
"""
# aqui se usa healpy para transformar de coordenadas galacticas a coordenadas
# ecuatoriales.
r = hp.Rotator(coord=['G','C'], deg=False)
theta_c, phi_c = r(gal_coords[1,:]*rad, gal_coords[0,:]*rad)
np.savetxt('cosa_d.txt',np.array([theta_c,phi_c]).T)
"""
data = np.loadtxt('70MHz.dat')
hp.mollview(np.log10(data), coord='G', rot=(gal_coords[0,0],(90.- gal_coords[1,0])))
hp.projplot(gal_coords[1,:]*rad,gal_coords[0,:]*rad,'.k',lonlat=False)

#hp.mollview(np.log10(data), coord=['G','C'], rot=(gal_coords[0,0],(90.-gal_coords[1,0])), flip='geo')
hp.mollview(np.log10(data), coord=['G','C'], rot=(theta_c[0]/rad, (90. - (phi_c[0]/rad))))#, flip='geo')

hp.projplot(theta_c, phi_c ,'.k',lonlat=False,coord='C')

plt.show()

#np.savetxt('cosa_d.txt',np.array([theta_c,phi_c]).T)
"""
