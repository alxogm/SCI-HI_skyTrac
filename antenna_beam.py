import numpy as np
import matplotlib.pyplot as plt
import datetime
from sidereal import *
from Transformations import *
import seaborn as sb
import healpy as hp
sb.set_style('white')
rad = np.pi/180.

#antenna pattern
#direct = '/Users/Oleg/Documents/SCI-Bayes/antenna_beam/'
data = np.loadtxt('70MHz.txt')

#this data is in the ISO coordinate system
theta = data[:,0]*rad #sim latitud [-90,90]
phi = data[:,1]*rad #longitude [-90,90]
temp = data[:,2]

#Isla Guadalupe coordinates in geo s.
glon = -118.3011 * rad
glat = 28.9733 * rad
e_theta =(90. - 28.9733) * rad #colatitude or polar angle

#check for some changes in the angles. THE REFERENCE SYSTEM.
#coordinate transformation from spherical (r=1,theta,phi) to cartesian(x,y,z)
x = np.sin(theta) * np.cos(phi)
y = np.sin(phi) * np.sin(theta)
z = np.cos(theta)


X = np.vstack((x,y,z)).T


#def euler matrices
B = np.matrix([[np.cos(glon),np.sin(glon),0.0],[-np.sin(glon),np.cos(glon),0.],\
[0.,0.,1.]]).A

C = np.matrix([[1.,0.,0.],[0.,np.cos(e_theta),np.sin(e_theta)],\
[0.,-np.sin(e_theta),np.cos(e_theta)]]).A

D = np.matrix([[1.,0.,0.,],[0.,1.,0.],[0.,0.,1.]]).A

EM = np.dot(np.dot(B,C),D)

#X_prime = np.dot(EM,X[0])

#transformation X' = EM*X.T[i]
rot_theta=[]
rot_phi=[]
for i in range(len(phi)):
    X_prime = np.dot(EM,X[i])
    rot_theta.append((np.pi/2.) - np.arccos(X_prime[2])) #latitude
    #rot_theta.append((np.pi/2.) - np.arccos(X_prime[2]))
    rot_phi.append(np.arctan2(X_prime[1],X_prime[0])) #maybe reverse!


dec = rot_theta
h = datetime.timedelta(minutes=int(1000))
ref_time = datetime.datetime(2013,6,1,23,0,0)



gal_x=[]
for i in range(len(rot_phi)):
    l_g, b_g = gal_sys(rot_phi[i],rot_theta[i], ref_time)
    b_g = 90. - b_g #colatitude (galactic)
    gal_x.append(np.array([l_g*rad, b_g*rad]))

gal_x = np.transpose(gal_x)

#reference point for the rotation
ref_time1 = datetime.datetime(2013,6,1,23,0,0)
phi_g, theta_g = gal_sys(glon, glat, ref_time1)
phi_g = phi_g*rad
theta_g = (theta_g)*rad
datos1, datos2 = np.loadtxt('cosa_d.txt',unpack=True)

nside = 512/16
pix = hp.ang2pix(nside, gal_x[1,:], gal_x[0,:])
#pix2 = hp.ang2pix(nside, RA, dec)

bmap = np.zeros(hp.nside2npix(nside), dtype=np.double)
bmap[pix] = temp
hp.mollview(bmap, coord=['G','C'], rot=(datos1[0]/rad, 90.-(datos2[0]/rad)))
hp.projplot(datos1, datos2,'.k',lonlat=False,coord='C')
plt.show()
