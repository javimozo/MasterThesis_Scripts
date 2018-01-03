# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 19:14:50 2016

"""

import glob
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

filename1 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Javier_ERAI_winds_icedrift_OSISAF/Marta_ssh/mssh_1993-2002.nc')          
filename4 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Sandbox/Data/nc/1978-1979_nc/NWP_nhKwok_aggr_1979010212-1979010412.nc')          

data1 = Dataset(filename1[0], mode='r')
data4 = Dataset(filename4[0], mode='r')

ssh = data1.variables['ssh'][:]              # extract specific variables from data
lat = data4.variables['lat1'][:]
lon = data4.variables['lon1'][:] 

lat_rad = np.deg2rad(lat)

Ug = np.zeros([47,45],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
Vg = np.zeros([47,45],dtype=np.float)    # 2D array to save LOCAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)

dx = 100000.0
dy = 100000.0
g = 9.81
#f = 2*((2*np.pi)/(24*60*60))*np.sin(np.pi/2)        # f-plane approximation

for i in range(1,46):
    for ii in range(1,44):
        f = 2*((2*np.pi)/(24*60*60))*np.sin(lat_rad[i,ii])         # exact latitude
        ug = - (g/f)*((ssh[i+1,ii]-ssh[i-1,ii])/(2*dy))
        Ug[i,ii] = ug
        vg = (g/f)*((ssh[i,ii+1]-ssh[i,ii-1])/(2*dx))
        Vg[i,ii] = vg

modulus = np.sqrt(Ug**2 + Vg**2)
U = Ug/modulus                      # normalized vectors
V = Vg/modulus
        
#############################
### Plotting with Basemap ###
#############################

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')
#m.drawparallels(np.arange(-80.,81.,10.))
#m.drawmeridians(np.arange(-180.,181.,20.))
#m.drawmapboundary(fill_color='aqua')

#lat=np.flipud(lat)
#lon=np.flipud(lon)
   
x,y = m(lon,lat,inverse=False)#,indexing='ij')

plt.figure(1)

#cs1 = m.contourf(x,y,ssh,100)#,np.arange(0.008,1.3,0.01))
#
#plt.hold(True)

Q = m.quiver(x,y,U,V,modulus,scale=40)#,headlength=5, headwidth=5, latlon=False,
            # units='xy', width=8, angles='xy', scale_units='xy', cmap=cm.seismic,
#qk = plt.quiverkey(Q, 0.15, 0.9, 1, '20 m/s', labelpos='W')
cb = m.colorbar()
cb.set_label('m/s')
plt.set_cmap('jet')
#cb.set_clim(np.nanmin(modulus),np.nanmax(modulus))
plt.title("Geostrophic sea currents / ssh _ 1993-2002")

