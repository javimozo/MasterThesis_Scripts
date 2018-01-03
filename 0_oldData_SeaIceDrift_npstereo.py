"""
Created on Mon May 23 11:17:43 2016

Plot sea ice drift on a north polar stereographic basemap.

The longitude lon_0 is at 6-o'clock, and the latitude circle boundinglat is 
tangent to the edge of the map at lon_0. 
Default value of lat_t (latitude of true scale) is pole
"""

import glob
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

b =range(1978,2004)                         # range of years to compute

UU=np.zeros([2115,],dtype=np.float)         # 2D array to save LOCAL (X) ICE VELOCITY values
VV=np.zeros([2115,],dtype=np.float)         # 2D array to save MERIDIONAL (X) ICE VELOCITY values both for each year averages

UU[:,] = np.nan
VV[:,] = np.nan

filenames1 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Sandbox/Data/nc/1978-1979_nc/NWP_nhKwok_aggr_1979010212-1979010412.nc')
    
datum1 = Dataset(filenames1[0], mode='r')   # lat and lon info from .nc file

print 'filename: ', filenames1[0]

lat1=datum1.variables['lat1'][:]
lon1=datum1.variables['lon1'][:] 

m = -1                                      # time counter m initial value

for i in b:                                 # loop over every year at a time
    pep = i
    pip = i+1
    PEP = str(pep) + '-' + str(pip)
    print 'PEP: ', PEP
   
    filenames2 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Sandbox/Data/vec/'+PEP+'_vec/*.vec')
    
    Iu=np.zeros([2115,len(filenames2)],dtype=np.float)          # 2D array to save LOCAL (X) ICE VELOCITY values
    Iv=np.zeros([2115,len(filenames2)],dtype=np.float)          # 2D array to save MERIDIONAL (X) ICE VELOCITY values both from opening .vec file
    
    n = -1                                                      # time counter n initial value

    m += 1                                                      # # time counter m
    
    for filename in filenames2:             # loops over every file for a specific year
        f = file(filename, 'r')
        datum2 = np.genfromtxt(f, skip_header=1, usecols=(4,5,6,7,8,9))     # creates array from tabular data in f
        f.close()
        
        print 'filename: ', filename    
        
        ui=datum2[:,0]                  # extract specific variables from data
        vi=datum2[:,1]
        
        for i in range(len(ui)):
            if datum2[i,2] == 0 and datum2[i,3] == 0 and datum2[i,4] == 0 and datum2[i,5] == 0:
                ui[i] = np.nan
                vi[i] = np.nan            # checks for grid points with Nans
    
        if len(ui) < 2115:              # checks for ice drift data that lacks grid's last row
           a = np.zeros([45,])            # creates row
           a[:] = np.nan                  # fills it w/ Nans
           ui = np.hstack((ui,a))         # stack it to the main block
           vi = np.hstack((vi,a))
           
        n += 1                       # time counter n
        
        Iu[:,n] = ui                # saves LOCAL (X) ICE VELOCITY values in columns 
        Iv[:,n] = vi                # saves MERIDIONAL (Y) ICE VELOCITY values in columns both from opening .vec file

    if m == 0:        
        IU = Iu
        IV = Iv        
    else:
        IU = np.hstack((IU,Iu))
        IV = np.hstack((IV,Iv))      
        
#        uu+=ui       # kms/2.days
#        vv+=vi   

#        uu=(uu/(len(filenames1)))*(1000.0/(2*86400))
#        vv=(vv/(len(filenames1)))*(1000.0/(2*86400))

for ii in range(2115):

    UI = IU[ii,:]            # extracts every row one by one 
    VI = IV[ii,:]            # which represents ice drift time series for each grid point
    k = np.isnan(UI)         # get rid on Nans
    UI = UI[~k]
    VI = VI[~k]
    
    if len(UI) > 15:
        UU[ii] = np.mean(UI)*(1000.0/(2*86400))
        VV[ii] = np.mean(VI)*(1000.0/(2*86400))
    
#UU=(UU/len(b))*(1000.0/(2*86400))     # to convert to m/s 
#VV=(VV/len(b))*(1000.0/(2*86400))     # 1000/(2*86400) = 0.005787037

modulus = np.sqrt(UU**2 + VV**2)
U = UU/modulus                      # normalized vectors
V = VV/modulus

#############################
### Plotting with Basemap ###
#############################

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')
#m.drawparallels(np.arange(-80.,81.,20.))
#m.drawmeridians(np.arange(-180.,181.,20.))
#m.drawmapboundary(fill_color='aqua')

lat1=np.flipud(lat1)
lon1=np.flipud(lon1)
   
x,y = m(lon1,lat1,inverse=False)#,indexing='ij')

Q=m.quiver(x,y,U,V,modulus,scale=40)#,units='xy',headlength=5,headwidth=5,latlon=False,
            #width=8,angles='xy',scale_units='xy',cmap=cm.seismic,
cb = m.colorbar()
cb.set_label('m/s')
#cb.set_clim(np.nanmin(modulus),np.nanmax(modulus))
#qk = plt.quiverkey(Q, 0.1, 0.95, 1.5, 'm/s', labelpos='W')
plt.title("SeaIceDrift _ " + str(b[0]) + '-' + str(b[-1]+1))

#plt.show()
