"""
Plot winds data for period 1978-2004
on a north polar stereographic basemap.

The longitude lon_0 is at 6-o'clock, and the latitude circle boundinglat is 
tangent to the edge of the map at lon_0. 
Default value of lat_t (latitude of true scale) is pole

"""

import glob
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

WU=np.zeros([47,45],dtype=np.float)         # 2D array to save LOCAL (X) WIND VELOCITY values   
WV=np.zeros([47,45],dtype=np.float)         # 2D array to save MERIDIONAL (Y) WIND VELOCITY values both for each year averages

b = range(1978,2004)                        # range of years to compute

for i in b:                                 # loop over every year at a time
    pep = i
    pip = i+1
    PEP = str(pep) + '-' + str(pip)
    print 'PEP: ', PEP

    filenames = glob.glob('C:/Users/admin/Desktop/Master_thesis/Sandbox/Data/nc/'+PEP+'_nc/*.nc')          # list of files to compute
      
    Wu=np.zeros([47,45],dtype=np.float)         # 2D array to save LOCAL (X) WIND VELOCITY values 
    Wv=np.zeros([47,45],dtype=np.float)         # 2D array to save MERIDIONAL (Y) WIND VELOCITY values both from opening .nc file
    
    for filename in filenames:                  # loops over every file for a specific year
        data = Dataset(filename, mode='r')      # creates array from netCDF data
        
        print 'filename: ', filename    
        
        uw=data.variables['u_wind_avg'][:]      # extract specific variables from data
        vw=data.variables['v_wind_avg'][:]
        lat1=data.variables['lat1'][:]
        lon1=data.variables['lon1'][:] 

#        uw=np.flipud(uw)
#        vw=np.flipud(vw)
       
        Wu+=uw                                  # cummulative sum of whole year
        Wv+=vw    
    
    Wu=Wu/(len(filenames))                      # yearly average wind 
    Wv=Wv/(len(filenames))

    WU += Wu                                    # cummulative sum of whole period
    WV += Wv

WU = WU/len(b)                                  # period average wind 
WV = WV/len(b)
    
modulus = np.sqrt(WU**2 + WV**2)                # normalized vectors
U = WU/modulus                      
V = WV/modulus

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

x,y = m(lon1,lat1)          # transforms from spherical til flat

Q = m.quiver(x,y,U,V,modulus,scale=40)#,units='xy',width=8,angles='xy',scale_units='xy',\
            #,headlength=5,headwidth=5,latlon=False)#,cmap=cm.seismic,
cb = m.colorbar()
plt.set_cmap('jet')
cb.set_label('m/s')
#cb.set_clim(np.nanmin(modulus),np.nanmax(modulus))
#qk = plt.quiverkey(Q, 0.1, 0.95, 1.5, 'm/s', labelpos='W')
plt.title("Winds _ " + str(b[0]) + '-' + str(b[-1]+1) + " average")

#plt.show()
