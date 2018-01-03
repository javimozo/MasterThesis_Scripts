# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 23:28:13 2017

@author: admin
"""
"""
Created on Wed Aug 31 22:00:00 2016

"""

import glob
import numpy as np
import numpy.ma as ma
from numpy.linalg import inv
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap

b =range(2012,2015)
B = 'Summer'
#B = 'Winter'
#B = ''

ResR=np.zeros([21063,len(b)],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
ResI=np.zeros([21063,len(b)],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
ResR[:,:] = np.nan
ResI[:,:] = np.nan

m = -1          # initial value of counter for years

for i in b:         # loop over year by year
    PIP = str(i)
    POP = str(i+1)
    PEP = PIP + '-' + POP
    print 'PEP: ', PEP

    m += 1          # time counter for years
    
    filenames1 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Javier_ERAI_winds_icedrift_OSISAF/SortedData_OSISAF_IceWind/WindSpeed_10m/' + PEP + '/' + B + '/*.nc')
    filenames2 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Javier_ERAI_winds_icedrift_OSISAF/SortedData_OSISAF_IceWind/SeaIceDriftSpeed/' + PEP + '/' + B + '/*.nc')
    
    Wu=np.zeros([21063,len(filenames1)],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    Wv=np.zeros([21063,len(filenames1)],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    Iu=ma.zeros([21063,len(filenames2)],dtype=np.float)    # 2D array to save LOCAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)
    Iv=ma.zeros([21063,len(filenames2)],dtype=np.float)    # 2D array to save MERIDIONAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)
        
    n = -1       # initial value of counter for wind files
    nn = -1        # initial value of counter for iceDrift files
    
    for filename in filenames1:                         # extracts data from each .nc file one by one
        datum1 = Dataset(filename, mode='r')
        
        print 'filename: ', filename        
        lat1=datum1.variables['lat1'][:]
        lon1=datum1.variables['lon1'][:]              
        uw=datum1.variables['u_wind_avg'][:]
        vw=datum1.variables['v_wind_avg'][:]
    
        n += 1
               
        Wu[:,n] = uw.flat           
        Wv[:,n] = vw.flat      

    for filename in filenames2:                         # extracts data from each .nc file one by one
        datum2 = Dataset(filename, mode='r')
            
        print 'filename: ', filename        
            
        lat=datum2.variables['lat'][:]
        lon=datum2.variables['lon'][:]
        ui = datum2.variables['dX'][:]
        try:
            vi = datum2.variables['dY_v1p4'][:]
        except KeyError:
            vi = datum2.variables['dY'][:]
        
        nn += 1
        
        ui = ui[0,:,:]*(1000.0/(2*86400))
        vi = vi[0,:,:]*(1000.0/(2*86400))
        
        Iu[:,nn] = ma.ravel(ui)             # saves LOCAL (X) ICE VELOCITY values in columns 
        Iv[:,nn] = ma.ravel(vi)             # saves MERIDIONAL (Y) ICE VELOCITY values in columns    
    
    for element in range(21063):         # loop over every gridpoint
        
        print 'element: ', element
        
        WU = Wu[element,:]              # extracts one by one every row (which represents wind values for each grid point)
        WV = Wv[element,:]        
        IU = Iu[element,:]        
        IV = Iv[element,:]              # extracts one by one every row (which represents ice drift values for each grid point)

        k = ma.getmask(IU)                # checks for nans

        WC = np.array(WU,dtype=np.complex)    
        WC.imag = (WV)
        IC = np.array(IU,dtype=np.complex)
        IC.imag = (IV)

        WC = WC[~k]                     # removes nans
        IC = IC[~k]

        if len(IC) > 5:                # minimum lengh of time series to ???
            ID = np.eye(len(WC))
            IC = IC[:,np.newaxis]                             # transforms the horizontal array into a vertical one
            WC = [WC,np.ones_like(WC,dtype=np.complex)]       # attaches an array of ones to the list of complex velocity values
            WC = np.mat(WC)                                   # converts a list of a 1D list and 1D array into a 2 row matrix
            WC = WC.T                                         # transposes the matrix into a 2 column one (1 x 2)
            WCH = WC.getH()                                   # transforms it into a Hermitian (conjugate transpose) matrix

            Res = (ID - (WC*(inv(WCH*WC))*WCH))*IC
#            Res = np.array(Res)        
#            Res = np.reshape(Res,(len(Res),))        
            ResR[element,m] = np.sqrt(np.mean(np.square(Res.real)))
            ResI[element,m] = np.sqrt(np.mean(np.square(Res.imag)))   

RESr = np.zeros([21063,],dtype=np.float)
RESi = np.zeros([21063,],dtype=np.float)

for ii in range(21063):
    RESr[ii] = np.nanmean(ResR[ii,:])
    RESi[ii] = np.nanmean(ResI[ii,:])

RESr = np.reshape(RESr,(177,119))     # reshapes the A parameter 1D array into 2D for plotting with contourf
RESi = np.reshape(RESi,(177,119))     # reshapes the A parameter 1D array into 2D for plotting with contourf

#############################
### Plotting with Basemap ###
#############################

m = Basemap(projection='npstere',boundinglat=61.5,lat_0=90,lat_ts=70,lon_0=-35,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')
#m.drawparallels(np.arange(-80.,81.,5.))
#m.drawmeridians(np.arange(-180.,181.,20.))
#m.drawmapboundary(fill_color='aqua')

x,y = m(lon1,lat1,inverse=False)#,indexing='ij')

plt.figure(1)
c2 = m.contourf(x,y,RESr,np.arange(0.0,0.11,0.0001))
cb = m.colorbar(ticks=(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11))
cb.set_label('')
plt.set_cmap('jet')
#cb.set_clim(np.nanmin(RESr),np.nanmax(RESr))
plt.title('ResReal _ ' + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=61.5,lat_0=90,lat_ts=70,lon_0=-35,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(2)
c2 = m.contourf(x,y,RESi,np.arange(0.0,0.11,0.0001))
cb = m.colorbar(ticks=(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11))
cb.set_label('')
plt.set_cmap('jet')
#cb.set_clim(np.nanmin(RESi),np.nanmax(RESi))
plt.title('ResImag _ ' + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=61.5,lat_0=90,lat_ts=70,lon_0=-35,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')
