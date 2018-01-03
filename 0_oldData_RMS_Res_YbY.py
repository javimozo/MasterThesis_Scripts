# -*- coding: utf-8 -*-
"""
Created on Tue Aug 02 13:13:59 2016

Computes average squared mean residuals for the inverse problem U = AG + c + e  
from wind and sea ice drift time series

Vectors converted to complex numbers
"""

import glob
import numpy as np
from numpy.linalg import inv
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

b =range(1978,2004)

RESr = np.zeros([2115,len(b)],dtype=np.complex)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
RESi = np.zeros([2115,len(b)],dtype=np.complex)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
RESc = np.zeros([2115,len(b)],dtype=np.complex)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
RESr[:,:] = np.nan
RESi[:,:] = np.nan
RESc[:,:] = np.nan

m = -1

for i in b:
    pep = i
    pip = i+1
    PEP = str(pep) + '-' + str(pip)
    print 'PEP: ', PEP

    m += 1    
    
    filenames1 = glob.glob('C:/Users/javier/Desktop/Master_thesis/Sandbox/Data/nc/' + PEP + '_nc/*.nc')
    filenames2 = glob.glob('C:/Users/javier/Desktop/Master_thesis/Sandbox/Data/vec/' + PEP + '_vec/*.vec')
    
    Wu=np.zeros([2115,len(filenames1)],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    Wv=np.zeros([2115,len(filenames1)],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    Iu=np.zeros([2115,len(filenames2)],dtype=np.float)    # 2D array to save LOCAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)
    Iv=np.zeros([2115,len(filenames2)],dtype=np.float)    # 2D array to save MERIDIONAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)
    
    n = -1        # initial value of counter for .nc files
    nn = -1       # initial value of counter for .vec files
    
    for filename in filenames1:                         # extracts data from each .nc file one by one
        datum1 = Dataset(filename, mode='r')
        
        print 'filename: ', filename        
        
        lat1=datum1.variables['lat1'][:]
        lon1=datum1.variables['lon1'][:] 

        uw=datum1.variables['u_wind_avg'][:]
        vw=datum1.variables['v_wind_avg'][:]

        uw=np.flipud(uw)
        vw=np.flipud(vw)
    
        n += 1
        Wu[:,n]=uw.flat         # saves LOCAL (X) WIND VELOCITY values in columns after flattening the 2D array from .nc file into 1D
        Wv[:,n]=vw.flat         # saves MERIDIONAL (Y) WIND VELOCITY values in columns after flattening the 2D array from .nc file into 1D
               
    for filename in filenames2:             # extracts data from each .vec file one by one
        f = file(filename, 'r')
        datum2 = np.genfromtxt(f, skip_header=1, usecols=(4,5,6,7,8,9))
        f.close()
            
        print 'filename: ', filename        
            
        ui=datum2[:,0]*(1000.0/(2*86400))       # converts to m/s from km/2days
        vi=datum2[:,1]*(1000.0/(2*86400))
        
        for i in range(len(ui)):            # checks for grid points with Nans
            if datum2[i,2] == 0 and datum2[i,3] == 0 and datum2[i,4] == 0 and datum2[i,5] == 0:
                ui[i]=np.nan
                vi[i]=np.nan
    
        if len(ui) < 2115:
           a=np.zeros([45,])
           a[:]=np.nan
           ui=np.hstack((ui,a))
           vi=np.hstack((vi,a))
    
        nn += 1
        Iu[:,nn]=ui         # saves LOCAL (X) ICE VELOCITY values in columns 
        Iv[:,nn]=vi         # saves MERIDIONAL (Y) ICE VELOCITY values in columns 
            
    Wc=np.array(Wu,dtype=np.complex)      # combines 2 real arrays (u,v) to create a complex array
    Wc.imag=(Wv)
    
    Ic=np.array(Iu,dtype=np.complex)
    Ic.imag=(Iv)
        
    for element in range(2115):
        
        print 'element: ', element
        
        WC = Wc[element,:]                                # extracts one by one every row (which represents wind values for each grid point)
        IC = Ic[element,:]                                # extracts one by one every row (which represents ice drift values for each grid point)
        k = np.isnan(IC)    
        WC = WC[~k]
        IC = IC[~k]

        if len(IC) > 5:
            ID = np.eye(len(WC))
            IC = IC[:,np.newaxis]                             # transforms the horizontal array into a vertical one
            WC = [WC,np.ones_like(WC,dtype=np.complex)]       # attaches an array of ones to the list of complex velocity values
            WC = np.mat(WC)                                   # converts a list of a 1D list and 1D array into a 2 row matrix
            WC = WC.T                                         # transposes the matrix into a 2 column one
            WCH = WC.getH()                                   # transforms it into a Hermitian (conjugate transpose) matrix
                 
            ResC = (ID - (WC*(inv(WCH*WC))*WCH))*IC
            
            ResR = ResC.real
            ResI = ResC.imag
            ResC_MSq = np.sqrt(np.mean(np.square(ResC)))
            ResR_MSq = np.sqrt(np.mean(np.square(ResR)))
            ResI_MSq = np.sqrt(np.mean(np.square(ResI)))
            RESr[element,m]=ResR_MSq         
            RESi[element,m]=ResI_MSq         
            RESc[element,m]=ResC_MSq         
#            Res = np.array(Res)        
#            Res = np.reshape(Res,(len(Res),))

ResMSqR = np.zeros([2115,],dtype=np.complex)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
ResMSqI = np.zeros([2115,],dtype=np.complex)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
ResMSqC = np.zeros([2115,],dtype=np.complex)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
ResMSqR[:,] = np.nan
ResMSqI[:,] = np.nan
ResMSqC[:,] = np.nan

for ii in range(len(ResMSqR)):
    ResMSqR[ii]=np.mean(RESr[ii,:])         # saves LOCAL (X) ICE VELOCITY values in columns           
    ResMSqI[ii]=np.mean(RESi[ii,:])         # saves LOCAL (X) ICE VELOCITY values in columns           
    ResMSqC[ii]=np.mean(RESc[ii,:])         # saves LOCAL (X) ICE VELOCITY values in columns           
    
ResReal = ResMSqC.real
ResImag = ResMSqC.imag    
    
ResMSqR = np.reshape(ResMSqR,(47,45))
ResMSqI = np.reshape(ResMSqI,(47,45))
ResReal = np.reshape(ResReal,(47,45))
ResImag = np.reshape(ResImag,(47,45))

#############################
### Plotting with Basemap ###
#############################

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')
#m.drawparallels(np.arange(-80.,81.,20.))
#m.drawmeridians(np.arange(-180.,181.,20.))
#m.drawmapboundary(fill_color='aqua')

lat1=np.flipud(lat1)
lon1=np.flipud(lon1)
   
x,y = m(lon1,lat1,inverse=False)#,indexing='ij')

plt.figure(1)
c3 = m.contourf(x,y,ResReal,100)
cb = m.colorbar()
#cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
cb.set_label('')
plt.title("mean squared complexRes.real _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(2)
c3 = m.contourf(x,y,ResMSqR,100)
cb = m.colorbar()
#cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
cb.set_label('')
plt.title("mean squared realRes _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(3)
c2 = m.contourf(x,y,ResMSqR,np.arange(0, 0.0035, 0.0001))
cb = m.colorbar()
#cb.set_clim(np.nanmin(ResR),np.nanmax(ResR))
cb.set_label('')
plt.title("mean squared realRes _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(4)
c3 = m.contourf(x,y,ResImag,100)
cb = m.colorbar()
#cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
cb.set_label('')
plt.title("mean squared complexRes.imag _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(5)
c3 = m.contourf(x,y,ResMSqI,100)
cb = m.colorbar()
#cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
cb.set_label('')
plt.title("mean squared imagRes _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(6)
c3 = m.contourf(x,y,ResMSqI,np.arange(0, 0.0035, 0.0001))
cb = m.colorbar()
#cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
cb.set_label('')
plt.title("mean squared imagRes _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')
