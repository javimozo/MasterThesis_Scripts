# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 23:44:29 2016

Computes C (sea currents) and A (modulus:|A| and angle:theta) coefficients  
from wind and sea ice drift as an inverse problem

Plots them on a north polar stereographic basemap.

Vectors converted to complex conjugated variables
"""

import glob
import numpy as np
import scipy.stats as st
from numpy.linalg import inv
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap

b =range(1978,2004)

SHR=np.zeros([2115,len(b)],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
SHI=np.zeros([2115,len(b)],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)

m = -1

for i in b:
    pep = i
    pip = i+1
    PEP = str(pep) + '-' + str(pip)
    print 'PEP: ', PEP

    m += 1    
    
    filenames1 = glob.glob('C:/Users/javier/Desktop/Master_thesis/Sandbox/Data/nc/'+PEP+'_nc/*.nc')
    filenames2 = glob.glob('C:/Users/javier/Desktop/Master_thesis/Sandbox/Data/vec/'+PEP+'_vec/*.vec')
    
    Wu=np.zeros([2115,len(filenames1)],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    Wv=np.zeros([2115,len(filenames1)],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    Iu=np.zeros([2115,len(filenames2)],dtype=np.float)    # 2D array to save LOCAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)
    Iv=np.zeros([2115,len(filenames2)],dtype=np.float)    # 2D array to save MERIDIONAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)
    
    datum1 = Dataset(filenames1[0], mode='r')
    
    print 'filename: ', filenames1[0]
    
    lat=datum1.variables['lat'][:]
    lon=datum1.variables['lon'][:]
    lat1=datum1.variables['lat1'][:]
    lon1=datum1.variables['lon1'][:] 
    
    n = -1        # initial value of counter for .nc files
    nn = -1       # initial value of counter for .vec files
    
    for filename in filenames1:                         # extracts data from each .nc file one by one
        datum1 = Dataset(filename, mode='r')
        
        print 'filename: ', filename        
             
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
    
    SHreal = np.zeros([2115,],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    SHimag = np.zeros([2115,],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    
    SHreal[:,] = np.nan
    SHimag[:,] = np.nan
    
    for element in range(2115):
        
        print 'element: ', element
        
        WC = Wc[element,:]                                # extracts one by one every row (which represents wind values for each grid point)
        IC = Ic[element,:]                                # extracts one by one every row (which represents ice drift values for each grid point)
        k = np.isnan(IC)    
        WC = WC[~k]
        IC = IC[~k]

        if len(IC) > 10:
            ID = np.eye(len(WC))
            IC = IC[:,np.newaxis]                             # transforms the horizontal array into a vertical one
            WC = [WC,np.ones_like(WC,dtype=np.complex)]       # attaches an array of ones to the list of complex velocity values
            WC = np.mat(WC)                                   # converts a list of a 1D list and 1D array into a 2 row matrix
            WC = WC.T                                         # transposes the matrix into a 2 column one
            WCH = WC.getH()                                   # transforms it into a Hermitian (conjugate transpose) matrix
                 
            Res = (ID - (WC*(inv(WCH*WC))*WCH))*IC
            Res = np.array(Res)        
            Res = np.reshape(Res,(len(Res),))        
    
            shReal = st.shapiro(Res.real)
            shImag = st.shapiro(Res.imag)
            
#            SHreal[element,] = shReal[0]
#            SHimag[element,] = shImag[0]
            
            if shReal[1] > 0.05:
                SHreal[element,] = shReal[0]       # appends the SEA CURRENTS VELOCITY values to a list as they're generated
            else:
                SHreal[element,] = 0.0
                
            if shImag[1] > 0.05:
                SHimag[element,] = shImag[0]       # appends the A parameter values to a list as they're generated
            else:
                SHimag[element,] = 0.0

    SHR[:,m] = SHreal         # saves LOCAL (X) ICE VELOCITY values in columns 
    SHI[:,m] = SHimag         # saves MERIDIONAL (Y) ICE VELOCITY values in columns 

shR = np.zeros([2115,],dtype=np.float)
shI = np.zeros([2115,],dtype=np.float)
   
for ii in range(2115):
    shR[ii] = np.nanmean(SHR[ii,:])
    shI[ii] = np.nanmean(SHI[ii,:])

shR = np.reshape(shR,(47,45))     # reshapes the A parameter 1D array into 2D for plotting with contourf
shI = np.reshape(shI,(47,45))     # reshapes the A parameter 1D array into 2D for plotting with contourf

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

plt.figure(11)
c1 = m.contourf(x,y,shR,100)
cb = m.colorbar()
cb.set_label('')
#cb.set_clim(np.nanmin(shR),np.nanmax(shR))
plt.title("S-W stat _ ResidualsReal _ (p_val>0.05) _ " + str(b[0]) + '-' + str(b[-1]+1) + "average (YbyY)")

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(21)
c2 = m.contourf(x,y,shI,100)
cb = m.colorbar()
cb.set_label('')
#cb.set_clim(np.nanmin(shI),np.nanmax(shI))
plt.title("S-W stat _ ResidualsImag _ (p_val>0.05) _ " + str(b[0]) + '-' + str(b[-1]+1) + "average (YbyY)")

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

#plt.show()


