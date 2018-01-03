# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 00:32:22 2016

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

b = range(1978,2004)

ADR = np.zeros([2115,len(b)],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
ADI = np.zeros([2115,len(b)],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
#ADC=np.zeros([2115,len(b)],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)

m = -1

for i in b:
    pep = i
    pip = i+1
    PEP = str(pep) + '-' + str(pip)
    print 'PEP: ', PEP

    m += 1    
    
    filenames1 = glob.glob('C:/Users/javier/Desktop/Master_thesis/Sandbox/Data/nc/'+PEP+'_nc/*.nc')
    filenames2 = glob.glob('C:/Users/javier/Desktop/Master_thesis/Sandbox/Data/vec/'+PEP+'_vec/*.vec')
    
    Wu = np.zeros([2115,len(filenames1)],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    Wv = np.zeros([2115,len(filenames1)],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    Iu = np.zeros([2115,len(filenames2)],dtype=np.float)    # 2D array to save LOCAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)
    Iv = np.zeros([2115,len(filenames2)],dtype=np.float)    # 2D array to save MERIDIONAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)
    
    datum1 = Dataset(filenames1[0], mode='r')
    
    print 'filename: ', filenames1[0]
    
#    lat=datum1.variables['lat'][:]
#    lon=datum1.variables['lon'][:]
    lat1 = datum1.variables['lat1'][:]
    lon1 = datum1.variables['lon1'][:] 
    
    n = -1        # initial value of counter for .nc files
    nn = -1       # initial value of counter for .vec files
    
    for filename in filenames1:                         # extracts data from each .nc file one by one
        datum1 = Dataset(filename, mode='r')
        
        print 'filename: ', filename        
             
        uw = datum1.variables['u_wind_avg'][:]
        vw = datum1.variables['v_wind_avg'][:]

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
    
    ADreal = np.zeros([2115,],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    ADimag = np.zeros([2115,],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
#    ADcomp = np.zeros([2115,],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    
    ADreal[:,] = np.nan
    ADimag[:,] = np.nan
#    ADcomp[:,] = np.nan
    
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
    
            adReal = st.anderson(Res.real)
            adImag = st.anderson(Res.imag)
            
#            AbsRes = abs(Res)
#            adComp = st.anderson(AbsRes)
            
#            if adReal[0] < adReal[1][2]:
#                ADreal[element,] = adReal[0]/adReal[1][2]       # appends the SEA CURRENTS VELOCITY values to a list as they're generated
#            else:
#                ADreal[element,] = np.nan
#                
#            if adImag[0] < adImag[1][2]:
#                ADimag[element,] = adImag[0]/adImag[1][2]       # appends the A parameter values to a list as they're generated
#            else:
#                ADimag[element,] = np.nan

            ADreal[element,] = adReal[0]/adReal[1][2]
            ADimag[element,] = adImag[0]/adImag[1][2]
#            ADcomp[element,] = adComp[0]/adComp[1][2]

    ADR[:,m]=ADreal         # saves LOCAL (X) ICE VELOCITY values in columns 
    ADI[:,m]=ADimag         # saves MERIDIONAL (Y) ICE VELOCITY values in columns 
#    ADC[:,m]=ADcomp         # saves MERIDIONAL (Y) ICE VELOCITY values in columns 



adR=np.zeros([2115,],dtype=np.float)
adI=np.zeros([2115,],dtype=np.float)
#adC=np.zeros([2115,],dtype=np.float)
   
for ii in range(2115):
    adR[ii] = np.nanmean(ADR[ii,:])
    
    adi = ADI[ii,:]
    kk = np.isinf(adi)
    adi = adi[~kk]
    
    adI[ii] = np.nanmean(adi)
#    adC[ii] = np.nanmean(ADC[ii,:])

adR=np.reshape(adR,(47,45))     # reshapes the A parameter 1D array into 2D for plotting with contourf
adI=np.reshape(adI,(47,45))     # reshapes the A parameter 1D array into 2D for plotting with contourf
#adC=np.reshape(adC,(47,45))     # reshapes the A parameter 1D array into 2D for plotting with contourf

#############################
### Plotting with Basemap ###
#############################

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')
#m.drawparallels(np.arange(-80.,81.,5.))
#m.drawmeridians(np.arange(-180.,181.,20.))
#m.drawmapboundary(fill_color='aqua')

lat1=np.flipud(lat1)
lon1=np.flipud(lon1)
   
x,y = m(lon1,lat1,inverse=False)#,indexing='ij')

plt.figure(1)
c1 = m.contourf(x,y,adR,np.arange(0.0,3.4,0.001))
cb = m.colorbar()
cb.set_label('')
#cb.set_clim(np.nanmin(adR),np.nanmax(adR))
plt.title("A-D stat _ ResReal _ (stat/critVal) _ " + str(b[0]) + '-' + str(b[-1]+1) +" average (YbyY)")

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(2)
c2 = m.contourf(x,y,adI,np.arange(0.0,3.4,0.001))
cb = m.colorbar()
cb.set_label('')
#cb.set_clim(np.nanmin(adI),np.nanmax(adI))
plt.title("A-D stat _ ResImag _ (stat/critVal) _ " + str(b[0]) + '-' + str(b[-1]+1) +" average (YbyY)")

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

#plt.figure(3)
#c2 = m.contourf(x,y,adC,100)
#cb = m.colorbar()
#cb.set_label('')
#cb.set_clim(np.nanmin(adC),np.nanmax(adC))
#plt.title("Anderson-Darlig stat for Res.comp _ YearByYear _ " + str(b[0]) + '-' + str(b[-1]+1))
#
#m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
#            resolution='l',round=False)
#m.drawcoastlines()
#m.fillcontinents(color='grey')#,lake_color='aqua')

#plt.show()

