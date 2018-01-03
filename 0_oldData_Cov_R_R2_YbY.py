# -*- coding: utf-8 -*-
"""
Created on Fri May 27 13:50:31 2016

Computes covariance and correlation coefficient matrices 
of wind and sea ice drift (magnitude and angle)

Plots them on a north polar stereographic basemap.

Vectors converted to complex conjugated variables
"""

import glob
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
#import matplotlib as mpl

b =range(1978,2004)

KOV=np.zeros([2115,len(b)],dtype=np.complex)    # 2D array to save COVARIANCE values as they're computed from each .nc file (colunmns)
ERRE=np.zeros([2115,len(b)],dtype=np.complex)    # 2D array to save CORRELATION COEFFICIENT values as they're extracted from each .nc file (colunmns)
ERRE2=np.zeros([2115,len(b)],dtype=np.complex)    # 2D array to save R2 values as they're extracted from each .nc file (colunmns)

m = -1

for i in b:
    pep = i
    pip = i+1
    PEP = str(pep) + '-' + str(pip)
    print 'PEP: ', PEP

    m += 1    

    filenames1 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Sandbox/Data/nc/'+PEP+'_nc/*.nc')
    filenames2 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Sandbox/Data/vec/'+PEP+'_vec/*.vec')
    
    Wu=np.zeros([2115,len(filenames1)])   # 2D array to save local wind velocity values as they're extracted from each file (colunmns)
    Wv=np.zeros([2115,len(filenames1)])   # 2D array to save meridional wind velocity values as they're extracted from each file (colunmns)
    Iu=np.zeros([2115,len(filenames2)])   # 2D array to save local ice velocity values as they're extracted from each file (colunmns)
    Iv=np.zeros([2115,len(filenames2)])   # 2D array to save meridional wind velocity values as they're extracted from each file (colunmns)
    
    datum1 = Dataset(filenames1[0], mode='r')
    
    print 'filename: ', filenames1[0]
    
#    lat=datum1.variables['lat'][:]
#    lon=datum1.variables['lon'][:]
    lat1=datum1.variables['lat1'][:]
    lon1=datum1.variables['lon1'][:] 
    
    n=-1        # initial value of counter for .nc files loop
    nn=-1       # initial value of counter for .vec files loop
    
    for filename in filenames1:
        datum1 = Dataset(filename, mode='r')
        
        print 'filename', filename
        
        uw=datum1.variables['u_wind_avg'][:]
        vw=datum1.variables['v_wind_avg'][:]

        uw=np.flipud(uw)
        vw=np.flipud(vw)
            
        n+=1                    # counter
        Wu[:,n]=uw.flat         # saves LOCAL (X) WIND VELOCITY values in columns after flattening the 2D array from .nc file into 1D
        Wv[:,n]=vw.flat         # saves MERIDIONAL (Y) WIND VELOCITY values in columns after flattening the 2D array from .nc file into 1D
                   
    for filename in filenames2:
        f = file(filename, 'r')
        datum2 = np.genfromtxt(f, skip_header=1, usecols=(4,5,6,7,8,9))
        f.close()
        
        print 'filename', filename
    
        ui=datum2[:,0]*(1000.0/(2*86400))       # extract data column file by file
        vi=datum2[:,1]*(1000.0/(2*86400))       # converts to m/s from km/2days
    
        for i in range(len(ui)):
            if datum2[i,2] == 0 and datum2[i,3] == 0 and datum2[i,4] == 0 and datum2[i,5] == 0:
           # if np.sum((datum2[i,2],datum2[i,3],datum2[i,4],datum2[i,5]) == 0:
                ui[i]=np.nan
                vi[i]=np.nan
                
        if len(ui) < 2115:
           a=np.zeros([45,])
           a[:]=np.nan
           ui=np.hstack((ui,a))
           vi=np.hstack((vi,a))
        
        nn+=1               # counter
        Iu[:,nn]=ui         # saves LOCAL (X) ICE VELOCITY values in columns 
        Iv[:,nn]=vi         # saves MERIDIONAL (Y) ICE VELOCITY values in columns 
        
    Wc=np.array(Wu,dtype=complex)      # creates a complex array from 2 real arrays
    Wc.imag=(Wv)
        
    Ic=np.array(Iu,dtype=complex)  
    Ic.imag=(Iv)

    Kov=np.zeros([2115,],dtype=np.complex)       # list to save MODULUS COVARIANCE values as they're generated
    R=np.zeros([2115,],dtype=np.complex)         # list to save MODULUS CORRELATION COEFFICIENT values as they're generated
    R2=np.zeros([2115,],dtype=np.complex)        # list to save MODULUS SQUARED CORRELATION COEFFICIENT values as they're generated
    
    Kov[:,].real = np.nan
    Kov[:,].imag = np.nan
    R[:,].real = np.nan
    R[:,].imag = np.nan
    R2[:,].real = np.nan
    R2[:,].imag = np.nan

    for element in range(2115):
        
        print 'element: ', element
        
        WC = Wc[element,:]        # extract values by rows which represent each grid point
        IC = Ic[element,:]
        
        k = np.isnan(IC)
        IC = IC[~k]
        WC = WC[~k]
                                               
        WIC = np.vstack((WC,IC))       # stack them on top of each other
                                       # row=variable ; columns=values
        
        ###############################################################################
        # I do not conjugate and/or transpose any of of the arrays 'cos I assume that #
        # functions cov, corrcoef already do it on one array when doing computations  #
        # When tried Hermitian I receive an error message: dimensions don't match     #
        # When just conjugate one array it plots funny results; with both arrays      #
        # conjugated results are same as with no conjugate and/or transposed          #
        # Check documentation: it does conjuagte and transpose                        #        
        ###############################################################################

        if len(IC) > 10:
            
            covar = np.cov(WIC,bias=0)        # computes TOTAL COVARIANCE at each grid point
            r = np.corrcoef(WIC)              # computes TOTAL CORRELATION COEFFICIENT at each grid point
            r2 = r**2                         # computes TOTAL SQUARED CORRELATION COEFFICIENT at each grid point
            
            Kov[element,] = covar[1,0]          # saves values obtained for TOTAL covariance at each grid point to list
            R[element,] = r[1,0]                # saves values obtained for TOTAL correlation coefficient at each grid point to list
            R2[element,] = r2[1,0]              # saves values obtained for TOTAL squared correlation coefficient at each grid point to list
    
    KOV[:,m] = Kov
    ERRE[:,m] = R
    ERRE2[:,m] = R2
    
Kovar = np.zeros([2115,],dtype=np.complex)
Erre = np.zeros([2115,],dtype=np.complex)
Erre2 = np.zeros([2115,],dtype=np.complex)
   
for ii in range(2115):
    Kovar[ii] = np.nanmean(KOV[ii,:])
    Erre[ii] = np.nanmean(ERRE[ii,:])
    Erre2[ii] = np.nanmean(ERRE2[ii,:])

absKov = abs(Kovar)
angKov = np.angle(Kovar,deg=True)
absR = abs(Erre)
angR = np.angle(Erre,deg=True)
absR2 = abs(Erre2)
angR2 = np.angle(Erre2,deg=True)

absKov = np.reshape(absKov,(47,45))           # reshapes them into 47x45 2D arrays
absR = np.reshape(absR,(47,45))
absR2 = np.reshape(absR2,(47,45))
angKov = np.reshape(angKov,(47,45))
angR = np.reshape(angR,(47,45))
angR2 = np.reshape(angR2,(47,45))

#############################
### Plotting with Basemap ###
#############################

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='gray')#,lake_color='aqua')

lat1 = np.flipud(lat1)        # flips upside down latitude data from netCDF file
lon1 = np.flipud(lon1)        # flips upside down longitude data from netCDF file

X,Y = m(lon1,lat1,inverse=False)

#plt.figure(11)
#cs1 = m.contourf(X,Y,absKov,np.arange(0.0,1.75,0.0001))
#cb = m.colorbar()
#cb.set_label('')
##cb.set_clim(np.nanmin(absKov),np.nanmax(absKov))
#plt.title("|Cov| _ Wind/Ice _ " + str(b[0]) + '-' + str(b[-1]+1) + " average (YbyY)")
#
#m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
#            resolution='l',round=False)
#m.drawcoastlines()
#m.fillcontinents(color='gray')#,lake_color='aqua')

#plt.figure(21)
#cs2 = m.contourf(X,Y,absR,np.arange(0.16,0.96,0.0001))
#cb = m.colorbar()
#cb.set_label('')
##cb.set_clim(np.nanmin(absR),np.nanmax(absR))
#plt.title("|R| _ Wind/Ice _ " + str(b[0]) + '-' + str(b[-1]+1) + " average (YbyY)")
#
#m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
#            resolution='l',round=False)
#m.drawcoastlines()
#m.fillcontinents(color='gray')#,lake_color='aqua')

plt.figure(31)
cs3 = m.contourf(X,Y,absR2,np.arange(0.0,1.0,0.001))#np.arange(0.0,1.0,0.01))
cb = m.colorbar(ticks=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
cb.set_label('')
#cb.set_clim(np.nanmin(absR2),np.nanmax(absR2))
plt.title("|R2| _ Wind/Ice _ " + str(b[0]) + '-' + str(b[-1]+1) + " average (YbyY)")

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='gray')#,lake_color='aqua')

#plt.figure(41)
#cs4 = m.contourf(X,Y,angKov,100)#np.arange(-40,5,0.5))
#cb = m.colorbar()
#cb.set_label('degrees')
##cb.set_clim(np.nanmin(angKov),np.nanmax(angKov))
#plt.title("angleCov _ Wind/Ice _ " + str(b[0]) + '-' + str(b[-1]+1) + " average (YbyY)")
#
#m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
#            resolution='l',round=False)
#m.drawcoastlines()
#m.fillcontinents(color='gray')#,lake_color='aqua')
#
#plt.figure(51)
#cs5 = m.contourf(X,Y,angR,100)#np.arange(-40,5,0.5))
#cb = m.colorbar()
#cb.set_label('degrees')
##cb.set_clim(np.nanmin(angR),np.nanmax(angR))
#plt.title("angleR _ Wind/Ice _ " + str(b[0]) + '-' + str(b[-1]+1) + " average (YbyY)")
#
#m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
#            resolution='l',round=False)
#m.drawcoastlines()
#m.fillcontinents(color='gray')#,lake_color='aqua')
#
#plt.figure(61)
#cs6 = m.contourf(X,Y,angR2)#,np.arange(-40,5,0.5))
#cb = m.colorbar()
#cb.set_label('degrees')
##cb.set_clim(np.nanmin(angR2),np.nanmax(angR2))
#plt.title("angleR2 _ Wind/Ice _ " + str(b[0]) + '-' + str(b[-1]+1) + " average (YbyY)")
#
#m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
#            resolution='l',round=False)
#m.drawcoastlines()
#m.fillcontinents(color='gray')#,lake_color='aqua')

#plt.show()

############################################################################

#KOV = m.transform_scalar(Kov,lon11,lat11,47,45,returnxy=False,masked=True) ???

#im = m.pcolormesh(lon1,lat1,Kov,shading='flat',cmap=plt.cm.jet,latlon=True) ???


