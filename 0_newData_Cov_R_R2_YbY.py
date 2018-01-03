# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 00:01:33 2016

Computes covariance and correlation coefficient matrices 
of wind and sea ice drift (magnitude and angle)

Plots them on a north polar stereographic basemap.

Vectors converted to complex conjugated variables
"""

import glob
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
#import matplotlib as mpl

b = range(2012,2015)
#B = 'Summer'
#B = 'Winter'
B = ''

KOV=np.zeros([21063,len(b)],dtype=np.complex)    # 2D array to save COVARIANCE values as they're computed from each .nc file (colunmns)
ERRE=np.zeros([21063,len(b)],dtype=np.complex)    # 2D array to save CORRELATION COEFFICIENT values as they're extracted from each .nc file (colunmns)
ERRE2=np.zeros([21063,len(b)],dtype=np.complex)    # 2D array to save R2 values as they're extracted from each .nc file (colunmns)

m = -1

for i in b:
    pep = i
    pip = i+1
    PEP = str(pep) + '-' + str(pip)
    print 'PEP: ', PEP

    m += 1    

    filenames1 = glob.glob('C:/Users/javier/Desktop/Master_thesis/Javier_ERAI_winds_icedrift_OSISAF/SortedData_OSISAF_IceWind/WindSpeed_10m/' + PEP + '/' + B + '/*.nc')
    filenames2 = glob.glob('C:/Users/javier/Desktop/Master_thesis/Javier_ERAI_winds_icedrift_OSISAF/SortedData_OSISAF_IceWind/SeaIceDriftSpeed/' + PEP + '/' + B + '/*.nc')

    Wu=np.zeros([21063,len(filenames1)],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    Wv=np.zeros([21063,len(filenames1)],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    Iu=ma.zeros([21063,len(filenames2)],dtype=np.float)    # 2D array to save LOCAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)
    Iv=ma.zeros([21063,len(filenames2)],dtype=np.float)    # 2D array to save MERIDIONAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)

#    datum1 = Dataset(filenames1[0], mode='r')
#    
#    print 'filename: ', filenames1[0]
        
    n=-1        # initial value of counter for .nc files loop
    nn=-1       # initial value of counter for .vec files loop
    
    for filename in filenames1:
        datum1 = Dataset(filename, mode='r')
        
        print 'filename', filename
        
#        lat=datum1.variables['lat'][:]
#        lon=datum1.variables['lon'][:]
        lat1=datum1.variables['lat1'][:]
        lon1=datum1.variables['lon1'][:] 

        uw=datum1.variables['u_wind_avg'][:]
        vw=datum1.variables['v_wind_avg'][:]

#        uw=np.flipud(uw)
#        vw=np.flipud(vw)
            
        n+=1                    # counter

        Wu[:,n] = uw.flat         # saves LOCAL (X) WIND VELOCITY values in columns after flattening the 2D array from .nc file into 1D
        Wv[:,n] = vw.flat         # saves MERIDIONAL (Y) WIND VELOCITY values in columns after flattening the 2D array from .nc file into 1D
                   
    for filename in filenames2:                         # extracts data from each .nc file one by one
        datum2 = Dataset(filename, mode='r')
            
        print 'filename: ', filename        
            
        ui=datum2.variables['dX'][:]
        vi=datum2.variables['dY'][:]
        
        nn += 1
        
        ui = ui[0,:,:]*(1000.0/(2*86400))
        vi = vi[0,:,:]*(1000.0/(2*86400))
        
        Iu[:,nn] = ma.ravel(ui)             # saves LOCAL (X) ICE VELOCITY values in columns 
        Iv[:,nn] = ma.ravel(vi)             # saves MERIDIONAL (Y) ICE VELOCITY values in columns 
        
    Kov=np.zeros([21063,],dtype=np.complex)       # list to save MODULUS COVARIANCE values as they're generated
    R=np.zeros([21063,],dtype=np.complex)         # list to save MODULUS CORRELATION COEFFICIENT values as they're generated
    R2=np.zeros([21063,],dtype=np.complex)        # list to save MODULUS SQUARED CORRELATION COEFFICIENT values as they're generated
    
    Kov[:,].real = np.nan
    Kov[:,].imag = np.nan
    R[:,].real = np.nan
    R[:,].imag = np.nan
    R2[:,].real = np.nan
    R2[:,].imag = np.nan

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
                                               
        WIC = np.vstack((WC,IC))       # stack them on top of each other
                                       # row=variable ; columns=values
        
        ###############################################################################
        # I do not conjugate and/or transpose any of of the arrays 'cos I assume that #
        # functions cov, coorcoef already do it on one array when doing computations  #
        # When tried Hermitian I receive an error message: dimensions don't match     #
        # When just conjugate one array it plots funny results; with both arrays      #
        # conjugated results are same as with no conjugate and/or transposed          #
        ###############################################################################

        if len(IC) > 10:
            
            covar=np.cov(WIC,bias=0)        # computes TOTAL COVARIANCE at each grid point
            r=np.corrcoef(WIC)              # computes TOTAL CORRELATION COEFFICIENT at each grid point
            r2=r**2                         # computes TOTAL SQUARED CORRELATION COEFFICIENT at each grid point
            
            Kov[element,] = covar[1,0]          # saves values obtained for TOTAL covariance at each grid point to list
            R[element,] = r[1,0]                # saves values obtained for TOTAL correlation coefficient at each grid point to list
            R2[element,] = r2[1,0]              # saves values obtained for TOTAL squared correlation coefficient at each grid point to list
    
    KOV[:,m] = Kov
    ERRE[:,m] = R
    ERRE2[:,m] = R2
    
Kovar=np.zeros([21063,],dtype=np.complex)
Erre=np.zeros([21063,],dtype=np.complex)
Erre2=np.zeros([21063,],dtype=np.complex)
   
for ii in range(21063):
    Kovar[ii] = np.nanmean(KOV[ii,:])
    Erre[ii] = np.nanmean(ERRE[ii,:])
    Erre2[ii] = np.nanmean(ERRE2[ii,:])

absKov=abs(Kovar)
angKov=np.angle(Kovar,deg=True)
absR=abs(Erre)
angR=np.angle(Erre,deg=True)
absR2=abs(Erre2)
angR2=np.angle(Erre2,deg=True)

absKov=np.reshape(absKov,(177,119))           # reshapes them into 47x45 2D arrays
absR=np.reshape(absR,(177,119))
absR2=np.reshape(absR2,(177,119))
angKov=np.reshape(angKov,(177,119))
angR=np.reshape(angR,(177,119))
angR2=np.reshape(angR2,(177,119))

#############################
### Plotting with Basemap ###
#############################

m = Basemap(projection='npstere',boundinglat=61.5,lat_0=90,lat_ts=70,lon_0=-35,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='gray')#,lake_color='aqua')
m.drawparallels(np.arange(-80.,81.,10.))
m.drawmeridians(np.arange(-180.,181.,20.))

#lat1=np.flipud(lat1)        # flips upside down latitude data from netCDF file
#lon1=np.flipud(lon1)        # flips upside down longitude data from netCDF file

X,Y = m(lon1,lat1,inverse=False)

plt.figure(1)
cs1 = m.contourf(X,Y,absKov,np.arange(0.008,1.3,0.01))
cb = m.colorbar(cs1)
cb.set_label('')
#cb.set_clim(np.nanmin(absKov),np.nanmax(absKov))
#plt.title("|Cov| _ Wind/Ice _ Oct " + str(b[0]) + ' - Sep ' + str(b[-1]+1) + " average")
#plt.title("|Cov| _ Wind/Ice _ Summer _ " + str(b[-1]+1) + " _ average")
plt.title("|Cov| _ Wind/Ice _ " + B + " _ " + str(b[0]) + '-' + str(b[-1]+1) + " average")
#plt.title("|Cov| _ Wind/Ice _ Winter _ " + str(b[0]) + '-' + str(b[-1]+1) + " average")

m = Basemap(projection='npstere',boundinglat=61.5,lat_0=90,lat_ts=70,lon_0=-35,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='gray')#,lake_color='aqua')
m.drawparallels(np.arange(-80.,81.,10.))
m.drawmeridians(np.arange(-180.,181.,20.))

plt.figure(2)
cs2 = m.contourf(X,Y,absR,100)#np.arange(0.0,1.0,0.01))
cb = m.colorbar()
cb.set_label('')
#cb.set_clim(np.nanmin(absR),np.nanmax(absR))
#plt.title("|R| _ Wind/Ice _ Oct " + str(b[0]) + ' - Sep ' + str(b[-1]+1) + " average")
#plt.title("|R| _ Wind/Ice _ Summer _ " + str(b[-1]+1) + " _ average")
plt.title("|R| _ Wind/Ice _ " + B + " _ " + str(b[0]) + '-' + str(b[-1]+1) + " average")
#plt.title("|R| _ Wind/Ice _ Winter _ " + str(b[0]) + '-' + str(b[-1]+1) + " average")

m = Basemap(projection='npstere',boundinglat=61.5,lat_0=90,lat_ts=70,lon_0=-35,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='gray')#,lake_color='aqua')
m.drawparallels(np.arange(-80.,81.,10.))
m.drawmeridians(np.arange(-180.,181.,20.))

plt.figure(3)
cs3 = m.contourf(X,Y,absR2,np.arange(0.0,1.0,0.01))
cb = m.colorbar()
cb.set_label('')
#cb.set_clim(np.nanmin(absR2),np.nanmax(absR2))
#plt.title("|R2| _ Wind/Ice _ Oct " + str(b[0]) + ' - Sep ' + str(b[-1]+1) + " average")
#plt.title("|R2| _ Wind/Ice _ Summer _ " + str(b[-1]+1) + " _ average")
plt.title("|R2| _ Wind/Ice _ " + B + " _ " + str(b[0]) + '-' + str(b[-1]+1) + " average")
#plt.title("|R2| _ Wind/Ice _ Winter _ " + str(b[0]) + '-' + str(b[-1]+1) + " average")

m = Basemap(projection='npstere',boundinglat=61.5,lat_0=90,lat_ts=70,lon_0=-35,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='gray')#,lake_color='aqua')
m.drawparallels(np.arange(-80.,81.,10.))
m.drawmeridians(np.arange(-180.,181.,20.))

plt.figure(4)
cs4 = m.contourf(X,Y,angKov,np.arange(-100.0,10.0,0.1))
cb = m.colorbar()
cb.set_label('degrees')
#cb.set_clim(np.nanmin(angKov),np.nanmax(angKov))
#plt.title("angleCov _ Wind/Ice _ Oct " + str(b[0]) + ' - Sep ' + str(b[-1]+1) + " average")
#plt.title("angleCov _ Wind/Ice _ Summer _ " + str(b[-1]+1) + " _ average")
plt.title("angleCov _ Wind/Ice _ " + B + " _ " + str(b[0]) + '-' + str(b[-1]+1) + " average")
#plt.title("angleCov _ Wind/Ice _ Winter _ " + str(b[0]) + '-' + str(b[-1]+1) + " average")

m = Basemap(projection='npstere',boundinglat=61.5,lat_0=90,lat_ts=70,lon_0=-35,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='gray')#,lake_color='aqua')
m.drawparallels(np.arange(-80.,81.,10.))
m.drawmeridians(np.arange(-180.,181.,20.))

plt.figure(5)
cs5 = m.contourf(X,Y,angR,np.arange(-100.0,10.0,0.1))
cb = m.colorbar()
cb.set_label('degrees')
#cb.set_clim(np.nanmin(angR),np.nanmax(angR))
#plt.title("angleR _ Wind/Ice _ Oct " + str(b[0]) + ' - Sep ' + str(b[-1]+1) + " average")
#plt.title("angleR _ Wind/Ice _ Summer _ " + str(b[-1]+1) + " _ average")
plt.title("angleR _ Wind/Ice _ " + B + " _ " + str(b[0]) + '-' + str(b[-1]+1) + " average")
#plt.title("angleR _ Wind/Ice _ Winter _ " + str(b[0]) + '-' + str(b[-1]+1) + " average")

m = Basemap(projection='npstere',boundinglat=61.5,lat_0=90,lat_ts=70,lon_0=-35,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='gray')#,lake_color='aqua')
m.drawparallels(np.arange(-80.,81.,10.))
m.drawmeridians(np.arange(-180.,181.,20.))

plt.figure(6)
cs6 = m.contourf(X,Y,angR2,np.arange(-100.0,10.0,0.1))
cb = m.colorbar()
cb.set_label('degrees')
#cb.set_clim(np.nanmin(angR2),np.nanmax(angR2))
#plt.title("angleR2 _ Wind/Ice _ Oct " + str(b[0]) + ' - Sep ' + str(b[-1]+1) + " average")
#plt.title("angleR2 _ Wind/Ice _ Summer _ " + str(b[-1]+1) + " _ average")
plt.title("angleR2 _ Wind/Ice _ " + B + " _ " + str(b[0]) + '-' + str(b[-1]+1) + " average")
#plt.title("angleR2 _ Wind/Ice _ Winter _ " + str(b[0]) + '-' + str(b[-1]+1) + " average")

m = Basemap(projection='npstere',boundinglat=61.5,lat_0=90,lat_ts=70,lon_0=-35,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='gray')#,lake_color='aqua')
m.drawparallels(np.arange(-80.,81.,10.))
m.drawmeridians(np.arange(-180.,181.,20.))

