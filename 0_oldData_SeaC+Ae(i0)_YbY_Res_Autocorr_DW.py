# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 11:19:37 2016

Computes the autocorrelation function of the residuals 
from the inverse problem u = AG + c + e  
of wind and sea ice drift data

Plots the results on a north polar stereographic map from Basemap.

Vectors are converted from cartesian to complex conjugated variables
"""

import glob
import numpy as np
from statsmodels.tsa.stattools import acf
from statsmodels.stats.stattools import durbin_watson
from statsmodels.stats.diagnostic import acorr_ljungbox
from numpy.linalg import inv
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap

b = range(1978,2004)

ACreal=np.zeros([2115,len(b)],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
ACimag=np.zeros([2115,len(b)],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
DWreal=np.zeros([2115,len(b)],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
DWimag=np.zeros([2115,len(b)],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
LBreal=np.zeros([2115,len(b)],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
LBimag=np.zeros([2115,len(b)],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)

m = -1

for i in b:
    pep = i
    pip = i+1
    PEP = str(pep) + '-' + str(pip)
    print 'PEP: ', PEP

    m += 1    
    
    filenames1 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Sandbox/Data/nc/'+PEP+'_nc/*.nc')
    filenames2 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Sandbox/Data/vec/'+PEP+'_vec/*.vec')
    
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
    
    acReal = np.zeros([2115,],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    acImag = np.zeros([2115,],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    dwReal = np.zeros([2115,],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    dwImag = np.zeros([2115,],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    lbReal = np.zeros([2115,],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    lbImag = np.zeros([2115,],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    
    acReal[:,] = np.nan
    acImag[:,] = np.nan
    dwReal[:,] = np.nan
    dwImag[:,] = np.nan
    lbReal[:,] = np.nan
    lbImag[:,] = np.nan
   
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
#            Res = np.array(Res)        
#            Res = np.reshape(Res,(len(Res),))
            
#            Res_m = np.mean(Res)
#            Res_sd = np.std(Res)
    
            acR = acf(Res.real,unbiased=False,nlags=10,qstat=True,fft=False)#,alpha=.05)
            acI = acf(Res.imag,unbiased=False,nlags=10,qstat=True,fft=False)#,alpha=.05)
            
            dwR = durbin_watson(Res.real)
            dwI = durbin_watson(Res.imag)
            
            lbR = acorr_ljungbox(Res.real,lags=10)
            lbI = acorr_ljungbox(Res.imag,lags=10)
            
            lag=1
            acReal[element,] = acR[0][lag]       # appends the SEA CURRENTS VELOCITY values to a list as they're generated
            acImag[element,] = acI[0][lag]
            dwReal[element,] = dwR[0]       # appends the SEA CURRENTS VELOCITY values to a list as they're generated
            dwImag[element,] = dwI[0]
            lbReal[element,] = lbR[0][lag]       # appends the SEA CURRENTS VELOCITY values to a list as they're generated
            lbImag[element,] = lbI[0][lag]
               
    ACreal[:,m]=acReal         # saves LOCAL (X) ICE VELOCITY values in columns 
    ACimag[:,m]=acImag         # saves MERIDIONAL (Y) ICE VELOCITY values in columns 
    DWreal[:,m]=dwReal         # saves LOCAL (X) ICE VELOCITY values in columns 
    DWimag[:,m]=dwImag         # saves MERIDIONAL (Y) ICE VELOCITY values in columns 
    LBreal[:,m]=lbReal         # saves LOCAL (X) ICE VELOCITY values in columns 
    LBimag[:,m]=lbImag         # saves MERIDIONAL (Y) ICE VELOCITY values in columns 

ACR=np.zeros([2115,],dtype=np.float)
ACI=np.zeros([2115,],dtype=np.float)
DWR=np.zeros([2115,],dtype=np.float)
DWI=np.zeros([2115,],dtype=np.float)
LBR=np.zeros([2115,],dtype=np.float)
LBI=np.zeros([2115,],dtype=np.float)
   
for ii in range(2115):
    ACR[ii] = np.nanmean(ACreal[ii,:])
    ACI[ii] = np.nanmean(ACimag[ii,:])
    DWR[ii] = np.nanmean(DWreal[ii,:])
    DWI[ii] = np.nanmean(DWimag[ii,:])
    LBR[ii] = np.nanmean(LBreal[ii,:])
    LBI[ii] = np.nanmean(LBimag[ii,:])

ACR=np.reshape(ACR,(47,45))     # reshapes the A parameter 1D array into 2D for plotting with contourf
ACI=np.reshape(ACI,(47,45))     # reshapes the A parameter 1D array into 2D for plotting with contourf
DWR=np.reshape(DWR,(47,45))     # reshapes the A parameter 1D array into 2D for plotting with contourf
DWI=np.reshape(DWI,(47,45))     # reshapes the A parameter 1D array into 2D for plotting with contourf
LBR=np.reshape(LBR,(47,45))     # reshapes the A parameter 1D array into 2D for plotting with contourf
LBI=np.reshape(LBI,(47,45))     # reshapes the A parameter 1D array into 2D for plotting with contourf

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

plt.figure(1)
c1 = m.contourf(x,y,ACR,100)
cb = m.colorbar()
cb.set_label('')
cb.set_clim(np.nanmin(ACR),np.nanmax(ACR))
plt.title("Autocorr _ Res.real _ # lags: " + str(lag) + " _ " + str(b[0]) + "-" + str(b[-1]+1) + " average (YbyY)")

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(2)
c2 = m.contourf(x,y,ACI,100)
cb = m.colorbar()
cb.set_label('')
cb.set_clim(np.nanmin(ACI),np.nanmax(ACI))
plt.title("Autocorr _ Res.imag _ # lags: " + str(lag) + " _ " + str(b[0]) + "-" + str(b[-1]+1) + " average (YbyY)")

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(3)
c1 = m.contourf(x,y,DWR,100)
cb = m.colorbar()
cb.set_label('')
cb.set_clim(np.nanmin(DWR),np.nanmax(DWR))
plt.title("DW stat for Res.real _ YearByYear _ # lags: " + str(lag) + " _ " + str(b[0]) + "-" + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(4)
c2 = m.contourf(x,y,DWI,100)
cb = m.colorbar()
cb.set_label('')
cb.set_clim(np.nanmin(DWI),np.nanmax(DWI))
plt.title("DW stat for Res.imag _ YearByYear _ # lags: " + str(lag) + " _ " + str(b[0]) + "-" + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(5)
c1 = m.contourf(x,y,LBR,100)
cb = m.colorbar()
cb.set_label('')
cb.set_clim(np.nanmin(LBR),np.nanmax(LBR))
plt.title("LB stat _ Res.real _ # lags: " + str(lag) + " _ " + str(b[0]) + "-" + str(b[-1]+1) + " average (YbyY)")

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(6)
c2 = m.contourf(x,y,LBI,100)
cb = m.colorbar()
cb.set_label('')
plt.set_cmap('jet')
cb.set_clim(np.nanmin(LBI),np.nanmax(LBI))
plt.title("LB stat _ Res.Imag _ # lags: " + str(lag) + " _ " + str(b[0]) + "-" + str(b[-1]+1) + " average (YbyY)")

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

