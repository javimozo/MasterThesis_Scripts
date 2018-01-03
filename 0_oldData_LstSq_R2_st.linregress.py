# -*- coding: utf-8 -*-
"""
Created on Wed Dec 06 11:30:52 2017

@author: admin
"""

import glob
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import scipy.stats as st

b =range(1978,2004)

#KOV=np.zeros([2115,len(b)],dtype=np.complex)    # 2D array to save COVARIANCE values as they're computed from each .nc file (colunmns)
#ERRE=np.zeros([2115,len(b)],dtype=np.complex)    # 2D array to save CORRELATION COEFFICIENT values as they're extracted from each .nc file (colunmns)
ERRE2=np.zeros([2115,len(b)],dtype=np.complex)    # 2D array to save R2 values as they're extracted from each .nc file (colunmns)
ERRE2[:,].real = np.nan         # changes zeros to nans
ERRE2[:,].imag = np.nan

m = -1

for i in b:
    pep = i
    pip = i+1
    PEP = str(pep) + '-' + str(pip)
    print ('PEP: ', PEP)

    m += 1    

    filenames1 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Sandbox/Data/nc/'+PEP+'_nc/*.nc')
    filenames2 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Sandbox/Data/vec/'+PEP+'_vec/*.vec')
    
    Wu=np.zeros([2115,len(filenames1)])   # 2D array to save local wind velocity values as they're extracted from each file (colunmns)
    Wv=np.zeros([2115,len(filenames1)])   # 2D array to save meridional wind velocity values as they're extracted from each file (colunmns)
    Iu=np.zeros([2115,len(filenames2)])   # 2D array to save local ice velocity values as they're extracted from each file (colunmns)
    Iv=np.zeros([2115,len(filenames2)])   # 2D array to save meridional wind velocity values as they're extracted from each file (colunmns)
    
    datum1 = Dataset(filenames1[0], mode='r')
    
    print ('filename: ', filenames1[0])
    
#    lat=datum1.variables['lat'][:]
#    lon=datum1.variables['lon'][:]
    lat1=datum1.variables['lat1'][:]
    lon1=datum1.variables['lon1'][:] 
    
    n=-1        # initial value of counter for .nc files loop
    nn=-1       # initial value of counter for .vec files loop
    
    for filename in filenames1:
        datum1 = Dataset(filename, mode='r')
        
        print ('filename', filename)
        
        uw=datum1.variables['u_wind_avg'][:]
        vw=datum1.variables['v_wind_avg'][:]

        uw=np.flipud(uw)
        vw=np.flipud(vw)
            
        n+=1                    # counter
        Wu[:,n]=uw.flat         # saves LOCAL (X) WIND VELOCITY values in columns after flattening the 2D array from .nc file into 1D
        Wv[:,n]=vw.flat         # saves MERIDIONAL (Y) WIND VELOCITY values in columns after flattening the 2D array from .nc file into 1D
                   
    for filename in filenames2:
        f = open(filename, 'r')
        datum2 = np.genfromtxt(f, skip_header=1, usecols=(4,5,6,7,8,9))
        f.close()
        
        print ('filename', filename)
    
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

    for element in range(2115):
        
        print ('element: ', element)
        
        WC = Wc[element,:]        # extract values by rows which represent each grid point
        IC = Ic[element,:]
        
        k = np.isnan(IC)
        IC = IC[~k]
        WC = WC[~k]
                                               
        WIC = np.vstack((WC,IC))       # stack them on top of each other
                                       # row=variable ; columns=values
        if len(IC) > 10:
            
#            covar = np.cov(WIC,bias=0)        # computes TOTAL COVARIANCE at each grid point
            r = np.corrcoef(WIC)              # computes TOTAL CORRELATION COEFFICIENT at each grid point
            r2 = r*r                         # computes TOTAL SQUARED CORRELATION COEFFICIENT at each grid point
            
#            KOV[element,m] = covar[1,0]          # saves values obtained for TOTAL covariance at each grid point to list
#            ERRE[element,m] = r[1,0]                # saves values obtained for TOTAL correlation coefficient at each grid point to list
            ERRE2[element,m] = r2[1,0]              # saves values obtained for TOTAL squared correlation coefficient at each grid point to list

ERRE2 = abs(ERRE2)
    
Slope_R2 = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
Slope_R2[:,] = np.nan                           # changes zeros to nans
#Int = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
#Int[:,] = np.nan                           # changes zeros to nans
#r_val_R2 = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
#r_val_R2[:,] = np.nan                           # changes zeros to nans
p_val_R2 = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
p_val_R2[:,] = np.nan                           # changes zeros to nans
#std_R2 = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
#std_R2[:,] = np.nan                           # changes zeros to nans

for element in range(2115):         # loop over every gridpoint
    
    print 'element: ', element
    
    R2 = ERRE2[element,:]              # extracts one by one every row (which represents ice drift values for each grid point)
    Years = np.arange(1,len(b)+1,1)        # creates an array that includes all years to consider
    k = np.isnan(R2)                    # checks for nans
    R2 = R2[~k]                      # removes nans
    Years = Years[~k]                      #    "      "  # DO I HAVE TO ???? YESSSS !!!!

    if len(R2) > 2:                     # at least 2 points to fit a line
       slope_R2, intercept_R2, r_value_R2, p_value_R2, std_err_R2 = st.linregress(Years,R2)

       Slope_R2[element,] = slope_R2       # saves the A parameter values to a list as they're generated
#       Int[element,] = intercept
#       r_val_R2[element,] = r_value_R2
       p_val_R2[element,] = p_value_R2
#       std_R2[element,] = std_err_R2
       
Slope_R2 = np.reshape(Slope_R2,(47,45))
#Int = np.reshape(Int,(47,45))
#r_val_R2 = np.reshape(r_val_R2,(47,45))
p_val_R2 = np.reshape(p_val_R2,(47,45))
#std_R2 = np.reshape(std_R2,(47,45))

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

#plt.figure(1)
#plt.set_cmap('bwr')
#c1 = m.contourf(x,y,r_val_R2,100)#np.arange(-0.0015, 0.0015, 0.00001))#np.arange(-0.003,0.003,0.00001))
#cb = m.colorbar()
##cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
#cb.set_label('')
#plt.title("r_value of linear regression on R2 _ " + str(b[0]) + '-' + str(b[-1]+1))
#
#m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
#            resolution='l',round=False)
##m.drawcoastlines()
#m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(2)
plt.set_cmap('jet')
c1 = m.contourf(x,y,p_val_R2,np.arange(0.0,1.0,0.0001))
cb = m.colorbar(ticks=[0.10, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
#cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
cb.set_label('')
plt.title("p_value of linear regression on R2 _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(3)
plt.set_cmap('bwr')
c1 = m.contourf(x,y,Slope_R2,np.arange(-0.028,0.028,0.0001))
cb = m.colorbar(ticks=[-0.025, -0.020, -0.015, -0.010, -0.005, 0.000, 0.005, 0.010, 0.015, 0.020, 0.025])
#cb.set_clim(-0.02,0.02)
cb.set_label('')
plt.title("LstSq_Slope of linear regression on R2 _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

#plt.figure(4)
#plt.set_cmap('jet')
#c2 = m.contourf(x,y,std_R2,np.arange(0.0,0.0005,0.00001))
#cb = m.colorbar()
##cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
#cb.set_label('')
#plt.title("Std of linear regression on Res_R _ " + str(b[0]) + '-' + str(b[-1]+1))
#
#m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
#            resolution='l',round=False)
##m.drawcoastlines()
#m.fillcontinents(color='grey')#,lake_color='aqua')
