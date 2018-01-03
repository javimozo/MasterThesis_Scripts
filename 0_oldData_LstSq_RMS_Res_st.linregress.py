# -*- coding: utf-8 -*-
"""
Created on Tue Dec 05 22:56:40 2017

@author: admin
"""

import glob
import numpy as np
from numpy.linalg import inv
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import scipy.stats as st

b = range(1978,2004)

RESr = np.zeros([2115,len(b)],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
RESi = np.zeros([2115,len(b)],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
RESr[:,:] = np.nan
RESi[:,:] = np.nan

m = -1

for i in b:
    pep = i
    pip = i+1
    PEP = str(pep) + '-' + str(pip)
    print 'PEP: ', PEP

    m += 1    
    
    filenames1 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Sandbox/Data/nc/' + PEP + '_nc/*.nc')
    filenames2 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Sandbox/Data/vec/' + PEP + '_vec/*.vec')
    
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

        if len(IC) > 10:
            ID = np.eye(len(WC))
            IC = IC[:,np.newaxis]                             # transforms the horizontal array into a vertical one
            WC = [WC,np.ones_like(WC,dtype=np.complex)]       # attaches an array of ones to the list of complex velocity values
            WC = np.mat(WC)                                   # converts a list of a 1D list and 1D array into a 2 row matrix
            WC = WC.T                                         # transposes the matrix into a 2 column one
            WCH = WC.getH()                                   # transforms it into a Hermitian (conjugate transpose) matrix
                 
            ResC = (ID - (WC*(inv(WCH*WC))*WCH))*IC
            
            RESr[element,m] = np.sqrt(np.mean(np.square(ResC.real)))
            RESi[element,m] = np.sqrt(np.mean(np.square(ResC.imag)))


Slope_R = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
Slope_R[:,] = np.nan                           # changes zeros to nans
#Int = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
#Int[:,] = np.nan                           # changes zeros to nans
#r_val_R = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
#r_val_R[:,] = np.nan                           # changes zeros to nans
p_val_R = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
p_val_R[:,] = np.nan                           # changes zeros to nans
#std_R = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
#std_R[:,] = np.nan                           # changes zeros to nans

for element in range(2115):         # loop over every gridpoint
    
    print 'element: ', element
    
    ResR = RESr[element,:]              # extracts one by one every row (which represents ice drift values for each grid point)
    Years = np.arange(1,len(b)+1,1)        # creates an array that includes all years to consider
    k = np.isnan(ResR)                    # checks for nans
    ResR = ResR[~k]                      # removes nans
    Years = Years[~k]                      #    "      "  # DO I HAVE TO ???? YESSSS !!!!

    if len(ResR) > 2:                     # at least 2 points to fit a line
       slope_R, intercept_R, r_value_R, p_value_R, std_err_R = st.linregress(Years,ResR)

       Slope_R[element,] = slope_R       # saves the A parameter values to a list as they're generated
#       Int[element,] = intercept
#       r_val_R[element,] = r_value_R
       p_val_R[element,] = p_value_R
#       std_R[element,] = std_err_R
       
Slope_R = np.reshape(Slope_R,(47,45))
#Int = np.reshape(Int,(47,45))
#r_val_R = np.reshape(r_val_R,(47,45))
p_val_R = np.reshape(p_val_R,(47,45))
#std_R = np.reshape(std_R,(47,45))


Slope_I = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
Slope_I[:,] = np.nan                           # changes zeros to nans
#Int = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
#Int[:,] = np.nan                           # changes zeros to nans
#r_val_I = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
#r_val_I[:,] = np.nan                           # changes zeros to nans
p_val_I = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
p_val_I[:,] = np.nan                           # changes zeros to nans
#std_I = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
#std_I[:,] = np.nan                           # changes zeros to nans

for element in range(2115):         # loop over every gridpoint
    
    print 'element: ', element
    
    ResI = RESi[element,:]              # extracts one by one every row (which represents ice drift values for each grid point)
    Years = np.arange(1,len(b)+1,1)        # creates an array that includes all years to consider
    k = np.isnan(ResI)                    # checks for nans
    ResI = ResI[~k]                      # removes nans
    Years = Years[~k]                      #    "      "  # DO I HAVE TO ???? YESSSS !!!!

    if len(ResI) > 2:                     # at least 2 points to fit a line
       slope_I, intercept_I, r_value_I, p_value_I, std_err_I = st.linregress(Years,ResI)

       Slope_I[element,] = slope_I       # saves the A parameter values to a list as they're generated
#       Int[element,] = intercept
#       r_val_I[element,] = r_value_I
       p_val_I[element,] = p_value_I
#       std_I[element,] = std_err_I
       
Slope_I = np.reshape(Slope_I,(47,45))
#Int = np.reshape(Int,(47,45))
#r_val_I = np.reshape(r_val_I,(47,45))
p_val_I = np.reshape(p_val_I,(47,45))
#std_I = np.reshape(std_I,(47,45))

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
#c1 = m.contourf(x,y,r_val_R,100)#np.arange(-0.0015, 0.0015, 0.00001))#np.arange(-0.003,0.003,0.00001))
#cb = m.colorbar()
##cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
#cb.set_label('')
#plt.title("r_value of linear regression on Res_R _ " + str(b[0]) + '-' + str(b[-1]+1))
#
#m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
#            resolution='l',round=False)
##m.drawcoastlines()
#m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(21)
plt.set_cmap('jet')
c1 = m.contourf(x,y,p_val_R,np.arange(0.0,1.0,0.0001))
cb = m.colorbar(ticks=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
#cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
cb.set_label('')
plt.title("p_value of linear regression on Res_R _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(31)
plt.set_cmap('bwr')
c1 = m.contourf(x,y,Slope_R,np.arange(-0.003,0.003,0.000001))
cb = m.colorbar(ticks=[-0.0025, -0.0020, -0.0015, -0.0010, -0.0005, 0.0000, 0.0005, 0.0010, 0.0015, 0.0020, 0.0025])
#cb.set_clim(-0.02,0.02)
cb.set_label('')
plt.title("LstSq_Slope of linear regression on Res_R _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

#plt.figure(4)
#plt.set_cmap('jet')
#c2 = m.contourf(x,y,std_R,np.arange(0.0,0.0005,0.000001))
#cb = m.colorbar()
##cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
#cb.set_label('')
#plt.title("Std of linear regression on Res_R _ " + str(b[0]) + '-' + str(b[-1]+1))
#
#m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
#            resolution='l',round=False)
##m.drawcoastlines()
#m.fillcontinents(color='grey')#,lake_color='aqua')

#plt.figure(5)
#plt.set_cmap('bwr')
#c1 = m.contourf(x,y,r_val_I,100)#np.arange(-0.0015, 0.0015, 0.00001))#np.arange(-0.003,0.003,0.00001))
#cb = m.colorbar()
##cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
#cb.set_label('')
#plt.title("r_value of linear regression on Res_I _ " + str(b[0]) + '-' + str(b[-1]+1))
#
#m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
#            resolution='l',round=False)
##m.drawcoastlines()
#m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(61)
plt.set_cmap('jet')
c1 = m.contourf(x,y,p_val_I,np.arange(0.0,1.0,0.0001))
cb = m.colorbar(ticks=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
#cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
cb.set_label('')
plt.title("p_value of linear regression on Res_I _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(71)
plt.set_cmap('bwr')
c1 = m.contourf(x,y,Slope_I,np.arange(-0.003,0.003,0.000001))
cb = m.colorbar(ticks=[-0.0025, -0.0020, -0.0015, -0.0010, -0.0005, 0.0000, 0.0005, 0.0010, 0.0015, 0.0020, 0.0025])
#cb.set_clim(-0.02,0.02)
cb.set_label('')
plt.title("LstSq_Slope of linear regression on Res_I _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

#plt.figure(8)
#plt.set_cmap('jet')
#c2 = m.contourf(x,y,std_I,np.arange(0.0,0.0005,0.000001))
#cb = m.colorbar()
##cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
#cb.set_label('')
#plt.title("Std of linear regression on Res_I _ " + str(b[0]) + '-' + str(b[-1]+1))
#
#m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
#            resolution='l',round=False)
##m.drawcoastlines()
#m.fillcontinents(color='grey')#,lake_color='aqua')

            