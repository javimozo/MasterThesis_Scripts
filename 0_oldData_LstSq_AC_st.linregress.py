# -*- coding: utf-8 -*-
"""
Created on Tue Dec 05 19:35:45 2017

@author: admin
"""

import glob
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import scipy.stats as st

b = range(1978,2004)

C = np.zeros([2115,len(b)],dtype=np.complex)    # 2D array to save A and C values 
A = np.zeros([2115,len(b)],dtype=np.complex)    # as they're computed for each grid point (rows) for every year (colunmns)
C[:,].real = np.nan         # changes zeros to nans
C[:,].imag = np.nan
A[:,].real = np.nan
A[:,].imag = np.nan

m = -1          # time counter initial value

for i in b:         # loop over year by year
    pep = i
    pip = i+1
    PEP = str(pep) + '-' + str(pip)
    print 'PEP: ', PEP

    m += 1          # time counter
    
    filenames1 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Sandbox/Data/nc/'+PEP+'_nc/*.nc')
    filenames2 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Sandbox/Data/vec/'+PEP+'_vec/*.vec')
    
    Wu = np.zeros([2115,len(filenames1)],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    Wv = np.zeros([2115,len(filenames1)],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    Iu = np.zeros([2115,len(filenames2)],dtype=np.float)    # 2D array to save LOCAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)
    Iv = np.zeros([2115,len(filenames2)],dtype=np.float)    # 2D array to save MERIDIONAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)
    
    datum1 = Dataset(filenames1[0], mode='r')
    
    print 'filename: ', filenames1[0]
    
#    lat = datum1.variables['lat'][:]
#    lon = datum1.variables['lon'][:]
    lat1 = datum1.variables['lat1'][:]
    lon1 = datum1.variables['lon1'][:] 
    
    n = -1        # initial value of counter for .nc files
    nn = -1       # initial value of counter for .vec files
    
    for filename in filenames1:                         # extracts data from each .nc file one by one
        datum1 = Dataset(filename, mode='r')
        
        print 'filename: ', filename        
             
        uw = datum1.variables['u_wind_avg'][:]
        vw = datum1.variables['v_wind_avg'][:]

        uw = np.flipud(uw)
        vw = np.flipud(vw)
    
        n += 1
        Wu[:,n] = uw.flat         # saves LOCAL (X) WIND VELOCITY values in columns after flattening the 2D array from .nc file into 1D
        Wv[:,n] = vw.flat         # saves MERIDIONAL (Y) WIND VELOCITY values in columns after flattening the 2D array from .nc file into 1D
               
    for filename in filenames2:             # extracts data from each .vec file one by one
        f = file(filename, 'r')
        datum2 = np.genfromtxt(f, skip_header=1, usecols=(4,5,6,7,8,9))
        f.close()
            
        print 'filename: ', filename        
            
        ui = datum2[:,0]*(1000.0/(2*86400))       # converts to m/s from km/2days
        vi = datum2[:,1]*(1000.0/(2*86400))
        
        for i in range(len(ui)):           
            if datum2[i,2] == 0 and datum2[i,3] == 0 and datum2[i,4] == 0 and datum2[i,5] == 0:
                ui[i] = np.nan
                vi[i] = np.nan            # checks for grid points with Nans
    
        if len(ui) < 2115:              # checks for ice drift data that lacks grid's last row
           a = np.zeros([45,])
           a[:] = np.nan
           ui = np.hstack((ui,a))
           vi = np.hstack((vi,a))
    
        nn += 1
        Iu[:,nn] = ui         # saves LOCAL (X) ICE VELOCITY values in columns 
        Iv[:,nn] = vi         # saves MERIDIONAL (Y) ICE VELOCITY values in columns 
            
    Wc = np.array(Wu,dtype=np.complex)      # combines 2 real arrays (u,v) to create a complex array
    Wc.imag=(Wv)
    
    Ic = np.array(Iu,dtype=np.complex)
    Ic.imag=(Iv)
       
    for element in range(2115):         # loop over every gridpoint
        
        print 'element: ', element
        
        WC = Wc[element,:]              # extracts one by one every row (which represents wind values for each grid point)
        IC = Ic[element,:]              # extracts one by one every row (which represents ice drift values for each grid point)
        k = np.isnan(IC)                # checks for nans
        WC = WC[~k]                     # removes wind grid point where there's no ice
        IC = IC[~k]                     # removes nans

        if len(IC) > 10:                # minimum lengh of time series to ???
            IC = IC[:,np.newaxis]                             # transforms the horizontal array into a vertical one
            WC = [WC,np.ones_like(WC,dtype=np.complex)]       # attaches an array of ones to the list of complex velocity values
            WC = np.mat(WC)                                   # converts a list of a 1D list and 1D array into a 2 row matrix
            WC = WC.T                                         # transposes the matrix into a 2 column one (1 x 2)
            WCH = WC.getH()                                   # transforms it into a Hermitian (conjugate transpose) matrix
                 
            AC = (inv(WCH*WC)*WCH*IC)         # computes the inversion // gives a 1x2 matrix with values of parameter A and sea currents C (in complex form)
                    
            C[element,m] = AC[1,0]       # appens the SEA CURRENTS VELOCITY values to a list as they're generated
            A[element,m] = AC[0,0]       # appens the A parameter values to a list as they're generated

Cabs = abs(C)
Cang = np.angle(C, deg=1)
Aabs = abs(A)
Aang = np.angle(A, deg=1)


Slope_A = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
Slope_A[:,] = np.nan                           # changes zeros to nans
#Int = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
#Int[:,] = np.nan                           # changes zeros to nans
r_val_A = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
r_val_A[:,] = np.nan                           # changes zeros to nans
p_val_A = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
p_val_A[:,] = np.nan                           # changes zeros to nans
std_A = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
std_A[:,] = np.nan                           # changes zeros to nans

for element in range(2115):         # loop over every gridpoint
    
    print 'element: ', element
    
    absA = Aabs[element,:]              # extracts one by one every row (which represents ice drift values for each grid point)
    Years = np.arange(1,len(b)+1,1)        # creates an array that includes all years to consider
    k = np.isnan(absA)                    # checks for nans
    absA = absA[~k]                      # removes nans
    Years = Years[~k]                      #    "      "  # DO I HAVE TO ???? YESSSS !!!!

    if len(absA) > 2:                     # at least 2 points to fit a line
       slope_A, intercept_A, r_value_A, p_value_A, std_err_A = st.linregress(Years,absA)

       Slope_A[element,] = slope_A       # saves the A parameter values to a list as they're generated
#       Int[element,] = intercept
       r_val_A[element,] = r_value_A
       p_val_A[element,] = p_value_A
       std_A[element,] = std_err_A
       
Slope_A = np.reshape(Slope_A,(47,45))
#Int = np.reshape(Int,(47,45))
r_val_A = np.reshape(r_val_A,(47,45))
p_val_A = np.reshape(p_val_A,(47,45))
std_A = np.reshape(std_A,(47,45))

Slope_angA = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
Slope_angA[:,] = np.nan                           # changes zeros to nans
#Int = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
#Int[:,] = np.nan                           # changes zeros to nans
r_val_angA = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
r_val_angA[:,] = np.nan                           # changes zeros to nans
p_val_angA = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
p_val_angA[:,] = np.nan                           # changes zeros to nans
std_angA = np.zeros([2115,],dtype=np.float)     # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
std_angA[:,] = np.nan                           # changes zeros to nans

for element in range(2115):         # loop over every gridpoint
    
    print 'element: ', element
    
    angA = Aang[element,:]              # extracts one by one every row (which represents ice drift values for each grid point)
    Years = np.arange(1,len(b)+1,1)        # creates an array that includes all years to consider
    k = np.isnan(angA)                    # checks for nans
    angA = angA[~k]                      # removes nans
    Years = Years[~k]                      #    "      "  # DO I HAVE TO ???? YESSSS !!!!

    if len(angA) > 2:                     # at least 2 points to fit a line
       slope_angA, intercept_angA, r_value_angA, p_value_angA, std_err_angA = st.linregress(Years,angA)

       Slope_angA[element,] = slope_angA       # saves the A parameter values to a list as they're generated
#       Int[element,] = intercept
       r_val_angA[element,] = r_value_angA
       p_val_angA[element,] = p_value_angA
       std_angA[element,] = std_err_angA
       
Slope_angA = np.reshape(Slope_angA,(47,45))
#Int = np.reshape(Int,(47,45))
r_val_angA = np.reshape(r_val_angA,(47,45))
p_val_angA = np.reshape(p_val_angA,(47,45))
std_angA = np.reshape(std_angA,(47,45))


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
plt.set_cmap('bwr')
c1 = m.contourf(x,y,r_val_A,100)#np.arange(-0.0015, 0.0015, 0.00001))#np.arange(-0.003,0.003,0.00001))
cb = m.colorbar()
#cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
cb.set_label('')
plt.title("r_value of linear regression on |A| _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(2)
plt.set_cmap('jet')
c1 = m.contourf(x,y,p_val_A,np.arange(0.0,1.0,0.0001))
cb = m.colorbar(ticks=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
#cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
cb.set_label('')
plt.title("p_value of linear regression on |A| _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(3)
plt.set_cmap('bwr')
c1 = m.contourf(x,y,Slope_A,np.arange(-0.0004,0.0004,0.000001))
cb = m.colorbar(ticks=[-0.0003, -0.0002, -0.0001, 0.0000, 0.0001, 0.0002, 0.0003])
#cb.set_clim(-0.02,0.02)
cb.set_label('')
plt.title("LstSq_Slope of linear regression on |A| _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(4)
plt.set_cmap('jet')
c2 = m.contourf(x,y,std_A,np.arange(0.0,0.0005,0.000001))
cb = m.colorbar()
#cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
cb.set_label('')
plt.title("Std of linear regression on |A| _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(5)
plt.set_cmap('bwr')
c1 = m.contourf(x,y,r_val_angA,100)#np.arange(-0.0015, 0.0015, 0.00001))#np.arange(-0.003,0.003,0.00001))
cb = m.colorbar()
#cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
cb.set_label('')
plt.title("r_value of linear regression on |A| _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(6)
plt.set_cmap('jet')
c1 = m.contourf(x,y,p_val_angA,np.arange(0.0,1.0,0.0001))
cb = m.colorbar()#ticks=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
#cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
cb.set_label('')
plt.title("p_value of linear regression on |A| _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(7)
plt.set_cmap('bwr')
c1 = m.contourf(x,y,Slope_angA,np.arange(-0.00004,0.00004,0.00001))
cb = m.colorbar()#ticks=[-0.0003, -0.0002, -0.0001, 0.0000, 0.0001, 0.0002, 0.0003])
#cb.set_clim(-0.02,0.02)
cb.set_label('')
plt.title("LstSq_Slope of linear regression on |A| _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(8)
plt.set_cmap('jet')
c2 = m.contourf(x,y,std_angA,100)#np.arange(0.0,0.0005,0.000001))
cb = m.colorbar()
#cb.set_clim(np.nanmin(ResMSqI),np.nanmax(ResMSqI))
cb.set_label('')
plt.title("Std of linear regression on |A| _ " + str(b[0]) + '-' + str(b[-1]+1))

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')
