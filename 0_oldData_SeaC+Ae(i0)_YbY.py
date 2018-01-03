# -*- coding: utf-8 -*-
"""
Created on Mon May 30 10:21:49 2016

Computes coefficients C (sea currents) and A (modulus:|A| and angle:theta)   
from wind and sea ice drift data as an inverse problem

Plots them on a north polar stereographic basemap.

Vectors converted to complex conjugated variables
"""

import glob
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap

b = range(1978,2004)

CC=np.zeros([2115,len(b)],dtype=np.complex)    # 2D array to save A and C values 
AA=np.zeros([2115,len(b)],dtype=np.complex)    # as they're computed for each grid point (rows) for every year (colunmns)

m = -1          # time counter initial value

for i in b:         # loop over year by year
    pep = i
    pip = i+1
    PEP = str(pep) + '-' + str(pip)
    print 'PEP: ', PEP

    m += 1          # time counter
    
    filenames1 = glob.glob('C:/Users/javier/Desktop/Master_thesis/Sandbox/Data/nc/'+PEP+'_nc/*.nc')
    filenames2 = glob.glob('C:/Users/javier/Desktop/Master_thesis/Sandbox/Data/vec/'+PEP+'_vec/*.vec')
    
    Wu=np.zeros([2115,len(filenames1)],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    Wv=np.zeros([2115,len(filenames1)],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    Iu=np.zeros([2115,len(filenames2)],dtype=np.float)    # 2D array to save LOCAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)
    Iv=np.zeros([2115,len(filenames2)],dtype=np.float)    # 2D array to save MERIDIONAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)
    
    datum1 = Dataset(filenames1[0], mode='r')
    
    print 'filename: ', filenames1[0]
    
#    lat=datum1.variables['lat'][:]
#    lon=datum1.variables['lon'][:]
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
        
        for i in range(len(ui)):           
            if datum2[i,2] == 0 and datum2[i,3] == 0 and datum2[i,4] == 0 and datum2[i,5] == 0:
                ui[i]=np.nan
                vi[i]=np.nan            # checks for grid points with Nans
    
        if len(ui) < 2115:              # checks for ice drift data that lacks grid's last row
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
    
    C=np.zeros([2115,],dtype=np.complex)    # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
    A=np.zeros([2115,],dtype=np.complex)    
    
    C[:,].real = np.nan         # changes zeros to nans
    C[:,].imag = np.nan
    A[:,].real = np.nan
    A[:,].imag = np.nan
    
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
                    
            C[element,] = AC[1,0]       # appens the SEA CURRENTS VELOCITY values to a list as they're generated
            A[element,] = AC[0,0]       # appens the A parameter values to a list as they're generated

    CC[:,m]=C        # saves LOCAL (X) ICE VELOCITY values in columns 
    AA[:,m]=A        # saves MERIDIONAL (Y) ICE VELOCITY values in columns 

CCC=np.zeros([2115,],dtype=complex)
AAA=np.zeros([2115,],dtype=complex)

for ii in range(2115):                  # loops over every grid point to compute the mean over the whole time series
    CCC[ii] = np.nanmean(CC[ii,:])
    AAA[ii] = np.nanmean(AA[ii,:])

Mean = np.nanmean(AAA)      
StD = np.nanstd(AAA)

MEANabs = abs(Mean)
MEANang = np.angle(Mean, deg=1)
STDabs = abs(StD)
STDang = np.angle(StD, deg=1)

AAA=np.reshape(AAA,(47,45))             # reshapes the A parameter 1D array into 2D for plotting with contourf
absA=abs(AAA)                           # |A| informs of coupling btw wind/ice (and internal ice stresses)
ThA=np.angle(AAA,deg=True)              # angle btw wind/ice motion vectors

modulus = abs(CCC)
#modulus = np.sqrt(C.real**2 + C.imag**2)*(1000/(2*86400))

c = CCC/modulus          # normalized vectors

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
Q = m.quiver(x,y,c.real,c.imag,modulus,scale=50)#,headlength=5, headwidth=5, latlon=False,
            # units='xy', width=8, angles='xy', scale_units='xy', cmap=cm.seismic,
#qk = plt.quiverkey(Q, 0.15, 0.9, 1, '20 m/s', labelpos='W')
cb = m.colorbar()
cb.set_label('m/s')
#cb.set_clim(0,0.15)
plt.title("C (SeaCurrents) _ " + str(b[0]) + '-' + str(b[-1]+1) + " average")

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(2)
c3 = m.contourf(x,y,absA,np.arange(0,0.028,0.0001))
cb = m.colorbar()
cb.set_label('')
#cb.set_clim(np.nanmin(absA),np.nanmax(absA))
plt.title("|A| _ Wind/Ice _ " + str(b[0]) + '-' + str(b[-1]+1) + " average (YbyY)")

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

plt.figure(3)
c4 = m.contourf(x,y,ThA,np.arange(-30,8,0.5))
cb = m.colorbar()
#cb.set_cmap('bwr')
cb.set_label('degrees')
#cb.set_clim(np.nanmin(ThA),np.nanmax(ThA))
plt.title("exp(i0) _ Wind/Ice _ " + str(b[0]) + '-' + str(b[-1]+1) + " average (YbyY)")

m = Basemap(projection='npstere',boundinglat=65,lat_0=90,lat_ts=70,lon_0=-45,\
            resolution='l',round=False)
m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')

#plt.show()

