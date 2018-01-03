# -*- coding: utf-8 -*-
"""
Created on Sun Aug 28 15:28:09 2016

Computes coefficients C (sea currents) and A (modulus:|A| and angle:theta)   
from wind and sea ice drift data as an inverse problem

Plots them on a north polar stereographic basemap.

Vectors converted to complex conjugated variables
"""

import glob
import numpy as np
import numpy.ma as ma
from numpy.linalg import inv
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap

b = range(2014,2015)
#B = 'Summer'
B = 'Winter'
#B = ''

CC = np.zeros([21063,len(b)],dtype=np.complex)    # 2D array to save A and C values 
AA = np.zeros([21063,len(b)],dtype=np.complex)    # as they're computed for each grid point (rows) for every year (colunmns)

m = -1          # initial value of counter for years

for i in b:         # loop over year by year
    PIP = str(i)
    POP = str(i+1)
    PEP = PIP + '-' + POP
    print ('PEP: ', PEP)

    m += 1          # time counter for years
    
    filenames1 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Javier_ERAI_winds_icedrift_OSISAF/SortedData_OSISAF_IceWind/WindSpeed_10m/' + PEP + '/' + B + '/*.nc')
    filenames2 = glob.glob('C:/Users/admin/Desktop/Master_thesis/Javier_ERAI_winds_icedrift_OSISAF/SortedData_OSISAF_IceWind/SeaIceDriftSpeed/' + PEP + '/' + B + '/*.nc')
    
    Wu = np.zeros([21063,len(filenames1)],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    Wv = np.zeros([21063,len(filenames1)],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values as they're extracted from each .nc file (colunmns)
    Iu = ma.zeros([21063,len(filenames2)],dtype=np.float)    # 2D array to save LOCAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)
    Iv = ma.zeros([21063,len(filenames2)],dtype=np.float)    # 2D array to save MERIDIONAL (X) ICE VELOCITY values as they're extracted from each .vec file (colunmns)
        
    n = -1       # initial value of counter for iceDrift files
    nn = -1        # initial value of counter for wind files
    
    for filename in filenames1:                         # extracts data from each .nc file one by one
        datum1 = Dataset(filename, mode='r')
        
        print ('filename: ', filename)
        
        lat1 = datum1.variables['lat1'][:]
        lon1 = datum1.variables['lon1'][:] 
        uw = datum1.variables['u_wind_avg'][:]
        vw = datum1.variables['v_wind_avg'][:]
    
        nn += 1
               
        Wu[:,nn] = uw.flat           
        Wv[:,nn] = vw.flat      

    for filename in filenames2:                         # extracts data from each .nc file one by one
        datum2 = Dataset(filename, mode='r')
            
        print ('filename: ', filename)
        
        lat = datum2.variables['lat'][:]       
        lon = datum2.variables['lon'][:]
        ui = datum2.variables['dX'][:]
        
        try:
            vi = datum2.variables['dY_v1p4'][:]
        except KeyError:
            vi = datum2.variables['dY'][:]
        
        n += 1
        
        ui = ui[0,:,:]*(1000.0/(2*86400))
        vi = vi[0,:,:]*(1000.0/(2*86400))
        
        Iu[:,n] = ma.ravel(ui)             # saves LOCAL (X) ICE VELOCITY values in columns 
        Iv[:,n] = ma.ravel(vi)             # saves MERIDIONAL (Y) ICE VELOCITY values in columns 

    C = np.zeros([21063,],dtype=np.complex)    # 2D array of zeros to save A nad C values as they're computed for each year (colunmns)
    A = np.zeros([21063,],dtype=np.complex)    
    
    C[:,].real = np.nan         # changes zeros to nans
    C[:,].imag = np.nan
    A[:,].real = np.nan
    A[:,].imag = np.nan
    
    for element in range(21063):         # loop over every gridpoint
        
        print ('element: ', element)
        
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

        if len(IC) > 10:                # minimum lengh of time series to ???
            IC = IC[:,np.newaxis]                             # transforms the horizontal array into a vertical one
            WC = [WC,np.ones_like(WC,dtype=np.complex)]       # attaches an array of ones to the list of complex velocity values
            WC = np.mat(WC)                                   # converts a list of a 1D list and 1D array into a 2 row matrix
            WC = WC.T                                         # transposes the matrix into a 2 column one (1 x 2)
            WCH = WC.getH()                                   # transforms it into a Hermitian (conjugate transpose) matrix
                 
            AC = (inv(WCH*WC)*WCH*IC)         # computes the inversion // gives a 1x2 matrix with values of parameter A and sea currents C (in complex form)
                    
            C[element,] = AC[1,0]       # appens the SEA CURRENTS VELOCITY values to a list as they're generated
            A[element,] = AC[0,0]       # appens the A parameter values to a list as they're generated

    CC[:,m] = C        # saves LOCAL (X) ICE VELOCITY values in columns 
    AA[:,m] = A        # saves MERIDIONAL (Y) ICE VELOCITY values in columns 

CCC = np.zeros([21063,],dtype=complex)
AAA = np.zeros([21063,],dtype=complex)

AAA[:,].real = np.nan
AAA[:,].imag = np.nan

for ii in range(21063):                  # loops over every grid point to compute the mean over the whole time series
    CCC[ii] = np.nanmean(CC[ii,:])
    AAA[ii] = np.nanmean(AA[ii,:])
    
Mean = np.nanmean(AAA)      
StD = np.nanstd(AAA)

MEANabs = abs(Mean)
MEANang = np.angle(Mean, deg=1)
STDabs = abs(StD)
STDang = np.angle(StD, deg=1)

AAA =np.reshape(AAA,(177,119))             # reshapes the A parameter 1D array into 2D for plotting with contourf
absA = abs(AAA)                           # |A| informs of coupling btw wind/ice (and internal ice stresses)
ThA = np.angle(AAA,deg=True)              # angle btw wind/ice motion vectors

modulus = abs(CCC)
#modulus = np.sqrt(C.real**2 + C.imag**2)*(1000/(2*86400))

c = CCC/modulus          # normalized vectors

#############################
### Plotting with Basemap ###
#############################

m = Basemap(projection='npstere',boundinglat=61.5,lat_0=90,lat_ts=70,lon_0=-35,\
            resolution='l',round=False)
#m.etopo()
#m.drawcoastlines()
m.fillcontinents(color='grey',lake_color='aqua')
m.drawparallels(np.arange(-80.,81.,10.))
m.drawmeridians(np.arange(-180.,181.,20.))
#m.drawmapboundary(fill_color='aqua')

#lat1=np.flipud(lat1)
#lon1=np.flipud(lon1)
   
x,y = m(lon1,lat1,inverse=False)#,indexing='ij')

plt.figure(1)
Q = m.quiver(x,y,c.real,c.imag,modulus,scale=50)#,headlength=5, headwidth=5, latlon=False,
            # units='xy', width=8, angles='xy', scale_units='xy', cmap=cm.seismic,
#qk = plt.quiverkey(Q, 0.15, 0.9, 1, '20 m/s', labelpos='W')
cb = m.colorbar()
cb.set_label('m/s')
cb.set_clim(0,0.135)
#plt.title("C (SeaCurrents) _ Oct " + str(b[0]) + ' - Sep ' + str(b[-1]+1) + " _ average")
plt.title("C (SeaCurrents) _ " + B + ' ' + str(b[0]) + '-' + str(b[-1]+1) + " _ average")
#plt.title("C (SeaCurrents) _ Summer " + str(b[0]+1) + " _ average")
#plt.title("C (SeaCurrents) _ Summer " + str(b[0]+1) + '-' + str(b[-1]+1) + " _ average")

m = Basemap(projection='npstere',boundinglat=61.5,lat_0=90,lat_ts=70,lon_0=-35,\
            resolution='l',round=False)
#m.etopo()
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')
m.drawparallels(np.arange(-80.,81.,10.))
m.drawmeridians(np.arange(-180.,181.,20.))

plt.figure(2)
c2 = m.contourf(x,y,absA,np.arange(0.0001, 0.025, 0.0001))
cb = m.colorbar()
cb.set_label('')
#cb.set_clim(np.nanmin(absA),np.nanmax(absA))
#plt.title("|A| _ Wind/Ice _ Oct " + str(b[0]) + ' - Sep ' + str(b[-1]+1) + " _ average")
plt.title("|A| _ Wind/Ice _ " + B + ' ' + str(b[0]) + '-' + str(b[-1]+1) + " _ average")
#plt.title("|A| _ Wind/Ice _ Summer " + str(b[0]+1) + " _ average")
#plt.title("|A| _ Wind/Ice _ Summer " + str(b[0]+1) + '-' + str(b[-1]+1) + " _ average")

m = Basemap(projection='npstere',boundinglat=61.5,lat_0=90,lat_ts=70,lon_0=-35,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')
m.drawparallels(np.arange(-80.,81.,10.))
m.drawmeridians(np.arange(-180.,181.,20.))

plt.figure(3)
c3 = m.contourf(x,y,ThA,np.arange(-45.0, 5.0, 0.05))
cb = m.colorbar()
cb.set_label('degrees')
#cb.set_clim(np.nanmin(ThA),np.nanmax(ThA))
#plt.title("exp(i0) _ Wind/Ice _ Oct " + str(b[0]) + ' - Sep ' + str(b[-1]+1) + " _ average")
plt.title("angA _ Wind/Ice _ " +  B + ' ' + str(b[0]) + '-' + str(b[-1]+1) + " _ average")
#plt.title("exp(i0) _ Wind/Ice _ Summer " + str(b[0]+1) + " _ average")
#plt.title("exp(i0) _ Wind/Ice _ Summer " + str(b[0]+1) + '-' + str(b[-1]+1) + " _ average")

m = Basemap(projection='npstere',boundinglat=61.5,lat_0=90,lat_ts=70,lon_0=-35,\
            resolution='l',round=False)
#m.drawcoastlines()
m.fillcontinents(color='grey')#,lake_color='aqua')
m.drawparallels(np.arange(-80.,81.,10.))
m.drawmeridians(np.arange(-180.,181.,20.))

#plt.show()

