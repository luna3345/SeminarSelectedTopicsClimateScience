# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 20:04:21 2022

@author: lule
"""
import numpy as np
import xarray as xr
import glob as glob
import pandas as pd
import csv

#https://stackoverflow.com/a/51715491
def checkerboard(shape):
    return np.indices(shape, dtype='object').sum(axis=0) % 2


# opens the netcdf file
ds = xr.open_dataset("data_ganzEuropa/ECMWF_IFS-ea-0001_fc_F32-T42_yggdrasil-0.7.2_altitude_2008-07-15T00_1.5x3.0_europe.nc")

# converts into a panda data frame
new = ds.to_dataframe()


#parameters of the grid
dlat = 1.5
dlong = 3
latmin = 37.8
latmax = 52.2
longmin = -7
longmax = 18


#definining the Richardson Grid
lat = np.arange(latmin,latmax, dlat)
long = np.arange(longmin,longmax+dlong,dlong)

latP = np.arange(latmin, latmax-dlat, 2*dlat) 
latM = np.arange(latmin+dlat,latmax+dlat,2*dlat) 

longP = np.arange(longmin,longmax+dlong,2*dlong)
longM = np.arange(longmin+dlong,longmax+dlong,2*dlong)


heights = np.array([0,2000,4200,7200,11800.0])
heights_dyn = np.array([11.543, 7.048, 4.113, 1.959, 0.0 ])*1000 #..m geodynamic height of standard levels
p_bar = np.array([200, 400, 600, 800, 1013])*100 #in Pa, std pressure levels
p_bar_uv = np.array([100,300,500,700,900])*100 # in Pa, for U&V btw pressure levls

#reducing the grid size to the size of the given grid to reduce the amount of variables
dsp = ((ds.sel(longitude=long,method = "nearest")).sel(latitude=lat,method = "nearest")).sel(altitude=heights, method = "nearest").isel(time=0)
new1 = dsp.to_dataframe() #only dataframes can be shown in spyder variable explorer


#                               Setting up input files for the fortran program
###################################################################################################
#read surface pressure and convert to 10xmmHg (value in Pa times 0.0075006156130264)
surf_p = dsp.pressure_2_meter
surf_p.data = surf_p.data*10*0.0075006156130264
surf_p.attrs['units'] = '10*mmHg'

    
#creating a checkerboard, where every second entry into the grid is a blank tab, e.g.
# 010                         340.0
# 101 translates to :    400.0     500.0
a = checkerboard([surf_p.sizes['latitude'],surf_p.sizes['longitude']])
for i in range(a.shape[0]):
    for j in range(a.shape[1]):
        if a[i,j]==0:
            a[i,j] = '\t'
        else:
            a[i,j]=str(surf_p.isel(latitude=i,longitude=j).data[()]) 
    

#save data to txt file for PS.dat
np.savetxt('PS.csv',a,fmt='%.4s',delimiter='',header='DDR: Surface pressure (mm Hg X 10)')



#######################################################################################################
#Z.dat gives geopotential height at 200mb,400mb,600mb,800mb
dsz = ((ds.sel(longitude=longP,method = "nearest")).sel(latitude=latP,method = "nearest")).isel(time=0) #reduce size to grid of P points
lp = ds.sel(longitude=longP,method = "nearest").longitude.data #save values for lat in this array
lap = ds.sel(latitude=latP,method = "nearest").latitude.data #save values for long in this array
dsz_df = dsz.to_dataframe() #convert to pandas dataframe for easier handling

#create array to hold geopotential values for matching pressure surfaces
geopot= np.zeros((4,dsz.sizes['latitude'], dsz.sizes['longitude']))

#go through each lat|long point and search for altitude where pressure=p_bar, then save geopotential height for that
#find nearest value in pandas dataframe (first answer):
#https://codereview.stackexchange.com/questions/204549/lookup-closest-value-in-pandas-dataframe
for myp in range(p_bar.size-1):
    for mylat in range(lap.size):
        for mylon in range(lp.size):
                sub = dsz_df.loc[lap[mylat],lp[mylon]]
                idx = sub['pressure'].sub(p_bar[myp]).abs().idxmin() #find index of current pressure surface
                geopot[myp,mylat,mylon] = sub.loc[idx].geopotential_height  #save that geopotential height to array

for i in range(geopot.shape[0]):
    #write to csv file
    #first open the csv on write access, to delete any previous contents
    # for i>0, open on append access, to append data block to csv file, without deleting the previous contents
    if i == 0:
        with open('Z.csv', 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f, delimiter =" ")
        
            # write the header
            header = 'Geopotential height at %.4s mb' %(p_bar[i]/100)
            writer.writerow([header])
            
            # write the data
            writer.writerows(np.int32(geopot[i,:,:]))
        
    else:
        with open('Z.csv', 'a', encoding='UTF8', newline='') as f:
            writer = csv.writer(f, delimiter =" ")
        
            # write the header
            header = 'Geopotential height at %.4s mb' %(p_bar[i]/100)
            writer.writerow([header])
            
            # write the data
            writer.writerows(np.int32(geopot[i,:,:]))


####################################################################################################
#wind speed and direction at 100mb,300mb,500mb,700mb, 900mb at M points!
#only latitude changes 
dsv = ((ds.sel(longitude=longM,method = "nearest")).sel(latitude=latM,method = "nearest")).isel(time=0) #reduce size to grid of P points
lp = ds.sel(longitude=longM,method = "nearest").longitude.data #save values for lat in this array
lap = ds.sel(latitude=latM,method = "nearest").latitude.data #save values for long in this array
dsv_df = dsv.to_dataframe() #convert to pandas dataframe for easier handling

#create array to hold uv values for matching pressure surfaces
u = np.zeros((5,dsv.sizes['latitude'], dsv.sizes['longitude']))
v = np.zeros((5,dsv.sizes['latitude'], dsv.sizes['longitude']))
dddff = np.zeros((5,dsv.sizes['latitude'], dsv.sizes['longitude'])) #wind angle(ddd) and speed(ff)


#go through each lat|long point and search for altitude where pressure=p_bar, then save windspeed at that altitude&pressure
for myp in range(p_bar_uv.size):
    for mylat in range(lap.size):
        for mylon in range(lp.size):
                sub = dsv_df.loc[lap[mylat],lp[mylon]]
                idx = sub['pressure'].sub(p_bar_uv[myp]).abs().idxmin()
                u[myp,mylat,mylon] = sub.loc[idx].eastward_wind
                v[myp,mylat,mylon] = sub.loc[idx].northward_wind
                
                #calculate wind speed and turn unit into knots
                wind_speed = (np.sqrt(u[myp,mylat,mylon]**2 + v[myp,mylat,mylon]**2))/0.514   #wind speed
                
                #! this is the wrong way !
                #https://stackoverflow.com/questions/21484558/how-to-calculate-wind-direction-from-u-and-v-wind-components-in-r
                #wind_dir = 270-(np.arctan2(u[myp,mylat,mylon], v[myp,mylat,mylon])*180/np.pi) #wind angle
                
                #correct wind direction (from https://confluence.ecmwf.int/pages/viewpage.action?pageId=133262398)
                wind_dir = np.mod(180+180/np.pi*np.arctan2(u[myp,mylat,mylon], v[myp,mylat,mylon]),360)
                
                #combine it into degrees in front of the . and speed after the .
                #which is the format used in the fortran program
                dddff[myp,mylat,mylon] = int(wind_dir)+wind_speed/100

for i in range(dddff.shape[0]):
    #write to csv file
    #first open the csv on write access, to delete any previous contents
    # for i>0, open on append access, to append data block to csv file, without deleting the previous contents

    if i == 0:
        #the first time open in write mode, so file will be cleared
        with open('UV.csv', 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f, delimiter =" ")
        
            # write the header
            header = 'DDR: Wind Speed and Direction at %.4s mb' %(p_bar_uv[i]/100)
            writer.writerow([header])
            
            # write the data
            writer.writerows(np.round(dddff[i,:,:],2))
    else: 
        with open('UV.csv', 'a', encoding='UTF8', newline='') as f:
            writer = csv.writer(f, delimiter =" ")
        
            # write the header
            header = 'DDR: Wind Speed and Direction at %.4s mb' %(p_bar_uv[i]/100)
            writer.writerow([header])
            
            # write the data
            writer.writerows(np.round(dddff[i,:,:],2))



########################################################################################################################
#height of 100mbar-200mbar surface, fortran program calculates stratospheric temperature from that              
dst = ((ds.sel(longitude=longP,method = "nearest")).sel(latitude=latP,method = "nearest")).isel(time=0) #reduce size to grid of P points
lp = ds.sel(longitude=longP,method = "nearest").longitude.data #save values for lat in this array
lap = ds.sel(latitude=latP,method = "nearest").latitude.data #save values for long in this array
dst_df = dst.to_dataframe() #convert to pandas dataframe for easier handling

#create array to hold thickness for pressure surfaces 100mb-200mb
t= np.zeros((dst.sizes['latitude'], dst.sizes['longitude']))

#go through each lat|long point and search for altitude where pressure=p_bar, then save geopotential height for that
#find nearest value in pandas dataframe (first answer):
#https://codereview.stackexchange.com/questions/204549/lookup-closest-value-in-pandas-dataframe
for mylat in range(lap.size):
    for mylon in range(lp.size):
            sub = dst_df.loc[lap[mylat],lp[mylon]]
            idx = sub['pressure'].sub(10000).abs().idxmin() #find index of 100mb pressure surface
            idx2 = sub['pressure'].sub(20000).abs().idxmin() #find index of 200mb pressure surface
            t[mylat,mylon] = sub.loc[idx].geopotential_height-sub.loc[idx2].geopotential_height

#save to csv
np.savetxt('myT.csv',t,fmt='%.4s',delimiter=' ',header='DDR: Thickness of 200-100mb layer, for Temperature.')           

#########################################################################################################################
##read output from Fortran program from csv file into panda dataframe
##and then save this dataframe as a csv, with latitude and longitude values
#rowarray = np.zeros(shape=(10,9))
#with open('pend.csv') as csv_file:
#    csv_reader = csv.reader(csv_file, delimiter='\t')
#    line_count = 0
#    for row in csv_reader:
#        rowarray[line_count] = row[0].split() #read values line by line and split each line into a list at whitespace
#        line_count += 1
#    print(f'Processed {line_count} lines.')
#
#surf_p_end = surf_p.to_dataframe()
#surf_p_end = surf_p_end.drop('time',axis=1) #since time is the same for every value (same day), drop this information
#for mylat in range(surf_p.sizes['latitude']):
#    for mylon in range(surf_p.sizes['longitude']):
#        #the following line replaces the data in the existing dataframe with the wanted data
#        #this is (seems) the easiest way to generate a correct dataset with lat/lon values etc.
#        surf_p_end.at[surf_p.longitude.data[mylon],surf_p.latitude.data[mylat]]=rowarray[mylat,mylon]
#
#surf_p_end.to_csv('p_end.csv',sep=';',float_format='%.2f') #write data out to a csv, with lat/lon values for each data entry

print('all done')     