# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 21:35:31 2022

@author: lule
"""

# Script to take in the output from the fortran program that has been saved into a seperate csv file
# and plot the sea level pressure onto a map of europe
# adapted partially from:

#https://unidata.github.io/python-gallery/examples/HILO_Symbol_Plot.html
#and
#https://www.unidata.ucar.edu/support/help/MailArchives/python/msg00197.html

###############################
# Imports

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
import xarray as xr
import csv
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

###############################
# Function for finding and plotting max/min points


def plot_maxmin_points(lon, lat, data, extrema, nsize, symbol, color='k',
                       plotValue=True, transform=None):
    """
    This function will find and plot relative maximum and minimum for a 2D grid. The function
    can be used to plot an H for maximum values (e.g., High pressure) and an L for minimum
    values (e.g., low pressue). It is best to used filetered data to obtain  a synoptic scale
    max/min value. The symbol text can be set to a string value and optionally the color of the
    symbol and any plotted value can be set with the parameter color
    lon = plotting longitude values (2D)
    lat = plotting latitude values (2D)
    data = 2D data that you wish to plot the max/min symbol placement
    extrema = Either a value of max for Maximum Values or min for Minimum Values
    nsize = Size of the grid box to filter the max and min values to plot a reasonable number
    symbol = String to be placed at location of max/min value
    color = String matplotlib colorname to plot the symbol (and numerica value, if plotted)
    plot_value = Boolean (True/False) of whether to plot the numeric value of max/min point
    The max/min symbol will be plotted on the current axes within the bounding frame
    (e.g., clip_on=True)
    """
    from scipy.ndimage.filters import maximum_filter, minimum_filter

    if (extrema == 'max'):
        data_ext = maximum_filter(data, nsize, mode='nearest')
    elif (extrema == 'min'):
        data_ext = minimum_filter(data, nsize, mode='nearest')
    else:
        raise ValueError('Value for hilo must be either max or min')

    mxy, mxx = np.where(data_ext == data)

    for i in range(len(mxy)):
        ax.text(lon[mxy[i], mxx[i]], lat[mxy[i], mxx[i]], symbol, color=color, size=24,
                clip_on=True, horizontalalignment='center', verticalalignment='center',
                transform=transform)
        ax.text(lon[mxy[i], mxx[i]], lat[mxy[i], mxx[i]],
                '\n' + str(int(data[mxy[i], mxx[i]])),
                color=color, size=12, clip_on=True, fontweight='bold',
                horizontalalignment='center', verticalalignment='top', transform=transform)


###############################
plt.close('all')
#get data (initial data and data from runs)
        
# opens the netcdf file
#ds = xr.open_dataset("ECMWF_IFS-ea-0001_fc_F32-T42_yggdrasil-0.7.2_altitude_2008-07-16T00_1.5x3.0_europe.nc")
ds = xr.open_dataset("ECMWF_IFS-ea-0001_an_F32-T42_yggdrasil-0.7.2_altitude_2008-07-15T00_1.5x3.0_europe.nc")
                       

# converts into a panda data frame
new = ds.to_dataframe()


#parameters of the grid
dlat = 1.5
dlong = 3
latmin = 37.8
latmax = 52.2
longmin = -7
longmax = 17


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



#read surface pressure and convert to 10xmmHg (value in Pa times 0.0075006156130264)
surf_p_ds = dsp.pressure_2_meter
surf_p_ds.data = surf_p_ds.data/100
surf_p_ds.attrs['units'] = 'hPa'
    
surf_p = surf_p_ds.to_dataframe()
surf_p = surf_p.drop('time',axis=1)

#mslp = dsp.pressure_2_meter.data
mslp = gaussian_filter(dsp.pressure_2_meter.data, sigma=1.0)
mslp = np.transpose(mslp)


##########################################################################################################################

#read output from Fortran program from csv file into panda dataframe
rowarray = np.zeros(shape=(10,9))
with open('pend2.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    line_count = 0
    for row in csv_reader:
        rowarray[line_count] = row[0].split()
        line_count += 1
    print(f'Processed {line_count} lines.')

surf_p_end = surf_p_ds.to_dataframe()
surf_p_end = surf_p_end.drop('time',axis=1) #since time is the same for every value (same day), drop this information
for lat in range(surf_p_ds.sizes['latitude']):
    for lon in range(surf_p_ds.sizes['longitude']):
        #the following line replaces the data in the existing dataframe with the wanted data
        #this is (seems) the easiest way to generate a correct dataset with lat/lon values etc.
        #dataframe.at is how you access a value to overwrite it, (opposed to .sel etc)
        surf_p_end.at[surf_p_ds.longitude.data[lon],surf_p_ds.latitude.data[lat]]=rowarray[lat,lon]

surf_p_end_ds = surf_p_end.to_xarray()
mslp2 = surf_p_end_ds.pressure_2_meter.data
mslp2 = gaussian_filter(mslp2, sigma=1.2) #smooth data to reduce count of contours in the plot
mslp2 = np.transpose(mslp2) #because of the HiLow function (plot_maxmin) the values have to be transposed (switches lat/lon as first index)
       
################################################### PLOTTING THE MAP(S) ################################################
# Set map and data projections for use in mapping

# Set projection of map display
mapproj = ccrs.LambertConformal(central_latitude=50., central_longitude=11.)

# Set projection of data
dataproj = ccrs.PlateCarree()

# Grab data for plotting state boundaries
states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')

######################################################################################################################
# Create figure and plot data
fig = plt.figure(1, figsize=(17., 11.))
ax = plt.subplot(111, projection=mapproj)

# Set extent and plot map lines
ax.set_extent([-10., 25, 35, 56.], ccrs.PlateCarree())
ax.coastlines('50m', edgecolor='black', linewidth=0.75)
ax.add_feature(states_provinces, edgecolor='black', linewidth=0.5)

# Plot thickness with multiple colors
kw_clabels = {'fontsize': 11, 'inline': True, 'inline_spacing': 5, 'fmt': '%i',
              'rightside_up': True, 'use_clabeltext': True}

# Plot MSLP
clevmslp = np.arange(900., 1100., 5)
XX,YY = np.meshgrid(dsp.longitude.data,dsp.latitude.data)
cs2 = ax.contour(XX,YY, mslp, clevmslp, colors='k', linewidths=1.25,
                 linestyles='solid',transform=dataproj)
plt.clabel(cs2, **kw_clabels)
ax.contourf(XX,YY, mslp,transform=dataproj)

# Use definition to plot H/L symbols
plot_maxmin_points(XX,YY, mslp, 'max', 8, symbol='H', color='b',  transform=dataproj)
plot_maxmin_points(XX, YY, mslp, 'min', 5, symbol='L', color='r', transform=dataproj)

# Put on some titles
plt.title('Surface pressure (hPa) with Highs and Lows', loc='left')
plt.title('VALID: {}'.format(ds.time.data), loc='right')


#add axis labels 
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.2, linestyle='--',x_inline=False, y_inline=False)
gl.xlabels_top = False
#gl.ylabels_left = False
gl.xlabel_style = {'rotation': 0}

#since gridlines is using xlabel, we have to add axis labels manually
ax.text(-0.07, 0.55, 'Latitude', va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax.transAxes)
ax.text(0.5, -0.05, 'Longitude', va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes)

plt.show()


# Create figure and plot data
fig = plt.figure(2)
ax = plt.subplot(111, projection=mapproj)

# Set extent and plot map lines
ax.set_extent([-10., 25, 35, 56.], ccrs.PlateCarree())
ax.coastlines('50m', edgecolor='black', linewidth=0.75)
ax.add_feature(states_provinces, edgecolor='black', linewidth=0.5)
# Plot MSLP
clevmslp = np.arange(900., 1100., 5)
XX,YY = np.meshgrid(dsp.longitude.data,dsp.latitude.data)
cs2 = ax.contour(XX,YY, mslp2, clevmslp, colors='k', linewidths=1.25,
                 linestyles='solid',transform=dataproj)
plt.clabel(cs2, **kw_clabels)
ax.contourf(XX,YY, mslp2,transform=dataproj)

# Use definition to plot H/L symbols
plot_maxmin_points(XX,YY, mslp2, 'max', 8, symbol='H', color='b',  transform=dataproj)
plot_maxmin_points(XX, YY, mslp2, 'min', 5, symbol='L', color='r', transform=dataproj)

# Put on some titles
plt.title('Simulated Surface pressure (hPa) with Highs and Lows', loc='left')
plt.title('VALID: 2008-07-16', loc='right')

#add axis labels 
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.2, linestyle='--',x_inline=False, y_inline=False)
gl.xlabels_top = False
#gl.ylabels_left = False
gl.xlabel_style = {'rotation': 0}

#since gridlines is using xlabel, we have to add axis labels manually
ax.text(-0.07, 0.55, 'Latitude', va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax.transAxes)
ax.text(0.5, -0.05, 'Longitude', va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes)


plt.show()

########################################################################################################
#plot differences in pressure distribution

# Create figure and plot data
fig = plt.figure()
ax = plt.subplot(111, projection=mapproj)

# Set extent and plot map lines
ax.set_extent([-10., 25, 35, 56.], ccrs.PlateCarree())
ax.coastlines('50m', edgecolor='black', linewidth=0.75)
ax.add_feature(states_provinces, edgecolor='black', linewidth=0.5)

# Plot thickness with multiple colors
kw_clabels = {'fontsize': 11, 'inline': True, 'inline_spacing': 5, 'fmt': '%i',
              'rightside_up': True, 'use_clabeltext': True}

# Plot MSLP
clevmslp = np.arange(-40., 20., 5)
XX,YY = np.meshgrid(dsp.longitude.data,dsp.latitude.data)
cs2 = ax.contour(XX,YY, (mslp-mslp2), clevmslp, colors='k', linewidths=1.25,
                 linestyles='solid',transform=dataproj)
plt.clabel(cs2, **kw_clabels)
ax.contourf(XX,YY, (mslp-mslp2),transform=dataproj)

# Put on some titles
plt.title('Surface pressure difference (hPa) between Reanalysis and Simulation', loc='left')


#add axis labels 
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.2, linestyle='--',x_inline=False, y_inline=False)
gl.xlabels_top = False
#gl.ylabels_left = False
gl.xlabel_style = {'rotation': 0}

#since gridlines is using xlabel, we have to add axis labels manually
ax.text(-0.07, 0.55, 'Latitude', va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax.transAxes)
ax.text(0.5, -0.05, 'Longitude', va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes)

plt.show()

#########################################################################################################################
#PLOT GRIDPOINTS USED FOR SIMULATION & PLOTTING
# Create figure and plot data
fig = plt.figure()
ax = plt.subplot(111, projection=mapproj)

# Set extent and plot map lines
ax.set_extent([-10., 25, 35, 56.], ccrs.PlateCarree())
ax.coastlines('50m', edgecolor='black', linewidth=0.75)
ax.add_feature(states_provinces, edgecolor='black', linewidth=0.5)

XXred = XX[0:-1:2]

# the transform=ccrs.PlateCarree() is very important!! it tells the plot to go from lat/lon to position in plot
ax.scatter(XX,YY,s = 10, c = 'b', edgecolor = 'black', transform=ccrs.PlateCarree())


# Put on some titles
plt.title('Grid points', loc='left')

#add axis labels 
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.2, linestyle='--',x_inline=False, y_inline=False)
gl.xlabels_top = False
#gl.ylabels_left = False
gl.xlabel_style = {'rotation': 0}

#since gridlines is using xlabel, we have to add axis labels manually
ax.text(-0.07, 0.55, 'Latitude', va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax.transAxes)
ax.text(0.5, -0.05, 'Longitude', va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes)


plt.show()