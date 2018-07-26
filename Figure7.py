#!/usr/bin/env python

from obspy.core import read, UTCDateTime, Stream
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import obspy.signal
import numpy as np
from obspy.signal.invsim import evalresp
debug = True
import matplotlib as mpl

# Importing and applying font
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)

# Seismic station
net1 = 'IU'
sta1 = 'ANMO'
loc1 = '00'
chan1 = 'LH'

# Weather station
net2 = 'IU'
sta2 = 'ANMO'
loc2 = '50'
chan2 = 'LWS'

# Establishing start and end times for the figure
stime = UTCDateTime('2018-176T07:00:00')
etime = UTCDateTime('2018-182T00:00:00')

# Set current time to be the start time and define st as a Stream
ctime = stime
st = Stream()

# Read in data for the selected time window
while ctime <= etime:
    path1 = '/tr1/telemetry_days/' + net1 + '_' + sta1 + '/' + str(ctime.year) + '/' + str(ctime.year) \
        + '_' + str(ctime.julday).zfill(3) + '/' + loc1 + '_' + chan1 + '*'
    path2 = '/tr1/telemetry_days/' + net2 + '_' + sta2 + '/' + str(ctime.year) + '/' + str(ctime.year) \
        + '_' + str(ctime.julday).zfill(3) + '/' + loc2 + '_' + chan2 + '*'

    st += read(path1)
    st += read(path2)
    ctime += 24.*60.*60.
    
    # Merge traces and if there's a data gap, make it zero
    st.merge(fill_value = 0)

if debug:
    print(st)

# Trim stream between start and end time
st.trim(stime, etime)

# Reading in responses - responses for the XX network isn't yet with IRIS
if net1 == 'XX':
	resppath = ''

else:
	resppath = '/APPS/metadata/RESPS/'

# Removing response from seismic data traces
for tr in st.select(channel = 'LH*'):
    dic = {'filename': resppath + 'RESP.' + tr.id, 'date': tr.stats.starttime, 'units': 'VEL'}
    tr.simulate(seedresp = dic)
    tr.filter('bandpass',freqmin=1./100., freqmax=1./50.)
    tr.taper(0.05)
    tr.data *= 10**9

# Removing wind speed response and converting trace to an array       
for tr in st.select(channel = 'LWS*'):
    tr.data = tr.data.astype(np.float64)
    tr.data *= 0.1

# Creating empty lists to add to of wind speed and horizontal and vertical seismic data
ws=[]
LHH=[]
LHZ=[]

# Calculating the RMS velocity of the wind speed and horizontal and vertical seismic data and adding it to the lists
for stW in st.slide(window_length=240., step =60.):
    ws.append(np.sqrt(np.sum((stW.select(channel="LWS")[0].data)**2)/float(stW.select(channel="LWS")[0].stats.npts)))
    temp = np.sum((stW.select(channel="LH1")[0].data)**2 + (stW.select(channel="LH2")[0].data)**2)
    temp *= 1./(float(stW.select(channel="LH1")[0].stats.npts + stW.select(channel="LH2")[0].stats.npts))
    temp = np.sqrt(temp)
    LHH.append(temp)
    LHZ.append(np.sqrt(np.sum((stW.select(channel="LHZ")[0].data)**2)/float(stW.select(channel="LHZ")[0].stats.npts)))

# Establishing times as the independent variable to plot the weather and seismic data against, and converting to an array        
times = range(len(LHZ))
times = np.asarray(times)*60. 
times *= 1./(60.*60.)

# Converting data lists to arrays
ws = np.asarray(ws)
LHH = np.asarray(LHH)
LHZ = np.asarray(LHZ)

# Setting the maxiumum velocity for seismic data to in the first subplot to 150 nm/s to remove the outliers from an earthquake in Mexico
LHH[(LHH >= 150.)] = 150.
LHZ[(LHZ >= 150.)] = 150.

# Creating figure and first subplot, which is the time series of RMS seismic and wind speed data
fig = plt.figure(1, figsize=(12,12))

# Plotting the time series of wind speed data and adding labels, adjusting appearance, etc.
ax1 = plt.subplot(2,1,1)
ax1.text(-0.07, 1.05, '(a)', transform=ax1.transAxes,
      fontsize=26, fontweight='bold', va='top', ha='right')
ax1.set_xlabel('Time (hr)', fontsize = 19)
plt.tick_params(axis = 'y', labelsize = 15)
plt.tick_params(axis = 'x', labelsize = 15)
ax1.set_ylabel('RMS Wind Velocity (m/s)', fontsize = 19)
ax1.plot(times, ws, color = 'C0', label = "Wind Speed", linewidth = 0.8)

# Plotting the time series of seismic data on a twinned axis to share the same plot with wind speed with different y-axes
ax2 = ax1.twinx()
ax2.set_ylabel('RMS Seismic Velocity (nm/s)', fontsize = 19, rotation = 270, va = 'bottom')
plt.tick_params(axis = 'y', labelsize = 15)
ax2.plot(times, LHZ, color = 'C1', label = "LHZ " + sta1 + " " + loc1, linewidth = 0.8)
ax2.plot(times, LHH, color = 'C2', label = "LHH " + sta1 + " " + loc1, linewidth = 0.8)
ax2.plot(-10., 1., color = 'C0', label = 'Wind Speed')
plt.legend(loc = 2, fontsize = 15)
plt.xlim((min(times),max(times)))
plt.xlabel('Time (hours)', fontsize = 15)
plt.ylabel('RMS Seismic Velocity (nm/s)', fontsize = 15)
# plt.title('Envelope Calculations for Seismic Data and Wind Speed', fontsize = 24, y = 1.03)

# For second subplot, limiting RMS velocity of seismic data to 80 to avoid extreme outliers when calculating fits
LHZ[(LHZ >= 80.)] = 80.
LHH[(LHH >= 80.)] = 80.

# Establishing function to calculate the linear fits between seismic and wind velocity
def function(x, b, c):
    val = x*b + c
    return val

poptH, pcovH = curve_fit(function, ws, LHH)
if debug:
    print(poptH)
    print(pcovH)
    
poptZ, pcovZ = curve_fit(function, ws, LHZ)
if debug:
    print(poptZ)
    print(pcovZ)

# Plotting RMS velocity data on subplot two as a scatter plot with lines of fit    
ax= plt.subplot(2,1,2)
ax.text(-0.07, 1.05, '(b)', transform=ax.transAxes,
      fontsize=26, fontweight='bold', va='top', ha='right')
plt.scatter(ws, LHZ, s = 5, color = 'C1', label = "Wind and Vertical")
plt.scatter(ws, LHH, s = 5, color = 'C2', label = 'Wind and Horizontal')
plt.tick_params(axis = 'y', labelsize = 15)
plt.tick_params(axis = 'x', labelsize = 15)
plt.plot(ws, function(ws, *poptH), color = 'C2', label = 'Fit: ' + str(round(poptH[0],2)) + ' m/s per nm/s')
plt.plot(ws, function(ws, *poptZ), color = 'C1', label = 'Fit: ' + str(round(poptZ[0],2)) + ' m/s per nm/s')

# Labeling and adjusting appearance of second subplot
plt.ylabel('RMS Seismic Velocity (nm/s)', fontsize = 19)
plt.xlabel('RMS Wind Velocity (m/s)', fontsize = 19)
plt.xlim((min(ws), max(ws)))
plt.legend(loc = 2, fontsize = 15)

# Saving or showing figure
plt.savefig('RMS_Velocities_' + net1 + sta1 + '_' + loc1 + '_' + chan1 + '_176-182.jpg', format = 'jpeg', dpi = 400)
plt.show()
