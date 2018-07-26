#!/usr/bin/env python

from obspy.core import read, UTCDateTime, Stream
import matplotlib.pyplot as plt
import matplotlib as mpl
import obspy.signal
import numpy as np
from obspy.signal.invsim import evalresp
debug = True

# Import/apply font
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)

# Weather station
net1 = 'XX'
sta1 = 'FBA1'
loc1= '50'
chan1 = 'LWS'

# Seismic station
net2 = 'GS'
sta2 = 'ASL8'
loc2 = '00'
chan2 = 'LH'

# Start and end time of data for figure
stime = UTCDateTime('2018-176T07:00:00')
etime = UTCDateTime('2018-182T00:00:00')

# Set current time to be the start time and define st as a Stream
ctime = stime
st = Stream()

# Read in data for the selected time window
while ctime <= etime:
    path1 = '/tr1/telemetry_days/' + net1 + '_' + sta1 + '/' + str(ctime.year) + '/' + str(ctime.year) + '_' + str(ctime.julday).zfill(3) + '/' + loc1 + '_' + chan1 + '*'
    path2 = '/tr1/telemetry_days/' + net2 + '_' + sta2 + '/' + str(ctime.year) + '/' + str(ctime.year) + '_' + str(ctime.julday).zfill(3) + '/' + loc2 + '_' + chan2 + '*'

    st += read(path1)
    st += read(path2)
    ctime += 24.*60.*60.

# Merge traces and if there's a data gap, make it zero. Then trim stream between start and end time
st.merge(fill_value = 0.)
st.trim(stime, etime)

# Reading in responses - responses for the XX network isn't yet with IRIS
if net2 == 'XX':
	resppath = ''

else:
	resppath = '/APPS/metadata/RESPS/'

# Removing response from seismic data traces
for tr in st.select(channel = 'LH*'):
    dic = {'filename': resppath + 'RESP.' + tr.id, 'date': tr.stats.starttime, 'units': 'VEL'}
    tr.simulate(seedresp = dic)
    tr.filter('bandpass',freqmin=1./100., freqmax=1./50.)
    tr.taper(0.05)
    tr.data *= 10**6
    tr.detrend('constant')

# Removing wind speed response and converting trace to an array   
for tr in st.select(channel = 'LWS*'):
    tr.data = tr.data.astype(np.float64)
    tr.data *= 0.1

# Creating figure and moving subplots to have no space between them
fig = plt.figure(1, figsize=(12,16))
plt.subplots_adjust(hspace=0.001)

# Define labels for the legend
labels = ['Wind\nVelocity (m/s)', 'N-S Seismic\nVelocity (' + r'$\mu$' + 'm/s)', 'E-W Seismic\nVelocity (' + r'$\mu$' + 'm/s)', 'Vertical Seismic\nVelocity (' + r'$\mu$' + 'm/s)']

# Plot each trace on its own subplot, labeling axes, etc.
for idx, tr in enumerate(st):
    if tr.stats.channel == 'LWS':
        plt.subplot(411)
    elif tr.stats.channel == 'LH1':
        plt.subplot(412)
    elif tr.stats.channel == 'LH2':
        plt.subplot(413)
    else:
        plt.subplot(414)
    times = range(tr.stats.npts)
    times = np.asarray(times)/(tr.stats.sampling_rate) 
    times *= 1./(60.*60.)
    plt.plot(times, tr.data, label = tr.id, color = 'gray', linewidth = 0.7)
    plt.ylabel(labels[idx], fontsize = 19)
    leg = plt.legend(loc = 1, fontsize = 12)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)
    plt.xlim(0., 137.)
    if tr.stats.channel == 'LWS':
        plt.ylim(0., 1.1 * max(tr.data))
        plt.yticks([5.,10.,15.])
    else:
        plt.ylim(-8. * np.std(tr.data), 8. * np.std(tr.data))
        plt.yticks([-0.04,-0.02,0.,0.02,0.04])
    if (tr.stats.channel == "LH1") or (tr.stats.channel == "LH2"):
        plt.ylim((-.25,.25))
        plt.yticks([-0.2,-0.1,0.,0.1,0.2])
    plt.tick_params(labelsize = 15)
    plt.xticks([])  
        
# Subplot labels    
ax1 = plt.subplot(411)
ax1.text(-0.07, 1.05, '(a)', transform=ax1.transAxes,
      fontsize=26, fontweight='bold', va='top', ha='right')
# plt.title('Seismic and Wind Velocity', fontsize = 24, y = 1.03)
ax2 = plt.subplot(412)
ax2.text(-0.07, 1.05, '(b)', transform=ax2.transAxes,
      fontsize=26, fontweight='bold', va='top', ha='right')
ax3 = plt.subplot(413)
ax3.text(-0.07, 1.05, '(c)', transform=ax3.transAxes,
      fontsize=26, fontweight='bold', va='top', ha='right')      
ax4 = plt.subplot(414)
ax4.text(-0.07, 1.05, '(d)', transform=ax4.transAxes,
      fontsize=26, fontweight='bold', va='top', ha='right')
plt.xticks([0.,20.,40.,60.,80.,100.,120.])

# Horizontal axis labeling
plt.xlabel('Time (hr)', fontsize = 19)  
plt.tick_params(labelsize = 15)  

# Save or show figure
plt.savefig('Timeseries_' + sta2 + '_' + loc2 + '.jpg', format="JPEG", dpi=400)
plt.show()
