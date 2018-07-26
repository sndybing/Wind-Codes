#!/usr/bin/env python

from obspy.core import read, UTCDateTime, Stream
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import obspy.signal
import numpy as np
from obspy.signal.invsim import evalresp
debug = True
import matplotlib as mpl
from scipy import signal
from itertools import combinations

# Importing and applying font
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)

# Establishing windowing information and the frequency band we're looking at
nfft = 2**9
fmin = 1./20.
fmax = 1./8.
ws = []
wl = 14400.
stepsize = 3600.

# Weather station
wnet = 'XX'
wsta = 'FBA1'
wloc = '50'
wchan = 'LWS'

# Listing seismic stations to be used
netstas = ['GS_ASL8_00', 'GS_ASL9_00', 'GS_ASL9_10', 'XX_ENG5_60']

# Establishing start and end times for the figure
stime = UTCDateTime('2018-173T18:00:00')
etime = UTCDateTime('2018-204T15:00:00')

# Set current time to be the start time and define Streams for seismic and wind data
ctime = stime
st = Stream()
wst = Stream()

# Creating figure and moving subplots together
fig = plt.figure(1, figsize=(12,12))
plt.subplots_adjust(hspace=0.001)

# Reading in wind speed data
while ctime <= etime:
    wpath = '/tr1/telemetry_days/' + wnet + '_' + wsta + '/' + str(ctime.year) + '/' + str(ctime.year) + '_' + str(ctime.julday).zfill(3) + '/' + wloc + '_' + wchan + '*'
    wst += read(wpath)
    ctime += 24.*60.*60.

# Merge traces and if there's a data gap, make it zero. Then trim stream between start and end time	        
wst.merge(fill_value = 0)
wst.sort()
wst.trim(stime, etime)

# Merge traces and if there's a data gap, make it zero. Then trim stream between start and end time	
for tr in wst:
    tr.data = tr.data.astype(np.float64)
    tr.data *= 0.1

if debug:
    print(wst)

# Adding wind data to list established above with windowing
for stW in wst.slide(window_length=wl, step=stepsize):
    ws.append(np.mean(stW.select(channel="LWS")[0].data))

if debug:
    print(len(ws))

# Turning wind data into an array and establishing the times (independent) variable to plot the wind
ws = np.asarray(ws)
timesws = range(len(ws))
timesws = np.asarray(timesws)*stepsize
timesws *= 1./(60.*60.)

# Creating empty lists of handles and labels to add to and eventually use in creating the legend
handles = []
labels = []

# Plotting seismic coherences of LHZ data on different subplots
for idx, chan in enumerate(['LHZ']):
    st = Stream()
    
    # Reading in LHZ data for all stations
    for netsta in netstas:
        net, sta, loc = netsta.split('_')
        ctime = stime
        while ctime <= etime:
            path = '/tr1/telemetry_days/' + net + '_' + sta + '/' + str(ctime.year) + '/' + str(ctime.year) \
                + '_' + str(ctime.julday).zfill(3) + '/' + loc + '_' + chan + '*'
            st += read(path)
            ctime += 24.*60.*60.
    
    # Merge traces and if there's a data gap, make it zero. Then trim stream between start and end time      
	st.merge(fill_value = 0)
	st.sort()
	st.trim(stime, etime)
    
    if debug:
        print(st)
    
    # Creating list of colors for the traces
    colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
    
    # Calculating coherence with windowing between each possible combination of two stations and plotting on subplots
    for idx2, comb in enumerate(combinations(range(len(netstas)),2)):
        ax1 = fig.add_subplot(6,1,idx2+1)
        
        # Creating empty list of coherences to add to
        coh = []
        for stW in st.slide(window_length=wl, step=stepsize):
            f, Cxy = signal.coherence(stW[comb[0]].data, stW[comb[1]].data, fs = 1., nperseg = nfft, noverlap = int(0.5*nfft))
            meanc = np.mean(Cxy[(f>=fmin)&(f<=fmax)])
            coh.append(meanc)
        
        # Establishing times (independent) variable to plot coherences against and converting coh to an array
        times = range(len(coh))
        times = np.asarray(times)*stepsize
        times *= 1./(60.*60.)
        coh= np.asarray(coh)
        
        # Creating labels and adding to subplots
        lab = st[comb[0]].stats.station + ' ' + st[comb[0]].stats.location + ' and '
        lab += st[comb[1]].stats.station + ' ' + st[comb[1]].stats.location
        ax1.plot(times, coh, label = lab, color = colors[idx2], linewidth = 1.3)
        plt.yticks([0.1, 0.5, 0.9])
        plt.tick_params(axis = 'y', labelsize = 15)
        
        # Adding x-ticks and label only on the bottom subplot
        if idx2 == 5:
            plt.tick_params(axis = 'x', labelsize = 15)
            plt.xticks([0.,100.,200.,300.,400.,500.,600.])
            plt.xlabel('Time (hr)', fontsize = 19) 
        
        # Twinning x-axis to plot wind data and coherences on different axes
        ax2 = ax1.twinx()
        ax2.plot(timesws, ws, color = 'gray', label = "Wind Speed", linewidth = 0.8)
        
        # Removing x-ticks on all plots but the last
        if idx2 != 5:
            plt.xticks([])
        
        # Adding other labels and modifying plot appearances    
        plt.yticks([2.,4.,6.])
        plt.xlim(0.,648.)
        plt.tick_params(axis = 'y', labelsize = 15)
        
        # Adding seismic data to handles and labels to make legend  
        hand, lab = ax1.get_legend_handles_labels()
        handles += hand
        labels += lab
    
    # Adding wind data to handles and labels to make legend     
    hand, lab = ax2.get_legend_handles_labels()
    handles += hand
    labels += lab
    
    # Creating centered x and y-axis labels to span all subplots
    fig.text(0.04, 0.5, 'Coherence ($\gamma^2$)', ha = 'center', va = 'center', rotation = 'vertical', fontsize = 20)
    fig.text(0.95, 0.5, 'Wind Speed (m/s)', ha = 'center', va = 'center', rotation = 270, fontsize = 20)
    
    # Creating legend centered below plots
    ax1 = plt.gca()
    leg = fig.legend(handles, labels, loc = 'lower center', ncol = 4, fontsize = 15)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)
    plt.subplots_adjust(bottom = 0.14)

# Adding subplot labels
xs =[.06, .06, .06, .06, .06, .06]
ys=[.86, .73, .61, .48, .36, .24]
letters =['a','b', 'c', 'd', 'e', 'f']
for triple in zip(xs, ys, letters):
    plt.text(triple[0], triple[1], '(' + triple[2] + ')', fontsize=24, transform=plt.gcf().transFigure)

# Saving or showing figure
plt.savefig('Coherence_TimeSeries.jpg', format="JPEG", dpi=400)
plt.show()
