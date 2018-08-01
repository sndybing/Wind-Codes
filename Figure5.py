#!usr/bin/env python

from obspy.core import read, UTCDateTime, Stream
from matplotlib.mlab import csd
import numpy as np
from obspy.signal.invsim import evalresp
import matplotlib.pyplot as plt
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
import glob
import sys
import matplotlib as mpl

debug = False

# Importing and applying font
mpl.rc('font', family = 'serif')
mpl.rc('font', serif = 'Times') 
mpl.rc('text', usetex = True)

# Listing seismic stations to be used
netstas = ['GS_ASL8_00', 'GS_ASL9_00', 'GS_ASL9_10', 'XX_ENG5_60', 'IU_ANMO_00']

# Establishing start and end times for the figure
stimes = [UTCDateTime('2018-177T02:24:00'), UTCDateTime('2018-176T09:27:00')]
etimes = [UTCDateTime('2018-177T08:37:00'), UTCDateTime('2018-176T15:34:00')]

nfft = 2**12
windlap = 2**6

# Plotting low wind data and high wind data on different subplots
for idx2, sandetime in enumerate(zip(stimes, etimes)): 
    stime = sandetime[0]
    etime = sandetime[1]
    
    # Reading in LHZ data for each station     
    for idx, chan in enumerate(['LHZ']):
        st = Stream()
        for netsta in netstas:
            net, sta, loc = netsta.split('_')
            ctime = stime
            while ctime <= etime:
                path = '/tr1/telemetry_days/' + net + '_' + sta + '/' + str(ctime.year) + '/' + str(ctime.year) + '_' + str(ctime.julday).zfill(3) + '/' + loc + '_' + chan + '*'
                st += read(path)
                ctime += 24.*60.*60.
        
        # Merge traces and if there's a data gap, make it zero. Then trim stream between start and end time
        st.merge(fill_value = 0)
        st.sort()
        st.trim(stime, etime)
    
    # Creating figure and list of colors for traces    
    fig = plt.figure(1, figsize = (12,12))
    colors = ['C0', 'C1', 'C2', 'C3', 'C4']
    
    # Removing responses and calculating PSDs for each station  
    for idx, tr in enumerate(st):
        if tr.stats.network == 'XX':
            resppath = 'RESP.'
        else:
            resppath = '/APPS/metadata/RESPS/RESP.'
        power, freq = csd(tr.data, tr.data, NFFT = nfft, noverlap = windlap, Fs = 1./tr.stats.delta, scale_by_freq = True)
        freq = freq[1:]
        power = power[1:]
        power = np.absolute(power)
        resppath += tr.id
        resp = evalresp(t_samp = tr.stats.delta, nfft = nfft, filename = resppath, date = tr.stats.starttime, station = tr.stats.station, 
                                channel = tr.stats.channel, locid = tr.stats.location, network = tr.stats.network, units = 'ACC')
        resp = resp[1:]
        powerR = 10.*np.log10(power/np.abs(resp)**2)
		
        period = 1./freq
        print(tr.id + ' ' + str(np.mean(powerR[(period>=10.)&(period<=15.)])))
						   
        pernlnm, nlnm = get_nlnm()
        pernhnm, nhnm = get_nhnm()
        
        # Plotting data on either the low or high wind subplot
        if idx2 == 0:
            plt.subplot(211)
            plt.semilogx(1./freq, powerR, label = tr.id, color = colors[idx], linewidth = 0.7)     
        else:
            plt.subplot(212)
            plt.semilogx(1./freq, powerR, label = tr.id, color = colors[idx], linewidth = 0.7)  

# Adjusting appearance of low wind subplot and adding labels
plt.subplots_adjust(hspace = 0)
ax = plt.subplot(211)
plt.semilogx(pernlnm, nlnm, 'k', label = 'NLNM/NHNM')
plt.semilogx(pernhnm, nhnm, 'k')
plt.xlim(8.,20.)
plt.ylim(-165.,-140.)
plt.yticks([-160.,-155.,-150.,-145.])
plt.tick_params(axis = 'y', labelsize = 17)
plt.tick_params(axis = 'x', which = 'both', bottom = False, labelbottom = False)
plt.text(8.15, -133., 'Low Wind', fontsize = 22, color = 'black')
ax.text(-0.03, 1., '(a)', transform=ax.transAxes,
      fontsize=26, fontweight='bold', va='top', ha='right')

# Adjusting appearance of high wind subplot and adding labels
ax2 = plt.subplot(212)
plt.semilogx(pernlnm, nlnm, 'k', label = 'NLNM/NHNM')
plt.semilogx(pernhnm, nhnm, 'k')
plt.xlim(8.,20.)
plt.ylim(-165.,-140.)
plt.yticks([-160.,-155.,-150.,-145.])
plt.tick_params(axis = 'x', which = 'both', bottom = True, labelbottom = True, labelsize = 17)
plt.tick_params(axis = 'y', labelsize = 17)
plt.xticks([10.,15.])
plt.xlabel('Period (s)', fontsize = 20)
plt.text(8.15, -133., 'High Wind', fontsize = 22, color = 'black')
ax2.text(-0.03, 1., '(b)', transform=ax2.transAxes,
      fontsize=26, fontweight='bold', va='top', ha='right')

# Creating a y-axis label to extend across both subplots
fig.text(0.04, 0.5, 'Power (dB) (rel. 1 $(m/s/s)^2$ /Hz)', ha = 'center', va = 'center', rotation = 'vertical', fontsize = 20)

# Creating legend of all traces and plotting centered below the plots
handles, labels = ax.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()
handles = handles[:]
labels = labels[:]
ax = plt.gca()
leg = fig.legend(handles, labels, loc = 'lower center', ncol = 3, fontsize = 15)
for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)
plt.subplots_adjust(bottom = 0.14)

# Saving or showing figure
plt.savefig('8-20s_Incoherence_PSD.jpg', format="JPEG", dpi=400)
plt.show()
