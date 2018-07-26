#!/usr/bin/env python

from obspy.core import read, UTCDateTime, Stream
import matplotlib.pyplot as plt
from scipy import signal
import matplotlib as mpl
from itertools import combinations
import numpy as np
debug = True

# Importing and applying font
mpl.rc('font', family = 'serif')
mpl.rc('font', serif = 'Times') 
mpl.rc('text', usetex = True)

# Creating figure
fig = plt.figure(1, figsize=(18,14))
#plt.suptitle('Coherences Between Seismic and Weather Sensors', y = 0.96, fontsize = 28)

nfft = 2**9

# Listing seismic stations to be used
netstas = ['GS_ASL8_00', 'GS_ASL9_00', 'GS_ASL9_10', 'XX_ENG5_60']

# Establishing start and end times for the low wind and high wind period chosen for the figure
stimes = [UTCDateTime('2018-177T02:24:00'), UTCDateTime('2018-176T09:27:00')]
etimes = [UTCDateTime('2018-177T08:37:00'), UTCDateTime('2018-176T15:34:00')]

# Plotting low wind data and high wind data in different columns of subplots
for idx2, sandetime in enumerate(zip(stimes, etimes)): 
    stime = sandetime[0]
    etime = sandetime[1]
    
    #### Seismic data ####
    
    # Plotting each channel of seismic data on different subplots    
    for idx, chan in enumerate(['LH1', 'LH2', 'LHZ']):
        st = Stream()
        
        # Reading in data for each station
        for netsta in netstas:
            net, sta, loc = netsta.split('_')
            ctime = stime
            while ctime <= etime:
                path = '/tr1/telemetry_days/' + net + '_' + sta + '/' + str(ctime.year) + '/' + str(ctime.year) + '_' + str(ctime.julday).zfill(3) + '/' + loc + '_LH*'
                st += read(path)
                ctime += 24.*60.*60.
        
        # Merge traces and if there's a data gap, make it zero. Then trim stream between start and end time	
        st.merge(fill_value = 0)
        st.sort()
        st.trim(stime, etime)
        
        # Correcting seismometer orientation from North for each instrument   
        for netsta in netstas:
            net, sta, loc = netsta.split("_")
            if netsta == "GS_ASL9_00":
                angle = -62.44
            elif netsta == "GS_ASL9_10":
                angle= -0.4263
            elif netsta == "XX_ENG5_60":
                angle = 79.56
            else:
                angle = 0.
            if angle != 0.:
                cosd = np.cos(np.deg2rad(angle))
                sind = np.sin(np.deg2rad(angle))
                dataLH1 = cosd * st.select(station = sta, location = loc, channel = 'LH1')[0].data + \
                    - sind * st.select(station = sta, location = loc, channel = 'LH2')[0].data
                
                dataLH2 = sind * st.select(station = sta, location = loc, channel = 'LH1')[0].data + \
                        cosd * st.select(station = sta, location = loc, channel = 'LH2')[0].data
                st.select(network=net, station=sta, location=loc, channel= "LH1")[0].data = dataLH1
                st.select(network=net, station=sta, location=loc, channel= "LH2")[0].data = dataLH2
        st = st.select(channel=chan)
        
        # Plotting in either the low wind or high wind column        
        if idx2 == 0:
            plt.subplot(4, 2, 2*(idx) + 1) 
        else:
            plt.subplot(4, 2, 2*(idx) + 2)
             
        # Calculating coherence between each possible combination of two stations and plotting on subplots
        for comb in combinations(range(len(netstas)),2):
            f, Cxy = signal.coherence(st[comb[0]].data, st[comb[1]].data, fs = st[0].stats.delta, nperseg = nfft, noverlap = int(0.5*nfft))
            per = 1./f
            if (idx2 == 0) and (chan == "LH1"):
                ax = plt.subplot(4,2,1)
                lab = st[comb[0]].stats.station + ' ' + st[comb[0]].stats.location + ' and '
                lab += st[comb[1]].stats.station + ' ' + st[comb[1]].stats.location 
                ax.semilogx(per, Cxy, label = lab)
            else:
                plt.semilogx(per, Cxy)
        
        # Adjusting appearance of plots and labels
        plt.xlim((4.,500.))
        plt.ylim((0.0,1.0))
        plt.yticks([0.2, 0.4, 0.6, 0.8])
        plt.tick_params(axis = 'y', labelsize = 17)
        plt.tick_params(axis = 'x', which = 'both', bottom = False, labelbottom = False)
        
        if chan == 'LH1':
            chanlabel = 'LHN'
        elif chan == 'LH2':
            chanlabel = 'LHE'
        else:
            chanlabel = chan
      
        plt.text(4.5, 0.5, chanlabel, fontsize = 22)
        
        plt.subplots_adjust(hspace = 0)

        plt.subplot(4, 2, 1)
        #plt.title('Calm', y = 1.07, fontsize = 20)
        plt.subplot(4, 2, 2)
        #plt.title('Windy', y = 1.07, fontsize = 20)
        
    ##### Weather data #####
    
    # Listing weather stations used
    wnetstas = ['IU_ANMO_50', 'XX_FBA1_50']
    
    # Plotting each channel of weather station information on the same subplot
    for idx, chan in enumerate(['LWS', 'LDO', 'LWD']):
        wst = Stream()
        
        # Reading in weather data for each station
        for wnetsta in wnetstas:
            wnet, wsta, wloc = wnetsta.split('_')
            ctime = stime
            while ctime <= etime:
                wpath = '/tr1/telemetry_days/' + wnet + '_' + wsta + '/' + str(ctime.year) + '/' + str(ctime.year) + '_' + str(ctime.julday).zfill(3) + '/' + wloc + '_' + chan + '*'
                wst += read(wpath)
                ctime += 24.*60.*60.
        
        # Merge traces and if there's a data gap, make it zero. Then trim stream between start and end time
        wst.merge(fill_value = 0)
        wst.sort()
        wst.trim(stime, etime)
		
		# Plotting in either the low wind or high wind column
        if idx2 == 0:
            ax2 = plt.subplot(4, 2, 7)
        else:
            plt.subplot(4, 2, 8)
        
        # Calculating coherence between each possible combination of two stations and plotting on subplots                        
        for comb in combinations(range(len(wnetstas)),2):
            wf, wCxy = signal.coherence(wst[comb[0]].data, wst[comb[1]].data, fs = wst[0].stats.delta, nperseg = nfft, noverlap = int(0.5*nfft))
            wper = 1./wf
            if idx2 == 0 and chan != 'LDO':
                lab = wst[comb[0]].stats.station + ' ' + chan + ' and ' + wst[comb[1]].stats.station + ' ' + chan
                ax2.semilogx(wper, wCxy, label = lab, color = str(float(3*idx)/4.))
            if chan == 'LDO':
                lab = wst[comb[0]].stats.station + ' ' + chan + ' and ' + wst[comb[1]].stats.station + ' ' + chan
                plt.semilogx(wper, wCxy, ':', color = 'k', label = lab)
            else:
                lab = wst[comb[0]].stats.station + ' ' + chan + ' and ' + wst[comb[1]].stats.station + ' ' + chan
                plt.semilogx(wper, wCxy, color = str(float(3*idx)/4.), label=lab)
            
        # Adjusting appearances of subplots and labels        
        plt.xlim((4.,500.))
        plt.xlabel('Period (s)', fontsize = 20)
        plt.ylim((0.0,1.0))
        plt.yticks([0.2, 0.4, 0.6, 0.8])
        plt.tick_params(axis = 'y', labelsize = 17)
        plt.tick_params(axis = 'x', labelsize = 17)
        plt.text(4.5, 0.38, 'Pressure - - -', fontsize = 18, color = 'black')
        plt.text(4.5, 0.58, 'Wind Speed', fontsize = 18, color = 'black')
        plt.text(4.5, 0.47, 'Wind Direction', fontsize = 18, color = 'gray')

# Grabbing handles and labels for each legend object
handles, labels = ax.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()
if debug:
    print(handles2)
    print(handles)

# Creating legend of all traces and centering below the subplots
handles = handles + [handles2[0]] + [handles2[2]] + [handles2[3]]
labels = labels + [labels2[0]] + [labels2[2]] + [labels2[3]]
ax = plt.gca()
fig.legend(handles, labels, loc = 'lower center', ncol = 3, fontsize = 17)

# Creating y-axis labels and centering on the left of each column of subplots
fig.text(0.08, 0.5, 'Coherence ($\gamma^2$)', ha = 'center', va = 'center', rotation = 'vertical', fontsize = 20)
fig.text(0.5, 0.5, 'Coherence ($\gamma^2$)', ha = 'center', va = 'center', rotation = 'vertical', fontsize = 20)

# Making more space for the legend at the bottom
plt.subplots_adjust(bottom = 0.17)

# Adding subplot labels
xs =[.09, .09, .09, .09, .51, .51, .51, .51]
ys=[.87, .69, .52, .34, .87, .69, .52, .34]
letters =['a','b', 'c', 'd', 'e', 'f', 'g', 'h']
for triple in zip(xs, ys, letters):
    plt.text(triple[0], triple[1], '(' + triple[2] + ')', fontsize=26, transform=plt.gcf().transFigure)

# Saving or showing figure
plt.savefig('Coherences_wdir_176-177_6hr.jpg', format="JPEG", dpi=400)
plt.show()
