#!/usr/bin/env python

from obspy.core import read, UTCDateTime, Stream
from matplotlib import mlab, rcParams
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from obspy.signal.invsim import evalresp
import matplotlib as mpl

# Importing and applying font
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)

debug = False
nfft = 2**12

# Weather station
sta = 'FBA1'
net = 'XX'

# Seismic station
sta2 = 'ASL9'
net2 = 'GS'
loc2 = '00'
chan2 = 'LHZ'

# Start and end time of data for figure
stime = UTCDateTime('2018-176T07:00:00')
etime = UTCDateTime('2018-182T00:00:00')

# Set current time to be the start time and define st as a Stream
ctime = stime
st = Stream()
stws = Stream()

if debug:
	print(etime)

# Read in data for the selected time window
while ctime <= etime:
    path = '/tr1/telemetry_days/' + net2 + '_' + sta2 + '/' + str(ctime.year) + '/' + str(ctime.year) + '_' + str(ctime.julday).zfill(3) + '/' + loc2 + '_' + chan2 + '*'
    pathws = '/tr1/telemetry_days/' + net + '_' + sta + '/' + str(ctime.year) + '/' + str(ctime.year) + '_' + str(ctime.julday).zfill(3) + '/50_LWS*'
    st += read(path)
    stws += read(pathws)
    ctime += 24.*60.*60.

if debug:
	print(pathws)
	print(st)

# Merge traces and if there's a data gap, make it zero. Then trim stream between start and end time	
st.merge(fill_value = 0.)
stws.merge(fill_value = 0.)
st.trim(stime,etime)
stws.trim(stime,etime)

# Removing wind speed response and converting trace to an array
for tr in stws:
    tr.data = tr.data.astype(np.float64)
    tr.data *= 0.1

# Creating spectrogram
specgram, freq, time = mlab.specgram(st[0].data, Fs = st[0].stats.sampling_rate, NFFT = nfft,
                                         pad_to = nfft * 2, noverlap = int(0.25 * nfft), scale_by_freq = True, mode='psd')	
specgram = specgram[1:,:]
freq = freq[1:]

# Reading in responses - responses for the XX network isn't yet with IRIS
if net2 == 'XX':
	resppath = ''

else:
	resppath = '/APPS/metadata/RESPS/'

# Removing response from seismic data traces	
resp = evalresp(t_samp = st[0].stats.delta, nfft = 2 * nfft, filename = resppath + 'RESP.' + net2 + '.' + \
                            sta2 + '.' + loc2 + '.' + chan2,  date = st[0].stats.starttime, station = sta2,
                            channel = chan2, network = net2, locid = loc2, units = 'ACC') 
resp = resp[1:]
specgram = (specgram.T/(np.abs(resp*np.conjugate(resp)))).T
specgram = 10. * np.log10(specgram)

# Decimating wind speed to have a data point every minute
tr.decimate(5)
tr.decimate(2)
tr.decimate(6)

# Establishing the time (independent) variable
time *= 1. / (60. * 60.)
timews = np.asarray(range(tr.stats.npts))/60.

# Creating the figure and labeling it, etc.
fig = plt.figure(1, figsize=(16,8))
ax1 = fig.add_subplot(111)
plt.xlabel('Time (hr)', fontsize = 19)
plt.xticks()
plt.tick_params(labelsize = 15)
att1 = ax1.pcolormesh(time, 1./freq, specgram,  vmin = -180., vmax = -140., cmap = 'bone')
ax1.set_yscale('log')
ax1.set_ylim(10., 1000.)

# Plotting the wind speed on a second y axis on the same figure
ax2 = ax1.twinx()
ax2.plot(timews, tr.data, color = 'aqua', linewidth = 0.7, alpha = 0.7, label = 'Wind Speed')
leg = plt.legend(fontsize = 18)
for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)
ax2.set_ylim((0., 1.2 * max(tr.data)))
plt.xlim(min(time),max(time)) 
#plt.title('Vertical Seismic Noise and Wind Speed', fontsize = 24, y = 1.02)
ax1.set_ylabel('Seismic Period (s)', fontsize = 19)
ax2.set_ylabel('Wind Speed (m/s)', fontsize = 19, rotation = 270, va = 'bottom')
ax2.tick_params(labelsize = 15)

# Creating power colorbar
cb = plt.colorbar(att1, orientation = 'horizontal')
cb.set_label('Seismic Power (dB rel. 1 ($(m/s^2)^2/Hz$))', fontsize = 17)
cb.ax.tick_params(labelsize = 15)

# Saving or showing figure
plt.savefig(sta + 'windspec_' + sta2 + '_' + loc2 + '_' + chan2 + '.jpg', format="JPEG", dpi=400)
# plt.show()
