#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from obspy.core import UTCDateTime

debug = True

#Set font parameters from matplotlib. 
import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)


sta = "TUC"
loc = "10"
year = 2017
BHH =[]
BHZ=[]
ws=[]
times = []

for day in range(1,366):
    hour = 0.
    try:
    #if True:
        f = open("DATA/" + sta +"_" + loc + "_" + str(year) + "_" + str(day).zfill(3) + ".csv")
        for line in f:
            BHH.append(float(line.split(",")[0]))
            BHZ.append(float(line.split(",")[1]))
            ws.append(float(line.split(",")[2])/10.)
            times.append(float(day) + hour/24.)
            hour += 1.
        f.close()
    except:
        BHH += [0]*24
        BHZ += [0]*24
        ws  += [0]*24
        for hour in range(0,24):
            times.append(day + float(hour)/24.)

times = np.asarray(times)

dows = []
for time in times:
    print(time)
    tstr = str(year) + '_' + str(int(round(time,0))).zfill(3) + "T" + str(int(24*(time % 1))).zfill(2) +":00:00"
    print(tstr)
    dows.append(float(UTCDateTime(tstr).weekday) + 24.*(time % 1)/24.)  

if debug:
    print(len(ws))

BHH = np.asarray(BHH)
BHZ = np.asarray(BHZ)
BHH = 20.*np.log10(BHH)
BHZ = 20.*np.log10(BHZ)
dows = np.asarray(dows)
#BHH[BHH >= 2.*np.std(BHH)] = 2.*np.std(BHH)    
#BHZ[BHZ >= 2.*np.std(BHZ)] = 2.*np.std(BHZ) 
fig = plt.figure(1, figsize=(16,16))
ax1 = plt.subplot(3,1,1)
#plt.title(sta + ' ' + loc + ' ' + str(year) + ' 5 Hz to 15 Hz')
ax1.plot(times, BHH, '.', label="BHH")
ax1.plot(times, BHZ, '.', label="BHZ")
plt.ylim((-150,-60))
plt.xlabel('Day of Year')
ax1.set_ylabel('RMS Power (dB)')
ax2= ax1.twinx()
ax2.plot(times, ws,label="WS", color="C2", alpha=.4)
ax2.set_ylabel("Wind Speed (m/s)", rotation=270, labelpad = 30.  )
plt.xlabel('Day of Year')
#plt.legend()
#plt.xlim((min(times), max(times)))
plt.xlim((1., max(times)))
plt.subplot(3,1,2)
plt.plot(ws, BHH, '.', label="BHH")
plt.plot(ws, BHZ, '.', label="BHZ")
#plt.legend()
plt.xlim((min(ws),max(ws)))
plt.ylim((-150,-60))
plt.ylabel('RMS Power (dB)')
plt.xlabel('Wind Speed (m/s)')
ax3 = plt.subplot(3,1,3)
for x in range(7):
    if x == 0:
        ax3.axvspan(x + 14./24.,x + 23./24.,alpha=.3, color='.5', label='Local 9 to 5 Hours')
    else:
        ax3.axvspan(x + 14./24.,x + 23./24.,alpha=.3, color='.5')
plt.plot(dows, BHH, '.')
plt.plot(dows, BHZ, '.')
for num in range(7):
    plt.plot([num,num], [-180., -40.],'k')
plt.ylim((-150,-60))
plt.xticks(np.asarray(range(7)) + 0.5 , ["Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"])
plt.ylim((-150., -60.))
plt.xlim((0.,7.))
plt.xlabel('Day of the Week (UTC)')
plt.ylabel('RMS Power (dB)')

#plt.subplot(2,1,1)
#plt.plot(times, BHH,'.', label="BHH")
#plt.plot(times, BHZ, '.', label="BHZ")
#plt.subplot(2,1,2)
#plt.plot(times, ws, '.', label='ws')
#plt.legend()


handles, labels = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()
handles3, labels3 = ax3.get_legend_handles_labels()
handles += handles2 + handles3
labels += labels2 + labels3
ax = plt.gca()
fig.legend(handles, labels, loc = 'lower center', ncol = 4, fontsize = 17)

xs =[.05, .05, .05]
ys=[.89 ,.61, .34 ]
letters =['a','b', 'c']

for triple in zip(xs, ys, letters):
    plt.text(triple[0], triple[1], '(' + triple[2] + ')', fontsize=22, transform=plt.gcf().transFigure)

plt.savefig('High_frequency_' + str(sta) + '_' + loc + '_' + str(year) + '.jpg', format="JPEG", dpi = 400) 




plt.show()
