#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from obspy.core import UTCDateTime

debug = True

sta = "COLA"
loc = "10"
year = 2017
BHH =[]
BHZ=[]
ws=[]
times =[]

for day in range(1,366):
    hour = 0.
    try:
    #if True:
        f = open("DATA/" + sta +"_" + loc + "_" + str(year) + "_" + str(day).zfill(3) + ".csv")
        for line in f:
            BHH.append(float(line.split(",")[0]))
            BHZ.append(float(line.split(",")[1]))
            ws.append(float(line.split(",")[2])/10.)
            times.append(day + hour/24.)
            hour += 1.
        f.close()
    except:
        BHH += [0]*24
        BHZ += [0]*24
        ws  += [0]*24

times = range(len(ws))
times = np.asarray(times) + 1
times = times.astype(np.float64)
times *= 1./24.        
if debug:
    print(len(ws))
    
times = times[:-1]
BHH = BHH[:-1]
BHZ = BHZ[:-1]
ws = ws[:-1]    
    
dows =[]
for time in times:
    print(time)
    print("Hour:" + str(24*(time % 1)))
    tstr = str(year) + '_' + str(int(round(time,0)+ 1.)).zfill(3) + "T" + str(int(24*(time % 1))).zfill(2) +":00:00"
    print(tstr)
    dows.append(float(UTCDateTime(tstr).weekday) + 24.*(time % 1)/24.)

BHH = np.asarray(BHH)
BHZ = np.asarray(BHZ)
BHH = 20.*np.log10(BHH)
BHZ = 20.*np.log10(BHZ)

print(BHH)
dows = np.asarray(dows)



#BHH[BHH >= 2.*np.std(BHH)] = 2.*np.std(BHH)    
#BHZ[BHZ >= 2.*np.std(BHZ)] = 2.*np.std(BHZ) 
fig = plt.figure(1)
plt.plot(dows, BHH,'.')
plt.plot(dows,BHZ,'.')
plt.show()







#ax1 = plt.subplot(2,1,1)
#plt.title(sta + ' ' + loc + ' 2017 5 Hz to 15 Hz')
#ax1.plot(times, BHH, '.', label="BHH")
#ax1.plot(times, BHZ, '.', label="BHZ")
#ax1.set_ylabel('RMS Power (dB)')
#ax2= ax1.twinx()
#ax2.plot(times, ws,label="WS", color="C3", alpha=.3)
#ax2.set_ylabel("Wind Speed (m/s)")
#plt.xlabel('Day of Year')
#plt.legend()
#plt.xlim((min(times), max(times)))
#plt.subplot(2,1,2)
#plt.plot(ws, BHH, '.', label="BHH")
#plt.plot(ws, BHZ, '.', label="BHZ")
#plt.xlim((min(ws),max(ws)))
#plt.ylim((-150,-60))
#plt.ylabel('RMS Power (dB)')
#plt.xlabel('Wind Speed (m/s)')
##plt.subplot(2,1,1)
##plt.plot(times, BHH,'.', label="BHH")
##plt.plot(times, BHZ, '.', label="BHZ")
##plt.subplot(2,1,2)
##plt.plot(times, ws, '.', label='ws')
#plt.legend()
#plt.show()


