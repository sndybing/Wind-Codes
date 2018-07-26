#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
debug = True

#Import font parameters from matplotlib
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

z = np.linspace(0., 100., 1000) # meters
freqs = [1./10, 1./50.]
# pairs for Young's modulus and Poisson's ratio
# Granite, loose rock
pairs =[[50.,0.2],[2., 0.2]]

def lam(E,v):
    return E*v/((1+v)*(1-2*v))

def mu(E,v):
   return E/(2*(1+v))

fig = plt.figure(1, figsize=(12,12))
fig.text(0.5, 0.06, 'Attenuation (dB)', ha = 'center', va = 'center', fontsize = 20)
plt.subplots_adjust(wspace=0.001)
plt.subplot(1,2,1)
for idx, pair in enumerate(pairs):
    clam = lam(pair[0],pair[1])
    cmu = mu(pair[0],pair[1])
    if debug:
        print(clam)
        print(cmu)
    c = 5.
    for idx2, freq in enumerate(freqs):
        colors2 = ['red', 'black']
        colors = ['xkcd:green', 'xkcd:azure']
        w0= 2.*np.pi*freq
        term= ((clam + cmu)/(clam+2*cmu))
        print('Here is term:' + str(term))
        G= (1 + term*((w0)/c)*z)*np.exp(-w0*z/c)
        if idx == 1:
            plt.plot(np.abs(10*np.log10(G**2)),z,":", label=str(int(1./freq)) + ' s E=' + str(int(pair[0])) + ' GPa', linewidth=6., color = colors[idx2])
        else:
            plt.plot(np.abs(10*np.log10(G**2)),z, label=str(int(1./freq)) + ' s E=' + str(int(pair[0])) + ' GPa', linewidth=2., color = colors2[idx2])

plt.ylim((min(z),max(z)))
ax = plt.gca()
ax.invert_yaxis()
plt.ylabel('Depth (m)')
#plt.xlabel('Attenuation (dB)')
plt.xlim((0., 50.))
plt.xticks([0., 10., 20., 30., 40.])
plt.legend(loc=4)
UzG = (1./(4.*np.pi*mu(pairs[0][0],pairs[0][1])))*((2.*(1-pairs[0][1])/z) + (1/z))
UzS = (1./(4.*np.pi*mu(pairs[1][0],pairs[1][1])))*((2.*(1-pairs[1][1])/z) + (1/z))

print(pairs[0][0])
plt.subplot(1,2,2)
plt.plot((np.absolute(20.*np.log10(UzG))),z, label='Noise E=' + str(int(pairs[0][0])) +  ' GPa', color = 'xkcd:violet', linewidth = 2.)
plt.plot((np.absolute(20.*np.log10(UzS))),z, label='Noise E=' + str(int(pairs[1][0])) + ' GPa', color = 'orange', linewidth = 2.)
#plt.ylim((min(z),max(z)))
plt.ylim((min(z),max(z)))
plt.yticks([0., 20., 40., 60., 80., 100.],[])
ax = plt.gca()
ax.invert_yaxis()
#plt.ylabel('Depth (m)')
#plt.xlabel('Attenuation (dB)')
#plt.xlim((0.1, 50.))

plt.xticks([0., 20., 40., 60., 80.])
plt.legend(loc=3)
xs =[.09, .51]
ys=[.91 ,.91 ]
letters =['a','b']

for triple in zip(xs, ys, letters):
    plt.text(triple[0], triple[1], '(' + triple[2] + ')', fontsize=22, transform=plt.gcf().transFigure)

plt.savefig('Noisemodel_edit.jpg',format='JPEG',dpi=400)
plt.show()
