#!/usr/bin/python

from numpy import *
from pylab import *
from matplotlib import rc
from matplotlib import rcParams

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern']})

A=23373131.742021047
B=1.27227214352889927E-20
r0=1.5e18
pi = 3.14159265
rho_ext=9.7573962e-23
alphaexp=1.5

pc2cm = 3.08568025e18
Msun = 1.989e33

numpoints=2000

x = logspace(0,19,numpoints)
rho = zeros(numpoints)
rho_e = zeros(numpoints)
M_r = zeros(numpoints)
for i in range(numpoints):
    if x[i] < r0:
        rho[i] = B
    else:
        rho[i] = A*x[i]**(-alphaexp)
    if i > 0:
        dx = x[i]-x[i-1]
        M_r[i]=M_r[i-1]+4*pi*x[i]**2*dx*rho[i]

print rho[0]
plot(x/pc2cm,rho)
xlabel('Radial distance (pc)')
ylabel('Density [g/cm$^3$]')
title("Thomas Peters' density profile")

savefig('peters_density_profile.pdf')
savefig('peters_density_profile.eps')

#plot(x/pc2cm,M_r/Msun)

show()
