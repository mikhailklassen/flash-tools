#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from toolbox import *

sinks_evol_file = 'sinks_evol_fixed.dat'

f = open(sinks_evol_file,'r')
header = f.readline().split()
f.close()

data = np.genfromtxt(sinks_evol_file,skip_header=1)
tags = data[:,header.index('part_tag')]
time = data[:,header.index('time')]
mass = data[:,header.index('mass')]

taglist = np.unique(tags)

plt.figure()
tag_idx = np.where(tags == taglist[0])
time_part,mass_part = time[tag_idx]/secyr,mass[tag_idx]/Msun
plt.plot(time_part,mass_part)
plt.hold(True)
for i in range(1,len(taglist)):
    tag_idx = np.where(tags == taglist[i])
    time_part,mass_part = time[tag_idx]/secyr,mass[tag_idx]/Msun
    plt.plot(time_part,mass_part)
plt.xlabel('Time (yrs)')
plt.ylabel(r'Mass ($M_{\odot}$)')
plt.title('Mass evolution')
plt.savefig('mass_accretion.png')
