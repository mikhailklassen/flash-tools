#!/usr/bin/env python
'''
This routine plots the cooling rate used by FLASH for molecular line cooling, based on the data by
Neufeld, Lepp and Melnick, (1995), ApJS, 100, 132

It is possible to create a plot for different temperatures by modifying the tplots variable.
'''
# Temperature at which to plot
tplots = [10.0,40.0,100.0,1000.0,2000.0]

# Import modules
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

# Use TeX
rc('text', usetex=True)
rc('font', family='serif')

# Read in the cooling table, skipping the header line
try:
   open('cool.dat')
except IOError as e:
   print 'ERROR: Could not locate cool.dat file in current working directory. Exiting.'
   sys.exit()

cool_table = np.genfromtxt('cool.dat',skip_header=1)

# Break up the table into columns
temp = 10**cool_table[1:,0] # in Kelvin
dens = 10**cool_table[1:,1] # in particles/cm^3 (I assume)
cooling_power = 10**cool_table[1:,5]

# Internal monologue
'''
Okay, so now it get's a little annoying because of how the data is structured.

cool.dat has the cooling power for different densities and temperature. For a given
density and temperature, the gas cools at a particular rate.

So to proceed, we can either plot the cooling rate as a function of density for a given temperature,
                          OR      the cooling rate as a function of temperature for a given density.

I'm going to hold temperature constant and plot cooling power versus density.
'''
# Find the unique values of temperature
templist = np.unique(temp)

# Define a function that gets the nearest value in an array
def find_nearest(array,value):
    idx=(np.abs(array-value)).argmin()
    return array[idx]

# Create plot
ax = plt.subplot(1,1,1)
plt.hold(True)

for i in range(len(tplots)):
    # Find the unique temperature value in the table nearest to tplots[i]
    nearest_temp = find_nearest(templist,tplots[i])

    # Find the indices at which the temperature is equal to nearest_temp
    idx = np.where(temp == nearest_temp)[0]

    # Plot the cooling power as a function of density, holding T at nearest_temp
    plt.loglog(dens[idx],cooling_power[idx])

plt.xlabel(r'n(H$_2$) [cm$^{-3}$]')
plt.ylabel(r'n(H$_2$) $\Lambda$ [erg s$^{-1}$ (H$_2$)$^{-1}$]')
#plt.axis([1e2,1e10,1e-29,1e-23])

# Create an array of strings from tplot for the legend, but first converting to integers
tplots_labels = map(str,map(int,tplots))
for i in range(len(tplots_labels)):
    tplots_labels[i] = tplots_labels[i]+' K'
plt.legend(tplots_labels,loc=2)

# Save file
filename = 'flash_cooling_table_array'
plt.savefig(filename+'.png')
plt.savefig(filename+'.eps')

# Show figure
plt.show()

