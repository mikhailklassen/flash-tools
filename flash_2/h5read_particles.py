#!/usr/bin/python
'''
This script uses the h5py module to read in the particle data from a single HDF5 file written by FLASH2.5.

It is an example script that can be adapted to do more complex tasks.

Author: Mikhail Klassen
Date:   July 4, 2011
Email:  klassm@mcmaster.ca
'''

import warnings
warnings.filterwarnings('ignore')
import numpy
import h5py
import sys
warnings.filterwarnings('always')

if len(sys.argv) > 2:
    sys.exit("Too many command line arguments!")

if sys.argv[1]=='-h':
    print ' '
    print ' h5read_particles.py'
    print ' ---------'
    print ' Reads in an HDF5 data file specified on the command line following the usage instructions below.'
    print ' Reads in the particle data from the specified file, then prints all the particle properties to the screen.'
    print ' '
    print ' Usage:'
    print '  ./h5read_particles.py datafile.hdf5'
    print ' '
    sys.exit(0)

filename = sys.argv[1]

# Open the specified file in read-only mode
f = h5py.File(filename,'r')

# Extract the particles groups
particles = f.get('particle tracers')

# Number of particles
nparticles = len(particles)

# Number of particle properties
nprops = len(particles[0])

# Loop through particles and list their properties
for j in range(nparticles):
    print 'Particle #',j+1
    print '--------------'
    for i in range(nprops):
        print '{0:25} : {1:15g}'.format(particles.dtype.names[i],particles.value[0][i])
        #print particles.dtype.names[i],' : ',particles.value[0][i]
    print ' '
    print ' '

