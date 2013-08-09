#!/usr/bin/env python

# Import required modules
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import sys
import os
import numpy as np
import matplotlib.pyplot as plt


sigma = 5.6705119000000E-05
Msun = 1.9889225E33
Lsun = 3.839E33
Rsun = 69550000000.0

# Function definitions
def strictly_increasing(x):
    dx = np.diff(x)
    return np.all(dx > 0)

if len(sys.argv) == 1:
    sys.exit("Must specify an input file. Call 'clean_sinks_evol.py -h' for help")

# Open datafile
if sys.argv[1] == '-h':
    # Print helpful explanation
    print '\n CLEAN_SINKS_EVOL.PY'
    print ' ==================='
    print '\n Due to restarts in the simulation, SINKS_EVOL.DAT data may contain jumps to earlier time'
    print ' t1 ... t2 ... t3 ... [CHECKPOINT: t4] ... t5 ... t6 ... <SIMULATION STOPPED>'
    print '                      [< RESTARTING >] ... t5 ... t6 ... t7 ...'
    print '\n SINKS_EVOL.DAT now contains duplicate t5 and t6 entries.'
    print ' This script will scrub out the earlier duplicates.'

    print '\n CALLING SEQUENCE:'
    print ' clean_sinks_evol.py sinks_evol.dat'
    print '\n You specify the sinks_evol.dat file as the argument to the script.'

    sys.exit()

# Determine SINKS_EVOL.dat from first command line argument
datafile = sys.argv[1]
relpath = os.path.dirname(datafile)
if len(relpath)==0:
    relpath='.'

creation_info_file = relpath+'/'+'sinks_creation_info.dat'

# Error checking: Do sinks_evol.dat and sinks_creation_info.dat exist?
try:
    open(datafile)
except IOError as e:
    print 'Could not find datafile: {0}'.format(datafile)

try:
    open(creation_info_file)
except IOError as e:
    print 'Could not find associated sinks_creation_info.dat. Make sure it resides in the same dir as sinks_evol.dat'

# Read data files
data = np.genfromtxt(datafile)
print 'SINKS_EVOL.DAT file read.'
creation_info = np.genfromtxt(creation_info_file)
taglist = creation_info[1:,0].astype(int)
print 'SINKS_CREATION_INFO.DAT file read.'

tags = data[1:,0].astype(int)
time = data[1:,1]

# Check if there are other sink particles that were present from the start.
# These would not have been listed in the file sinks_creation_info.dat
extras = np.nonzero(~np.setmember1d(tags,taglist))[0] # Finds unique indices where 
                                                      # entries in 'tags' do not 
                                                      # match entries in 'taglist'
for i in extras:
    taglist = np.concatenate(([tags[i]],taglist))

nparticles = len(taglist)
print '\nThere are {0} sink particles in our simulation'.format(nparticles)
print 'Sink particle tags:'
print taglist

# Prepare sinks_evol_fixed.dat for writing
f = open(relpath+'/'+'sinks_evol_fixed.dat','w')
g = open(datafile,'r')
header = g.readline()
hwords = str.split(header)
g.close()
if 'L_int' in hwords:
    f.write(header)
else:
    header = header[:-1]+'            L_int'+'      ZAMS Radius\n'
    f.write(header)
    zams = np.genfromtxt(relpath+'/'+'ZAMS.dat')
    massarr = zams[:,0]            # In solar masses
    lumarr = 10**zams[:,1] * Lsun  # In erg/s
    temparr = 10**zams[:,2]        # In Kelvin


def calc_vals(mass):
    sc = -1 
    for j in range(len(massarr)):
        if mass/Msun>=massarr[j]:
            sc = j 
    if sc>=0:
        lint = lumarr[sc]
        T = temparr[sc]
        radius = (lint / (4 * np.pi * sigma * T**4 ))**(0.5)
    else:
        lint = 0.0
        T = 3000.0
        radius = 0.0
    return lint,T,radius



# Perform loop over particles
for i in range(nparticles):
    print '\nReading particle {0:2d}, Tag: {1:10d}...'.format(i,taglist[i])
    indices = np.where(tags == taglist[i])[0]
    timearr = time[indices]

    print 'Checking monotonicity of time array'

    if not strictly_increasing(timearr):
        print 'Time array not strictly increasing'
        # Find the shortest timestep in our simulation
        dts = timearr[:-1]-timearr[1:]   # Array of timesteps
        dtmin = min(abs(dts))            # Smallest timestep

        # Use lexsort to sort timearr, but in the order they occur (idx)
        idx = range(len(timearr))
        idlex = np.lexsort((idx,timearr))
        
        # Create tsorted array, sorted according to lexsort
        tsorted = timearr[idlex]
        
        # Take difference between tsorted and tsorted shifted by 1
        # Find the indices where the difference is smaller than smallest timestep/2
        ind = np.where(abs(tsorted[1:]-tsorted[:-1]) < dtmin/2)[0]
        
        # Deleting array entries at those indices results in duplicates being eliminated
        tsifted = np.delete(tsorted,ind)
        
        # Now do this for the all the current sink particles data
        data_fixed = np.delete(data[1:,:][indices,:][idlex],ind,axis=0)

        numlines = np.shape(data_fixed)[0]
        numcols = np.shape(data_fixed)[1]

        for i in range(numlines):
            # First column is tag number
            line = ' %(tag)16i' % {"tag": data_fixed[i,0].astype(int)}
            for j in range(1,numcols):
                if hwords[j]=='stage':
                    # If the current column is the protostellar stage, format as integer
                    line = line + ' %(number)16i' % {"number": data_fixed[i,j].astype(int)}
                else:
                    # Otherwise format as float
                    line = line + ' %(number)16.9E' % {"number": data_fixed[i,j]}
            if 'L_int' not in hwords:
                lint,temp,radius = calc_vals(data_fixed[i,14])
                line = line + ' %(number)16.9E' % {"number": lint}
                line = line + ' %(number)16.9E' % {"number": radius}

            # Add terminal line break
            line = line + '\n'
            f.write(line)

    else:
        # Time array monotonically increasing
        # Write sink particle data to file and go to next sink particle
        print 'Time array monotonically increasing.'
        print 'Writing data to file.'

        data_fixed = data[1:,:][indices,:]
        
        numlines = np.shape(data_fixed)[0]
        numcols = np.shape(data_fixed)[1]

        for i in range(numlines):
            # First column is tag number
            line = ' %(tag)16i' % {"tag": data_fixed[i,0].astype(int)}
            for j in range(1,numcols):
                if hwords[j]=='stage':
                    # If the current column is the protostellar stage, format as integer
                    line = line + ' %(number)16i' % {"number": data_fixed[i,j].astype(int)}
                else:
                    # Otherwise format as float
                    line = line + ' %(number)16.9E' % {"number": data_fixed[i,j]}
            if 'L_int' not in hwords:
                lint,temp,radius = calc_vals(data_fixed[i,14])
                line = line + ' %(number)16.9E' % {"number": lint}
                line = line + ' %(number)16.9E' % {"number": radius}
            # Add terminal line break
            line = line + '\n'
            f.write(line)
        
        print 'Continuing...'


f.close()


