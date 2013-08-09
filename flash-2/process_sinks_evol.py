import sys
import numpy as np
from toolbox import *

def strictly_increasing(x):
    dx = np.diff(x)
    return np.all(dx > 0)

def calc_vals(mass):
    Lint = ZAMS_luminosity(mass)
    radius = ZAMS_radius(mass)
    T = (Lint / (4.0 * np.pi * radius**2 *sb))**(0.25)
    return lint,T,radius

def process_sinks_evol(datafile,outdir):
    try:
        open(datafile)
    except IOError as e:
        print 'Could not find datafile: {0}'.format(datafile)
        sys.exit('Could not open sinks_evol.dat.')

    # Prepare sinks_evol_fixed.dat for writing
    f = open(outdir+'/'+'sinks_evol_fixed.dat','w')
    g = open(datafile,'r')
    h = open(outdir+'/'+'sinks_evol.dat','w')
    header = g.readline()
    hwords = str.split(header)
    h.write(header)
    line = g.readline()
    # Copy the data file to output directory, removing surplus headers
    while line != '':
        if 'part_tag' not in line:
            h.write(line)
        line = g.readline()
    h.close()
    g.close()

    # Read data files
    data = np.genfromtxt(outdir+'/sinks_evol.dat',skip_header=1)
    print 'SINKS_EVOL.DAT file read.'

    tags = data[:,0].astype(int)
    time = data[:,1]
    taglist = np.unique(tags)
    nparticles = len(taglist)
    print '\nThere are {0} sink particles in our simulation'.format(nparticles)
    print 'Sink particle tags:'
    print taglist

    if 'L_int' in hwords:
        f.write(header)
    else:
        header = header[:-1]+'            L_int'+'      ZAMS Radius\n'
        f.write(header)

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
            data_fixed = np.delete(data[:,:][indices,:][idlex],ind,axis=0)

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

            data_fixed = data[:,:][indices,:]
            
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
    final_masses = []
    data = np.genfromtxt(outdir+'/sinks_evol_fixed.dat', skip_header=1)
    tags = data[:,0]
    taglist = np.unique(tags)
    massarr = data[:,hwords.index('mass')]
    for i in range(nparticles):
        indices = np.where(tags == taglist[i])
        final_masses.append(massarr[indices][-1])
    final_masses = np.array(final_masses)
    final_masses = final_masses[np.where(final_masses >= 0.001*Msun)[0]]
    print final_masses
    sinks_info = {}
    sinks_info['nparticles'] = nparticles
    sinks_info['final_masses'] = final_masses
    sinks_info['mean_mass'] = np.mean(final_masses)
    sinks_info['median_mass'] = np.median(final_masses)
    sinks_info['max_mass'] = max(final_masses)
    sinks_info['min_mass'] = min(final_masses)
    return sinks_info
