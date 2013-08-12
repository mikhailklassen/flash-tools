#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
try:
    from yt.mods import *
    yt_exists = True
except:
    print 'yt not found, or not currently in PYTHONPATH.'
    yt_exists = False
    sys.exit()
import os
import cPickle as pickle
import glob
import numpy as np
from flash_generic.toolbox import *
import subprocess as sub

# Determine directory paths
outdirname = 'diagnostic_output'
srcdir = os.getcwd()
print 'Analysis source directory: ',srcdir
datadir = os.path.dirname(srcdir)
print 'Data directory: ',datadir
outdir = datadir+'/'+outdirname
print 'Diagnostic output directory: ',outdir

# Perform checks
print 'Performing checks.'
if os.path.exists(datadir+'/flash.par'):
    parfile_exists = True
    print 'Found file flash.par.'
else:
    parfile_exists = False
    print 'Parameter file flash.par not found.'

if os.path.exists(datadir+'/flash.dat'):
    datfile_exists = True
    print 'Found file flash.dat.'
else:
    datfile_exists = False
    print 'flash.dat file not found.'

if os.path.exists(datadir+'/sinks_evol.dat'):
    datfile_exists = True
    print 'Found file sink_evol.dat.'
else:
    datfile_exists = False
    print 'sink_evol.dat file not found.'

# Get checkpoint files and plotfiles
print 'Agglomerating checkpoint files.'
chkfiles = glob.glob(datadir+'/*hdf5_chk_*')
chkfiles.sort()
print 'Found {0} checkpoint files.\n'.format(len(chkfiles))

print 'Agglomerating plot files.'
pltfiles = glob.glob(datadir+'/*hdf5_plt_cnt_*')
pltfiles.sort()
print 'Found {0} plot files.\n'.format(len(pltfiles))

# Determine a "first" and "last" chk/plt-file
if len(chkfiles) != 0:
    first_file = chkfiles[0]
    pf = load(first_file)
    last_file = chkfiles[-1]
    pflast = load(last_file)
elif len(pltfiles) != 0:
    first_file = pltfiles[0]
    pf = load(first_file)
else:
    sys.exit('Could not locate any checkpoint files or plotfiles.')
if len(pltfiles) != 0:
    last_file = pltfiles[-1]
    pflast = load(last_file)
print 'Using file {0} to estimate "initial" state of the simulation.'.format(first_file)
print 'Using file {0} to estimate "final" state of the simulation.\n'.format(last_file)

# Analyse parameter file
print 'Extracting simulation parameters from {0}.'.format(first_file)
flashparms = {}
if parfile_exists:
    parfile = open(datadir+'/flash.par','r')
    parlines = parfile.readlines()
    for line in parlines:
        line = line.replace('"', '').strip()
        isParm = False
        if not line: continue
        if len(line.split()) == 3 and line.split()[1] == '=': isParm = True
        if isParm:
            parm, val = line.split()[0], line.split()[2]
            flashparms[parm] = val
parameters = pf.parameters
print 'Found {0} simulation parameters.\n'.format(len(parameters))

# Analyse turbulent velocity field info file
if yt_exists:
    try: 
        vel_binfilename  = pf.parameters["velfile"]
        velfile_exists = True
        vel_infofilename = vel_binfilename.replace('.bin','.info')
        print 'Extracting information about the turbulent velocity field from {0}.'.format(vel_infofilename)
        vel_infofile     = open(datadir+'/'+vel_infofilename,'r')
        vel_info         = vel_infofile.readlines()
        vel_infofile.close()
        for line in vel_info:
            sol_decomp = False
            comp_decomp = False
            turb_type = None
            splitline = line.split(':')
            if splitline[0].strip().find('power spec indx hi') == 0:
                turb_spec_index = splitline[1].strip()
            if splitline[0].strip().find('k-decomp. solenoid') == 1:
                sol_decomp = True # Decomposition into solenoidal modes
            if splitline[0].strip().find('k-decomp. compress') == 1:
                comp_decomp = True # Decomposition into compressive modes

    except:
        velfile_exists = False
        print 'No velocity info file associated with this simulation.\n'

# Create output directory
print 'Creating output directory.'
parameters['outdir'] = outdir
if not os.path.exists(outdir):
    os.mkdir(outdir)
print 'Created {0}.\n'.format(outdir)

# Copy output files
summary_template = srcdir+'/templates/summary_template.tex'
summary_filename = 'summary.tex'
summary_file = outdir+'/'+summary_filename
print 'Copying from '+summary_template+' to '+summary_file
cp = sub.Popen(["cp",summary_template,summary_file],stdout=sub.PIPE,stderr=sub.PIPE)
output, errors = cp.communicate()
print errors
cp.stdout.close()

# Define function that replaces keywords with values
def sed_replace(keyword,value):
    if type(value) == float or type(value) == np.float64:
        value = '{0:12.6g}'.format(value)
    substitution = 's/'+keyword+'/'+str(value)+'/'
    sed = sub.Popen(['sed','-i',substitution,summary_file],stdout=sub.PIPE,stderr=sub.PIPE)
    output,errors = sed.communicate()
    if len(errors) != 0:
        sys.exit(errors)
    sed.stdout.close()
    print '{0}: {1}'.format(keyword,value)

## PHYSICAL SIMULATION PARAMETERS
print 'Processing physical simulation parameters.\n'
R0 = np.mean(pf.domain_width)/2.0 # For now, the radius is just half the box size
cloud_radiusCM = R0
sed_replace('@CLOUDRADCM',cloud_radiusCM)
cloud_radiusAU = R0/AU
sed_replace('@CLOUDRADAU',cloud_radiusAU)
cloud_radiusPC = R0/pc
sed_replace('@CLOUDRADPC',cloud_radiusPC)
total_cloud_mass = pf.h.all_data().quantities["TotalMass"]()*Msun 
sed_replace('@TOTCLOUDMASS',total_cloud_mass/Msun)
mean_mass_density = total_cloud_mass / np.product(pf.domain_width)
sed_replace('@MEANMASSDENS',mean_mass_density )
mean_number_density = mean_mass_density/pf.parameters["mu_mol"]/proton_mass
sed_replace('@MEANNUMDENS',mean_number_density)
mu_mol = pf.parameters["mu_mol"]
sed_replace('@MUMOL',mu_mol )
#temperature = pf.h.all_data().quantities["WeightedAverageQuantity"]("Temperature","CellMassMsun")   # Mass-weighted average temperature
temperature = pf.parameters["isothermal_temp"]
sed_replace('@TEMPERATURE',temperature)
sound_speed = np.sqrt(kb*temperature/mu_mol/proton_mass)
sed_replace('@SOUNDSPEED',sound_speed)
rms_Mach_number = pf.parameters["v_unit"]/sound_speed
sed_replace('@MACHNUMBER',rms_Mach_number)
turbulent_power_index = -int(turb_spec_index)
sed_replace('@TURBPOWER',turbulent_power_index)
#turbulent_field_type = solenoidal/compressive/mix/rand
#sed_replace('@TURBTYPE',turbulent_field_type)
mean_freefall_time = np.sqrt(3.0*np.pi / (32.0 * G * mean_mass_density))
sed_replace('@FREEFALLTIME',mean_freefall_time/secyr)
sound_crossing_time = np.mean(pf.domain_width)/sound_speed
sed_replace('@SOUNDCROSSTIME',sound_crossing_time/secyr)
turbulent_crossing_time = np.mean(pf.domain_width)/sound_speed/rms_Mach_number
sed_replace('@TURBCROSSTIME',turbulent_crossing_time/secyr)
Jeans_length =  np.sqrt( 15.0 * kb * temperature / (4.0 * np.pi * G * mu_mol * proton_mass * mean_mass_density))
sed_replace('@JEANSLENGTH',Jeans_length/AU)
Jeans_volume = 4.0/3.0 * np.pi * Jeans_length**3
sed_replace('@JEANSVOLUME',Jeans_volume/AU**3)
Jeans_mass = mean_mass_density * (Jeans_length)**3
sed_replace('@JEANSMASS',Jeans_mass/Msun)
magnetic_field_flux_density = pf.parameters["magnetic_field_val"]
sed_replace('@MAGFIELD',magnetic_field_flux_density*1.0e6)
Alfven_speed = magnetic_field_flux_density/np.sqrt(mean_mass_density) # magnetic permeability is 1 in CGS units
sed_replace('@ALFVENSPEED',Alfven_speed )
Alfven_time = np.mean(pf.domain_width)/Alfven_speed 
sed_replace('@ALFVENTIME',Alfven_time/secyr )
# Mass-to-flux / Critical mass-to-flux (via Seifried et al. 2011, Mouschovias & Spitzer 1976, ApJ, 210, 326)
mass_to_flux_ratio = total_cloud_mass /(np.pi * cloud_radiusCM**2 * magnetic_field_flux_density) * np.sqrt(G)/0.13
sed_replace('@MASSTOFLUX',mass_to_flux_ratio)
rigid_rotation_omega, is_rotating = find_rigid_rotation_omega(parameters)
sed_replace('@RIGIDOMEGA',rigid_rotation_omega )
rotational_energy_fraction = find_rotation_beta(total_cloud_mass,R0,rigid_rotation_omega)
sed_replace('@ROTENERGYFRAC',rotational_energy_fraction*100 )


## NUMERICAL SIMULATION PARAMETERS
print 'Processing numerical simulation parameters.\n'
simulation_box_sizeCM = np.mean(pf.domain_width)
sed_replace('@BOXSIZECM',simulation_box_sizeCM)
simulation_box_sizeAU = np.mean(pf.domain_width)/AU
sed_replace('@BOXSIZEAU',simulation_box_sizeAU)
simulation_box_sizePC = np.mean(pf.domain_width)/pc
sed_replace('@BOXSIZEPC',simulation_box_sizePC)
simulation_box_volume = np.product(pf.domain_width/pc)
sed_replace('@BOXVOLUME',simulation_box_volume)
smallest_cell_size = pf.h.get_smallest_dx()/AU
sed_replace('@MINGRIDSIZE',smallest_cell_size)
minimum_refinement_level = pf.parameters["lrefine_min"]
sed_replace('@LREFMIN',minimum_refinement_level)
maximum_refinement_level = pf.parameters["lrefine_max"]
sed_replace('@LREFMAX',maximum_refinement_level)
max_gas_density = pf.parameters["max_gas_density"]
sed_replace('@MAXGASDENS',max_gas_density)
max_number_density = max_gas_density/mu_mol/proton_mass
sed_replace('@MAXNUMDENS',max_number_density)
sink_particle_accretion_radius = pf.parameters["r_accretion"]/AU
sed_replace('@SINKRADACCR',sink_particle_accretion_radius)
cfl = pf.parameters["cfl"]
sed_replace('@CFL',cfl)


## SIMULATION OUTCOMES
print 'Processing log file.'
from flash_2.process_log_file import *
logfile = datadir+'/'+pf.parameters["log_file"]
loginfo = process_log_file(logfile)
print 'Done.\n'

print 'Processing sinks_evol.dat'
from flash_2.process_sinks_evol import *
sinks_info = process_sinks_evol(datadir+'/sinks_evol.dat',outdir)
print 'Done.\n'

print 'Processing simulation outcomes.'
machine_info = loginfo["sysinfo"]
machine_info = machine_info.replace('_','\\\_')
print machine_info
sed_replace('@MACHINEINFO',machine_info)
num_chk_files = len(chkfiles)
sed_replace('@NUMCHKFILES',num_chk_files)
num_plt_files = len(pltfiles)
sed_replace('@NUMPLTFILES',num_plt_files)
num_t_steps = loginfo["nsteps"] 
sed_replace('@NUMSTEPS',num_t_steps)
total_wallclock_time = loginfo["wallclock"]
sed_replace('@REALTIME',total_wallclock_time/24.0)
total_cpuhours = loginfo["cpuhours"]
sed_replace('@CPUHOURS',total_cpuhours)
final_simulation_time = pflast.current_time/secyr
sed_replace('@FINALTIME',final_simulation_time)
final_simulation_timestep = pflast["timestep"]/secyr
sed_replace('@DELTATFIN',final_simulation_timestep)
num_sink_particles = sinks_info["nparticles"] 
sed_replace('@NUMSINKS',num_sink_particles)
max_sink_mass = sinks_info['max_mass']
sed_replace('@MAXSINKMASS',max_sink_mass/Msun)
min_sink_mass = sinks_info['min_mass']
sed_replace('@MINSINKMASS',min_sink_mass/Msun)
mean_sink_mass = sinks_info['mean_mass']
sed_replace('@MEANSINKMASS',mean_sink_mass/Msun)
median_sink_mass = sinks_info['median_mass']
sed_replace('@MEDIANSINKMASS',median_sink_mass/Msun)

## Make a plot of the evolution of the density profile
from flash_generic import profiles
profiles.density_profile_1D_evolution(pltfiles[0:2],outdir)

## COMPILE SUMMARY.TEX
print 'Compiling diagnostic summary file.'
os.chdir('../'+outdirname)
pdflatex = sub.Popen(['pdflatex',summary_filename],stdout=sub.PIPE,stderr=sub.PIPE)
output, errors = pdflatex.communicate()
print errors
pdflatex.stdout.close()


# Get the radial profile
#xsize, R0, rho_c = radial_profile(parameters,chkfiles)
