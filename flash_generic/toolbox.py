try:
    from yt.mods import *
except:
    import sys
    sys.exit("yt not found or not in PYTHONPATH")
import numpy as np
import matplotlib.pyplot as plt

# List of constants in CGS units
pc              = 3.08568025e18       # Parsec in cm
AU              = 1.49598e13          # AU in cm
Msun            = 1.98892e33          # Solar mass in g
Rsun            = 6.955e10            # Radius of the Sun in cm
Lsun            = 3.839e33            # Solar luminosity in erg
secyr           = 31556926.0          # Seconds in a year
G               = 6.6725985e-8        # Gravitational constant in cm^3 g^-1 s^-2
sb              = 5.6705119e-5        # Stefan-Boltzmann constant in erg cm^-2 s^-1 K^-4
hconst          = 6.6260755400000E-27 # Planck's constant in erg Hz^-1
c               = 2.9979245800000E+10 # Speed of light in cm s^-1
kb              = 1.3806581200000E-16 # Boltzmann's constant in erg K^-1
electron_charge = 4.8032068150000E-10 # Charge of the electron in esu
electron_mass   = 9.1093897540000E-28 # Mass of the electron in g
proton_mass     = 1.6726231100000E-24 # Mass of the proton in g
fine_structure  = 7.2973530764000E-03 # Fine structure constant
avogadro        = 6.0221367360000E+23 # Avogadro's constant
gas_constant    = 8.3145119843000E+07 # Ideal gas constant in cm^2 s^-2 K^-1
wien            = 2.8977562400000E-01 # Wien's constant in cm K

def select_scale(length):
    '''
    Chooses an "optimal" scale for plotting purposes, i.e. should the distance axis
    be in units of cm, solar radii, AU, pc, kpc, Mpc, etc.
    '''
    scale = 'cm'
    if length/Rsun > 0.001:
        scale = 'Rsun'
    if length/Rsun > 1000.0:
        scale = 'AU'
    if length/AU > 10000.0:
        scale = 'pc'
    if length/pc > 1000.0:
        scale = 'kpc'
    if length/pc > 1.e6:
        scale = 'Mpc'
    return scale

def get_times(files):
    '''
    Returns an array of times in the base units of the simulation from an input
    array of plot files or checkpoint files.
    '''
    ts = TimeSeriesData.from_filenames(files)
    times = np.zeros(len(ts)) 
    for i, pf in enumerate(ts):
        times[i] = pf.current_time
    return times

def index_containing_substring(the_list, substring):
    for i, s in enumerate(the_list):
        if substring in s:
              return i
    return -1

def find_rigid_rotation_omega(parameters):
    is_rotating = False
    keys = parameters.keys()
    try:
        omega_idx = index_containing_substring(keys,'omega')
        print 'Found rotation parameter {0}.'.format(keys[omega_idx])
        omega = parameters[keys[omega_idx]]
        is_rotating = True
    except:
        omega = 0.0
    return omega, is_rotating

def find_rotation_beta(M,R,omega):
    # Calculate the ratio of rotational energy to gravitational binding
    # energy, crudely approximating the system as a uniform density
    # sphere.
    I = 2.0/5.0 * M * R**2
    Krot = 0.5 * I * omega**2
    U = 3.0/5.0 * G * M**2 / R
    frac = Krot/U
    return frac

def radial_profile(parameters,chkfiles):
    xmin, xmax = float(parameters['xmin']), float(parameters['xmax']) 
    ymin, ymax = float(parameters['ymin']), float(parameters['ymax']) 
    zmin, zmax = float(parameters['zmin']), float(parameters['zmax']) 
    xsize = xmax-xmin
    xradius = xsize/2.0
    ysize = ymax-ymin
    yradius = ysize/2.0
    zsize = zmax-zmin
    zradius = zsize/2.0
    try:
        mu_mol = float(parameters['mu_mol'])
    except:
        mu_mol = 2.14
 
    pf = load(chkfiles[0])
    sphere = pf.h.sphere(pf.domain_center, (1., "pc"))
    rad_profile = BinnedProfile1D(sphere, 300, "Radiuspc", 0.0, 1., log_space=False)
    rad_profile.add_fields("Density")
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.axvline(x=1.595e17/pc,color='#aaaaaa')
    ax1.semilogy(rad_profile["Radiuspc"], rad_profile["Density"])
    ax1.set_ylim(min(rad_profile["Density"]),1.1*max(rad_profile["Density"]))
    ax1.grid()
    ax1.set_title('Density Profile')
    ax1.set_xlabel('Radius (pc)')
    ax1.set_ylabel(r'Mass Density (g cm$^{-3}$)')
    ax2 = ax1.twinx()
    ax2.set_yscale('log')
    y1, y2 = ax1.get_ylim()
    ax2.set_ylim(y1/mu_mol/proton_mass,y2/mu_mol/proton_mass)
    ax2.set_ylabel(r'Number Density (cm$^{-3}$)')
    fig.savefig(parameters['outpath']+'/density_profile.png') 
    # Same plot, but log-log axes 
    plt.clf()
    ax1 = fig.add_subplot(111)
    ax1.axvline(x=1.595e17/pc,color='#aaaaaa')
    ax1.loglog(rad_profile["Radiuspc"], rad_profile["Density"])
    ax1.grid()
    ax1.set_title('Density Profile')
    ax1.set_xlabel('Radius (pc)')
    ax1.set_ylabel(r'Mass Density (g cm$^{-3}$)')
    ax2 = ax1.twinx()
    ax2.set_yscale('log')
    y1, y2 = ax1.get_ylim()
    ax2.set_ylim(y1/mu_mol/proton_mass,y2/mu_mol/proton_mass)
    ax2.set_ylabel(r'Number Density (cm$^{-3}$)')
    fig.savefig(parameters['outpath']+'/density_profile_log.png') 
    
    rho_c = max(rad_profile["Density"])
    return xsize, xradius, rho_c

#def cloud_mass(parameters):
#    '''
#    Takes a dict "parameters" and determines the cloud mass.
#    '''
#    return

#def mean_density(parameters):
#    return rho_mean, n_mean

#def freefall_time(parameters):
#    return
#
#def rms_Mach_number(parameters):
#    return
#
#def sound_cross_time(parameters):
#    return
#
#
#def turbulent_cross_time(parameters):
#    return
#
#
#def Jeans_length(parameters):
#    return
#
#
#def Jeans_volume(parameters):
#    return
#
#
#def Jeans_mass(parameters):
#    return
#
#
#def simulation_box_size(parameters):
#    return
#
#
#def smallest_cell_size(parameters):
#    return
#
#
#def max_gas_density(parameters):
#    return
#
#
#def max_number_density(parameters):
#    return
#
#
#def sink_r_accr(parameters):
#    return
#
#
#def number_sinks(pltfiles):
#    return
#

def ZAMS_radius(mass):
    '''
    Based on Tout et al. 1996
    Calculate the ZAMS radius of a star for a given mass (expressed in solar units).
    Returns the radius in solar units. Accurate for stars of solar metallicity.
    '''
    theta = 1.71535900
    iota = 6.59778800
    kappa = 10.08855000
    llambda = 1.01249500
    mu = 0.07490166
    nu = 0.01077422
    xi = 3.08223400
    omicron = 17.84778000
    ppi = 0.00022582

    mm = mass   # solar masses
    Rms = (theta*mm**(2.5) + iota*mm**(6.5) + kappa*mm**(11) + llambda*mm**(19) + mu*mm**(19.5))
    Rms = Rms / (nu + xi*mm**(2) + omicron*mm**(8.5) + mm**(18.5) + ppi*mm**(19.5))

    return Rms

def ZAMS_luminosity(mass):
    '''
    Based on Tout et al. 1996
    Calculate the ZAMS luminosity of a star for a given mass (expressed in solar units).
    Returns the luminosity in solar units. Accurate for stars of solar metallicity.
    '''
    alpha = 0.39704170
    beta = 8.52762600
    gamm = 0.00025546
    delta = 5.43288900
    epsil = 5.56357900
    zeta = 0.78866060
    eta = 0.00586685
     
    mm = mass   # solar masses
    Lms = (alpha*mm**(5.5)+beta*mm**(11))
    Lms = Lms / (gamm + mm**(3) + delta*mm**(5) + epsil*mm**(7) + zeta*mm**(8) + eta*mm**(9.5))

    return Lms

def planck_phot(fr,T):
    '''
    Planck function, with the integration over solid angle already taken,
    divided by the photon energy. When integrated, we get the total number
    of photons.
    Returns the number of photons emitted at a given frequency, for a
    blackbody at a given temperature.
    '''
    tpic2 = 2*pi/(c**2)
    rfr = hconst/kb/T
    nphot = tpic2 * fr**2 / ( exp(fr*rfr) - 1.0 )

    return nphot

def plot_profiles(parameters):
    xmin, xmax = float(parameters['xmin']), float(parameters['xmax']) 
    ymin, ymax = float(parameters['ymin']), float(parameters['ymax']) 
    zmin, zmax = float(parameters['zmin']), float(parameters['zmax']) 
    xsize = xmax-xmin
    xradius = xsize/2.0
    ysize = ymax-ymin
    yradius = ysize/2.0
    zsize = zmax-zmin
    zradius = zsize/2.0
    try:
        mu_mol = float(parameters['mu_mol'])
    except:
        mu_mol = 2.14
    try:
        profile = int(parameters['density_profile'])
    except:
        sys.exit('density_profile parameter not defined. Check the flash.par parameter file.')
    if profile == 1: # Power-law profile
        print 'Power-law profile'
        dens_inner_radius = float(parameters['dens_inner_radius'])
        dens_outer_radius = float(parameters['dens_outer_radius'])
        M_total = float(parameters['M_total'])
        dens_power_law = float(parameters['dens_power_law'])
        density_contrast = float(parameters['density_contrast'])
 
        fita = dens_power_law * (3.0 - dens_power_law) * M_total / \
               (8.0*np.pi) * dens_outer_radius**(dens_power_law-3.0) / \
               dens_inner_radius**(dens_power_law+2.0)
        rho_c = (3.0 - dens_power_law)*M_total / (4.0*np.pi)* \
                dens_outer_radius**(dens_power_law - 3.0) / \
                dens_inner_radius**(dens_power_law) * \
                (1.0 + dens_power_law / 2.0)
        rho0 = (3.0 - dens_power_law)*M_total/(4.0*np.pi*dens_outer_radius**3.0)
        
        radius = xradius 
        radius = np.linspace(0,radius,300)        
        dens = np.zeros(len(radius))
        for i in range(len(dens)):
            if (radius[i] <= dens_inner_radius):
                dens[i] = -fita*radius[i]**2.0 + rho_c
            elif (radius[i] <= dens_outer_radius) and (radius[i] > dens_inner_radius):
                dens[i] = (3.0 - dens_power_law)*M_total / \
                       (4.0*np.pi*dens_outer_radius**(3.0-dens_power_law)* \
                       radius[i]**dens_power_law)
            else:
                dens[i] = density_contrast * rho0
        # Plot the density profile
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.axvline(x=2.440868032664758e-4,color='#aaaaaa')
        ax1.axvline(x=dens_inner_radius/pc,color='#aaaaaa')
        ax1.axvline(x=dens_outer_radius/pc,color='#aaaaaa')
        ax1.semilogy(radius/pc,dens)
        ax1.set_xlim(radius[0]/pc,1.1*radius[-1]/pc)
        ax1.set_ylim(min(dens),1.1*max(dens))
        ax1.grid()
        ax1.set_title('Density Profile')
        ax1.set_xlabel('Radius (pc)')
        ax1.set_ylabel(r'Mass Density (g cm$^{-3}$)')
        ax2 = ax1.twinx()
        ax2.set_yscale('log')
        y1, y2 = ax1.get_ylim()
        ax2.set_ylim(y1/mu_mol/proton_mass,y2/mu_mol/proton_mass)
        ax2.set_ylabel(r'Number Density (cm$^{-3}$)')
        plt.savefig(parameters['outpath']+'/density_profile.png')
        # Same plot, but log-log axes 
        plt.clf()
        ax1 = fig.add_subplot(111)
        ax1.axvline(x=dens_inner_radius/pc,color='#aaaaaa')
        ax1.axvline(x=dens_outer_radius/pc,color='#aaaaaa')
        ax1.loglog(radius/pc,dens)
        ax1.grid()
        ax1.set_xlim(radius[0]/pc,1.1*radius[-1]/pc)
        ax1.set_ylim(min(dens),1.1*max(dens))
        ax1.set_title('Density Profile')
        ax1.set_xlabel('Radius (pc)')
        ax1.set_ylabel(r'Mass Density (g cm$^{-3}$)')
        ax2 = ax1.twinx()
        ax2.set_yscale('log')
        y1, y2 = ax1.get_ylim()
        ax2.set_ylim(y1/mu_mol/proton_mass,y2/mu_mol/proton_mass)
        ax2.set_ylabel(r'Number Density (cm$^{-3}$)')
        plt.savefig(parameters['outpath']+'/density_profile_log.png')
