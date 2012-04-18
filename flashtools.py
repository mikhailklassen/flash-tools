# Some Python functions useful in working with data from FLASH

def read_sinks_evol(sinks_evol_file):
    '''Read in a sinks_evol.dat file from a FLASH2.5 simulation and return a Python dictionary object where the keys are the column headers.'''
    from numpy import genfromtxt
    f = open(sinks_evol_file)
    headers = f.readline()
    headers = headers.split()
    ncols = len(headers)
    data = genfromtxt(sinks_evol_file)

    a = {}

    for i in range(ncols):
        k = headers[i]
        a[k] = data[1:,i]

    return a

def hosokawa_radius(mass,accretionrate):
    '''Calculates the best fit protostellar radius by interpolating to the Hosokawa & Omukai (2008) protostellar tracks
    using the current mass and accretion rate, expressed in solar masses and solar masses per year, respectively.
    Returns the protostellar radius in units of solar radius.'''
    from numpy import log10

    mass_log = log10(mass)
    accrrate_log = log10(accretionrate)

    msline = 10.0**((log10(10.0) - log10(2.0)) / (log10(70.0) - log10(2.0)) * (mass_log - log10(2.0)) + log10(2.0))

    anfM_log = log10(0.1)
    anfR_log_int = (1.05 - 0.15) / (-3.0 + 6.0) * (accrrate_log + 6.0) + 0.15
    anfbumpM_log_int = (0.85 - 0.25) / 3.0 * (accrrate_log + 6.0) + 0.25
    anfbumpR_log_int = (1.45 - 0.25) / 3.0 * (accrrate_log + 6.0) + 0.25
    maxbumpM_log_int = (1.05 - 0.30) / 3.0 * (accrrate_log + 6.0) + 0.30
    maxbumpR_log_int = (2.00 - 0.50) / 3.0 * (accrrate_log + 6.0) + 0.50
    endbumpM_log_int = (1.3 - 0.7) / 3.0 * (accrrate_log + 6.0) + 0.70
    endbumpR_log_int = (0.75 - 0.15) / 3.0 * (accrrate_log + 6.0) + 0.15
    endM_log = log10(99.0)
    endR_log = log10(10.2)

    anftobump_log_slope = (anfbumpR_log_int - anfR_log_int) / (anfbumpM_log_int - anfM_log)
    bumpup_log_slope = (maxbumpR_log_int - anfbumpR_log_int) / (maxbumpM_log_int - anfbumpM_log_int)
    bumpdown_log_slope = (endbumpR_log_int - maxbumpR_log_int) / (endbumpM_log_int - maxbumpM_log_int)
    end_log_slope = (endR_log - endbumpR_log_int) / (endM_log - endbumpM_log_int)

    if mass_log <= anfbumpM_log_int:
        R = 10.0**(anftobump_log_slope * (mass_log - anfM_log) + anfR_log_int)
    elif mass_log <= maxbumpM_log_int:
        R = 10.0**(bumpup_log_slope * (mass_log - anfbumpM_log_int) + anfbumpR_log_int)
    elif mass_log <= endbumpM_log_int:
        R = 10.0**(bumpdown_log_slope * (mass_log - maxbumpM_log_int) + maxbumpR_log_int)
    elif mass_log <= endM_log:
        R = 10.0**(end_log_slope * (mass_log - endbumpM_log_int) + endbumpR_log_int)

    radius = max(R, msline)

    return radius

