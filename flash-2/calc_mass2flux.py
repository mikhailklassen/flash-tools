#!/usr/bin/env python

import numpy as np
import sys

# Newton's gravitational constant
G = 6.6725985e-8                    # Newton's constant in CGS units
Msun = 1.98892e33                   # Mass of the sun in g
pc = 3.08568025e18                  # Parsec in cm
AU = 1.49598e13                     # AU in cm


parfile = 'flash.par'

parfile = open(parfile)

mag_field_val = -999
M_total       = -999
xmin          = -999
xmax          = -999

for line in parfile:
	if 'magnetic_field_val' in line:
		mag_field_val = float(line.split()[-1])
	if 'M_total' in line:
		M_total = float(line.split()[-1])
	if 'xmin' in line:
		xmin = float(line.split()[-1])
	if 'xmax' in line:
		xmax = float(line.split()[-1])

if mag_field_val == -999:
	sys.exit('Parameter file flash.par does not define any magnetic field strength.')
else:
	print '\nMagnetic field strength: {0:10.6g} Gauss'.format(mag_field_val)

if M_total == -999:
	print '\nUnable to determine total system mass.'
	try:
		M_total=float(raw_input('Please enter total system mass (in g): '))
	except ValueError:
		print "Invalid entry. Aborting."
		sys.exit()

if xmin != -999 and xmax != -999:
	R = (xmax-xmin)/2
	print '\nSystem radius: {0:9.6g} cm'.format(R)
	print '           or: {0:9.6g} AU'.format(R/AU)
	print '           or: {0:9.6g} pc'.format(R/pc)
else:
	sys.exit('Error calculating system radius. Could not identify either xmin or xmax in flash.par')

print '\nSystem mass:   {0:9.6g} g'.format(M_total)
print '         or, {0:9.6g} Msun'.format(M_total/Msun)
	

mu = M_total/(np.pi * R**2 * mag_field_val) * np.sqrt(G)/0.13

print '\nThe mass-to-flux ratio, normalized to the critical mass-to-flux ratio, is:'.format(mu)
print '\n            {0}'.format(mu)

