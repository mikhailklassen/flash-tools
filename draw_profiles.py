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

# Get checkpoint files and plotfiles
print 'Agglomerating checkpoint files.'
chkfiles = glob.glob(datadir+'/*hdf5_chk_*')
chkfiles.sort()
print 'Found {0} checkpoint files.\n'.format(len(chkfiles))

print 'Agglomerating plot files.'
pltfiles = glob.glob(datadir+'/*hdf5_plt_cnt_*')
pltfiles.sort()
print 'Found {0} plot files.\n'.format(len(pltfiles))

# Create output directory
print 'Creating output directory.'
if not os.path.exists(outdir):
    os.mkdir(outdir)
print 'Created {0}.\n'.format(outdir)

## Make a plot of the evolution of the density profile
from flash_generic import profiles
profiles.density_profile_1D_evolution(pltfiles,outdir)

