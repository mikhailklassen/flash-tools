#!/usr/bin/env python
"""
Example script that takes a plot file as input and samples the
data values along a ray. This is a LineOut query in VisIt, but
much simpler to perform. The ray data structure is a Python
dictionary with the grid variables as keys. The 't' key is 
used to step along the ray. This is more useful when the ray
is not directly along one of the coordinate axes.
"""
import numpy as np
import matplotlib.pyplot as plt
from yt.mods import *
import sys

pltfile = sys.argv[1]

pf = load(pltfile)

ray = pf.h.ray([0,0,0],[1,0,0]/pf['unitary'])

radius = ray['x']
t = ray['t']
dens = ray['Dens']
temp = ray['Temp']
