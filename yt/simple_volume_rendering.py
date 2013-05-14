'''
Based on the Simple Volume Rendering script from yt-project.org.
http://yt-project.org/doc/cookbook/simple_plots.html#simple-volume-rendering

Usage:
    python simple_volume_rendering.py PLOTFILE
'''
from yt.mods import *
import sys

# Load the dataset.
plotfile = sys.argv[1]
pf = load(plotfile)

# Create a data container (like a sphere or region) that
# represents the entire domain.
dd = pf.h.all_data()

# Get the minimum and maximum densities.
mi, ma = dd.quantities["Extrema"]("Density")[0]

# Parameters
numlevels = 5                           # Number of Gaussians in the transfer function
widths = 0.02                           # Width of the Gaussians
cmap = "spectral"                       # Color map
center = [0.5, 0.5, 0.5]/pf["unitary"]  # Center of the volume rendering
look = [0.5, 0.2, 0.7]/pf["unitary"]    # Look direction for the camera
width = 1.0/pf["unitary"]               # Width of the image. Use this to zoom in/out
resolution = 512                        # Number of pixels in each dimension of the final image

# Create a transfer function to map field values to colors.
# We bump up our minimum to cut out some of the background fluid
tf = ColorTransferFunction((np.log10(mi)+1, np.log10(ma)))

# Add several Gausians to the transfer function, evenly spaced
# between the min and max specified above with widths w
tf.add_layers(numlevels, w=widths, colormap=cmap)

# Create the camera object
cam = pf.h.camera(center, look, width, resolution, tf)

# Create a snapshot.
# The return value of this function could also be accepted, modified (or saved
# for later manipulation) and then put written out using write_bitmap.
# clip_ratio applies a maximum to the function, which is set to that value
# times the .std() of the array.
cam.snapshot("%s_volume_rendered.png" % pf, clip_ratio=8.0)
