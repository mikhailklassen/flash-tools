import matplotlib.pyplot as plt
import skimage
#from skimage.feature import hog
from skimage import data, color, exposure
from filament_tools import hro

fname = 'hydro_filaments.jpg'
#fname = 'test6.bmp'
image = skimage.io.imread(fname,as_grey=True)

fd, hog_image = hro(image, orientations=8, pixels_per_cell=(10, 10),
                    cells_per_block=(1, 1), visualise=True)

plt.figure(figsize=(8, 4))

plt.subplot(121).set_axis_off()
plt.imshow(image, cmap=plt.cm.gray)
plt.title('Input image')

# Rescale histogram for better display
hog_image_rescaled = exposure.rescale_intensity(hog_image, in_range=(0, 0.02))

plt.subplot(122).set_axis_off()
plt.imshow(hog_image_rescaled, cmap=plt.cm.gray)
plt.title('Histogram of Oriented Gradients')
plt.show()
