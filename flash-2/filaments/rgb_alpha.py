def rgb_alpha_w(colors,alpha):
    '''
    Assuming a white background, compute the new RGB values for an input RGB color array
    of the form [r,g,b] and an alpha transparency factor. Note, Python handles color
    values as ranging [0-1] and similarly, alpha ranges [0-1]. 
    '''
    import numpy as np
    colors = np.array(colors)
    blended_rgb = alpha * colors + (1 - alpha) * np.array([1,1,1])

    return blended_rgb

def rgb_alpha_k(colors,alpha):
    '''
    Assuming a black background, compute the new RGB values for an input RGB color array
    of the form [r,g,b] and an alpha transparency factor. Note, Python handles color
    values as ranging [0-1] and similarly, alpha ranges [0-1]. 
    '''
    import numpy as np
    colors = np.array(colors)
    blended_rgb = alpha * colors + (1 - alpha) * np.array([0,0,0])

    return blended_rgb

def blend(color1,color2,blend_factor):
    '''
    Blends two colors, specified as RGB arrays with values ranging [0-1]. The blending
    factor controls the degree to which color1 is blended with color2.
    '''
    import numpy as np
    color1 = np.array(color1)
    color2 = np.array(color2)
    blended_rgb = blend_factor * color1 + (1 - blend_factor) * color2

    return blended_rgb
