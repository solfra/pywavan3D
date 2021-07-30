import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

def mom0(filename, affiche = 1, sauv = 0, sauvName = '') :
    """
    Calcule the moment 0 of a cube whithout using spectralcube

    input : 
    - filname : name of the fits cube 
    - affiche (option 1 by default) : show (1) or not (0) the moment 0 map
    - sauv (optional, default 0) : save parameter of the animation, if 1 : save, if 0 : no save
    - sauvName (optional, default '') : add detail at the name of the animation. 
            By default the name is animation_mom1_vFix_center.gif. Add name precision befor the .gif. 
            Attention, do not use space in the name !
    """
    HDU = fits.open(filename)
    cube = HDU[0].data
    header = HDU[0].header
    mom0 = np.sum(cube,axis=0)

    if affiche :
        plt.imshow(mom0,origin="lower")
        plt.title('moment 0 original cube') 
        plt.show()

    if sauv :
        del header['NAXIS3']
        del header['PC3_1']
        del header['PC3_2']
        del header['PC1_3']
        del header['PC2_3']
        del header['PC3_3']
        del header['CTYPE3']
        del header['CRVAL3']
        del header['CDELT3']
        del header['CRPIX3']
        del header['CUNIT3']
        header['NAXIS']=2
        fits.writeto("mom0{}.fits".format(sauvName),mom0, header,overwrite= True) 