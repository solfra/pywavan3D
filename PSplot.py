from PIL.Image import NONE
import matplotlib.pyplot as plt
import numpy as np
import aplpy
from astropy.io import fits
from spectral_cube import SpectralCube
import astropy.units as u
import subprocess
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic

def PSplot(n, wk="wave_k.npy",total="s1at.npy",gaussian="s1ag.npy",nonGauss="s1ang.npy") :
    """
    Fonction for plot and show a power spectrum
    Load the element return by fan_trans3D()

    input : 
    - n : image postion in velocity field for witch make power spectrum
    - wk,total,gaussian,nonGauss (option, by defult the same as the fan_trans return) : npy files to loads for make the plot. 

    return :
    Nothing, just showthe plot 
    """
    t = np.load(total)
    g = np.load(gaussian)
    ng = np.load(nonGauss)
    w_k = np.load(wk)

    plt.plot(w_k, t[n],label="total")
    plt.plot(w_k, g[n],"o",label="gaussian")
    plt.plot(w_k, ng[n],"x",label="non gaussian")
    plt.xscale("log")
    plt.yscale("log")
    plt.ylabel('P(k)')
    plt.xlabel('k (arcmin$^{-1}$)')
    plt.title('canal number {}'.format(n))
    plt.legend()
    plt.show()

def animPSplot(wk="wave_k.npy",total="s1at.npy",gaussian="s1ag.npy",nonGauss="s1ang.npy", lim=[10**-10,1], sauv = 0, sauvName = '',**kwargs ):
    """
    Fonction for plot and show power spectrum animation
    permit to show also the non gaussian image or the gaussian image or the two at the same time
    Load the element return by fan_trans3D()

    input : 
    - wk,total,gaussian,nonGauss (option, by defult the same as the fan_trans return) : npy files to loads for make the plot. 
    - lim : list, ylim for fix the scale during the animation 
    - sauv (option, default 0) : save parameter of the animation, if 1 : save, if 0 : no save
    - sauvName (option, default '') : add detail at the name of the animation. 
        By default the name is animPS.gif. Add name precision befor the .gif
        Attention, do not use space in the name !

    keyword sup :
     - ngImg : str, place of the non-gaussian fits cube
     - gImg : str, place of the gaussian fits cube

    return :
    Nothing, just showthe plot 
    """
    print(kwargs)
    t = np.load(total)
    g = np.load(gaussian)
    ng = np.load(nonGauss)
    w_k = np.load(wk)
    fig_all = plt.figure(1, figsize=(20,10))

    N = len(t) #determine the len of the animation

    for i in range(N) :
        if 'ngImg' or 'gImg' in kwargs :
            plt.subplot(121)
        plt.plot(w_k, t[i],label="total")
        plt.plot(w_k, g[i],"o",label="gaussian")
        plt.plot(w_k, ng[i],"x",label="non gaussian")
        plt.xscale("log")
        plt.yscale("log")
        plt.ylim(lim[0],lim[1]) 
        plt.title("canal number {}".format(i))
        plt.legend()

        if 'ngImg' in kwargs and 'gImg' in kwargs:
            filname = kwargs.get('ngImg')
            filname2 = kwargs.get('gImg')
            fig = aplpy.FITSFigure(filname, figure=fig_all,slices=[i],subplot=(2,2,4))
            fig.show_colorscale(cmap='gray')
            fig.add_colorbar()
            fig.set_title("coherent part number {}".format(i))

            fig = aplpy.FITSFigure(filname2, figure=fig_all,slices=[i], subplot=(2,2,2))
            fig.show_colorscale(cmap='gray')
            fig.add_colorbar()
            fig.set_title("gaussian part number {}".format(i))

        elif 'ngImg' in kwargs :
            filname = kwargs.get('ngImg')
            print('ng')
            fig = aplpy.FITSFigure(filname, figure=fig_all,slices=[i],subplot=(1, 2, 2))
            fig.show_colorscale(cmap='gray')
            fig.add_colorbar()
            fig.set_title("coherent part number {}".format(i))

        elif 'gImg' in kwargs :
            filname = kwargs.get('gImg')
            print('g')
            fig = aplpy.FITSFigure(filname, figure=fig_all,slices=[i],subplot=(1, 2, 2))
            fig.show_colorscale(cmap='gray')
            fig.add_colorbar()
            fig.set_title("gaussian part number {}".format(i))


        plt.pause(0.04)
        if sauv :
            plt.savefig('image_' + str(i).zfill(5) + '.png', dpi=100) # sauve fichier image (pour ensuite gif anime). dpi: precision de l'image
        plt.clf()
        plt.cla()
    plt.show()

    if sauv :
        subprocess.getoutput('convert image_*.png GIF:- | gifsicle --delay=10 --loop --optimize=2 --colors=256 --multifile - > animPS{}.gif'.format(sauvName)) # avec gifsicle (plus efficace).
        subprocess.getoutput('rm image_*.png') # efface les fichiers images temporaires

