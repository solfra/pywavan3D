import matplotlib.pyplot as plt
import numpy as np
import aplpy
from astropy.io import fits
from numpy.core.defchararray import center
from spectral_cube import SpectralCube
import astropy.units as u
import subprocess
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic


def mom1_anim_center(filename, regionName ,sauv=0, sauvName = ''):
    """
    Fonction who make an animation of the moment 1 map for a data cube in fits format
    The animation is performed from the center of the cube to the end
    Assume the cube is a Freq cube in radio
    Expresses the speed in kilometer per second

    input : 
    - filname : name of the fits cube 
    - regionName : the name of the region for make map title
    - sauv : save parameter of the animation, if 1 : save, if 0 : no save. 0 by default
    - sauvName : add detail at the name of the animation. 
            By default the name is animation_mom1_vFix_center.gif. Add name precision befor the .gif. 
            Attention, do not use space in the name !

    Return :
    Nothing, make and show plot
    """
    cube = SpectralCube.read(filename)  

    HDU = fits.open(filename)
    header = HDU[0].header
    w = wcs.WCS(header)

    N = header['NAXIS3']
    fig_all = plt.figure(1, figsize=(20,10))

    if N%2 == 1 :
        L = N/2+1

    else :
        L = N/2

    nii_cube = cube.with_spectral_unit(u.km/u.s,velocity_convention='radio', rest_value=header['RESTFRQ']*u.Hz)
    moment_1 = nii_cube.moment(order=1)
    vmin = np.nanmin(moment_1.hdu.data)-1
    vmax = np.nanmax(moment_1.hdu.data)+1

    for i in range(int(L)) :
        if N%2 == 1 :
            c = int(N/2)
            wmin = w.pixel_to_world(0,0,c-i)[1]
            wmax = w.pixel_to_world(0,0,c+i)[1]
        if N%2 == 0 :
            c = int(N/2)
            wmin = w.pixel_to_world(0,0,c-1-i)[1]
            wmax = w.pixel_to_world(0,0,c+i)[1]
        nii_subcube = cube.spectral_slab(wmin,wmax)
        nii_cube = nii_subcube.with_spectral_unit(u.km/u.s,velocity_convention='radio', rest_value=header['RESTFRQ']*u.Hz)
        moment_1 = nii_cube.moment(order=1)

        ax = plt.gca()
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
        f = aplpy.FITSFigure(moment_1.hdu,figure=fig_all)  
        f.show_colorscale(cmap = 'jet', vmin = vmin, vmax =vmax)
        f.add_colorbar()
        f.set_title("{} \n from {} to {} ".format(regionName,wmin, wmax))
        plt.pause(0.1)

        if sauv :
            plt.savefig('image_' + str(i).zfill(5) + '.png', dpi=100) # sauve fichier image (pour ensuite gif anime). dpi: precision de l'image
        
        plt.clf()
        plt.cla()

    plt.show()

    if sauv :
        subprocess.getoutput('convert image_*.png GIF:- | gifsicle --delay=10 --loop --optimize=2 --colors=256 --multifile - > animation_mom1_vFix_center{}.gif'.format(sauvName)) # avec gifsicle (plus efficace).
        subprocess.getoutput('rm image_*.png') # efface les fichiers images temporaires

def mom1_anim(filename, regionName = '', zero =0, vfix = True ,sauv=0, sauvName = ''):
    """
    Fonction who make an animation of the moment 1 map for a data cube in fits format
    The animation is performed rom one edge of the cube to the other
    Assume the cube is a Freq cube in radio
    Expresses the speed in kilometer per second

    input : 
    - filname : name of the fits cube 
    - regionName (option) : the name of the region for make map title
    - zero (option) : define where the cube start. Alows to reject the first chanel if it's juste noise 0 by default
    - vfix (option) : define if the color scale is fixe or not. True by default 
    - sauv (option) : save parameter of the animation, if 1 : save, if 0 : no save. 0 by default
    - sauvName (option) : add detail at the name of the animation. 
            By default the name is animation_mom1_vFix_center.gif. Add name precision befor the .gif
            Attention, do not use space in the name !

    Return :
    Nothing, make and show plot
    """
    cube = SpectralCube.read(filename)  

    HDU = fits.open(filename)
    header = HDU[0].header
    w = wcs.WCS(header)

    N = header['NAXIS3']
    fig_all = plt.figure(1, figsize=(20,10))

    if vfix :
        nii_cube = cube.with_spectral_unit(u.km/u.s,velocity_convention='radio', rest_value=header['RESTFRQ']*u.Hz)
        moment_1 = nii_cube.moment(order=1)
        vmin = np.nanmin(moment_1.hdu.data)-1
        vmax = np.nanmax(moment_1.hdu.data)+1

    for i in range(zero,N) :
        wmin = w.pixel_to_world(0,0,zero)[1]
        wmax = w.pixel_to_world(0,0,i)[1]
        nii_subcube = cube.spectral_slab(wmin,wmax)
        nii_cube = nii_subcube.with_spectral_unit(u.km/u.s,velocity_convention='radio', rest_value=header['RESTFRQ']*u.Hz)
        moment_1 = nii_cube.moment(order=1)
        
        ax = plt.gca()
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
        f = aplpy.FITSFigure(moment_1.hdu,figure=fig_all)
        if vfix :  
            f.show_colorscale(cmap = 'jet', vmin = vmin, vmax =vmax)
        else :
            f.show_colorscale(cmap = 'jet')
        f.add_colorbar()
        f.set_title("{} \n from {} to {} ".format(regionName,wmin, wmax))
        plt.pause(0.1)

        if sauv :
            plt.savefig('image_' + str(i).zfill(5) + '.png', dpi=100) # sauve fichier image (pour ensuite gif anime). dpi: precision de l'image
        
        plt.clf()
        plt.cla()

    plt.show()

    if sauv :
        subprocess.getoutput('convert image_*.png GIF:- | gifsicle --delay=10 --loop --optimize=2 --colors=256 --multifile - > animation_mom1_vFix{}.gif'.format(sauvName)) # avec gifsicle (plus efficace).
        subprocess.getoutput('rm image_*.png') # efface les fichiers images temporaires

def moment1Cube(filename, regionName = '', zero = 0, save = 0, sauvName ='' ):
    """
    Fonction who make a moment 1 map for a data cube in fits format
    The animation is performed rom one edge of the cube to the other
    Assume the cube is a Freq cube in radio
    Expresses the speed in kilometer per second

    input : 
    - filname : name of the fits cube 
    - regionName (option) : the name of the region for make map title
    - zero (option) : define where the cube start. Alows to reject the first chanel if it's juste noise 0 by default
    - sauv (option) : save parameter of moment map in fi fits, 1 : save, 0 : no save. 0 by default
    - sauvName (option) : add detail at the name of the animation. 
            By default the name is moment1.fits. Add name precision befor the .fits
            Attention, do not use space in the name !

    Return :
    Nothing, make and show plot
    """

    cube = SpectralCube.read(filename)  
    HDU = fits.open(filename)
    header = HDU[0].header 
    w = wcs.WCS(header)
    N = header['NAXIS3']

    wmin = w.pixel_to_world(0,0,zero)[1]
    wmax = w.pixel_to_world(0,0,N)[1]
    nii_subcube = cube.spectral_slab(wmin,wmax)
    nii_cube = nii_subcube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=header['RESTFRQ']*u.Hz)
    moment_1 = nii_cube.moment(order=1)

    f = aplpy.FITSFigure(moment_1.hdu)  
    f.show_colorscale(cmap='jet')
    f.add_colorbar()
    f.set_title("{} \n from {} to {} ".format(regionName,wmin, wmax))
    plt.show()

    if save == 1 :    
        moment_1.write('moment1{}.fits'.format(sauvName))
