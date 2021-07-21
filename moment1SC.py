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


def mom1_anim_center(filename, regionName='' ,sauv=0, sauvName = '', **kwargs):
    """
    Fonction who make an animation of the moment 1 map for a data cube in fits format
    The animation is performed from the center of the cube to the end
    Assume the cube is a Freq cube in radio
    Expresses the speed in kilometer per second

    input : 
    - filname : name of the fits cube 
    - regionName (optional, default ''): the name of the region for make map title
    - sauv (optional, default 0) : save parameter of the animation, if 1 : save, if 0 : no save
    - sauvName (optional, default '') : add detail at the name of the animation. 
            By default the name is animation_mom1_vFix_center.gif. Add name precision befor the .gif. 
            Attention, do not use space in the name !

    keyword :
    - cores : list of position of the core in pixel coordinate. Use world_to_pixel for exemple for have pixel position
                    list type [(x1,y1), (x2,y2)]

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
        if 'cores' in kwargs :
            lcore = kwargs.get('cores')
            for c in lcore :
                f.show_markers(c[0],c[1],coords_frame='pixel',facecolor='k')

        plt.pause(0.1)

        if sauv :
            plt.savefig('image_' + str(i).zfill(5) + '.png', dpi=100) # sauve fichier image (pour ensuite gif anime). dpi: precision de l'image
        
        plt.clf()
        plt.cla()

    plt.show()

    if sauv :
        subprocess.getoutput('convert image_*.png GIF:- | gifsicle --delay=10 --loop --optimize=2 --colors=256 --multifile - > animation_mom1_vFix_center{}.gif'.format(sauvName)) # avec gifsicle (plus efficace).
        subprocess.getoutput('rm image_*.png') # efface les fichiers images temporaires

def mom1_anim(filename, regionName = '', zero =0, vfix = True ,sauv=0, sauvName = '', **kwargs):
    """
    Fonction who make an animation of the moment 1 map for a data cube in fits format
    The animation is performed rom one edge of the cube to the other
    Assume the cube is a Freq cube in radio
    Expresses the speed in kilometer per second

    input : 
    - filname : name of the fits cube 
    - regionName (option, default '') : the name of the region for make map title
    - zero (option, default 0) : define where the cube start. Alows to reject the first chanel if it's juste noise
    - vfix (option, default True) : define if the color scale is fixe or not
    - sauv (option, default 0) : save parameter of the animation, if 1 : save, if 0 : no save
    - sauvName (option, default '') : add detail at the name of the animation. 
            By default the name is animation_mom1.gif. Add name precision befor the .gif
            Attention, do not use space in the name !

    keyword :
    - cores : list of position of the core in pixel coordinate. Use world_to_pixel for exemple for have pixel position
                    list type [(x1,y1), (x2,y2)]        

    Return :
    Nothing, make and show plot
    """
    cube = SpectralCube.read(filename)  

    HDU = fits.open(filename)
    header = HDU[0].header
    w = wcs.WCS(header)
    N = header['NAXIS3'] #len of the animation

    fig_all = plt.figure(1, figsize=(20,10))

    if vfix :
        # calculate total moment 1 for eterminate max and min
        nii_cube = cube.with_spectral_unit(u.km/u.s,velocity_convention='radio', rest_value=header['RESTFRQ']*u.Hz)
        moment_1 = nii_cube.moment(order=1)
        vmin = np.nanmin(moment_1.hdu.data)-1
        vmax = np.nanmax(moment_1.hdu.data)+1

    for i in range(zero,N) :
        wmin = w.pixel_to_world(0,0,zero)[1] #freq min
        wmax = w.pixel_to_world(0,0,i)[1] #freq max

        nii_subcube = cube.spectral_slab(wmin,wmax)
        nii_cube = nii_subcube.with_spectral_unit(u.km/u.s,velocity_convention='radio', rest_value=header['RESTFRQ']*u.Hz)
        moment_1 = nii_cube.moment(order=1)

        ax = plt.gca() #sert axes for plot
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)

        f = aplpy.FITSFigure(moment_1.hdu,figure=fig_all)
        if vfix :  
            f.show_colorscale(cmap = 'jet', vmin = vmin, vmax =vmax)
        else :
            f.show_colorscale(cmap = 'jet')
        f.add_colorbar()
        f.set_title("{} \n from {} to {} ".format(regionName,wmin, wmax))

        if 'cores' in kwargs :
            lcore = kwargs.get('cores')
            for c in lcore :
                f.show_markers(c[0],c[1],coords_frame='pixel',facecolor='k')

        plt.pause(0.1)

        if sauv :
            plt.savefig('image_' + str(i).zfill(5) + '.png', dpi=100) # sauve fichier image (pour ensuite gif anime). dpi: precision de l'image
        
        plt.clf() #remove figure in the plot
        plt.cla()

    plt.show()

    if sauv :
        subprocess.getoutput('convert image_*.png GIF:- | gifsicle --delay=10 --loop --optimize=2 --colors=256 --multifile - > animation_mom1{}.gif'.format(sauvName)) # avec gifsicle (plus efficace).
        subprocess.getoutput('rm image_*.png') # efface les fichiers images temporaires

def moment1Cube(filename, regionName = '', zero = 0, save = 0, sauvName ='' ,**kwargs):
    """
    Fonction who make a moment 1 map for a data cube in fits format
    The animation is performed rom one edge of the cube to the other
    Assume the cube is a Freq cube in radio
    Expresses the speed in kilometer per second

    input : 
    - filname : name of the fits cube 
    - regionName (option, default '') : the name of the region for make map title
    - zero (option, default 0) : define where the cube start. Alows to reject the first chanel if it's juste noise
    - sauv (option, default 0) : save parameter of moment map in fi fits, 1 : save, 0 : no save
    - sauvName (option, default '') : add detail at the name of the animation. 
            By default the name is moment1.fits. Add name precision befor the .fits
            Attention, do not use space in the name !

    keyword :
    - cores : list of position of the core in pixel coordinate. Use world_to_pixel for exemple for have pixel position
                    list type [(x1,y1), (x2,y2)]

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
    if 'cores' in kwargs :
        lcore = kwargs.get('cores')
        for c in lcore :
            f.show_markers(c[0],c[1],coords_frame='pixel',facecolor='k')
    plt.show()

    if save == 1 :    
        moment_1.write('moment1{}.fits'.format(sauvName))

def imgCompar(filname1, filname2, title ='', trigger = 0.2):
    """
    Fonction for compare 2 fits image
    Developed for compare moment 1 map coherent and total

    input :
    filname1, filname2 : the two filename of the image to compare
    title (optional) : title of the plot
    trigger (optional, default 0.2) : float, trigger to remove value close to 0 when make the image difference
    """

    HDU = fits.open(filname1)
    cube1 = HDU[0].data
    cube1[np.isnan(cube1)]=0

    HDU = fits.open(filname2)
    cube2 = HDU[0].data
    cube2[np.isnan(cube2)]=0

    diff = cube1-cube2
    diff[np.abs(diff)<=trigger]=np.nan

    plt.figure(figsize=(20,20))
    plt.imshow(diff,origin="lower",cmap='jet')
    plt.colorbar()
    plt.title("{}".format(title))
    plt.show()