import numpy as np
from pywavan import fan_trans,  nb_scale
from astropy.io import fits

def fan_trans3D(filename,nx, ny, fitsExp = False):
    """
    3D version of fan_trans from pywavan
    For this versio, fan_trans using the argument apodize and arrdim
    q is set to 2.0 and not using qdyn parmeter

    input :
    filename : fits data cube name
    nx, ny sizes including the zero value pixel padding
    fitsExp (option , false by default) : option for export coherent and gaussian part in fits
    """

    HDU = fits.open(filename)
    cube = HDU[0].data
    header = HDU[0].header
    reso = header ['CDELT2']*60
    M = nb_scale((nx,ny))
    N = header['NAXIS3']

    coherent_tot = []
    gaussian_tot = []
    s1a_tot_tot = []
    s1a_tot_gau = []
    s1a_tot_ng = []

    for i in range(N) :
        q = []
        q= [2.0]*M
        print("data number",i)
        wt, S11a, wave_k, S1a, q =  fan_trans(cube[i,:,:], reso=reso, angular=False,q=q,apodize = 0.98, arrdim = np.array([nx,ny]))

        coherent = np.sum(wt[M:2*M,:,:],axis=0)
        Gaussian = np.sum(wt[2*M:3*M,:,:],axis=0) 

        coherent_tot.append( coherent )
        gaussian_tot.append( Gaussian )
        s1a_tot_tot.append( S1a[0,:] )
        s1a_tot_ng.append(S1a[1,:])
        s1a_tot_gau.append(S1a[2,:])

    np.save("coh.npy", coherent_tot)
    np.save("gau.npy", gaussian_tot)
    np.save("s1at.npy",s1a_tot_tot)
    np.save("s1ag.npy",s1a_tot_gau)
    np.save("s1ang.npy",s1a_tot_ng)
    np.save("wave_k.npy",wave_k)

    if fitsExp :
        fits.writeto("nonGaussian.fits",coherent_tot.real, header, overwrite = True) 
        fits.writeto("Gaussian.fits",gaussian_tot.real, header, overwrite = True)

def partClean(part, cube) :
    img = np.load("/user/workdir/soldanof/data/w43_7_12_iso_cnts/cohP.npy")
    HDU = fits.open("/user/workdir/soldanof/ALMA/W43-MM1/W43-MM1_B3_spw0_7M12M_n2hp.image-isolated-contsub-crop_cut.fits")
    cube = HDU[0].data
    header = HDU[0].header
    N = header['NAXIS3']

    for i in range(N):
        img[i,:,:] += np.mean(cube[i,:,:])

        im_rmv = np.zeros((img.shape[1],img.shape[2]))
        img_2d = img[i,:,:]
        im_rmv[img_2d<=np.abs(cube[i,:,:]).mean()] = np.nan
        img[i,:,:] = img[i,:,:] + im_rmv