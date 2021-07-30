# pywavan3D
This code are devolped during my M1 stage whith Jean-François Robitaille and Frederique Motte at IPAG.    
It's a 3D version of MnGSeg developed by Jean-François Robitaille. This code make also new visualisation of a 3D data cube, including moment 0 and moment 1 map for the coherent part and animation of this moment map.   

This repository is organised as follows : 
- code : all useful fonction 
- g353 and w43 : all notebook for this region
- g353_feathered_isolated and w43_7+12m isolated : notebook and annimation make for this region
- 3d_test : test of a 3d wavelet, not workink

## Necesary package :
pywavan : https://github.com/jfrob27/pywavan           
spectral cube : https://spectral-cube.readthedocs.io/en/latest/           
astropy, numpy, mathplotlib, subprocess

Animationare make with gifsicle

## overview fonction :
### Pywavan3D :
 - fan_trans3D() : 3D version of fan_trans from pywavan
 - partClean() : clean the cube, to be used before make moment 1 map  

### Power Spectrum plot :
 - PSplot() : plot simply a power specctrum
 - animPSplot() : make animation of the power spectrum plot, whith or whithout gaussian and/or non-gaussian component

### moment 1 spectral cube : 
 - moment1Cube() : make the moment 1 map of a cube
 - mom1_anim_center() and mom1_anim(), tow way for animate moment 1 map
