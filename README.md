# pywavan3D
This code are devolped during my M1 stage whith Jean-François Robitaill and Frederique Motte at IPAG.    
It's a 3D version of MnGSeg developed by Jean-François Robitaill. This code make also new visualisation of a 3D data cube, including moment 0 and moment 1 map for the coherent part.

## Necesary package :
pywavan : https://github.com/jfrob27/pywavan           
spectral cube : https://spectral-cube.readthedocs.io/en/latest/           
astropy, numpy, mathplotlib, subprocess

Animationare make with gifsicle

## overview fonction :
### Pywavan3D :
 - fan_trans3D() : 3D version of fan_trans from pywavan

### Power Spectrum plot :
 - PSplot() : plot simply a power specctrum
 - animPSplot() : make animation of the power spectrum plot, whith or whithout gaussian and/or non-gaussian component

### moment 1 : 
 - moment1Cube() : make the moment 1 map of a cube
 - mom1_anim_center() and mom1_anim(), tow way for animate moment 1 map
