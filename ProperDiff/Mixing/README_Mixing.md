# Documenation for sinusoidal mixing homework

The main file that does the computations is the file 'mixing.c'.
Plotting is done with a Python file named 'plot_mixF.py'

## Specifics of file 'mixing.c'
Firstly, the main part contains two functions that handle the diffusive sinusoidal mixing: 'p_diffuse' and 'iterate'. Essentially, 'iterate' calls the 'p_diffuse' a N amount of times until we have proper diffusion. I found 11 to be quite sufficient. Otherwise the images get very pixelated.  
---
Secondly, there is also great effort put into doing Fourier transform in this 'mixing.c' file as well. Main goal was to do the same thing as numpy.fft.fftshift() and numpy.fft.fft2() do in Python. For this, it is necessary to download the FFTW library and link it with the C compiler. In Ubuntu linux, it can be done by the command 'sudo apt-get install libfftw3-dev'. So in the C code, it consists of three main functions to handle the Fourier transform: fftshift, fft2 and fft1. 
