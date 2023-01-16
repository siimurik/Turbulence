# Documenation for sinusoidal mixing homework

The main file that does the computations is the file 'mixing.c'.
Plotting is done with a Python file named 'plot_mixF.py'.

## How to run the codes
I have written a simple bash .sh file that can be run with the command:

    bash run_mixing.sh

This bash file compiles the C code and plots it with Python.
To compile the C file, you have to use the command

    gcc mixing.c -o mix -lm -lfftw3

and then plot the data with the command

    python3 plot_mixF.py

## Images
There are two sets of images, main difference comes from what is used to take the Fourier transform: Python's NumPy commands or functions in the C code. 
Python's plot looks better, but the one generated with the C code gives a better overview what is done.

## Specifics of file 'mixing.c'
Firstly, the main part contains two functions that handle the diffusive sinusoidal mixing: 'p_diffuse' and 'iterate'. Essentially, 'iterate' calls the 'p_diffuse' a N amount of times until we have proper diffusion. I found 11 to be quite sufficient. Otherwise the images get very pixelated.  

The output of the mixing computations is written in the file 'colors.csv'.

Secondly, there is also great effort put into doing Fourier transform in this 'mixing.c' file as well. Main goal was to do the same thing as numpy.fft.fftshift() and numpy.fft.fft2() do in Python. For this, it is necessary to download the FFTW library and link it with the C compiler. In Ubuntu linux, it can be done by the command 

    sudo apt-get install libfftw3-dev

So in the C code, it consists of three main functions to handle the Fourier transform: fftshift, fft2 and fft1. Further documenation is written in the 'mixing.c' file.

The output of the Fourier transform computations is written in the file 'fourier.csv'.
