import csv
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
#====================================================================
data = pd.read_csv("dataMAIN.csv",header = None)
#print('\nAlgandmed\n',data.head())
origin = 'lower'
Tmax = len(data)
T = data.drop(data.columns[[Tmax]], axis=1)
N = len(T)
#print('\nAlgandmed\n',T.head())
# Implementation of matplotlib function
Z = T.values
Zf = np.fft.fft(T.values)
#print(Z)
#====================================================================
# Basically everu pixel is stored in a separate cell
# for this plotting method.
fig = plt.figure(figsize=(10,8))
colbar = plt.imshow(Zf.real, cmap='viridis', interpolation='none', 
                    extent=[0, N, 0, N], origin='lower', aspect='auto')
plt.clim(-1, 1) #fixes limits of color bar
fig.colorbar(colbar)
plt.show()
# For saving the images
fig.savefig('pixel_fourier.png', format='png', dpi=600)