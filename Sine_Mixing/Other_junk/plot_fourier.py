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
Zf = np.fft.fft(Z)
n = Z.size
Zfreq = np.fft.fftfreq(n)

#====================================================================
"""
fig, ax = plt.subplots(figsize=(12, 10))
# plots filled contour plot
ax.plot(Zfreq, Zf.real.values)
plt.show()
"""
"""
yf = np.fft.fft(Z,1,1,"forward")
#print(yf)
xf = np.fft.fftfreq(N)
plt.imshow(Zfreq,Zf.real)
plt.show()
"""
xf = np.fft.fftfreq(N)
nx = len(T)
#print(nx)
feature_x = np.linspace(0, 1, nx)
feature_y = np.linspace(0, 1, nx)
#print(feature_x)

# Creating 2-D grid of features
[X, Y] = np.meshgrid(xf, xf)
#[X, Y] = np.meshgrid(x_feat, y_feat)

fig, ax = plt.subplots(figsize=(12, 10))
# plots filled contour plot
ax.contourf(X, Y, Zf.real)
  
ax.set_title('Filled Contour Plot')
ax.set_xlabel('feature_x')
ax.set_ylabel('feature_y')
ax.set_aspect('equal', adjustable='box')
plt.show()