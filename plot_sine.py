import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker, cm

data = pd.read_csv("dataMAIN.csv",header = None)
xydf = pd.read_csv("x_y.csv", header = None)
#print('\nAlgandmed\n',data.head())
origin = 'lower'
Tmax = len(data)
T = data.drop(data.columns[[Tmax]], axis=1)
N = len(T)
#print('\nAlgandmed\n',T.head())
x_feat = xydf.drop(data.columns[[1]], axis=1)
y_feat = xydf.drop(data.columns[[0]], axis=1)
#print(x_feat.head())
# Implementation of matplotlib function
Z = T.values
#print(Z)
fig = plt.figure(figsize=(12,8))
colbar = plt.imshow(Z, cmap='viridis', interpolation='none', 
                    extent=[0, N, 0, N], origin='lower', aspect='auto')
plt.clim(-1, 1) #fixes limits of color bar
#cmap variante: Greys, hot, terrain, rainbow, cividis, viridis. viimane loengu näitele kõige sarnasem
fig.colorbar(colbar)

plt.show()
"""
nx = len(T)
#print(nx)
feature_x = np.linspace(0, 1, nx)
feature_y = np.linspace(0, 1, nx)
#print(feature_x)

# Creating 2-D grid of features
[X, Y] = np.meshgrid(feature_x, feature_y)
#[X, Y] = np.meshgrid(x_feat, y_feat)

fig, ax = plt.subplots(figsize=(12, 10))
# plots filled contour plot
ax.contourf(X, Y, T)
  
ax.set_title('Filled Contour Plot')
ax.set_xlabel('feature_x')
ax.set_ylabel('feature_y')
ax.set_aspect('equal', adjustable='box')
plt.show()
"""
