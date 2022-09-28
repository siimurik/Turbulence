import numpy as np
import matplotlib.pyplot as plt

N = 128
a = 25.0
c = np.zeros((N,N))
c_new = c

for i in range(N):
    for j in range(N):  
        c[i][j] = np.sin(i * np.pi * 2/N)

for i in range(N):
    for j in range(N):
        c_new[i][j] = c[(i + int(a*np.sin(j*2*np.pi/N))) & N-1][j] 

for i in range(N):
    for j in range(N):
        c[i][j] = c_new[i][(j + int(a*np.sin(i*2*np.pi/N))) & N-1]

#c = c_new.copy()


feature_x = np.linspace(0, 1, N)
feature_y = np.linspace(0, 1, N)
#print(feature_x)

# Creating 2-D grid of features
[X, Y] = np.meshgrid(feature_x, feature_y)

fig, ax = plt.subplots(1, 1)
# plots filled contour plot
ax.contourf(X, Y, c)
  
ax.set_title('Filled Contour Plot')
ax.set_xlabel('feature_x')
ax.set_ylabel('feature_y')
  
plt.show()