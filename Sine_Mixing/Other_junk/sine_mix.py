import numpy as np
import matplotlib.pyplot as plt

N = 100
a = 1.0
c = np.zeros((N,N))

for i in range(N):
    for j in range(N):
        d = int(a*np.sin(j*2*np.pi/N))
        if (i + d) < N:
            c_new[i+d][j] = c[i][j]
        elif (i + d) >= N:
            c_new[i+d-N][j] = c[i][j]
        else:
            c_new[i+d-N][j] = c[i][j]
        c[i][j] = c[j][i]

c_new = c.copy()
feature_x = np.linspace(0, 1, N)
feature_y = np.linspace(0, 1, N)
#print(feature_x)

# Creating 2-D grid of features
[X, Y] = np.meshgrid(feature_x, feature_y)

fig, ax = plt.subplots(1, 1)
# plots filled contour plot
ax.contourf(X, Y, c_new)
  
ax.set_title('Filled Contour Plot')
ax.set_xlabel('feature_x')
ax.set_ylabel('feature_y')
  
plt.show()