import numpy as np
import matplotlib.pyplot  as plt

N=512
#a=100.0
#i = (i+int(n*np.sin(j*2*np.pi/N)))


def mixit(a):
    c = np.zeros((N,N))
    c_new = np.zeros((N,N))
    #introducing the sin wave
    for i in range(N):
        for j in range(N):
            c[i][j] = np.sin(i*np.pi*2/N)

    #vertical mixing
    for i in range(N):
        for j in range(N):
            c_new[i][j] = c[(i+int(a*np.sin(j*2*np.pi/N))) & N-1][j]
            
    c = c_new.copy()

        #horizontal mixing
    for i in range(N):
        for j in range(N):
            c_new[i][j] = c[i][(j+int(a*np.sin(i*2*np.pi/N))) & N-1]

    c = c_new.copy()
    return c

#print(mixit(100))
#plotting
xpoints = np.linspace(0,1,N)
ypoints = np.linspace(0,1,N)
#plt.plot(xpoints)

[X,Y]=np.meshgrid(xpoints, ypoints)

fig = plt.figure(figsize=(8, 8))

for i in np.arange(1,7):
    plt.subplot(3, 2, i)
    plt.contourf(X,Y,mixit(100*i))
    plt.xlabel('amp = %d'%(100*i))

#ax.set_title('Filled Contour Plot')
plt.suptitle('Different Line Plots')
# Auto adjust
plt.tight_layout()
#plt.sub_xlabel('xpoints')
#ax.set_ylabel('ypoints')
#ax.set_aspect('equal', adjustable='box')
plt.show()


#==========================================================================
def nparray_tail(x: np.array, n:int):
    """
    Returns tail N elements of array.
    :param x: Numpy array.
    :param n: N elements to return on end.
    :return: Last N elements of array.
    """
    if n == 0:
        return x[0:0]  # Corner case: x[-0:] will return the entire array but tail(0) should return an empty array.
    else:
        return x[-n:]  # Normal case: last N elements of array.
#print(nparray_tail(c,N))
