import numpy as np
import matplotlib.pyplot as plt
import math

N = 2**2 #pildi lahutus
a = 500 #transformatsiooni tugevus
amount = 0.7 #diffusion amount if constant weights are used, 0 = none, 1 = too much
c1 = np.zeros((N,N))

#three lines below are for initial sine
for i in range(N): #sama kood, mis loengus oli: algne siinus
        for j in range(N):
            c1[i][j] = np.sin(j*np.pi*2/N) #siinus i-teljel   
#diffuse(amount, c1, True)
print(c1)
c2 = np.zeros((N,N)) #tegin kõik transformatsioonid eraldi: ei leidnud esialgu mõtet eraldi fuktsiooni kirjutama hakata
for i in range(N): #esimene transformatsioon koosinusega, j-teljel
    for j in range(N):
        c2[i][j] = c1[i][(j+int(a*np.cos(i*np.pi*2/N))) & N-1] #(math.floor teeb täpselt sama asja mis int, kui pos arvud)
matrix  = c2.copy()
print(matrix)
for i in range(N): #diffuses from left to right i think
    for j in range(N):       
        #if in_x_direction: #if transformation was made horizontally
        w1 = a*np.cos(j*np.pi*2/N) - math.floor(a*np.cos(j*np.pi*2/N)) #let us do only cos
        #print(' eqv w1 = ',w1)
        w2 = 1-w1
        matrix[i][j] = w1*matrix[i][j-1] + w2*matrix[i][j]
        """else: #else => transformation was made vertically
            w1 = a*np.cos(i*np.pi*2/N) - math.floor(a*np.cos(i*np.pi*2/N))
            w2 = 1-w1
            matrix[i][j] = w1*matrix[i-1][j] + w2*matrix[i][j]"""
#print(matrix)        
#if in_x_direction:
matrix = matrix.T[::-1].T #reverses array
print(matrix)
#else:
#    matrix = matrix[::-1]

for i in range(N): #diffuses from right to left (?)
    for j in range(N):
        #============ 2nd loop!
        #if in_x_direction:
        w1 = a*np.cos(j*np.pi*2/N) - math.floor(a*np.cos(j*np.pi*2/N)) #let us do only cos
        #print(' eqv w1 v2= ',w1)
        w2 = 1-w1
        matrix[i][j] = w1*matrix[i][j-1] + w2*matrix[i][j]
        """else:
            w1 = a*np.cos(i*np.pi*2/N) - math.floor(a*np.cos(i*np.pi*2/N))
            w2 = 1-w1
            matrix[i][j] = w1*matrix[i-1][j] + w2*matrix[i][j]"""

print(matrix)

fig = plt.figure(figsize=(8,8))
colbar = plt.imshow(matrix, cmap='viridis', interpolation='none', extent=[0, N, 0, N], origin='lower', aspect='auto')
plt.clim(-1, 1) #fixes limits of color bar
#cmap variante: Greys, hot, terrain, rainbow, cividis, viridis. viimane loengu näitele kõige sarnasem
fig.colorbar(colbar)

#plt.show()
