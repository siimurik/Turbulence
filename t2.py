import math as m
import numpy as np

N = 4
a = 500.0
c1 = np.zeros((N,N))
c2  = np.zeros((N,N))
#three lines below are for initial sine
for i in range(N): #sama kood, mis loengus oli: algne siinus
        for j in range(N):
            c1[i][j] = np.sin(j*np.pi*2/N)

for i in range(N): #esimene transformatsioon koosinusega, j-teljel
    for j in range(N):
        l = int(a*np.cos(i*np.pi*2/N))
        #print("l =", l)
        c2[i][j] = c1[i][(j+l) & N-1]
        #print((j+l),'\t', N-1, '\t', ((j+l) & N-1))
#print(c2)
for i in range(N):    
    for j in range(N):
        w1 = a*np.cos(j*np.pi*2/N) - m.floor(a*np.cos(j*np.pi*2/N))
        w2 = 1-w1
        #print("w1 =", w1,"\tw2 =", w2)
        #print("w1_uus = ", w2*c2[i][j])
        c2[i][j] = w1*c2[i][j-1] + w2*c2[i][j]

print(c2)