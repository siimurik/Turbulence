import numpy as np
# sampling rate
sr = 2000
# sampling interval
ts = 1.0/sr
print(ts)
t = np.arange(0,1,ts)

freq = 1.
x = 3*np.sin(2*np.pi*freq*t)
print(len(x))
print(np.fft.fft(x))
N = len(t)
N_half = int(N/2)
matrix = []
for i in range(N_half):
    for j in range(N_half):
        matrix[i][j] = x[N_half*i+j]
print(matrix)
freq = 4
x += np.sin(2*np.pi*freq*t)

freq = 7   
x += 0.5* np.sin(2*np.pi*freq*t)

