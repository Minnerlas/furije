#!/bin/python

from scipy.fft import fft, ifft
import matplotlib.pyplot as plt
from math import log, sin, cos, pi
import numpy as np

N = 65536
fsamp = 44100
dt = 1/fsamp

ni = 200
omega = 2*pi*ni
ni2 = 5000
omega2 = 2*pi*ni2

l = np.ndarray(N)

for i in range(N):
    l[i] = 2*sin(i*dt*omega) + sin(i*dt*omega2)

print(l)

plt.plot(l)
plt.show()


y = np.absolute(fft(l))
x = np.linspace(0.0, 1.0/(2.0*dt), N//2)

print(y)



plt.plot(x, 2/N * y[:N//2])
plt.show()
