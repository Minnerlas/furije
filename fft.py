#!/bin/python

from scipy.fft import fft, ifft
import matplotlib.pyplot as plt
from math import log, sin, cos, pi
import numpy as np

GRAFICI = 0

with open("izlaz.dat") as f:
    ul = f.read().split("\n")
    N, fsamp, cplx = map(int, ul[0].split(" "))
    ul = ul[1:]
    l = np.ndarray(N, dtype=complex) if cplx else np.ndarray(N)
    for i in range(N):
        r, im = map(float, ul[i].split(" "))
        l[i] = complex(r, im) if cplx else r

dt = 1/fsamp

print(N, fsamp)
print(l)

if GRAFICI:
    x = np.linspace(0.0, N*dt, N)
    plt.plot(x, l)
    plt.show()

lf = fft(l)
if GRAFICI:
    yf = 2/N * np.abs(lf[:N//2])
    #argf = np.angle(lf[:N//2])
    xf = np.linspace(0.0, 1.0/(2.0*dt), N//2)

    plt.plot(xf, yf)
    #plt.plot(xf, argf)
    plt.show()

with open("fft_py.dat", "w") as f:
    f.write(f"{N} {fsamp} 1")
    for i in lf:
        f.write(f"\n{i.real} {i.imag}")

