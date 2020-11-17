#!/bin/python

import matplotlib.pyplot as plt
from math import log, sin, cos, pi
import numpy as np
from sys import argv

f1, f2 = argv[1:3]

print(f1, f2)

with open(f1) as f:
    ul = f.read().split("\n")
    N1, fsamp1, cplx1 = map(int, ul[0].split(" "))
    ul = ul[1:]
    l1 = np.ndarray(N1, dtype=complex) if cplx1 else np.ndarray(N1)
    for i in range(N1):
        r, im = map(float, ul[i].split(" "))
        l1[i] = complex(r, im) if cplx1 else r

dt1 = 1/fsamp1

print("\n\n", f1 ,N1, fsamp1)
print(l1)

with open(f2) as f:
    ul = f.read().split("\n")
    N2, fsamp2, cplx2 = map(int, ul[0].split(" "))
    ul = ul[1:]
    l2 = np.ndarray(N2, dtype=complex) if cplx2 else np.ndarray(N2)
    for i in range(N2):
        r, im = map(float, ul[i].split(" "))
        l2[i] = complex(r, im) if cplx2 else r

dt2 = 1/fsamp2

print("\n\n", f2 ,N2, fsamp2)
print(l2)

xf1 = np.linspace(0.0, 1.0/(2.0*dt1), N1)
xf2 = np.linspace(0.0, 1.0/(2.0*dt2), N2)

plt.plot(xf1, np.abs(l1[:N1]))
plt.plot(xf2, np.abs(l2[:N2]))
plt.show()
