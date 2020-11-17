#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
from math import log, sin, cos, pi

cplx = 0
N = 4096*2
fsamp = 44100
dt = 1/fsamp

ni = 200
omega = 2*pi*ni
ni2 = 5000
omega2 = 2*pi*ni2

l = np.ndarray(N)

for i in range(N):
    l[i] = 2*sin(i*dt*omega) + sin(i*dt*omega2)

plt.plot(l)
plt.show()

izlaz = "izlaz.dat"

with open(izlaz, "w") as f:
    f.write(f"{N} {fsamp} {cplx}")
    for i in l:
        f.write(f"\n{i.real} {i.imag}")

