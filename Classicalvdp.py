# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 14:35:07 2022

@author: wigge
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint



mu=1
omega =1
beta=1

def dU_dx(U, t):
    # Here U is a vector such that y=U[0] and z=U[1]. This function should return [y', z']
    return [U[1], mu*(1-beta*U[0]**2)*U[1] - omega**2 *U[0]]

U0 = [0.01, 0]

T=40
h=0.1
t = np.linspace(0, T,int(T/h))
Us = odeint(dU_dx, U0, t)
xs=Us[:,0]
ys = Us[:,1]


plt.xlabel("x")
plt.ylabel("y")
plt.title("Van der Pol oscillator- Phase diagram")
plt.plot(xs,ys)
plt.xlim(-3,3)
plt.ylim(-3,3)
plt.show()

plt.plot(t,xs)
plt.xlabel("t")
plt.ylabel("x")
plt.title("Van der Pol oscillator")