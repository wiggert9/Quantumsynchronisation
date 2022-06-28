# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 17:25:08 2022

@author: wigge
"""

import matplotlib.pyplot as plt
import numpy as np
fig, ax = plt.subplots() 

ax.set_title("");
ax.set_xlabel("$\Delta\omega_{i,j}$")
ax.set_ylabel("$V$")


'''
p=1
k=2*p
posdelt = np.linspace(0, 50*k/4, 100)
negdelt = np.linspace(-50*k/4,0,100)
ax.plot(posdelt, 4/k *posdelt, 'r');
#ax.plot(negdelt, -4/k *negdelt,'r');

p=0.88
k=2*p
posdelt = np.linspace(0, 50*k/4, 100)
negdelt = np.linspace(-50*k/4,0,100)
ax.plot(posdelt, 4/k *posdelt, 'b');
ax.plot(negdelt, -4/k *negdelt,'b');

p=0.5
k=2*p
posdelt = np.linspace(0, 50*k/4, 100)
negdelt = np.linspace(-50*k/4,0,100)
ax.plot(posdelt, 4/k *posdelt, 'g');
ax.plot(negdelt, -4/k *negdelt,'g');
'''

p=0.84
k=2*p
posdelt = np.linspace(0, 50*k/4, 100)
negdelt = np.linspace(-50*k/4,0,100)
ax.plot(posdelt, 4/k *posdelt, 'r',label='All-to-all');
ax.plot(negdelt, -4/k *negdelt,'r');
p=0.6
k=2*p
posdelt = np.linspace(0, 50*k/4, 100)
negdelt = np.linspace(-50*k/4,0,100)
ax.plot(posdelt, 4/k *posdelt, 'b',label='Spade');
ax.plot(negdelt, -4/k *negdelt,'b');

p=0.45
k=2*p
posdelt = np.linspace(0, 50*k/4, 100)
negdelt = np.linspace(-50*k/4,0,100)
ax.plot(posdelt, 4/k *posdelt, 'g',label='Loop');
ax.plot(negdelt, -4/k *negdelt,'g');
p=0.4
k=2*p
posdelt = np.linspace(0, 50*k/4, 100)
negdelt = np.linspace(-50*k/4,0,100)
ax.plot(posdelt, 4/k *posdelt, 'pink',label='Tree and Tower');
ax.plot(negdelt, -4/k *negdelt,'pink');


p=0.25
k=2*p
posdelt = np.linspace(0, 50*k/4, 100)
negdelt = np.linspace(-50*k/4,0,100)
ax.plot(posdelt, 4/k *posdelt, 'grey',label='Chain');
ax.plot(negdelt, -4/k *negdelt,'grey');

plt.title('Synchronisation regime')
plt.legend()
plt.show()