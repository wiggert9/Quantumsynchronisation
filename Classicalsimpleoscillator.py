# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 15:31:56 2022

@author: wigge
"""

import numpy as np
import matplotlib.pyplot as plt

def classic():
    theta=0
    phi=np.pi
    dt=0.1
    thetalst=[]
    philst=[]
    r1=1
    r2=0.5
    loops=3
    d=3
    omega1=1
    omega2=2
    time=0
    timelst=[]
    while time<10:
        thetalst.append(theta)
        philst.append(phi)
        theta=(theta+dt*omega1 )%(2*np.pi)
        phi=(phi+dt*(omega2 +np.sin((theta-phi)))) %(2*np.pi)
        timelst.append(time)
        time=time+dt
    for j in range(len(thetalst)):
        fig,ax = plt.subplots()
        t1=thetalst[j]
        thetacirc=plt.Circle((d*np.cos(t1),d*np.sin(t1)),r1,color='b',label="Stimulus")
        t2=philst[j]
        phicirc=plt.Circle((d*np.cos(t2),d*np.sin(t2)),r2,color='r',label="Oscillator")
        ax.add_patch(thetacirc)
        ax.add_patch(phicirc)
        plt.axis('off')
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.legend(loc=10)
        plt.title("t="+str(round(timelst[j],3)))
        plt.savefig( str(j) +'circle4.png')      
        
classic()

