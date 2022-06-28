# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 15:36:06 2022

@author: wigge
"""
import numpy as np

import matplotlib.pyplot as plt

def initialcircles(d):
    global fig,ax
    n=d+1
    thetaarr=np.linspace(0,2*np.pi,n)
    distance = np.sqrt(2-2*np.cos(2*np.pi/n))
    r=distance/4
    fig,ax = plt.subplots()
    for i in range(n):
        t=thetaarr[i]
        circle=plt.Circle((np.cos(t),np.sin(t)), r,fill=False)
        ax.add_patch(circle)
    ax.set_xlim(-1-2*r,1+2*r)
    ax.set_ylim(-1-2*r,1+2*r)
    
def circlecolour(d,angle,colour,anglecheck):
    global fig,ax
    
    n=d+1
    distance = np.sqrt(2-2*np.cos(2*np.pi/n))
    r=distance/4
    t=angle*(2*np.pi/(d))
    if anglecheck==True:
        thetacirc=plt.Circle((np.cos(t),np.sin(t)),r,color=colour,label="Aanstuurder")
    else:
        thetacirc=plt.Circle((np.cos(t),np.sin(t)),0.5*r,color=colour,label="Oscillator")
    ax.add_patch(thetacirc)
    plt.axis('off')
    plt.legend(loc=10)
    
    

def discrete():
    diff=6
    d=16
    k=2
    omega1=1
    omega2=2
    startingdiff=diff
    thetalst=[]
    philst=[]
    deltalst=[]
    theta=startingdiff
    phi=0
    tau = omega1-omega2
    synctime=0
    delta = theta-phi
    def G(k,delta):
        if abs(delta)<=k or abs(d-delta)<=k:
            return delta
        else:
            return 0
    check=False
    steps=0
    
    thetalst.append(theta)
    philst.append(phi)
    deltalst.append(theta-phi)
    while check <5 and steps <20:
        
        Gt = G(k,delta)

        print()
        print(steps)
        print(theta)
        print(phi)
        print(Gt)
        print(delta)

        theta = (theta+omega1)%d
        phi = (phi+omega2 + Gt)%d
        delta=(delta+omega1-omega2-Gt)%d
        
        
        thetalst.append(theta)
        philst.append(phi)
        if Gt==tau or Gt-d==tau:
            check=check+1
            if synctime==100:
                synctime=steps
        else:
            check = 0 
            synctime=100
        steps=steps+1
    print("Steps until synchronisation is" , synctime , " steps")
    
    thetaarr=np.array(thetalst)
    phiarr=np.array(philst)
    
    ####Initial circles:
    for i in range(len(thetaarr)):
        string=""
        if i>=synctime:
            string= "Syncronised"
        initialcircles(d)
        circlecolour(d,thetaarr[i],'b',True)
        circlecolour(d,phiarr[i],'r',False)
        plt.title("t="+str(i))
        plt.savefig( str(i) +'circle5.png')
        
discrete()
