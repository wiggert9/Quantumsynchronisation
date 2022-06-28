# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 14:38:57 2022

@author: wigge
"""


import numpy as np
import matplotlib.pyplot as plt
import time
from toqito.channels import partial_trace
d=40

omega1=5
omega2=2
K=5

###initialising theta and phi

theta=np.zeros(d)
theta[0]=1
phi=np.zeros(2)
phi[0]=1

zerostate=np.array([1,0])
onestate=np.array([0,1])
#### Initialising stimulus rotation matrix

stimrot=np.zeros((d,d))

for i in range(d):
    stimrot[(omega1+i)%d,i]=1
    
# Rotation matrix for qubit

    
oscrot = np.zeros((2,2))
oscrot[0,0]=np.cos(np.pi/d)
oscrot[0,1]=-np.sin(np.pi/d)
oscrot[1,0]=np.sin(np.pi/d)
oscrot[1,1]=np.cos(np.pi/d)


####


def Wk(theta,phi):
    thetastate = np.argwhere(theta)[0][0]    
    
    ### first considering state zero
    if np.absolute(thetastate-d/2 *0 )>K:
        first = np.kron(np.kron(theta,zerostate,zerostate))
    else:
        phi_new= np.linalg.matrix_power(oscrot,thetastate).dot(zerostate)
        first = np.kron(np.kron(theta,phi_new),onestate)
        
    ####Now considering state one
    if np.absolute(thetastate-d/2 *1 )>K:
        second = np.kron(np.kron(theta,onestate),zerostate)
    else:
        phi_new= np.linalg.matrix_power(oscrot,thetastate).dot(zerostate)
        second = np.kron(np.kron(theta,phi_new),onestate)
        
    state= phi[0]*first+phi[1]*second
    print(state)
    rho=np.outer(state,np.transpose(state))
    return rho
        
        



def Pi(k):
    xi=np.linalg.matrix_power(oscrot,k).dot(zerostate)
    return np.outer(xi,np.transpose(xi))

def transformation2(omegat,sigma,i):
    pt = sigma[0,0]
    ptconj=sigma[1,1]
    if 40<i<81:
        R = np.linalg.matrix_power(oscrot,omega2)
        sigmanew = np.matmul(R,np.matmul(sigma,np.transpose(R)))
    elif np.abs(omegat) <=K:
        sigmanew= pt*Pi(omegat+omega2)+ptconj*Pi(d//2 +omega2)
    elif np.abs(omegat-d/2) <= K:
        sigmanew=pt * Pi(omega2)+ ptconj*Pi(omegat+omega2)
    else:
        R = np.linalg.matrix_power(oscrot,omega2)
        sigmanew = np.matmul(R,np.matmul(sigma,np.transpose(R)))
        
    return sigmanew



def plot1():
    global n,d,omega1,omega2,K
    n=400
    omega1=3
    omega2= 2
    K=2
    pt=np.zeros(n)
    theta=np.zeros(d)
    theta[0]=1
    phi=np.array([1,0])
    sigma=np.outer(phi,np.transpose(phi))
    for i in range(n):
        sigma =transformation2(omega1*i%d,sigma,0)
        pt[i]=sigma[0,0]
        
    plt.plot(pt)
    plt.show()

def plot3():
    global n,d,omega1,omega2,K
    n=120
    omega1=5
    omega2= 2
    K=5
    pt=np.zeros(n)
    theta=np.zeros(d)
    theta[0]=1
    phi=np.array([1,0])
    sigma=np.outer(phi,np.transpose(phi))
    for i in range(n):
        
        sigma =transformation2(omega1*i%d,sigma,i)
        pt[i]=sigma[0,0]
        
    plt.plot(pt)
    plt.ylim(0,1)
    plt.axvline(x=40,color="r")
    plt.axvline(x=80,color='r')
    plt.show()
    
def plot4():
    global n,d,omega1,omega2,K
    n=100
    omega1=5
    omega2= 2
    K=3
    pt=np.zeros(n)
    theta=np.zeros(d)
    theta[0]=1
    phi=np.array([1,0])
    sigma=np.outer(phi,np.transpose(phi))
    for i in range(n):
        sigma =transformation2(omega1*i%d,sigma,0)
        pt[i]=sigma[0,0]
        
    plt.plot(pt)
    plt.ylim(0,1)
    plt.show()

plot1()
plot3()
plot4()
