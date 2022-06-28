# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 16:09:46 2022

@author: wigge
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 15:57:15 2022

@author: wigge
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


N=11
n=((N-5))//2
est = 0.8
err=0.01
leftend = np.linspace(-(1+err)*est,-(1-err)*est,n)
rightend=np.linspace((1-err)*est,(1+err)*est,n)
middlebit = np.linspace(-(1-err)*est,(1-err)*est,5)
#a_arr=np.linspace(-L,L,N)
a_arr=np.concatenate((leftend,middlebit,rightend))
b_arr=np.copy(a_arr)
c_arr=np.copy(a_arr)
values=np.zeros((N,N,N))

def eq1(x,y,z,a):
    eq1=-2*np.sin(x)+np.sin(y-x)-np.sin(z)-np.sin(y)   -a
    return eq1
def eq2(x,y,z,b):
    eq2=-2*np.sin(y)-np.sin(y-x)-np.sin(x)-np.sin(z)-b
    return eq2
def eq3(x,y,z,c):
    eq3=-2*np.sin(z)-np.sin(y)-np.sin(x) -c
    return eq3
for i in range(N):
    print(i)
    for j in range(N):
        for k in range(N):
            a=a_arr[i]
            b=b_arr[j]
            c=c_arr[k]
            def equation(vars):
                x,y,z=vars
                eq1=-2*np.sin(x)+np.sin(y-x)-np.sin(z)-np.sin(y)   -a
                eq2=-2*np.sin(y)-np.sin(y-x)-np.sin(x)-np.sin(z)-b
                eq3=-2*np.sin(z)-np.sin(y)-np.sin(x) -c
                return [eq1,eq2,eq3]
            
            
            x,y,z=fsolve(equation,(0,0,0))
            
            if abs(eq1(x,y,z,a))>0.0001 or abs(eq2(x,y,z,b))>0.0001 or abs(eq3(x,y,z,c))>0.0001:
                values[i,j,k]=0
    
            else:
                values[i,j,k]=1
            


def pv(i):
    k=(N-1)//2
    
    #print(values[k+i,k+i,k+i])
    #print(values[k+i,k-i,k+i])
    #print(values[k-i,k-i,k+i])
    #print()
    
    return (values[k+i,k+i,k+i],values[k+i,k-i,k+i],values[k-i,k-i,k+i])
    
    
    
    
def findmaxlength(values):
    i=0
    found=False
    while found==False:
        s,r,t = pv(i)
        if s==0 or r==0 or t==0:
            print(a_arr[(N-1)//2 +i]/2)
            print(a_arr[(N-1)//2 +i-1]/2)
            found=True
        i=i+1
            
findmaxlength(values)