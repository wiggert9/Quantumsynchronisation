

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


N=11
n=((N-5))//2
est = 1.2
err=0.005
leftend = np.linspace(-(1+err)*est,-(1-err)*est,n)
rightend=np.linspace((1-err)*est,(1+err)*est,n)
middlebit = np.linspace(-(1-err)*est,(1-err)*est,5)
#L=2
#a_arr=np.linspace(-L,L,N)
a_arr=np.concatenate((leftend,middlebit,rightend))
b_arr=np.copy(a_arr)
c_arr=np.copy(a_arr)
values=np.zeros((N,N,N))

def eq1(x,y,z,a):
    eq1= -2*np.sin(x)+np.sin(y)-np.sin(y+x)-np.sin(x+y+z) -a
    return eq1
def eq2(x,y,z,b):
    eq2= -2*np.sin(y) +np.sin(z)+np.sin(x)-np.sin(y+x) -b
    return eq2
def eq3(x,y,z,c):
    eq3= -2*np.sin(z)+np.sin(y)+np.sin(y+x)-np.sin(z+y+x)-c
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
                eq1= -2*np.sin(x)+np.sin(y)-np.sin(y+x)-np.sin(x+y+z) -a
                eq2= -2*np.sin(y) +np.sin(z)+np.sin(x)-np.sin(y+x) -b
                eq3= -2*np.sin(z)+np.sin(y)+np.sin(y+x)-np.sin(z+y+x)-c
                return [eq1,eq2,eq3]
            
            
            x,y,z=fsolve(equation,(0,0,0))
            
            if abs(eq1(x,y,z,a))>0.0001 or abs(eq2(x,y,z,b))>0.0001 or abs(eq3(x,y,z,c))>0.0001:
                values[i,j,k]=0
    
            else:
                values[i,j,k]=1
            

def findsquare(values):
    square=np.ones((N,N,N))
    for i in range(N):
        if values[i,i,i]>1:
            square[i,:,:]=0
            square[:,i,:]=0
            square[:,:,i]=0
            maxlength = abs(a_arr[i])
        if values[i,N-i-1,i]>1:
            square[i,:]=0
            square[:,N-i-1]=0
        
    #print(maxlength)
    return square


sq=findsquare(values)


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