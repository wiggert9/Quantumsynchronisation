# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 12:55:35 2022

@author: wigge
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 12:38:12 2022

@author: wigge
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from toqito.channels import partial_trace
import cmath
import time

n=3
posdelt = np.linspace(0, 22, 2501)
negdelt = np.linspace(-22,0,2501)
posdeltfull = np.linspace(0, 25, 2501)
negdeltfull = np.linspace(-25,0,2501)


lst1 = [[np.sqrt(0.25)],[np.sqrt(0.75)]]
lst2 = [[np.sqrt(0.05)],[np.sqrt(0.95)]]

for i in range(n-2):
    lst1=lst1+[[0]]
    lst2= lst2+[[0]]

gamup =0.01
gamdown=10000

a1=tensor(tensor(tensor(destroy(n),qeye(n)),qeye(n)),qeye(n))
a2=tensor(tensor(tensor(qeye(n),destroy(n)),qeye(n)),qeye(n))
a3=tensor(tensor(tensor(qeye(n),qeye(n)),destroy(n)),qeye(n))
a4=tensor(tensor(tensor(qeye(n),qeye(n)),qeye(n)),destroy(n))

aa1=a1**2
aa2=a2**2
aa3=a3**2
aa4=a4**2

a1d=a1.dag()
a2d=a2.dag()
a3d=a3.dag()
a4d=a4.dag()

Domgsp = np.linspace(-25,25,25)
Vspplot = np.linspace(0,50,25)
Vsp = Vspplot*gamup
VClistlist =[]

c_opsstandard = [np.sqrt(gamup)*a1d,np.sqrt(gamup)*a2d , np.sqrt(gamup)*a3d,np.sqrt(gamup)*a4d, np.sqrt(gamdown)*aa1,np.sqrt(gamdown)*aa2,np.sqrt(gamdown)*aa3,np.sqrt(gamdown)*aa4]
omega1=8*np.pi
i=0

plt.plot()
plt.show()


###Chain
print("Chain")
i=0

VClistlist2=[]
VClistlist3=[]
VClistlist4=[]
VClistlist23=[]
VClistlist24=[]
VClistlist34=[]
for V in Vsp:
    print(i)
    i=i+1
    VClist2=[]
    VClist3=[]
    VClist4=[]
    VClist23=[]
    VClist24=[]
    VClist34=[]
    c_ops=c_opsstandard+[np.sqrt(V)*(a2-a1),np.sqrt(V)*(a3-a2),np.sqrt(V)*(a4-a3)]
    for Domg in Domgsp:
        
        ###Adjust which frequences you want, linear or constant
        
        omega2=omega1+1 *Domg*gamup
        #omega2=omega1+Domg*gamup
        omega3 =omega1+1* Domg*gamup
        #omega3=omega2+ Domg*gamup
        omega4=omega1+1 *Domg*gamup
        #omega4=omega3+Domg*gamup
        
        ####
        Hvdp = omega1*a1d*a1 + omega2*a2d*a2 +omega3*a3d*a3+omega4*a4d*a4
        firsttime=time.time()
        rho_ss = steadystate(Hvdp,c_ops)
        above2=(rho_ss*a1.dag()*a2).tr()
        below2=np.sqrt((rho_ss*a1.dag()*a1).tr() * (rho_ss*a2.dag()*a2).tr())
        C2=above2/below2
        VClist2.append(np.abs(C2))
        above3=(rho_ss*a1.dag()*a3).tr()
        below3=np.sqrt((rho_ss*a1.dag()*a1).tr() * (rho_ss*a3.dag()*a3).tr())
        C3=above3/below3
        VClist3.append(np.abs(C3))
        above4=(rho_ss*a1.dag()*a4).tr()
        below4=np.sqrt((rho_ss*a1.dag()*a1).tr() * (rho_ss*a4.dag()*a4).tr())
        C4=above4/below4
        VClist4.append(np.abs(C4))
        
        above23=(rho_ss*a2.dag()*a3).tr()
        below23=np.sqrt((rho_ss*a3.dag()*a3).tr() * (rho_ss*a2.dag()*a2).tr())
        C23=above23/below23
        VClist23.append(np.abs(C23))
        above24=(rho_ss*a2.dag()*a4).tr()
        below24=np.sqrt((rho_ss*a2.dag()*a2).tr() * (rho_ss*a4.dag()*a4).tr())
        C24=above24/below24
        VClist24.append(np.abs(C24))
        above34=(rho_ss*a3.dag()*a4).tr()
        below34=np.sqrt((rho_ss*a3.dag()*a3).tr() * (rho_ss*a4.dag()*a4).tr())
        C34=above34/below34
        VClist34.append(np.abs(C34))
    VClistlist2.append(VClist2)
    VClistlist3.append(VClist3)
    VClistlist4.append(VClist4)
    VClistlist23.append(VClist23)
    VClistlist24.append(VClist24)
    VClistlist34.append(VClist34)
    

Cvals2 = np.array(VClistlist2)
fig, ax = plt.subplots() 
cont0 = ax.contourf(Domgsp,Vspplot,Cvals2,levels=np.linspace(0,1,101),cmap=plt.cm.get_cmap('PuOr', 20))
ax.set_aspect(1)
ax.set_title("4 Chain -Quantum Arnold tongue"+ "$ |C_{\psi_{1,2}}|$");
ax.set_xlabel("$\Delta\omega_{i,1}$"+" (orders of "+r"$\gamma_{\uparrow}$"+")");
ax.set_ylabel("$V$"+" (orders of "+r"$\gamma_{\uparrow}$"+")")
cb =fig.colorbar(cont0)
def classiclines(k):
    global fig,ax
    posdelt = np.linspace(0, 50*k/4, 100)
    negdelt = np.linspace(-50*k/4,0,100)
    ax.plot(posdelt, 4/k *posdelt, 'k--');
    ax.plot(negdelt, -4/k *negdelt,'k--');

classiclines(0.5)
plt.show()

Cvals3 = np.array(VClistlist3)
fig, ax = plt.subplots() 
cont0 = ax.contourf(Domgsp,Vspplot,Cvals3,levels=np.linspace(0,1,101),cmap=plt.cm.get_cmap('PuOr', 20))
ax.set_aspect(1)
ax.set_title("4 Chain -Quantum Arnold tongue "+ "$ |C_{\psi_{1,3}}|$");
ax.set_xlabel("$\Delta\omega_{i,1}$"+" (orders of "+r"$\gamma_{\uparrow}$"+")");
ax.set_ylabel("$V$"+" (orders of "+r"$\gamma_{\uparrow}$"+")")
cb =fig.colorbar(cont0)
def classiclines(k):
    global fig,ax
    posdelt = np.linspace(0, 50*k/4, 100)
    negdelt = np.linspace(-50*k/4,0,100)
    ax.plot(posdelt, 4/k *posdelt, 'k--');
    ax.plot(negdelt, -4/k *negdelt,'k--');

classiclines(0.5)
plt.show()

Cvals4 = np.array(VClistlist4)
fig, ax = plt.subplots() 
cont0 = ax.contourf(Domgsp,Vspplot,Cvals4,levels=np.linspace(0,1,101),cmap=plt.cm.get_cmap('PuOr', 20))
ax.set_aspect(1)
ax.set_title("4 Chain -Quantum Arnold tongue "+ "$ |C_{\psi_{1,4}}|$");
ax.set_xlabel("$\Delta\omega_{i,1}$"+" (orders of "+r"$\gamma_{\uparrow}$"+")");
ax.set_ylabel("$V$"+" (orders of "+r"$\gamma_{\uparrow}$"+")")
cb =fig.colorbar(cont0)
def classiclines(k):
    global fig,ax
    posdelt = np.linspace(0, 50*k/4, 100)
    negdelt = np.linspace(-50*k/4,0,100)
    ax.plot(posdelt, 4/k *posdelt, 'k--');
    ax.plot(negdelt, -4/k *negdelt,'k--');

classiclines(0.5)
plt.show()

Cvals23 = np.array(VClistlist23)
fig, ax = plt.subplots() 
cont0 = ax.contourf(Domgsp,Vspplot,Cvals23,levels=np.linspace(0,1,101),cmap=plt.cm.get_cmap('PuOr', 20))
ax.set_aspect(1)
ax.set_title("4 Chain -Quantum Arnold tongue "+ "$ |C_{\psi_{2,3}}|$");
ax.set_xlabel("$\Delta\omega_{i,1}$"+" (orders of "+r"$\gamma_{\uparrow}$"+")");
ax.set_ylabel("$V$"+" (orders of "+r"$\gamma_{\uparrow}$"+")")
cb =fig.colorbar(cont0)
def classiclines(k):
    global fig,ax
    posdelt = np.linspace(0, 50*k/4, 100)
    negdelt = np.linspace(-50*k/4,0,100)
    ax.plot(posdelt, 4/k *posdelt, 'k--');
    ax.plot(negdelt, -4/k *negdelt,'k--');

classiclines(0.5)
plt.show()

Cvals24 = np.array(VClistlist24)
fig, ax = plt.subplots() 
cont0 = ax.contourf(Domgsp,Vspplot,Cvals24,levels=np.linspace(0,1,101),cmap=plt.cm.get_cmap('PuOr', 20))
ax.set_aspect(1)
ax.set_title("4 Chain -Quantum Arnold tongue "+ "$ |C_{\psi_{2,4}}|$");
ax.set_xlabel("$\Delta\omega_{i,1}$"+" (orders of "+r"$\gamma_{\uparrow}$"+")");
ax.set_ylabel("$V$"+" (orders of "+r"$\gamma_{\uparrow}$"+")")
cb =fig.colorbar(cont0)
def classiclines(k):
    global fig,ax
    posdelt = np.linspace(0, 50*k/4, 100)
    negdelt = np.linspace(-50*k/4,0,100)
    ax.plot(posdelt, 4/k *posdelt, 'k--');
    ax.plot(negdelt, -4/k *negdelt,'k--');

classiclines(0.5)
plt.show()

Cvals34 = np.array(VClistlist34)
fig, ax = plt.subplots() 
cont0 = ax.contourf(Domgsp,Vspplot,Cvals34,levels=np.linspace(0,1,101),cmap=plt.cm.get_cmap('PuOr', 20))
ax.set_aspect(1)
ax.set_title("4 Chain -Quantum Arnold tongue "+ "$ |C_{\psi_{3,4}}|$");
ax.set_xlabel("$\Delta\omega_{i,1}$"+" (orders of "+r"$\gamma_{\uparrow}$"+")");
ax.set_ylabel("$V$"+" (orders of "+r"$\gamma_{\uparrow}$"+")")
cb =fig.colorbar(cont0)
def classiclines(k):
    global fig,ax
    posdelt = np.linspace(0, 50*k/4, 100)
    negdelt = np.linspace(-50*k/4,0,100)
    ax.plot(posdelt, 4/k *posdelt, 'k--');
    ax.plot(negdelt, -4/k *negdelt,'k--');

classiclines(0.5)
plt.show()

