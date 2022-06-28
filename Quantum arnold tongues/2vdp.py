# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 11:25:31 2022

@author: wigge
"""

# -*- coding: utf-8 -*-
"""
Created on Fri May 13 17:12:05 2022

@author: wigge
"""

import numpy as np
import matplotlib.pyplot as plt
from toqito.channels import partial_trace
import cmath
from qutip import *
n=4
posdelt = np.linspace(0, 22, 2501)
negdelt = np.linspace(-22,0,2501)
posdeltfull = np.linspace(0, 25, 2501)
negdeltfull = np.linspace(-25,0,2501)

gamup =0.01
gamdown=10000
omega1=8*np.pi
a1=tensor(destroy(n),qeye(n))
a2=tensor(qeye(n),destroy(n))




Domgsp = np.linspace(-25,25,25)
Vspplot = np.linspace(0,50,25)
Vsp = Vspplot*gamup
VClistlist =[]
i=0
for V in Vsp:
    i=i+1
    VClist =[]
    c_ops = [np.sqrt(gamup)*a1.dag(),np.sqrt(gamup)*a2.dag() , np.sqrt(gamdown)*a1**2 , np.sqrt(gamdown)*a2**2, np.sqrt(V/2)*(a2-a1),np.sqrt(V/2)*(a1-a2)]
    for Domg in Domgsp:
        omega2=omega1+Domg*gamup
        Hvdp = omega1*a1.dag()*a1 + omega2 * a2.dag()*a2
        rho_ss = steadystate(Hvdp,c_ops)
        
        ####Enable these, and comment out the entropy if you want the complex correlator plots
        
        #above=(rho_ss*a1.dag()*a2).tr()
        #below=np.sqrt((rho_ss*a1.dag()*a1).tr() * (rho_ss*a2.dag()*a2).tr())
        #C=above/below
        #VClist.append(np.abs(C))
        
        ####
        rho_arr=np.array(rho_ss)
        diag= np.diagonal(rho_arr)
        rho_diag = Qobj(np.diag(diag))

        Slim=entropy_vn(rho_ss)
        Sdiag=entropy_vn(rho_diag)
        VClist.append(Sdiag-Slim)
    VClistlist.append(VClist)
Cvals = np.array(VClistlist)
fig, ax = plt.subplots() 
cont0 = ax.contourf(Domgsp,Vspplot,Cvals,levels=101,cmap=plt.cm.get_cmap('PuOr', 20))
ax.set_aspect(1)
ax.set_title("2 coupled oscillators - Relative entropy");
ax.set_xlabel("$\Delta\omega_{2,1}$"+" (orders of "+r"$\gamma_{\uparrow}$"+")");
ax.set_ylabel("$V$"+" (orders of "+r"$\gamma_{\uparrow}$"+")")
cb =fig.colorbar(cont0)


def classiclines(k):
    global fig,ax
    posdelt = np.linspace(0, 50*k/4, 100)
    negdelt = np.linspace(-50*k/4,0,100)
    ax.plot(posdelt, 4/k *posdelt, 'k--');
    ax.plot(negdelt, -4/k *negdelt,'k--');
    
    
classiclines(2)
plt.show()
