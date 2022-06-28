# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 16:18:55 2022

@author: wigge
"""
from qutip import *
import numpy as np
import matplotlib.pyplot as plt




n=3
posdelt = np.linspace(0, 22, 2501)
negdelt = np.linspace(-22,0,2501)
posdeltfull = np.linspace(0, 25, 2501)
negdeltfull = np.linspace(-25,0,2501)

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
def classiclines(k):
    global fig,ax
    posdelt = np.linspace(0, 50*k/4, 100)
    negdelt = np.linspace(-50*k/4,0,100)
    ax.plot(posdelt, 4/k *posdelt, 'k--');
    ax.plot(negdelt, -4/k *negdelt,'k--');


def extra_c_ops(system):
    if system == "Loop":
        return [np.sqrt(V)*(a2-a1),np.sqrt(V)*(a3-a2),np.sqrt(V)*(a4-a3),np.sqrt(V)*(a1-a4)]
    elif system == "Tree":
        return [np.sqrt(V)*(a2-a1),np.sqrt(V)*(a3-a1),np.sqrt(V)*(a4-a1)]
    elif system == "All":
        return [np.sqrt(V)*(a2-a1),np.sqrt(V)*(a3-a2),np.sqrt(V)*(a4-a3),np.sqrt(V)*(a1-a4),np.sqrt(V)*(a1-a3),np.sqrt(V)*(a2-a4)]
    elif system == "Chain":
        return [np.sqrt(V)*(a2-a1),np.sqrt(V)*(a3-a2),np.sqrt(V)*(a4-a3)]
    elif system == "Tower":
        return [np.sqrt(V)*(a2-a1),np.sqrt(V)*(a3-a1),np.sqrt(V)*(a4-a1),np.sqrt(V)*(a2-a3)]
    elif system == "Spade":
        return [np.sqrt(V)*(a2-a1),np.sqrt(V)*(a3-a2),np.sqrt(V)*(a4-a3),np.sqrt(V)*(a1-a4),np.sqrt(V)*(a1-a3)]
        
def classicalat(system):
    if system == "Loop":
        return 0.9
    elif system == "Tree":
        return 0.8
    elif system == "All":
        return 1.68
    elif system == "Chain":
        return 0.5
    elif system == "Tower":
        return 0.8
    elif system == "Spade":
        return 1.2

def qat(system,entropy):
    global V,fig,ax
    print(system)
    VClistlist =[]
    i=0
    for V in Vsp:
        print(i)
        i=i+1
        VClist =[]
        c_ops=c_opsstandard+extra_c_ops(system)
        for Domg in Domgsp:
            omega2=omega1+Domg*gamup
            omega3 =omega1+Domg*gamup
            omega4=omega1+Domg*gamup
            Hvdp = omega1*a1d*a1 + omega2*a2d*a2 +omega3*a3d*a3+omega4*a4d*a4
            firsttime=time.time()
            rho_ss = steadystate(Hvdp,c_ops)
            if entropy ==False:
                
                above=(rho_ss*a1.dag()*a2).tr()
                below=np.sqrt((rho_ss*a1.dag()*a1).tr() * (rho_ss*a2.dag()*a2).tr())
                C=above/below
                VClist.append(np.abs(C))
            else:
                rho_arr=np.array(rho_ss)
                diag= np.diagonal(rho_arr)
                rho_diag = Qobj(np.diag(diag))
                Slim=entropy_vn(rho_ss)
                Sdiag=entropy_vn(rho_diag)
                VClist.append(Sdiag-Slim)
        VClistlist.append(VClist)
    Cvals = np.array(VClistlist)
    fig, ax = plt.subplots() 
    if entropy ==False:
        top =1
    else:
        top=0.22
    cont0 = ax.contourf(Domgsp,Vspplot,Cvals,levels=np.linspace(0,top,101),cmap=plt.cm.get_cmap('PuOr', 20))
    ax.set_aspect(1)
    if entropy == False:
        ax.set_title("4" +system + "-Quantum Arnold tongue");
    else:
        ax.set_title("4" +system + "- Relative entropy")
    ax.set_xlabel("$\Delta\omega_{2,1}$"+" (orders of "+r"$\gamma_{\uparrow}$"+")");
    ax.set_ylabel("$V$"+" (orders of "+r"$\gamma_{\uparrow}$"+")")
    cb =fig.colorbar(cont0)
    classiclines(classicalat(system))
    plt.show()


qat("Chain",False)
qat("Tree",False)
qat("Loop",False)
qat("Tower",False)
qat("Spade",False)
qat("All",False)