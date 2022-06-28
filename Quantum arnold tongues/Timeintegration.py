# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 19:59:02 2022

@author: wigge
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 13:58:43 2022

@author: wigge
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 10:40:36 2022

@author: wigge
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May 31 12:25:47 2022

@author: wigge
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May 31 11:30:58 2022

@author: wigge
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from toqito.channels import partial_trace
import cmath
import time

n=3

    
lst1 = [[np.sqrt(1)],[np.sqrt(0)]]
lst2 = [[np.sqrt(0)],[np.sqrt(1)]]
lst3=  [[np.sqrt(1)],[np.sqrt(0)]]
lst4= [[np.sqrt(0)],[np.sqrt(1)]]
for i in range(n-2):
    lst1=lst1+[[0]]
    lst2= lst2+[[0]]
    lst3=lst3+[[0]]
    lst4=lst4+[[0]]

q1,q2,q3,q4=Qobj(lst1),Qobj(lst2),Qobj(lst3),Qobj(lst4)
psi=tensor(tensor(tensor(q1,q2),q3),q4)
rho=psi*psi.dag()

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

V=45*gamup
c_opsstandard = [np.sqrt(gamup)*a1d,np.sqrt(gamup)*a2d , np.sqrt(gamup)*a3d,np.sqrt(gamup)*a4d, np.sqrt(gamdown)*aa1,np.sqrt(gamdown)*aa2,np.sqrt(gamdown)*aa3,np.sqrt(gamdown)*aa4]
omega1=2*np.pi
c_ops=c_opsstandard+[np.sqrt(V)*(a2-a1),np.sqrt(V)*(a3-a2),np.sqrt(V)*(a4-a3)]
omega2=omega1+1*gamup
omega3 =omega2+1*gamup
omega4=omega3+1*gamup


tlst=np.linspace(0,100,1000)
Hvdp = omega1*a1d*a1 + omega2*a2d*a2 +omega3*a3d*a3+omega4*a4d*a4

firsttime=time.time()
result3=mesolve(Hvdp,rho,tlst,c_ops,options=Options(nsteps=100000))
N=600

Cabslst2=np.zeros(N)
Canglst2=np.zeros(N)

Cabslst3=np.zeros(N)
Canglst3=np.zeros(N)

Cabslst4=np.zeros(N)
Canglst4=np.zeros(N)

Cabslst23=np.zeros(N)
Canglst23=np.zeros(N)

Cabslst24=np.zeros(N)
Canglst24=np.zeros(N)

Cabslst34=np.zeros(N)
Canglst34=np.zeros(N)

entropylst=np.zeros(N)

tlist=np.linspace(0,60,N)
for i in range(N):
    #print(i)
    rho_ss=result3.states[i]
    above2=(rho_ss*a1.dag()*a2).tr()
    below2=np.sqrt((rho_ss*a1.dag()*a1).tr() * (rho_ss*a2.dag()*a2).tr())
    C2=above2/below2
    
    above3=(rho_ss*a1.dag()*a3).tr()
    below3=np.sqrt((rho_ss*a1.dag()*a1).tr() * (rho_ss*a3.dag()*a3).tr())
    C3=above3/below3

    above4=(rho_ss*a1.dag()*a4).tr()
    below4=np.sqrt((rho_ss*a1.dag()*a1).tr() * (rho_ss*a4.dag()*a4).tr())
    C4=above4/below4
    
    
    above23=(rho_ss*a2.dag()*a3).tr()
    below23=np.sqrt((rho_ss*a3.dag()*a3).tr() * (rho_ss*a2.dag()*a2).tr())
    C23=above23/below23
   
    above24=(rho_ss*a2.dag()*a4).tr()
    below24=np.sqrt((rho_ss*a2.dag()*a2).tr() * (rho_ss*a4.dag()*a4).tr())
    C24=above24/below24
    
    above34=(rho_ss*a3.dag()*a4).tr()
    below34=np.sqrt((rho_ss*a3.dag()*a3).tr() * (rho_ss*a4.dag()*a4).tr())
    C34=above34/below34
     
    Cabslst2[i]=np.abs(C2)
    Canglst2[i]=np.angle(C2)
    
    Cabslst3[i]=np.abs(C3)
    Canglst3[i]=np.angle(C3)
    
    Cabslst4[i]=np.abs(C4)
    Canglst4[i]=np.angle(C4)
    
    Cabslst23[i]=np.abs(C23)
    Canglst23[i]=np.angle(C23)
    
    Cabslst24[i]=np.abs(C24)
    Canglst24[i]=np.angle(C24)
    
    Cabslst34[i]=np.abs(C34)
    Canglst34[i]=np.angle(C34)

    rho_arr=np.array(rho_ss)
    diag= np.diagonal(rho_arr)
    rho_diag = Qobj(np.diag(diag))

    Slim=entropy_vn(rho_ss)
    Sdiag=entropy_vn(rho_diag)
    entropylst[i] = Sdiag-Slim
plt.plot(tlist,Cabslst2,label=r'$1,2$')
plt.plot(tlist,Cabslst3,label=r'$1,3$')
plt.plot(tlist,Cabslst4,label=r'$1,4$')
plt.plot(tlist,Cabslst23,label=r'$2,3$')
plt.plot(tlist,Cabslst24,label=r'$2,4$')
plt.plot(tlist,Cabslst34,label=r'$3,4$')
plt.legend(ncol=2)
plt.xlabel("Time(s)")
plt.ylabel("$|C_{i,j}|$")
plt.title("Correlation strength")
plt.show()

plt.plot(tlist,Canglst2,label=r'1,2')
plt.plot(tlist,Canglst3,label="1,3")
plt.plot(tlist,Canglst4,label="1,4")
plt.plot(tlist,Canglst23,label="2,3")
plt.plot(tlist,Canglst24,label="2,4")
plt.plot(tlist,Canglst34,label="3,4")
plt.xlabel("Time(s)")
plt.ylabel("$\Delta \phi _{i,j}$")
plt.legend(ncol=2)
plt.ylim(-0.5,0.5)
plt.title("Phase difference")
plt.show()



plt.plot(tlist,entropylst)
plt.xlabel("Time(s)")
plt.ylabel(r"$S(\rho_{diag})-S(\rho)$")
plt.title("Relative entropy")