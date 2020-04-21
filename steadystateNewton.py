#Steady state newton method
import matplotlib.pyplot as plt
import numpy as np
import random as rand
from cmath import *
from math import *

# NLS Newton 1D Method to solve:
# i u_t = -0.5u_xx +g(|u|^2)u
# We separate time dependence via u(x,t)-> u(x)*e^(i*mu*t)

############################

def con(a,b,axis):
    return np.concatenate((a,b),axis)

def sech(Y):
    return np.array([1/cosh(y) for y in Y])
    
N=201
g=-1
L=-10
R=10
x=np.linspace(L,R,num=N)
omega=0.1
V=0
dx=abs(x[1]-x[0])
mu=0.5
A=sqrt(2*mu)
pert=0.5
u0=A*sech(A*(x-(R+L)/2))
ONE=np.ones(N-1)
ONES=np.ones(N)
D2=np.diag(-2*ONES)+np.diag(ONE,-1)+np.diag(ONE,1)
D2[0,N-1]=1
D2[N-1,0]=1
D2=D2/(dx**2)

re=pert*(np.random.rand(N,1)-0.5)
up=u0+re #perturb IC
upi=np.zeros(N)
U=np.array(list(up)+list(upi))
it=0
err=1

while err>1e-5:
    it+=1
    Ur=U[0:N]
    Ui=U[N:2*N]
    J11=-0.5*D2 + np.diag(g*(3*(Ur**2)+Ui**2)+V+mu)
    J22=-0.5*D2+np.diag(g*(Ur**2 + 3*(Ui**2))+mu+V)
    J12=g*np.diag(2*Ur*Ui)
    J=con(con(J11,J12,1),con(J12,J22,1),0)
    U2=Ur**2 + Ui**2
    Fr=-0.5*np.dot(D2,Ur) +np.dot(g*U2 + V +mu,Ur) #same as below
    Fi=-0.5*np.dot(D2,Ui) +np.dot(g*U2 + V +mu,Ui) #.T missing after mu bracket
    F=con(Fr,Fi,0)
    DU=-np.dot(np.linalg.inv(J),F)
    U1=U+DU
    err=np.linalg.norm(U-U1)
    
    fig = plt.figure()
    plt.scatter(x,U[0:N],color="purple",marker=".")
    plt.scatter(x,U[N:2*N],color="green",marker=".")
    plt.plot(x,U1[0:N],color="r")
    plt.plot(x,U1[N:2*N],color="b")
    plt.legend(("Re current","Im current","Re previous","Im previous"))
    plt.title("iteration: "+str(it))
    plt.draw()
    plt.waitforbuttonpress(0)
    plt.close(fig)
    U=U1

u=[complex(U[i],U[N+i]) for i in range(0,N)] #solution to error confidence.


