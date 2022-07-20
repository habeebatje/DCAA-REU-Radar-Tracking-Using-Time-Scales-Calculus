# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 15:40:19 2022

@author: davis
"""

import numpy as np
#import timescalecalculus as tsc
from numpy import transpose as tr
from numpy import array as mat
from numpy import linalg
import matplotlib.pyplot as plt
import scipy.stats as stats

"""
Constant Mu KF
"""
#Initializing Matrices and T.S.
#By utilizing functions we can define the matrices in terms of variables and
#use varying mu's
def A(mu):
    return mat([[0,mu,0,0],[0,0,0,0],[0,0,0,mu],[0,0,0,0]])
def H(mu):
    return mat([[1,0,0,0],[0,0,1,0]])
def Covw(mu):
    return mat([[((mu**3)/3),((mu**2)/2),0,0],[((mu**2)/2),mu,0,0],[0,0,((mu**3)/3),((mu**2)/2)],[0,0,((mu**2)/2),mu]])
def R(mu):
    return (1/mu)*mat([[.1,0],[0,1]])
I = mat([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])

#Desired Values
n=100

tslist = [0]
g = []
for i in range (0,n-1):
    mu = stats.expon.rvs(loc=1,scale=.3)
    g.append(mu)
    tslist.append(tslist[i]+mu)
    i = i+1
#ts = tsc.timescale(tslist)


#Simulation
xvals = [mat([[0],[0],[0],[0]])]
zvals = [mat([[0],[0]])]
rsim = [0]
tsim = [0]
zt = [0]
zr=[0]
for k in range (0,n-1):
    mu = g[k]
    #system sim
    w1 = stats.norm.rvs(loc=0,scale=linalg.norm(Covw(g[k])))
    w2 = stats.norm.rvs(loc=0,scale=linalg.norm(Covw(g[k])))
    w3 = stats.norm.rvs(loc=0,scale=linalg.norm(Covw(g[k])))
    w4 = stats.norm.rvs(loc=0,scale=linalg.norm(Covw(g[k])))
    w = mat([[w1],[w2],[w3],[w4]])
    xnew = ((I + A(mu)).dot(xvals[k])) + w
    rsim.append(xnew[0,0])
    tsim.append(xnew[2,0])
    xvals.append(xnew)
for k in range (0,n-1):
    mu = g[k]
    v1 = stats.norm.rvs(loc=0,scale=linalg.norm(R(g[k])))
    v2 = stats.norm.rvs(loc=0,scale=linalg.norm(R(g[k])))
    v = mat([[v1],[v2]])
    znew = H(mu).dot(xvals[k+1]) + v
    zr.append(znew[0,0])
    zt.append(znew[1,0])
    zvals.append(znew)

#Initialize
P = I
Pminus = I
x = mat([[0],[0],[0],[0]])
xminus = mat([[0],[0],[0],[0]])
etar = [0] #Residual
etat = [0]
radius = [x[0,0]]
theta = [x[2,0]]

#Filter
for k in range (0,n-1):  
    mu = g[k]
    #Error Propogation
    Pminus = (I + (A(mu).dot(P))).dot(tr(I + (A(mu)))) + Covw(mu)
    P = linalg.inv(linalg.inv(Pminus) + (tr(H(mu)).dot(linalg.inv(R(mu)))).dot(H(mu)))
    #Estimation
    xminus = (I+A(mu)).dot(x)
    err = zvals[k] - H(mu).dot(xminus)
    x = xminus + ((P.dot(tr(H(mu)))).dot(linalg.inv(R(mu)))).dot(err)
    r = x[0,0]
    radius.append(r)
    th = x[2,0]
    theta.append(th)
    etar.append(err[0,0])
    etat.append(err[1,0])
    
#Plotting
plt.figure(1)
plt.plot(tslist,radius,'b',label='Estimation')
plt.plot(tslist,rsim,'r',label='True')
plt.plot(tslist,zr,'g',label='Measured')
#Blue is Estimation, Red is Sim
plt.xlabel('iterations')
plt.ylabel('radius values')
plt.title('Radius Actual v. Estimation')
plt.grid(True)
plt.legend()
plt.show()

plt.figure(2)
plt.plot(tslist,theta,'b',label='Estimation')
plt.plot(tslist,tsim,'r',label='True')
plt.plot(tslist,zt,'g',label='Measured')
#Blue is Estimation, Red is Sim
plt.xlabel('iterations')
plt.ylabel('angle values')
plt.title('Angle Actual v. Estimation')
plt.grid(True)
plt.legend()
plt.show()

plt.figure(3)
plt.plot(tslist,etar,'y-')
plt.xlabel('iterations')
plt.ylabel('Error')
plt.title('Error for Radius')
plt.grid(True)
plt.show()

plt.figure(4)
plt.plot(tslist,etat,'g-')
plt.xlabel('iterations')
plt.ylabel('Error')
plt.title('Error for Angle')
plt.grid(True)
plt.show()
    