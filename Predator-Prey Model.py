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
a = .4
c = .7
alpha = .5
gamma = .25
preyinitial = 600
predatorinitial = 100
def A(mu):
    return mat([[(-.5*a*c*(mu**2)), 
                 (-alpha*c*mu/gamma)],
                [(gamma*a*mu/alpha), 
                 (-.5*a*c*(mu**2))]])
def H(mu):
    return mat([[1,0],[0,1]])
def q1(mu):
    return ((alpha**2)*mu)+((alpha**2)*(c**2)*(mu**3)/3)
def q2(mu):
    return (alpha*gamma*(a-c)*(mu**2)/2)
def q3(mu):
    return (alpha*gamma*(a-c)*(mu**2)/2)
def q4(mu):
    return ((gamma**2)*mu)+((gamma**2)*(a**2)*(mu**3)/3)
def Covw(mu):
    return mat([[q1(mu), q2(mu)],[q3(mu), q4(mu)]])
def R(mu):
    return (1/mu)*mat([[.1,0],[0,1]])
I = mat([[1,0],[0,1]])
print(Covw(.5))
#Desired Values
n=300

tslist = [0]
g = []
for i in range (0,n-1):
    mu = stats.norm.rvs(loc=.5,scale=.05)
    g.append(mu)
    tslist.append(tslist[i]+mu)
    i = i+1
#ts = tsc.timescale(tslist)


#Simulation
xvals = [mat([[preyinitial],[predatorinitial]])]
zvals = [mat([[preyinitial],[predatorinitial]])]
preysim = [preyinitial+800]
predatorsim = [predatorinitial+500]
zpredator = [predatorinitial+500]
zprey=[preyinitial+800]

for k in range (0,n-1):
    mu = g[k]
    #system sim
    w1 = stats.norm.rvs(loc=0,scale=linalg.norm(Covw(g[k])))
    w2 = stats.norm.rvs(loc=0,scale=linalg.norm(Covw(g[k])))
    w = 5*mat([[w1],[w2]])
    xnew = ((I + A(mu)).dot(xvals[k])) + 10*(mat([[-alpha,0],[0,gamma]])).dot(w)
    preysim.append(xnew[0,0]+800)
    predatorsim.append(xnew[1,0]+500)
    xvals.append(xnew)
for k in range (0,n-1):
    mu = g[k]
    v1 = stats.norm.rvs(loc=0,scale=linalg.norm(R(g[k])))
    v2 = stats.norm.rvs(loc=0,scale=linalg.norm(R(g[k])))
    v = mat([[v1],[v2]])
    znew = H(mu).dot(xvals[k+1]) + v
    zprey.append(znew[0,0])
    zpredator.append(znew[1,0])
    zvals.append(znew)

#Initialize
P = I
Pminus = I
x = mat([[preyinitial],[predatorinitial]])
xminus = mat([[preyinitial],[predatorinitial]])
etaprey = [0] #Residual
etapredator = [0]
prey = [x[0,0]+800]
predator = [x[1,0]+500]

#Filter
for k in range (0,n-1):  
    mu = g[k]
    #Error Propogation
    Pminus = (I + ((A(mu)).dot(P))).dot(tr(I + A(mu)) + Covw(mu))
    P = linalg.inv(linalg.inv(Pminus) + (tr(H(mu)).dot(linalg.inv(R(mu)))).dot(H(mu)))
    #Estimation
    xminus = (I+(A(mu))).dot(x)
    err = zvals[k] - H(mu).dot(xminus)
    x = xminus + ((P.dot(tr(H(mu)))).dot(linalg.inv(R(mu)))).dot(err)
    s = x[0,0]
    prey.append(s+800)
    i = x[1,0]
    predator.append(i+500)
    etaprey.append(err[0,0])
    etapredator.append(err[1,0])
    
#Plotting
plt.figure(1)
plt.plot(tslist,prey,'b',label='Estimation')
plt.plot(tslist,preysim,'r',label='True')
#plt.plot(tslist,zr,'g',label='Measured')
#Blue is Estimation, Red is Sim
plt.xlabel('Time')
plt.ylabel('Population X')
plt.title('True Prey v. Estimation')
plt.grid(True)
plt.legend()
plt.show()

plt.figure(2)
plt.plot(tslist,predator,'b',label='Estimation')
plt.plot(tslist,predatorsim,'r',label='True')
#plt.plot(tslist,zt,'g',label='Measured')
#Blue is Estimation, Red is Sim
plt.xlabel('Time')
plt.ylabel('Population Y')
plt.title('True Preadtors v. Estimation')
plt.grid(True)
plt.legend()
plt.show()

plt.figure(3)
plt.plot(tslist,etaprey,'y-')
plt.xlabel('Time')
plt.ylabel('Error')
plt.title('Error for Prey')
plt.grid(True)
plt.show()

plt.figure(4)
plt.plot(tslist,etapredator,'g-')
plt.xlabel('Time')
plt.ylabel('Error')
plt.title('Error for Predators')
plt.grid(True)
plt.show()

plt.figure(5)
plt.plot(tslist,prey,'b',label='Prey')
plt.plot(tslist,predator,'g',label='Predator')
plt.xlabel('Time')
plt.ylabel('Populations')
plt.title('Predator and Prey')
plt.grid(True)
plt.legend()
plt.show()
    
    