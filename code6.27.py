# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 11:05:50 2022

@author: davis
"""

import numpy as np
#import timescalecalculus as tsc
from numpy import transpose as tr
from numpy import array as mat
from numpy import linalg
import matplotlib.pyplot as plt
import scipy.stats as stats
import random

"""
Transcribing Mathematica
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
    return (1/mu)*mat([[.001,0],[0,.1]])
I = mat([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])

tslist = [0]
g = []
for i in range (0,999):
    mu = random.uniform(.3,.9)
    g.append(mu)
    tslist.append(tslist[i]+mu)
    i = i+1
#ts = tsc.timescale(tslist)

#Simulation
#How to make a matrix of matrices?
#Talk to Dr. Cuchta

#Initialize
P = I
Pminus = I
x = mat([[0],[0],[0],[0]])
xminus = mat([[0],[0],[0],[0]])
eta = mat([[0]]) #Residual

#Filter
for k in range (0,999):  
    mu = g[k]
    #Error Propogation
    Pminus = (I + (A(mu).dot(P))).dot(tr(I + (A(mu)))) + Covw(mu)
    P = linalg.inv(linalg.inv(Pminus) + (tr(H(mu)).dot(linalg.inv(R(mu)))).dot(H(mu)))
    print(P)
    #Estimation
    xminus = (I+A(mu)).dot(x)
    #In order to update err and x, the simulation must work
    #err = zvals - H.dot(xminus)
    #x = xminus + ((P.dot(tr(H))).dot(linalg.inv(R(mu)))).dot(err)
    #iterator
    k = k+1    
    
"""
Finishing Scalar LQR Example
"""
S = []
S.insert(100,5)
s = 5
K = []
k = 99

while k >= 0:
    Knew = (1.05*.01*s)/(1+(s*(.01)**2))
    K.insert(0,Knew)
    s = 1 + ((1*s*(1.05)**2)/(1+(s*(.01)**2)))
    S.insert(0,s)
    k = k-1
    
x = [10]
u = []
for i in range (0,99):
    u.append(-K[i]*x[i])
    x.append((1.05*x[i])+(.01*u[i]))
    
domain = []
for i in range (1,101):
    domain.append(i)
domain2 = []
for i in range (1,102):
    domain2.append(i)

#plt.scatter(domain,x)
plt.scatter(domain2,S)