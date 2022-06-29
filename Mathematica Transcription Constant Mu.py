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
    mu = 1
    g.append(mu)
    tslist.append(tslist[i]+mu)
    i = i+1
#ts = tsc.timescale(tslist)

#Simulation
xvals = [mat([[0],[0],[0],[0]])]
zvals = [mat([[0],[0]])]
for k in range (0,999):
    mu = g[k]
    #system sim
    w1 = int(stats.norm.rvs(loc=0,scale=linalg.norm(Covw(g[k])),size=1))
    w2 = int(stats.norm.rvs(loc=0,scale=linalg.norm(Covw(g[k])),size=1))
    w3 = int(stats.norm.rvs(loc=0,scale=linalg.norm(Covw(g[k])),size=1))
    w4 = int(stats.norm.rvs(loc=0,scale=linalg.norm(Covw(g[k])),size=1))
    w = mat([[w1],[w2],[w3],[w4]])
    xnew = ((I + A(mu)).dot(xvals[k])) + w
    xvals.append(xnew)
for k in range (0,998):
    mu = g[k]
    v1 = int(stats.norm.rvs(loc=0,scale=linalg.norm(R(g[k])),size=1))
    v2 = int(stats.norm.rvs(loc=0,scale=linalg.norm(R(g[k])),size=1))
    v = mat([[v1],[v2]])
    znew = H(mu).dot(xvals[k+1]) + tr(v)
    zvals.append(znew)

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
    #Estimation
    xminus = (I+A(mu)).dot(x)
    err = zvals[k] - H(mu).dot(xminus)
    x = xminus + ((P.dot(tr(H(mu)))).dot(linalg.inv(R(mu)))).dot(err)