# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 10:24:34 2022

@author: habee
"""
import numpy as np
from numpy import transpose as tr
from numpy import array as mat
from numpy import mod as mod
from numpy import linalg
import matplotlib.pyplot as plt
import random

#Matrices
mu=1
A= mat([[0, mu, 0, 0], [0, 0, 0, 0], [0, 0, 0, mu], [0, 0, 0, 0]])
H= mat([[1,0,0,0],[0,0,1,0]])

Covw= mat([[((mu**3)/3),((mu**2)/2),0,0], [((mu**2)/2),mu,0,0], [0,0,((mu**3)/3),((mu**2)/2)],[0,0,((mu**2)/2),mu]])
# R= (1/mu)(mat([[0.001,0],[0,0.1]]))


                                 #*Set Up Time Scale*

tsl=1000
#timescale limit
MU= []
TS= [0]

xvals=mat([[[0],[0],[0],[0]]])
yvals=mat([[[0],[0]]])


i=1

#For[i = 1, i <= tsl, i++, AppendTo[mu, If[Mod[i, 3] == 0, 1, If[Mod[i, 3] == 1, 2, 3]]]];
while(i<=tsl):
    if(i%3==0):
        MU.append(1)
    elif(i%3==1):
        MU.append(2)
    elif(i%3==2):
        MU.append(3)
    i+=1
    

#For[i = 1, i <= tsl, i++, AppendTo[TS, Last[TS] + mu[[i]]]];
i=1
while(i<=tsl): 
    TS.append(TS[-1]+MU[i-1])   
    i+=1
    
                # (*Simulate the "True" state x and observations y*)
                
#For[i = 1, i <= tsl, i++, AppendTo[xvals, (IdentityMatrix[4] + A[mu[[i]]]) . Last[xvals] +
#Transpose[{RandomReal[MultinormalDistribution[{0, 0, 0, 0}, Covw[mu[[i]]]]]}]]]
i=1
identityFour=np.zeros((4,4))
randomMultinormalDist=np.random.multivariate_normal(mat[0,0,0,0], Covw[i])
while(i<=tsl):
    xvals.append((identityFour+A.dot(MU[i]).xvals[-1])+tr[randomMultinormalDist])
    i+=1
#figure this out at some point    
    
print(xvals)
    
