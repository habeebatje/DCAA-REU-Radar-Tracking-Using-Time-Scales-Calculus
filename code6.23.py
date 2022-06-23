# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 09:53:38 2022

@author: davis
"""
#Kalman Filter Discrete (It works!)
import numpy as np
from numpy import transpose as tr
from numpy import array as mat
from numpy import linalg

A =mat([[1, 1], [0, 1]]);
H = mat([[1, 0]]);
G = mat([[0, 1]]);
B = mat([[.5, 1]]);
u = mat([[9.8,9.8]])
R = mat([[2]]);
Rinv = linalg.inv(R)
Q = mat([[1.333, .5],[.5, 1]])
x = mat([[5],[0]])
P = mat([[1, 0],[0, 1]])
iterator = 0;
#This 
while iterator < 100:
    #Covariance Side
    Pp = (A.dot(P)).dot(tr(A)) + (G.dot(Q)).dot(tr(G))
    Ppinv = linalg.inv(Pp)
    Pnoninv = Ppinv + ((tr(H)).dot(Rinv)).dot(H)
    P = linalg.inv(Pnoninv)
    #State Side
    iterator = iterator + 1
    
#Optimal Feedback Control

    
    