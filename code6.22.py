# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 09:10:41 2022

@author: davis
"""

import numpy as np
from numpy import transpose as tr
from numpy import linalg
import timescalecalculus as tsc
import random
import matplotlib.pyplot as plt

"""
Matrix Experiments
"""

#A = np.array([[1,2],[3,4]])
#B = np.array([[5,6],[7,8]])
#C = np.array([[2,4],[6,8]])
#Ainv = linalg.inv(A)

#print(A.dot(B.dot(C)))
#print(A.dot(Ainv))

#AT = np.matrix.transpose(A)
#print(AT)
#print(A.dot(B.dot(AT)))

"""
Messing with generating a random ts
"""
#tslist = [0]
#i = 0
#while i<31:
 #   mu = random.uniform(.1,1.6)
  #  tslist.insert(i+1, tslist[i]+mu)
   # i = i+1

#def g(t):
 #   return t**2

#ts = tsc.timescale(tslist)
#print(ts)

#ys = [g(x) for x in ts.ts] #!!! works with list, not timescale
#print(ys)

#plt.scatter(tslist,ys)

"""
Discrete Kalman Test
"""
#Defining Matrices
A = np.array([[1, 1],[0, 1]])
#tr test-------
B = tr(A)
print(B)
#--------------
G = np.array([[0],[1]])
Q = 1                 
R = 1
H = np.array([1,0])
P = np.array([[2,0],[0,3]])
x = np.array([[0],[10]])
z = np.array([9.5])
#Filter
iterator = 0
while iterator<100:
    xp = A.dot(x)
    Pp = A.dot(P.dot(tr(A))) + G.dot(np.dot(Q,tr(G)))
    P = linalg.inv((linalg.inv(Pp))+(tr(H).dot(np.dot(R,H))))
    x = xp + P.dot(tr(H).dot(np.dot(R,(z-(H.dot(xp))))))
    iterator = iterator + 1
    print(iterator)
print(z-(H.dot(xp)))
print("soup")
u="bannana"
v=1+5
print(v)