"""
Created on Thu Jun 23 09:53:38 2022

@author: davis
"""

import numpy as np
from numpy import transpose as tr
from numpy import array as mat
from numpy import linalg
import matplotlib.pyplot as plt
import random

"""
Kalman Filter Discrete, Newtonian Example, Error Covariance Propogation
"""

A = mat([[1, 1], [0, 1]]);
H = mat([[1, 0]]);
G = mat([[0, 1]]);
B = mat([[.5, 1]]);
u = mat([[9.8]])
R = mat([[2]]);
Rinv = linalg.inv(R)
Q = mat([[1.333, .5],[.5, 1]])
x = mat([[5],[0]])
P = mat([[1, 0],[0, 1]])
iterator = 0;
#This 
while iterator < 10:
    #Covariance Side
    Pp = (A.dot(P)).dot(tr(A)) + (G.dot(Q)).dot(tr(G))
    Ppinv = linalg.inv(Pp)
    Pnoninv = Ppinv + ((tr(H)).dot(Rinv)).dot(H)
    P = linalg.inv(Pnoninv)
    print(P)
    #State Side
    iterator = iterator + 1
    
"""
Plotting
"""
#Defining Matrices
A = np.array([[1, 1],[0, 1]])
G = np.array([[0],[1]])
Q = np.array([[1]])
R = np.array([[1]])
Rinv = linalg.inv(R)
H = np.array([[1,0]])
P = np.array([[2,0],[0,3]])
x = np.array([[0],[10]])
z = 9.5 + random.uniform(0,1)
xs = [0]
ds = [0]
vs = [0]
#Filter
i = 0
n = 30        
while i<n:
    #Covariance Generator
    Pp = (A.dot(P)).dot(tr(A)) + (G.dot(Q)).dot(tr(G))
    Ppinv = linalg.inv(Pp)
    Pnoninv = Ppinv + ((tr(H)).dot(Rinv)).dot(H)
    P = linalg.inv(Pnoninv)
    #Measurement Update
    xp = A.dot(x)
    zr = z - H.dot(xp) #You can add float or int with a 1x1 array
    x = xp + ((P.dot(tr(H))).dot(Rinv)).dot(zr)
    #Plotting
    xs.insert(i,i)
    ds.insert(i,x[0,0])
    vs.insert(i,x[1,0])
    #Iteration and Z Updates
    z = z + 9.5 + random.uniform(0,1)
    i = i + 1
print(x)
xs.pop()
ds.pop()
vs.pop()
#plt.scatter(xs,ds)
plt.scatter(xs,vs)
    
    
