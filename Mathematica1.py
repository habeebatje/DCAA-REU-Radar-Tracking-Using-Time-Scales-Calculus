# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 10:24:34 2022

@author: habee
"""
import numpy as np
from numpy import transpose as tr
from numpy import array as mat
from numpy import linalg
import matplotlib.pyplot as plt
import random

#Matrices
mu=1
A(mu)= mat([[0, mu, 0, 0], [0, 0, 0, 0], [0, 0, 0, mu], [0, 0, 0, 0]])
H= mat([[1,0,0,0],[0,0,1,0]])

Covw(mu)= mat([[((mu**3)/3),((mu**2)/2),0,0], [((mu**2)/2),mu,0,0], [0,0,((mu**3)/3),((mu**2)/2)],[0,0,((mu**2)/2),mu]])
R(mu)= (1/mu)(mat([[0.001,0],[0,0.1]]))

#*Set Up Time Scale*

tsl=1000
#timescale limit
MU= 0
##FIX LATER!!!!!!!!!!!!!!!!!!!!!!!1

xvals=mat([[0],[0],[0],[0]])
yvals=mat([[0],[0]])

while