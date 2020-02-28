#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 12:25:48 2020

@author: Axel
"""
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
# This must come before the next imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))
from project.numerical.distribution_equations import DistributionEquations
from project.numerical.reaction_forces import AileronGeometry

E = DistributionEquations()
G = AileronGeometry()

def T(x):
    return E.torsion_x(x)



sample_steps = 50
min_x, max_x = min(G.station_x_coords()), max(G.station_x_coords())
x_steps = np.linspace(min_x, max_x, sample_steps)



C = 0.515     #[m]
h = 0.248     #[m]
t = 0.0011    #[m]
t_sp = 0.0022 #[m]
Gs = 28*10**9  #[N/m^2]
a = 0.4102    #[m]
A_1 = 0.0241  #[m^2]
A_2 = 0.0485  #[m^2]

A = np.zeros([3,3])

A[0,0] = 2*A_1
A[0,1] = 2*A_2
A[0,2] = 0
A[1,0] = 1/(2*A_1)*(np.pi*h/(Gs*t*2)+h/(Gs*t_sp))
A[1,1] = -1/(2*A_1)*h/(Gs*t_sp)
A[1,2] = -1
A[2,0] = -1/(2*A_2)*h/(Gs*t_sp)
A[2,1] = 1/(2*A_2)*(h/(Gs*t_sp)+2*a/(Gs*t))
A[2,2] = -1

F = np.array([0,0,0])
q = np.zeros([50,7])
 #index q : 0=x 1=Torque(x) 2=q1 3=q_2 4=tau_circ 5=tau_spar 6=tau_13/23 All signs cc positive

for i in range(0,50):
    F = np.array([T(x_steps[i]),0,0])
    R = np.linalg.solve(A, F)
    q[i,0] = x_steps[i]
    q[i,1] = T(x_steps[i])
    q[i,2] = R[0]
    q[i,3] = R[1]
    q[i,4] = R[0]/t
    q[i,5] = (R[0]-R[1])/t_sp
    q[i,6] = R[1]/t

plt.plot(q[0:50,0],q[0:50,4])
plt.plot(q[0:50,0],q[0:50,5])
plt.plot(q[0:50,0],q[0:50,6])

plt.show()












