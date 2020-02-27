#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 12:25:48 2020

@author: Axel
"""
import numpy as np


C = 0.515     #[m]
h = 0.248     #[m]
t = 0.0011    #[m]
t_sp = 0.0022 #[m]
G = 28*10**9  #[N/m^2]
a = 0.4102    #[m]
A_1 = 0.0241  #[m^2]
A_2 = 0.0485  #[m^2]

A = np.zeros([3,3])

A[0,0] = 2*A_1
A[0,1] = 2*A_2
A[0,2] = 0
A[1,0] = 1/(2*A_1)*(np.pi*h/(G*t*2)+h/(G*t_sp))
A[1,1] = -1/(2*A_1)*h/(G*t_sp)
A[1,2] = -1
A[2,0] = -1/(2*A_2)*h/(G*t_sp)
A[2,1] = 1/(2*A_2)*(h/(G*t_sp)+2*a/(G*t))
A[2,2] = -1
