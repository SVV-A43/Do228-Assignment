#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
title: deflection
project: 
created: 25/02/2020
author: lmaio
"""
# Main Structural Analysis model we need to develop

# Imports
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports
from project.numerical.Loading.aileron_geometry import AileronGeometry
from project.numerical.Loading.integration import def_integral
from project.numerical.reaction_forces import reaction_forces


# CODE...


# Equations of motion


def deflection_y(x, reaction_forces, q_x):
    G = AileronGeometry()
    r= reaction_forces

    def mac(x):
        return np.maximum(0, x)

    R1y = r[0]
    R2y = r[1]
    R3y = r[2]
    F   = r[6]
    C1  = r[7]
    C2  = r[8]

    terms = np.zeros((6,))

    terms[1] = -1* def_integral(q_x, 0, x, num_var_integrals=4, num_bins=20)

    terms[0] = R1y * mac(x-G.x1)**3 / 6
    terms[2] = F * np.sin(G.theta) / 6 * mac(x-G.x_a_1)**3
    terms[3] = R2y / 6 * mac(x-G.x2)**3
    terms[4] = -1* G.P * np.sin(G.theta) / 6 * mac(x-G.x_a_2)**3
    terms[5] = R3y / 6 * mac(x-G.x3)**3

    v = (-1/(G.E * G.I_zz) * terms).sum() + C1*x + C2

    return v

def plot_deflection_x():
    G = AileronGeometry()
    q_x = G.q_tilde()
    r, _ = reaction_forces()

    deflect_data = np.zeros((G.num_span_stations, 2))

    deflect_data[:, 0] = G.station_x_coords()

    for i, x in enumerate(deflect_data[:,0]):
        deflect_data[i, 1] = deflection_y(x, r, q_x)

    plt.plot(deflect_data[:,0], deflect_data[:,1])
    plt.show()




if __name__ == '__main__':
    plot_deflection_x()




