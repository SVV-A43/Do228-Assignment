#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
title: reaction_forces
project: Do228-Assignment
date: 2/23/2020
author: lmaio
"""

import os
import sys
import numpy as np
from tqdm import tqdm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports
from project.numerical.Loading.interpolation import InterpolateRBF
from project.numerical.Loading.integration import def_integral, def_integral
from project.numerical.Loading.aileron_geometry import AileronGeometry


# Known input variables:
def q_tilde(**kwargs):
    num_bins = kwargs.pop('num_bins', 100)
    aileron = AileronGeometry()
    q_x = []
    x_coords = []

    for station in range(aileron.num_span_stations):
        x, z, p = aileron.station_data(station)
        int_fn = InterpolateRBF(z, p)
        station_load = def_integral(int_fn.interpolate, min(z), max(z), num_bins=num_bins)
        q_x.append(station_load)
        x_coords.append(x[0])
    q_x = np.array(q_x)
    x_coords = np.array(x_coords)

    # Create interpolation function
    q_tilde_x = InterpolateRBF(x_coords, q_x)

    return q_tilde_x

def tau_tilde(**kwargs):
    num_bins = kwargs.pop('num_bins', 100)
    G = AileronGeometry()
    tau_x = []
    x_coords = []

    for station in range(G.num_span_stations):
        x, z, p = G.station_data(station)
        int_fn = InterpolateRBF(z, p)

        def tau_inner_fn(z):
            return np.multiply(int_fn.interpolate(z), np.squeeze(z - G.z_tilde))

        station_load = def_integral(tau_inner_fn, min(z), max(z), num_bins=num_bins)
        tau_x.append(station_load)
        x_coords.append(x[0])
    tau_x = np.array(tau_x)
    x_coords = np.array(x_coords)

    # Create interpolation function
    tau_tilde_x = InterpolateRBF(x_coords, tau_x)

    return tau_tilde_x


def equilibrium_eq_coefficients():
    # Geometry
    G = AileronGeometry()

    a = np.zeros((11,11))

    a[0, 0] = -(G.l_a - G.x1)
    a[0, 1] = -(G.l_a - G.x2)
    a[0, 2] = -(G.l_a - G.x3)
    a[0, 6] = -1*np.sin(G.theta)*(G.l_a - G.x_a_1)

    a[1, 3] = -(G.l_a - G.x1)
    a[1, 4] = -(G.l_a - G.x2)
    a[1, 5] = -(G.l_a - G.x3)
    a[1, 6] = -1*np.cos(G.theta)*(G.l_a - G.x_a_1)

    a[2, 0] = (G.z_h - G.z_tilde)
    a[2, 1] = (G.z_h - G.z_tilde)
    a[2, 2] = (G.z_h - G.z_tilde)
    a[2, 6] = np.sin(G.theta)*(G.z_h - G.z_tilde) - np.cos(G.theta)*(G.y_p)

    a[3, 0] = 1
    a[3, 1] = 1
    a[3, 2] = 1
    a[3, 6] = np.sin(G.theta)

    a[4, 3] = 1
    a[4, 4] = 1
    a[4, 5] = 1
    a[4, 6] = np.cos(G.theta)

    a[5, 0] = -1*(G.x2 - G.x1)**3 / (6 * G.E * G.I_zz)
    a[5, 6] = -1*(np.sin(G.theta) * (G.x2 - G.x_a_1)**3) / (6 * G.E * G.I_zz)
    a[5, 7] = G.x2
    a[5, 8] = 1

    a[6, 3] = -1*(G.x2 - G.x1)**3 / (6*G.E*G.I_yy)
    a[6, 6] = -1*(np.cos(G.theta) * (G.x2 - G.x_a_1)**3) / (6 * G.E * G.I_yy)
    a[6, 9] = G.x2
    a[6, 10] = 1

    a[7, 9] = G.x1
    a[7, 10] = 1

    a[8, 3] = -1*(G.x3 - G.x1)**3 / (6*G.E*G.I_yy)
    a[8, 4] = -1*(G.x3 - G.x2)**3 / (6*G.E*G.I_yy)
    a[8, 6] = -1*(np.cos(G.theta) * (G.x3 - G.x_a_1)**3) / (6 * G.E * G.I_yy)
    a[8, 9] = G.x3
    a[8, 10] = 1

    a[9, 7] = G.x1
    a[9, 8] = 1

    a[10, 0] = -1*(G.x3 - G.x1)**3 / (6*G.E*G.I_zz)
    a[10, 1] = -1*(G.x3 - G.x2)**3 / (6*G.E*G.I_zz)
    a[10, 6] = -1*(np.sin(G.theta) * (G.x3 - G.x_a_1)**3) / (6 * G.E * G.I_zz)
    a[10, 7] = G.x3
    a[10, 8] = 1

    return a

def equilibrium_eq_resultants():
    G = AileronGeometry()

    b = np.zeros((11,1))
    q_x = q_tilde()
    t_x = tau_tilde()

    b[0, 0] = def_integral(q_x.interpolate, 0, G.l_a) + (G.P * np.sin(G.theta) * (G.l_a - G.x_a_2))
    b[1, 0] = G.P * np.cos(G.theta) * (G.l_a - G.x_a_2)
    b[2, 0] = -1*def_integral(t_x.interpolate, 0, G.l_a) + (G.P * np.cos(G.theta) * G.y_p) - (G.P * np.sin(G.theta) * (G.z_h - G.z_tilde))
    b[3, 0] = -1*np.sin(G.theta) - def_integral(q_x.interpolate, 0, G.l_a)
    b[4, 0] = -1*G.P * np.cos(G.theta)
    b[5, 0] = def_integral(q_x.interpolate, 0, G.x2, num_var_integrals=4, num_bins=20) / (G.E * G.I_zz) + G.d_1*np.sin(G.theta)
    b[8, 0] = G.P * np.cos(G.theta) * (G.x3 - G.x_a_2)**3 / (6*G.E * G.I_yy) + G.d_3*np.sin(G.theta)
    b[9, 0] = def_integral(q_x.interpolate, 0, G.x1, num_var_integrals=4, num_bins=20) / (G.E * G.I_zz) - G.d_1*np.cos(G.theta)
    b[10, 0] = def_integral(q_x.interpolate, 0, G.x3, num_var_integrals=4, num_bins=20) / (G.E * G.I_zz) + G.P * np.sin(G.theta) / (6 * G.E * G.I_zz) - G.d_3*np.cos(G.theta)

    return -1 * b

def solve_reaction_forces():
    A = equilibrium_eq_coefficients()
    b = equilibrium_eq_resultants()
    r_force_names = ['R 1,y',
                     'R 2,y',
                     'R 3,y',
                     'R 1,z',
                     'R 2,z',
                     'R 3,z',
                     '    F',
                     '  C 1',
                     '  C 2',
                     '  C 3',
                     '  C 4']
    x = np.linalg.solve(A,b)
    return np.squeeze(x), r_force_names


if __name__ == '__main__':
    # A = equilibrium_eq_coefficients()
    # b = equilibrium_eq_resultants()

    x, names= solve_reaction_forces()
    for val, name in zip(x, names):
        print(f'{name}: {val}')
