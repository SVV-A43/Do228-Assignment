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
from mpl_toolkits import mplot3d

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports
from project.numerical.Loading.aileron_geometry import AileronGeometry
from project.numerical.Loading.integration import def_integral
from project.numerical.reaction_forces import reaction_forces


# CODE...


# Equations of motion

def mac(x): # Macauly step
    return np.maximum(0, x)


def deflection_y(x, reaction_forces, q_x):
    G = AileronGeometry()
    r = reaction_forces

    R1y = r[0]
    R2y = r[1]
    R3y = r[2]
    F   = r[6]
    C1  = r[7]
    C2  = r[8]

    terms = np.zeros((6,))

    terms[0] = R1y /6 * mac(x - G.x1)** 3
    terms[1] = -1* def_integral(q_x, 0, x, num_var_integrals=4, num_bins=20)
    terms[2] = F * np.sin(G.theta) / 6 * mac(x-G.x_a_1)**3
    terms[3] = R2y / 6 * mac(x-G.x2)**3
    terms[4] = -1* G.P * np.sin(G.theta) / 6 * mac(x-G.x_a_2)**3
    terms[5] = R3y / 6 * mac(x-G.x3)**3

    v = (-1/(G.E * G.I_zz) * terms).sum() + C1*x + C2

    return v

def deflection_z(x, reaction_forces):
    G = AileronGeometry()
    r = reaction_forces

    R1z = r[3]
    R2z = r[4]
    R3z = r[5]
    F = r[6]
    C3 = r[9]
    C4 = r[10]

    terms = np.zeros((5,))

    terms[0] = -1* R1z/ 6 * mac(x - G.x1)**3
    terms[1] = -1*F * np.cos(G.theta) / 6 * mac(x - G.x_a_1)**3
    terms[2] = -1*R2z / 6 * mac(x - G.x2)**3
    terms[3] = G.P*np.cos(G.theta) / 6 * mac(x - G.x_a_2)**3
    terms[4] = -1*R3z / 6 * mac(x - G.x3)**3

    w = (-1/(G.E * G.I_yy) * terms).sum() + C3*x + C4

    return w


def angle_of_twist(x, reaction_forces, t_x):
    G = AileronGeometry()
    r = reaction_forces

    # Terms of the eq excluding the +C5 value
    def variable_terms(x, r, t_x):
        R1y = r[0]
        R2y = r[1]
        R3y = r[2]
        F = r[6]


        terms = np.zeros((8,))

        terms[0] = -1*def_integral(t_x, 0, x, num_var_integrals=2)
        terms[1] = R1y*(G.z_h - G.z_tilde)*mac(x-G.x1)
        terms[2] = R2y*(G.z_h - G.z_tilde)*mac(x-G.x2)
        terms[3] = R3y*(G.z_h - G.z_tilde)*mac(x-G.x3)
        terms[4] = -1*F*np.cos(G.theta)*G.y_p*mac(x-G.x_a_1)
        terms[5] = F*np.sin(G.theta)*(0-G.z_tilde)*mac(x-G.x_a_1)
        terms[6] = G.P*np.cos(G.theta)*G.y_p*mac(x-G.x_a_2)
        terms[7] = -1*G.P*np.sin(G.theta)*(0-G.z_tilde)*mac(x-G.x_a_2)

        return (1/(G.G*G.J) * terms).sum()

    C5 = variable_terms(G.x_a_1, r, t_x)

    theta_rad = variable_terms(x, r, t_x) + -C5
    return np.degrees(theta_rad)






def plot_deflection_xz(steps=50):
    G = AileronGeometry()
    q_x = G.q_tilde()
    t_x = G.tau_tilde()
    r, _ = reaction_forces()

    min_x, max_x = min(G.station_x_coords()), max(G.station_x_coords())

    deflect_data = np.zeros((steps, 4))
    deflect_data[:, 0] = np.linspace(min_x, max_x, steps)

    for i, x in enumerate(deflect_data[:,0]):
        deflect_data[i, 1] = deflection_y(x, r, q_x)
        deflect_data[i, 2] = deflection_z(x, r)
        deflect_data[i, 3] = angle_of_twist(x, r, t_x)


    ax = plt.axes(projection='3d')
    ax.scatter3D(deflect_data[:,0], deflect_data[:,2], deflect_data[:,1], c=deflect_data[:,2])
    # ax.plot_trisurf(x_arr, z_arr, P_arr, cmap='viridis', edgecolor='none')
    ax.set_xlabel('Spanwise coordinate x [m]')
    ax.set_zlabel('Vertical Deflection y [m]')
    ax.set_ylabel('Lateral Deflection z [m]')
    ax.view_init(10, 45)
    ax.set_title('Deflection')
    ax.invert_xaxis()
    plt.show()

def plot_angle_of_twist(steps=50):
    G = AileronGeometry()
    t_x = G.tau_tilde()
    r, _ = reaction_forces()

    min_x, max_x = min(G.station_x_coords()), max(G.station_x_coords())

    deflect_data = np.zeros((steps, 2))
    deflect_data[:, 0] = np.linspace(min_x, max_x, steps)

    for i, x in enumerate(deflect_data[:, 0]):
        deflect_data[i, 1] = angle_of_twist(x, r, t_x)


    plt.scatter(deflect_data[:,0], deflect_data[:,1], color='red')
    plt.ylabel('Angle of twist [degrees]')
    plt.xlabel('x-coordinate (along span) [m]')
    plt.title(f'Angle of twist along span')
    # plt.legend()

    plt.show()


if __name__ == '__main__':
    # steps = int(input('Number of steps for y_deflection: '))
    # plot_deflection_xz(steps=50)
    plot_angle_of_twist(steps=500)



