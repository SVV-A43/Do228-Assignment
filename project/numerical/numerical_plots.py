#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
title: numerical_plots
project: 
created: 26/02/2020
author: lmaio
"""

import numpy as np
import os
import sys
import matplotlib.pyplot as plt


sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports
from project.numerical.Loading.aileron_geometry import AileronGeometry
from project.numerical.distribution_equations import DistributionEquations



#### Parts leading to deflection in Y direction

# ------------------------ LOCAL REFERENCE FRAME -----------------------

def plots_y_distribution(steps=50):
    G = AileronGeometry()
    E = DistributionEquations()

    min_x, max_x = min(G.station_x_coords()), max(G.station_x_coords())

    data = np.zeros((steps, 4))
    data[:, 0] = np.linspace(min_x, max_x, steps)
    data[:, 1] = E.deformations_in_global_coords(steps)[:, 1]
    for i, x in enumerate(data[:, 0]):
        data[i, 1] = E.deflection_y(x)
        data[i, 2] = E.moment_about_z(x)
        data[i, 3] = E.shear_y(x)


    fig, axs = plt.subplots(2, 2, sharex=True)
    axs[0, 0].plot(data[:, 0], data[:, 1])
    axs[0, 0].set(ylabel=r"$v(x)$ [m]")
    axs[0, 0].set_title('Deflection in y')

    axs[1, 0].plot(data[:, 0], data[:, 2])
    axs[1, 0].set(ylabel=r"$M_z (x)$ [Nm]")
    axs[1, 0].set_title("Bending moment about z")

    axs[1, 1].plot(data[:, 0], data[:, 3])
    axs[1, 1].set(ylabel=r"$S_y (x)$ [N]")
    axs[1, 1].set_title('Shear in y')

    for ax in axs.flat:
        ax.set(xlabel="x")

    # for ax in axs.flat:
    #     ax.label_outer()

    # plt.tight_layout()
    plt.show()

def plots_z_distribution(steps=50):
    G = AileronGeometry()
    E = DistributionEquations()

    min_x, max_x = min(G.station_x_coords()), max(G.station_x_coords())

    data = np.zeros((steps, 4))
    data[:, 0] = np.linspace(min_x, max_x, steps)
    data[:, 1] = E.deformations_in_global_coords(steps)[:, 2]
    for i, x in enumerate(data[:, 0]):
        data[i, 1] = E.deflection_z(x)
        data[i, 2] = E.moment_about_y(x)
        data[i, 3] = E.shear_z(x)


    fig, axs = plt.subplots(2, 2, sharex=True)
    axs[0, 0].plot(data[:, 0], data[:, 1])
    axs[0, 0].set(ylabel=r"$v(x)$ [m]")
    axs[0, 0].set_title('Deflection in z')

    axs[1, 0].plot(data[:, 0], data[:, 2])
    axs[1, 0].set(ylabel=r"$M_z (x)$ [Nm]")
    axs[1, 0].set_title("Bending moment about y")

    axs[1, 1].plot(data[:, 0], data[:, 3])
    axs[1, 1].set(ylabel=r"$S_z (x)$ [N]")
    axs[1, 1].set_title('Shear in z')

    for ax in axs.flat:
        ax.set(xlabel="x")

    plt.show()



def plots_torque_distribution(steps=50):
    G = AileronGeometry()
    E = DistributionEquations()

    min_x, max_x = min(G.station_x_coords()), max(G.station_x_coords())

    data = np.zeros((steps, 4))
    data[:, 0] = np.linspace(min_x, max_x, steps)

    for i, x in enumerate(data[:, 0]):
        data[i, 1] = E.angle_of_twist(x)
        data[i, 2] = E.torsion_x(x)

    gs = gridspec.GridSpec(2, 2)

    pl.figure()
    ax = pl.subplot(gs[0, 0])  # row 0, col 0
    pl.plot([0, 1])

    ax = pl.subplot(gs[0, 1])  # row 0, col 1
    pl.plot([0, 1])

    ax = pl.subplot(gs[1, :])  # row 1, span all columns
    pl.plot([0, 1])

    fig, axs = plt.subplots(2, sharex=True)
    axs[0].plot(data[:, 0], -1*data[:, 1])
    axs[0].set(ylabel=r"$\theta(x)$ [deg/m]")
    axs[0].set_title('Angle of Twist')

    axs[1].plot(data[:, 0], -1*data[:, 2])
    axs[1].set(ylabel=r"$T(x)$ [Nm]")
    axs[1].set_title("Torque")

    for ax in axs.flat:
        ax.set(xlabel="x")
        ax.label_outer()

    plt.show()


### -------------- Plots for Validation -----------------
def plot_lateral_deflection(steps=50, coord_sys='global', plot_3d=False):
    E = DistributionEquations()
    G = AileronGeometry()

    if coord_sys=='global':
        deflect_data = E.deformations_in_global_coords(steps=steps)
    else:
        deflect_data = E.deformations_in_local_coords(steps=steps)

    if plot_3d:
        ax3d = plt.axes(projection='3d')
        ax3d.scatter3D(deflect_data[:,0], deflect_data[:,2], deflect_data[:,1], c=deflect_data[:,2])
        # ax3d.plot_trisurf(x_arr, z_arr, P_arr, cmap='viridis', edgecolor='none')
        ax3d.set_xlabel('Spanwise coordinate x [m]')
        ax3d.set_zlabel('Vertical Deflection y [m]')
        ax3d.set_ylabel('Lateral Deflection z [m]')
        ax3d.view_init(10, 45)
        ax3d.set_title('Deflection')
        ax3d.invert_xaxis()
        plt.show()

    else:
        plt.plot(deflect_data[:, 0], deflect_data[:, 1]*1000)
        plt.title(f'Lateral Deflections; {coord_sys} reference frame; {G.aircraft}')
        plt.ylabel('v(x) [mm]')
        plt.xlabel('Spanwise x-location [m]')
        plt.show()



if __name__ == '__main__':
    # steps = int(input('Number of steps for y_deflection: '))
    # plot_lateral_deflection(steps=50, coord_sys='local')
    # plot_angle_of_twist(steps=50)

    plots_y_distribution()
    plots_z_distribution()
    # plots_torque_distribution()

    # plot_lateral_deflection()