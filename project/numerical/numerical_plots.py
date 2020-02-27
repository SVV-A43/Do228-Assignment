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
from mpl_toolkits import mplot3d

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports
from project.numerical.Loading.aileron_geometry import AileronGeometry
from project.numerical.reaction_forces import reaction_forces
from project.numerical.distribution_equations import DistributionEquations



#### Parts leading to deflection in Y direction

# ------------------------ LOCAL REFERENCE FRAME -----------------------

def plots_y_deformation(steps=50):
    G = AileronGeometry()
    E = DistributionEquations()

    min_x, max_x = min(G.station_x_coords()), max(G.station_x_coords())

    data = np.zeros((steps, 4))
    data[:, 0] = np.linspace(min_x, max_x, steps)

    for i, x in enumerate(data[:, 0]):
        data[i, 1] = E.deflection_y(x)
        data[i, 2] = E.moment_about_z(x)
        data[i, 3] = E.shear_y(x)


    fig, axs = plt.subplots(2, 2)
    axs[0, 0].plot(data[:, 0], data[:, 1])
    axs[0, 0].set(ylabel="v(x)")
    axs[0, 0].set_title('Deflections in local y')

    axs[1, 0].plot(data[:, 0], data[:, 2])
    axs[1, 0].set(ylabel="M_z(x)")
    axs[1, 0].set_title("Moment about local z'")

    axs[1, 1].plot(data[:, 0], data[:, 3])
    axs[1, 1].set(ylabel="Sy(x)")
    axs[1, 1].set_title('Shear in local y')

    for ax in axs.flat:
        ax.set(xlabel="x'")

    for ax in axs.flat:
        ax.label_outer()

    # plt.tight_layout()
    plt.show()





## -------------- TORSION ---------------------

def plot_angle_of_twist(steps=50):
    G = AileronGeometry()
    E = DistributionEquations()
    min_x, max_x = min(G.station_x_coords()), max(G.station_x_coords())

    deflect_data = np.zeros((steps, 2))
    deflect_data[:, 0] = np.linspace(min_x, max_x, steps)

    for i, x in enumerate(deflect_data[:, 0]):
        deflect_data[i, 1] = E.angle_of_twist(x)


    plt.plot(deflect_data[:,0], deflect_data[:,1], color='red')
    plt.ylabel('Angle of twist [degrees]')
    plt.xlabel('x-coordinate (along span) [m]')
    plt.title(f'Angle of twist along span')
    # plt.legend()

    plt.show()


### -------------- GENERAL PLOTS -----------------
def plot_lateral_deflection(steps=50, coord_sys='global', plot_3d=False):
    E = DistributionEquations()

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
        fig, axs = plt.subplots(2, sharex=True)
        fig.suptitle(f'Lateral deflections along span; {coord_sys} coordinate system')
        fig.subplots_adjust(top=0.9)

        axs[0].plot(deflect_data[:, 0], deflect_data[:, 1])
        axs[0].set(ylabel='Deflections in y-direction')

        axs[1].plot(deflect_data[:, 0], deflect_data[:, 2])
        axs[1].set(ylabel='Deflections in z-direction')
        axs[1].invert_yaxis()

        for ax in axs.flat:
            ax.set(xlabel='x coordinate along span')

        # for ax in axs.flat:
        #     ax.label_outer()

        # plt.tight_layout()
        plt.show()





if __name__ == '__main__':
    # steps = int(input('Number of steps for y_deflection: '))
    # plot_lateral_deflection(steps=50, coord_sys='local')
    # plot_angle_of_twist(steps=50)

    plots_y_deformation()