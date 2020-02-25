#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
title: reaction_forces
project: Do228-Assignment
date: 2/17/2020
author: lmaio
"""


# Imports
import sys
import os
import numpy as np
from definitions import AERO_LOADING_DATA

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports
from project.numerical.Loading.integration import def_integral
from project.numerical.Loading.interpolation import InterpolateRBF


# CODE...


class AileronGeometry():
    def __init__(self, filename=AERO_LOADING_DATA):
        self.C_a  = 0.515       # [m]
        self.l_a  = 2.691       # [m]
        self.x1   = 0.174       # [m]
        self.x2   = 1.051       # [m]
        self.x3   = 2.512       # [m]
        self.x_a  = 30.0 / 100  # [m]
        self.x_a_1 = self.x2 - self.x_a / 2     # [m]
        self.x_a_2 = self.x2 + self.x_a / 2     # [m]
        self.z_tilde = -0.35375    # [m] distance to shear center #FIXME: THIS IS FROM VERIFICATION MODEL
        self.h    = 24.8 / 100  # [m]
        self.d_1  = 1.034 / 100 # [m]
        self.d_3  = 2.066 / 100 # [m]
        self.z_h = -1* self.h / 2   # [m] #TODO: Verify this is true
        self.y_p = self.h / 2   # [m]
        self.theta = np.deg2rad(25)         # [rad]

        self.t_sk = 1.1         # [mm]
        self.t_sp = 2.2         # [mm]
        self.t_st = 1.2         # [mm]
        self.h_st = 1.5 * 10    # [mm]
        self.w_st = 3.0 * 10    # [mm]
        self.n_st = 11          # [-]
        self.I_xx = None
        self.I_yy = 5.643650631210155e-05  # [m^4]
        self.I_zz = 1.4221372629975417e-05 # [m^4]

        ### Constants:
        self.E = 73.1 * 10**9       # [Pa]
        self.P = 20.6 * 10 ** 3     # [N]

        self.__pressure = np.genfromtxt(filename, delimiter=',')
        self.num_span_stations = len(self.__pressure[0, :])
        self.num_chord_stations = len(self.__pressure[:, 0])

        self.__x_coords = []
        self.__z_coords = []


        # Calculate x-point_coords of stations along span
        for i in range(self.num_span_stations):
            self.__x_coords.append(self.__xi(i+1))

        # Calculate z-point coords of station along chord
        for j in range(self.num_chord_stations):
            self.__z_coords.append(self.__zi(j+1))

        num_data_pts = self.num_span_stations * self.num_chord_stations
        self.load_data = np.zeros((num_data_pts, 3))

        n = 0
        for i, x in enumerate(self.__x_coords):
            for j, z in enumerate(self.__z_coords):
                self.load_data[n, 0] = x
                self.load_data[n, 1] = z
                self.load_data[n, 2] = self.__pressure[j,i]
                n += 1


        # Create Array of the (z,x coordinates where the load data is found) LOCAL COORD SYSTEM AS DEFINED IN READER fig3
        # The order [z, x] here is to match that of the loading data, which is also z-coord per row
        self.__load_coords = np.empty((self.num_chord_stations, self.num_span_stations), dtype=object)
        for z in range(self.num_chord_stations):
            for x in range(self.num_span_stations):
                self.__load_coords[z, x] = (self.__z_coords[z], self.__x_coords[x])

    # Functions to calculate geometry
    def __th_xi(self, i):
        return (i - 1) / self.num_span_stations * np.pi

    def __xi(self, i):
        return .5 * (self.l_a / 2 * (1 - np.cos(self.__th_xi(i))) + self.l_a / 2 * (
                    1 - np.cos(self.__th_xi(i + 1))))

    def __th_zi(self, i):
        return (i - 1) / self.num_chord_stations * np.pi

    def __zi(self, i):
        return -.5 * (self.C_a / 2 * (1 - np.cos(self.__th_zi(i))) + self.C_a / 2 * (
                    1 - np.cos(self.__th_zi(i + 1))))

    # Methods to use:
    def station_x_coords(self):
        return self.__x_coords

    def station_z_coords(self):
        return self.__z_coords

    def loading_data_coords(self):
        return self.__load_coords

    def data_x_z_p(self):
        return self.load_data[:,0], self.load_data[:,1], self.load_data[:,2]

    def station_data(self, id):
        s = id * 81 #start idx
        e = (id + 1) * 81 #end idx
        if e > self.load_data.shape[0]:
            e = -1

        return self.load_data[s:e, 0], self.load_data[s:e, 1], self.load_data[s:e, 2]

    def q_tilde(self, **kwargs):
        num_bins = kwargs.pop('num_bins', 100)
        q_x = []
        x_coords = []

        for station in range(self.num_span_stations):
            x, z, p = self.station_data(station)
            int_fn = InterpolateRBF(z, p)
            station_load = def_integral(int_fn.interpolate, min(z), max(z), num_bins=num_bins)
            q_x.append(station_load)
            x_coords.append(x[0])
        q_x = np.array(q_x)
        x_coords = np.array(x_coords)

        # Create interpolation function
        q_tilde_x = InterpolateRBF(x_coords, q_x)

        return q_tilde_x

    def tau_tilde(self, **kwargs):
        num_bins = kwargs.pop('num_bins', 100)
        tau_x = []
        x_coords = []

        for station in range(self.num_span_stations):
            x, z, p = self.station_data(station)
            int_fn = InterpolateRBF(z, p)

            def tau_inner_fn(z):
                return np.multiply(int_fn.interpolate(z), np.squeeze(z - self.z_tilde))

            station_load = def_integral(tau_inner_fn, min(z), max(z), num_bins=num_bins)
            tau_x.append(station_load)
            x_coords.append(x[0])
        tau_x = np.array(tau_x)
        x_coords = np.array(x_coords)

        # Create interpolation function
        tau_tilde_x = InterpolateRBF(x_coords, tau_x)

        return tau_tilde_x



if __name__ == '__main__':
    aileron = AileronGeometry()
    x = aileron.station_x_coords()
    z = aileron.station_z_coords()

    # All data, coordinates repeated
    x_a, z_a, p_a = aileron.data_x_z_p()

    print(aileron.station_x_coords())
    print(aileron.station_z_coords())

