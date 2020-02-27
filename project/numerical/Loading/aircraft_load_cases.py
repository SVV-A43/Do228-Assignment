#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
title: aircraft_load_cases
project: 
created: 27/02/2020
author: lmaio
"""
# Imports
import sys
import os
import numpy as np

from definitions import AERO_LOADING_DATA

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports


class DornierDo228():
    def __init__(self, filename=AERO_LOADING_DATA):
        self.C_a = 0.515  # [m]
        self.l_a = 2.691  # [m]
        self.x1 = 0.174  # [m]
        self.x2 = 1.051  # [m]
        self.x3 = 2.512  # [m]
        self.x_a = 30.0 / 100  # [m]
        self.x_a_1 = self.x2 - self.x_a / 2  # [m]
        self.x_a_2 = self.x2 + self.x_a / 2  # [m]
        self.z_tilde = -0.35375  # [m] distance to shear center #FIXME: THIS IS FROM VERIFICATION MODEL
        self.h = 24.8 / 100  # [m]
        self.d_1 = 1.034 / 100  # [m]
        self.d_3 = 2.066 / 100  # [m]
        self.z_h = -1 * self.h / 2  # [m] #TODO: Verify this is true
        self.y_p = self.h / 2  # [m]
        self.theta = np.deg2rad(25)  # [rad]

        self.t_sk = 1.1  # [mm]
        self.t_sp = 2.2  # [mm]
        self.t_st = 1.2  # [mm]
        self.h_st = 1.5 * 10  # [mm]
        self.w_st = 3.0 * 10  # [mm]
        self.n_st = 11  # [-]

        self.I_xx = None
        self.I_yy = 5.643650631210155e-05  # [m^4]
        self.I_zz = 1.4221372629975417e-05  # [m^4]

        self.J = 0.00020500801555445113
        self.J_polar = self.I_yy + self.I_zz  # [m^4]

        ### Constants:
        self.E = 73.1 * 10 ** 9  # [Pa]
        self.G = 28 * 10 ** 9  # [Pa]

        ### Loads:
        self.P = 20.6 * 10 ** 3  # [N]
        self.__pressure = np.genfromtxt(filename, delimiter=',') * 10 ** 3

        ### Boundary Conditions
        self.bound_conds = np.zeros((11, 1))
        self.bound_conds[7, 0] = -1 * self.d_1 * np.sin(self.theta)
        self.bound_conds[8, 0] = -1 * self.d_3 * np.sin(self.theta)
        self.bound_conds[9, 0] = self.d_1 * np.cos(self.theta)
        self.bound_conds[10, 0] = self.d_3 * np.cos(self.theta)

        ### Class Attributes

        self.num_span_stations = len(self.__pressure[0, :])
        self.num_chord_stations = len(self.__pressure[:, 0])

        self.__x_coords = []
        self.__z_coords = []

        # Calculate x-point_coords of stations along span
        for i in range(self.num_span_stations):
            self.__x_coords.append(self.__xi(i + 1))

        # Calculate z-point coords of station along chord
        for j in range(self.num_chord_stations):
            self.__z_coords.append(self.__zi(j + 1))

        num_data_pts = self.num_span_stations * self.num_chord_stations
        self.load_data = np.zeros((num_data_pts, 3))

        n = 0
        for i, x in enumerate(self.__x_coords):
            for j, z in enumerate(self.__z_coords):
                self.load_data[n, 0] = x
                self.load_data[n, 1] = z
                self.load_data[n, 2] = self.__pressure[j, i]
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




# Used for validation
class Boeing737():
    def __init__(self, filename=AERO_LOADING_DATA):
        self.C_a = 0.515  # [m]
        self.l_a = 2.691  # [m]
        self.x1 = 0.174  # [m]
        self.x2 = 1.051  # [m]
        self.x3 = 2.512  # [m]
        self.x_a = 30.0 / 100  # [m]
        self.x_a_1 = self.x2 - self.x_a / 2  # [m]
        self.x_a_2 = self.x2 + self.x_a / 2  # [m]
        self.z_tilde = -0.35375  # [m] distance to shear center #FIXME: THIS IS FROM VERIFICATION MODEL
        self.h = 24.8 / 100  # [m]
        self.d_1 = 1.034 / 100  # [m]
        self.d_3 = 2.066 / 100  # [m]
        self.z_h = -1 * self.h / 2  # [m] #TODO: Verify this is true
        self.y_p = self.h / 2  # [m]
        self.theta = np.deg2rad(25)  # [rad]

        self.t_sk = 1.1  # [mm]
        self.t_sp = 2.2  # [mm]
        self.t_st = 1.2  # [mm]
        self.h_st = 1.5 * 10  # [mm]
        self.w_st = 3.0 * 10  # [mm]
        self.n_st = 11  # [-]

        self.I_xx = None
        self.I_yy = 5.643650631210155e-05  # [m^4]
        self.I_zz = 1.4221372629975417e-05  # [m^4]

        self.J = 0.00020500801555445113
        self.J_polar = self.I_yy + self.I_zz  # [m^4]

        ### Constants:
        self.E = 73.1 * 10 ** 9  # [Pa]
        self.G = 28 * 10 ** 9  # [Pa]

        ### Loads:
        self.P = 20.6 * 10 ** 3  # [N]
        self.__pressure = np.genfromtxt(filename, delimiter=',') * 10 ** 3

        ### Boundary Conditions
        self.bound_conds = np.zeros((11, 1))
        self.bound_conds[7, 0] = -1 * self.d_1 * np.sin(self.theta)
        self.bound_conds[8, 0] = -1 * self.d_3 * np.sin(self.theta)
        self.bound_conds[9, 0] = self.d_1 * np.cos(self.theta)
        self.bound_conds[10, 0] = self.d_3 * np.cos(self.theta)

        ### Class Attributes

        self.num_span_stations = len(self.__pressure[0, :])
        self.num_chord_stations = len(self.__pressure[:, 0])

        self.__x_coords = []
        self.__z_coords = []

        # Calculate x-point_coords of stations along span
        for i in range(self.num_span_stations):
            self.__x_coords.append(self.__xi(i + 1))

        # Calculate z-point coords of station along chord
        for j in range(self.num_chord_stations):
            self.__z_coords.append(self.__zi(j + 1))

        num_data_pts = self.num_span_stations * self.num_chord_stations
        self.load_data = np.zeros((num_data_pts, 3))

        n = 0
        for i, x in enumerate(self.__x_coords):
            for j, z in enumerate(self.__z_coords):
                self.load_data[n, 0] = x
                self.load_data[n, 1] = z
                self.load_data[n, 2] = self.__pressure[j, i]
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