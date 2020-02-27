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

from tqdm import tqdm
from definitions import AERO_LOADING_DATA_Do228, AERO_LOADING_DATA_B737

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports
from project.numerical.Loading.integration import def_integral
from project.numerical.Loading.interpolation import InterpolateRBF
from project.numerical.Loading.aircraft_configs import DornierDo228, Boeing737





### USE THIS SWITCH TO SET WHICH AIRCRAFT DATA TO USE ###
VALIDATION_MODE = True                          #
#########################################################
if VALIDATION_MODE:
    Aircraft = Boeing737
else:
    Aircraft = DornierDo228

class AileronGeometry(Aircraft):
    def __init__(self, filename=AERO_LOADING_DATA_Do228):
        super(AileronGeometry, self).__init__()


        ## SAME FOR ALL AIRCRAFT
        ### Boundary Conditions
        self.bound_conds = np.zeros((11, 1))
        self.bound_conds[7, 0] = -1 * self.d_1 * np.sin(self.theta)
        self.bound_conds[8, 0] = -1 * self.d_3 * np.sin(self.theta)
        self.bound_conds[9, 0] = self.d_1 * np.cos(self.theta)
        self.bound_conds[10, 0] = self.d_3 * np.cos(self.theta)

        ### Class Attributes
        if self.aircraft == 'boeing_737':
            self.pressure = np.expand_dims(self.pressure, axis=0)

        self.num_span_stations = len(self.pressure[0, :])
        self.num_chord_stations = len(self.pressure[:, 0])

        self.__x_coords = []
        self.__z_coords = []


        if self.aircraft == 'boeing_737':
            # Calculate x-point_coords of stations along span
            self.__x_coords = np.linspace(0, self.l_a, self.num_span_stations)

            # Calculate z-point coords of station along chord
            for j in range(self.num_chord_stations):
                self.__z_coords.append(-0.25 * self.C_a)

        else:
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
                self.load_data[n, 2] = self.pressure[j, i]
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
        '''
        :param id: station id between 0 and 40
        :return: x, z, p coords array
        '''
        s = id * 81 #start idx
        e = (id + 1) * 81 #end idx
        if e > self.load_data.shape[0]:
            e = -1

        # Returns
        return self.load_data[s:e, 0], self.load_data[s:e, 1], self.load_data[s:e, 2]

    def q_tilde(self, **kwargs):
        num_bins = kwargs.pop('num_bins', 100)
        self.q_x = []
        x_coords = []

        if self.aircraft == 'boeing_737':
            def q_tilde(x):
                return np.ones_like(x) * self.pressure[0, 0]
            return q_tilde

        else:
            for station in range(self.num_span_stations):
                x, z, p = self.station_data(station)
                int_fn = InterpolateRBF(z, p)
                station_load = def_integral(int_fn.interpolate, min(z), max(z), num_bins=num_bins)
                self.q_x.append(station_load)
                x_coords.append(x[0])
            self.q_x = np.array(self.q_x)
            x_coords = np.array(x_coords)

        # Create interpolation function
        q_tilde_x = InterpolateRBF(x_coords, self.q_x)

        return q_tilde_x.interpolate


    def tau_tilde(self, **kwargs):
        num_bins = kwargs.pop('num_bins', 100)
        tau_x = []
        x_coords = []

        if self.aircraft == 'boeing_737':
            moment_arm = -0.25*self.C_a - self.z_tilde
            def q_tilde(x):
                tau = self.pressure[0, 0] * moment_arm
                return np.ones_like(x) * tau
            return q_tilde



        for station in range(self.num_span_stations):
            x, z, p = self.station_data(station)

            t = p*(z) # This causes a negative moment

            int_fn = InterpolateRBF(z, t)

            # def tau_inner_fn(z):
            #     return np.multiply(int_fn.interpolate(z), np.squeeze(z - self.z_tilde))
            #     return int_fn.interpolate(z)


            station_load = def_integral(int_fn.interpolate, min(z), max(z), num_bins=num_bins)
            tau_x.append(station_load)
            x_coords.append(x[0])
        self.tau_x = np.array(tau_x)
        x_coords = np.array(x_coords)

        # Create interpolation function
        tau_tilde_x = InterpolateRBF(x_coords, self.tau_x)

        return tau_tilde_x.interpolate



if __name__ == '__main__':
    aileron = AileronGeometry()
    x = aileron.station_x_coords()
    z = aileron.station_z_coords()

    # All data, coordinates repeated
    x_a, z_a, p_a = aileron.data_x_z_p()

    print(aileron.station_x_coords())
    print(aileron.station_z_coords())

