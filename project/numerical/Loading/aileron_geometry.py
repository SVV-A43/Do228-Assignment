#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
title: reaction_forces
project: Do228-Assignment
date: 2/17/2020
author: lmaio
"""


# Imports
import numpy as np


# CODE...


class AileronGeometry():
    # Functions to calculate geometry
    def __th_xi(self, i):
        return (i - 1) / self.span_stations * np.pi

    def __xi(self, i):
        return .5 * (self.__la / 2 * (1 - np.cos(self.__th_xi(i))) + self.__la / 2 * (1 - np.cos(self.__th_xi(i + 1))))

    def __th_zi(self, i):
        return (i - 1) / self.chord_stations * np.pi

    def __zi(self, i):
        return -.5 * (self.__Ca / 2 * (1 - np.cos(self.__th_zi(i))) + self.__Ca / 2 * (1 - np.cos(self.__th_zi(i + 1))))





    def __init__(self, filename='./aero_loading_data/aerodynamicloaddo228.dat'):
        self.__Ca = 0.515
        self.__la = 2.961

        self.__pressure = np.genfromtxt(filename, delimiter=',')
        self.span_stations = len(self.__pressure[0, :])
        self.chord_stations = len(self.__pressure[:, 0])

        self.__x_coords = []
        self.__z_coords = []


        # Calculate x-point_coords of stations along span
        for i in range(self.span_stations):
            self.__x_coords.append(self.__xi(i+1))

        # Calculate z-point coords of station along chord
        for j in range(self.chord_stations):
            self.__z_coords.append(self.__zi(j+1))

        num_data_pts = self.span_stations * self.chord_stations
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
        self.__load_coords = np.empty((self.chord_stations, self.span_stations), dtype=object)
        for z in range(self.chord_stations):
            for x in range(self.span_stations):
                self.__load_coords[z, x] = (self.__z_coords[z], self.__x_coords[x])

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



if __name__ == '__main__':
    aileron = AileronGeometry()
    x = aileron.station_x_coords()
    z = aileron.station_z_coords()
    print(aileron.station_x_coords())
    print(aileron.station_z_coords())
    import matplotlib.pyplot as plt
    plt.plot(x,z)
    plt.show()
