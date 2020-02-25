#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
title: reaction_forces
project: Do228-Assignment
date: 2/20/2020
author: lmaio
"""

import os
import sys

import numpy as np

# This must come before the next imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))
from project.numerical.Loading.aileron_geometry import AileronGeometry


class InterpolateRBF():
    def __init__(self, *data_arrays, basis='linear', coeff_path=None,
                 save_path=None):
        '''
        :param X: 1D array of x-coordinates
        :param Z: 1D array of z-coordinates
        :param F: 1D array of data at known data points
        :param basis: 'linear'
        '''
        self._basis = basis
        self.__save_path = save_path
        self.phi_x = None

        # Save initial data
        coords = []
        for ar in data_arrays[:-1]:
            coords.append(np.asarray(ar))
        self._known_coords = np.array(coords)
        self._known_data = np.asarray(data_arrays[-1])

        # Compute coefficients:
        if coeff_path is None:
            self._compute_coeffs()
        else:
            self.coefficients = np.genfromtxt(coeff_path)


    def _dist_r(self, p1, p2):
        '''
        :param x1: point 1 [x1, z1]
        :param x2: point 2 [x2, z2]
        :return:
        '''
        return np.sqrt(((p2 - p1)**2).sum())


    def _phi(self, r):
        '''
        Radial basis of euclidean distance
        :param r: euclidean distance
        :param basis: 'linear'
        :return:
        '''
        if self._basis == 'linear':
            return r
        # Only linear RBF implemented for now
        raise NotImplementedError

    def _compute_coeffs(self):
        # self._known_coords = np.array([np.asarray(X), np.asarray(Z)])
        dim_n = len(self._known_coords[0, :])

        mat_A = np.zeros((dim_n, dim_n))

        for i in range(dim_n):
            for j in range(dim_n):
                r = self._dist_r(self._known_coords[:, i],
                                 self._known_coords[:, j])
                mat_A[i][j] = self._phi(r)

        self.coefficients = np.array([np.linalg.solve(mat_A, self._known_data)])
        if self.__save_path:
            np.savetxt(self.__save_path, self.coefficients, delimiter=',')


    def interpolate(self, input_coords):
        '''ONE VARIABLE ONLY'''
        if isinstance(input_coords, (float, int)):
            input_coords = [input_coords]
        coords = np.asarray(input_coords)

        if coords.ndim == 1:
            coords = np.array([coords]).T # Input coords along axis0

        def euclidean_norm_1d(x): # Calculate distance to known coords, using only first axis
            return np.sqrt((x.T - self._known_coords.T)**2)

        # Calc matrix phi, using only first axis (most recent added by integral fn)
        phi_r = np.apply_along_axis(euclidean_norm_1d, 0, coords)

        # Calculate dot product using the first axis of both matrices
        interp_pts = np.tensordot(np.squeeze(self.coefficients), np.squeeze(phi_r), axes=([0],[0]))

        return interp_pts


def select_station(station_id):
    '''
    :param station_id: station number
    :return:
    '''
    aileron = AileronGeometry()
    x, z, p = aileron.data_x_z_p()

    start = station_id * 81
    end = (station_id + 1) * 81
    if end > p.shape[0]:
        end = -1

    return x[start:end], z[start:end], p[start:end]


def main():
    x, d = np.random.rand(2, 50)
    my_inter = InterpolateRBF(x, d)

    xi = 0.5
    di = my_inter.interpolate(xi)
    print(di)


if __name__ == '__main__':
    main()
