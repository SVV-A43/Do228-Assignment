#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
title: verification_tests
project: Do228-Assignment
date: 2/17/2020
author: lmaio
"""

import os
import sys
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import Rbf

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports
from project.numerical.Loading.interpolation import InterpolateRBF
from project.numerical.Loading.integration import def_integral, indef_integral, indef_integral_v2


class LoadingTests(unittest.TestCase):

    @staticmethod
    def calc_error(real, approx):
        error = real - approx
        return np.abs(error) / real

    def test_uni_variate_linear_RBF(self):
        # Create random data to generate interpolation functions on
        x, d = np.random.rand(2, 50)

        # Create instance of our interpolation class
        my_interpolant = InterpolateRBF(x, d)

        # Create reference radial basis function
        ref_interpolant = Rbf(x, d, function='linear')

        # Generate points to test
        xi = np.linspace(0, 1, 200)
        # xi = zi = [0.5]

        di = np.array([ref_interpolant(xi)]).T # interpolated values
        fi = my_interpolant.interpolate(xi)

        assert (fi == di).all()

        ### Cleanup
        del x, d, my_interpolant, ref_interpolant, xi, di, fi


    # Our RBF interpolant is currently developed to handle multiple 1D datasets,
    #       rather than bi-variate data. This would be a potential future update
    #       This test is written for a previous version of the function where that did work
    @unittest.expectedFailure
    def test_bi_variate_linear_RBF(self):
        # Create random data to generate interpolation functions on
        x, z, d = np.random.rand(3, 50)

        # Create instance of our interpolation class
        my_i = InterpolateRBF(x, z, d)

        # Create reference radial basis function
        rbfi = Rbf(x, z, d, function='linear')

        # Generate points to test
        xi = zi = np.linspace(0, 1, 200)
        # xi = zi = [0.5]

        di = rbfi(xi, zi)  # interpolated values
        fi = my_i.interpolate(xi, zi)

        assert (fi == di).all()

        ### Cleanup
        del x, z, d, my_i, rbfi, xi, zi, di, fi


    def test_definite_integral(self):
        ### Test Setup
        # Create rbf interpolation function of random data points
        z, d = np.random.rand(2, 50)
        rbfi = Rbf(z, d, function='linear')

        # Run reference integration
        ref_integral = quad(rbfi, min(z), max(z))[0]

        # Our integral
        num_integral = def_integral(rbfi, min(z), max(z), num_bins=100)

        error = self.calc_error(ref_integral, num_integral[0])
        assert error < 0.01 # Error must be less than 1 %

        ### Cleanup
        del z, d, rbfi, ref_integral, num_integral, error

    def test_poly_definite_integral(self):
        ### Test Setup
        def fn_test1(x):
            return x**2 + 3*x + 1

        start = 0
        end = 2


        # Run reference integration
        ref_integral = quad(fn_test1, start, end)[0]
        ref_sol = 10.666667
        # Our integral
        num_integral = def_integral(fn_test1, start, end, num_bins=100)

        error = self.calc_error(ref_integral, num_integral[0])
        error2 = self.calc_error(ref_sol, num_integral[0])
        assert error < 0.001  # Error must be less than 0.1 %

        ### Cleanup



    def test_double_integral(self):
        ### Test Setup
        # Create single-variable functions
        def fn_test1(x):
            return 3*x

        def fn_test2(x):
            return x**2 + 3*x + 1

        start = 0
        end = 2

        manual_solution1 = 4
        manual_solution2 = 7.3333333333333333

        num_doub_int1 = indef_integral(fn_test1, start, end, num_bins=1000)
        num_doub_int2 = indef_integral(fn_test2, start, end, num_bins=1000)

        error1 = self.calc_error(manual_solution1, num_doub_int1)
        error2 = self.calc_error(manual_solution2, num_doub_int2)

        # Error must be less than 1 %
        assert error1 < 0.01
        assert error2 < 0.01

        # Cleanup
        del start, end, manual_solution2, manual_solution1, num_doub_int1, num_doub_int2, \
            error1, error2


    def test_n_indef_integrals(self):
        def fn_test(x):
            return 3*x

        start = 0
        end = 2

        # n=1 means once variable integral
        ref_sol_n_1 = 4
        ref_sol_n_2 = 2
        ref_sol_n_3 = 0.8

        num_sol_n_1 = indef_integral_v2(fn_test, start, end, num_var_integrals=1, num_bins=100)
        num_sol_n_2 = indef_integral_v2(fn_test, start, end, num_var_integrals=2, num_bins=100)
        num_sol_n_3 = indef_integral_v2(fn_test, start, end, num_var_integrals=3, num_bins=100)


        error1 = self.calc_error(ref_sol_n_1, num_sol_n_1)
        error2 = self.calc_error(ref_sol_n_2, num_sol_n_2)
        error3 = self.calc_error(ref_sol_n_3, num_sol_n_3)

        # Error must be less than 1 %
        assert error1 < 0.01
        assert error2 < 0.01
        assert error3 < 0.01

        # Cleanup
        del start, end, ref_sol_n_1, ref_sol_n_2, ref_sol_n_3,\
            num_sol_n_1, num_sol_n_2, num_sol_n_3, error1, error2, error3

    def test_n_indef_integrals_loops(self):
        def fn_test(x):
            return 3*x

        start = 0
        end = 2

        # n=1 means once variable integral
        ref_sol_n_1 = 4
        ref_sol_n_2 = 2
        ref_sol_n_3 = 0.8

        num_sol_n_1 = indef_integral(fn_test, start, end, num_var_integrals=1, num_bins=100)
        num_sol_n_2 = indef_integral(fn_test, start, end, num_var_integrals=2, num_bins=100)
        num_sol_n_3 = indef_integral(fn_test, start, end, num_var_integrals=3, num_bins=100)


        error1 = self.calc_error(ref_sol_n_1, num_sol_n_1)
        error2 = self.calc_error(ref_sol_n_2, num_sol_n_2)
        error3 = self.calc_error(ref_sol_n_3, num_sol_n_3)

        # Error must be less than 1 %
        assert error1 < 0.01
        assert error2 < 0.01
        assert error3 < 0.01

        # Cleanup
        del start, end, ref_sol_n_1, ref_sol_n_2, ref_sol_n_3,\
            num_sol_n_1, num_sol_n_2, num_sol_n_3, error1, error2, error3


