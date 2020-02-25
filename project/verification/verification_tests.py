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
from project.numerical.Loading.integration import def_integral, def_integral
from project.numerical.reaction_forces import equilibrium_eq_coefficients, equilibrium_eq_resultants
from project.numerical.Loading.aileron_geometry import AileronGeometry


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

        # Create reference radial basis functiond
        ref_interpolant = Rbf(x, d, function='linear')

        # Generate points to test
        xi = np.linspace(0, 1, 200)
        # xi = zi = [0.5]

        di = ref_interpolant(xi) # interpolated values
        fi = my_interpolant.interpolate(xi)

        error = self.calc_error(di, fi)
        assert max(error) < 0.01

        ### Cleanup
        del x, d, my_interpolant, ref_interpolant, xi, di, fi


    # Our RBF Interpolant is currently developed to handle ONLY multiple 1D datasets,
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
        interpolator = Rbf(z, d, function='linear')

        # Run reference integration
        ref_integral = quad(interpolator, min(z), max(z))[0]

        # Our integral
        num_integral = def_integral(interpolator, min(z), max(z), num_bins=1000)

        error = self.calc_error(ref_integral, num_integral)
        assert error < 0.01 # Error must be less than 1 %

        ### Cleanup
        del z, d, interpolator, ref_integral, num_integral, error

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

        error = self.calc_error(ref_integral, num_integral)
        error2 = self.calc_error(ref_sol, num_integral)
        assert error < 0.001  # Error must be less than 0.1 %
        assert error2 < 0.001

        ### Cleanup
        del start, end, ref_integral, ref_sol, num_integral, error, error2



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

        num_doub_int1 = def_integral(fn_test1, start, end, num_var_integrals=2, num_bins=100)
        num_doub_int2 = def_integral(fn_test2, start, end, num_var_integrals=2, num_bins=100)

        error1 = self.calc_error(manual_solution1, num_doub_int1)
        error2 = self.calc_error(manual_solution2, num_doub_int2)

        # Error must be less than 1 %
        assert error1 < 0.01
        assert error2 < 0.01

        # Cleanup
        del start, end, manual_solution2, manual_solution1, num_doub_int1, num_doub_int2, \
            error1, error2


    def test_n_variable_integrals(self):
        def fn_test(x):
            return 3*x

        start = 0
        end = 2

        # n=1 means once variable integral
        ref_sol_n_1 = 6     # 3/2 * x**2
        ref_sol_n_2 = 4     # 1/2 * x**3
        ref_sol_n_3 = 2     # 1/8 * x**4
        ref_sol_n_4 = 0.8   # 1/40 * x**5
        ref_sol_n_5 = 0.2666666      # 1/240 * x**6

        num_sol_n_1 = def_integral(fn_test, start, end, num_var_integrals=1, num_bins=100)
        num_sol_n_2 = def_integral(fn_test, start, end, num_var_integrals=2, num_bins=100)
        num_sol_n_3 = def_integral(fn_test, start, end, num_var_integrals=3, num_bins=100)
        num_sol_n_4 = def_integral(fn_test, start, end, num_var_integrals=4, num_bins=40)
        # Lower resolution due to memory restrictions (LONG RUNTIME > 1min):
        num_sol_n_5 = def_integral(fn_test, start, end, num_var_integrals=5, num_bins=40)


        error1 = self.calc_error(ref_sol_n_1, num_sol_n_1)
        error2 = self.calc_error(ref_sol_n_2, num_sol_n_2)
        error3 = self.calc_error(ref_sol_n_3, num_sol_n_3)
        error4 = self.calc_error(ref_sol_n_4, num_sol_n_4)
        error5 = self.calc_error(ref_sol_n_5, num_sol_n_5)

        # Error must be less than 1 %
        assert error1 < 0.01
        assert error2 < 0.01
        assert error3 < 0.01
        assert error4 < 0.01
        assert error5 < 0.01

        # Cleanup
        del start, end, ref_sol_n_1, ref_sol_n_2, ref_sol_n_3, ref_sol_n_4, ref_sol_n_5, \
            num_sol_n_1, num_sol_n_2, num_sol_n_3, num_sol_n_4, num_sol_n_5, \
            error1, error2, error3, error4, error5

    def test_equilibrium_coefficients(self):
        A = equilibrium_eq_coefficients()

        ref_A_0_1 = -1.64
        ref_A_10_1 = -4.99966 * 10**-7

        A_0_1 = A[0, 1]
        A_10_1 = A[10, 1]

        error1 = self.calc_error(ref_A_0_1, A_0_1)
        error2 = self.calc_error(ref_A_10_1, A_10_1)

        assert error1 < 0.01
        assert error2 < 0.01


    def test_equilibrium_resultants(self):
        G = AileronGeometry()
        b = equilibrium_eq_resultants()

        ref_b_4 = 18669
        ref_b_8 = -1* G.P*np.cos(G.theta)*(G.x3 - G.x_a_2)**3 / (6*G.E*G.I_yy)
        ref_b_8_manual = -0.001180

        b_4 = b[4, 0]
        b_8 = b[8, 0]

        error1 = self.calc_error(ref_b_4, b_4)
        error2 = self.calc_error(ref_b_8, b_8)

        assert error1 < 0.01
        assert error2 < 0.01