import sys
import os
import unittest
import numpy as np
from scipy.interpolate import Rbf
from scipy.integrate import quad, dblquad

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports
from project.numerical.Loading.interpolation import InterpolateRBF
from project.numerical.Loading.integration import def_integral, indef_integral

# Group A43 tests

class LoadingTests(unittest.TestCase):

    @staticmethod
    def calc_error(real, approx):
        error = real - approx
        return np.abs(error) / real


    def test_linear_RBF(self):


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
        num_integral = def_integral(rbfi, min(z), max(z))

        error = self.calc_error(ref_integral, num_integral)
        assert ( error < 0.001) # Error must be less than 0.1 %

        ### Cleanup
        del z, d, rbfi, ref_integral, num_integral, error

    def test_double_integral(self):

        ### Test Setup
        # Create rbf interpolation function of random data points

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