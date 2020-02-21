import sys
import os
import unittest
import numpy as np
from scipy.interpolate import Rbf
from scipy.integrate import quad

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports
from project.numerical.Loading.interpolation import InterpolateRBF
from project.numerical.Loading.integration import def_integral

# Group A43 tests

class LoadingTests(unittest.TestCase):

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

    def test_definite_integral(self):

        ### Test Setup
        # Create rbf interpolation function of random data points
        z, d = np.random.rand(2, 50)
        rbfi = Rbf(z, d, function='linear')

        # Run reference integration
        ref_integral = quad(rbfi, min(z), max(z))[0]

        # Our integral
        num_integral = def_integral(rbfi, min(z), max(z))

        error = ref_integral-num_integral
        er_percent = np.abs(error) / ref_integral * 100

        assert (er_percent < 0.1) # Error must be less than 0.1 %
