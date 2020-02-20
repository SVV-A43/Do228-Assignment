import sys
import os
import unittest
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports
from project.numerical.radial_basis import InterpolateRBF

# Group A43 tests

class InterpolationTests(unittest.TestCase):

    def test_linear_RBF(self):
        from scipy.interpolate import Rbf

        # Create random data to generate interpolation functions on
        x, z, d = np.random.rand(3, 50)

        # Create instance of our interpolation class
        my_i = InterpolateRBF(x, z, d)

        # Create reference radial basis function
        rbfi = Rbf(x, z, d, function='linear')

        # Generate points to test
        xi = zi = np.linspace(0, 1, 20)
        # xi = zi = [0.5]

        di = rbfi(xi, zi)  # interpolated values
        fi = my_i.interpolate(xi, zi)

        assert (fi == di).all()
