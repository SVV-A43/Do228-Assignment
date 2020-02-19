import os
import sys
import numpy as np
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next two imports
from project.numerical.interpolation import InterpolationRBF


class InterpolationTests(unittest.TestCase):

    def fn_test(self, x, y):
        return x**1.1 + (-1**y)*5 * y

    def gen_test_data(self, n):
        data = []
        for x in range(n):
            for y in range(n):
                f = self.fn_test(x, y)
                data.append((x, y, f))

        return data

    def test_rbf(self):
        n = 5
        sample_data = self.gen_test_data(n)
        pt = (n/3, n/2)

        gaussRBF = InterpolationRBF(sample_data, basis='linear')

        f_actual = self.fn_test(pt[0], pt[1])
        f_interp = gaussRBF.interpolate(pt)

        error = np.abs(f_actual-f_interp)

        assert error < 0.001

if __name__ == '__main__':
    inter_test = InterpolationTests()
    inter_test.test_rbf()