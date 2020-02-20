import pytest
import numpy as np
import math as m
import unittest
from project.numerical.Section_Properties.cross_section import variable_integrator, value_integrate
class Section_Tests(unittest.TestCase):

    # def test_perimeter(self):
    #     analytical_perimeter = m.pi*(h/2) + 2*m.sqrt((h/2)**2+(C-h/2)**2)
    #     numerical_perimeter = perimeter_check
    #     assert

    def test_variable_int(self):
        func = [[2,4],[1,2]]
        result_my = variable_integrator(func)
        result_analytical = [[1,4/3],[2,3]]
        assert result_my == result_analytical

    def test_value_int(self):
        result_my = value_integrate(0,m.pi/2,[[2],['c']])
        result_analytical = 2
        assert result_my == result_analytical
