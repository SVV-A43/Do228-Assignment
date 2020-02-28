#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
title: verification_tests_internal_stress
project:
created: 28/02/2020
author: kevin
"""

import os
import sys
import unittest
import math as m

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports
from project.numerical.Section_Properties.cross_section_properties import Cross_section_properties
from project.numerical.Loading.integration import def_integral


class InternalStress( unittest.TestCase ):
    global G
    G = Cross_section_properties()
    @staticmethod
    def calc_error(real, approx):
        error = real - approx
        rel_error = np.abs(error) / real
        if isinstance(rel_error, (float, int)):
            if rel_error == float('inf'):
                rel_error = 0
        elif rel_error.any() == float('inf'):
            rel_error = 0
        return rel_error

    # Tests of cross section properties
    def test_perimeter(self):
        """The perimeter calculated analytically should equal the perimeter calculated by the numerical model."""
        analytical_perimeter = m.pi*(G.h/2) + 2*m.sqrt((G.h/2) ** 2 + (G.C_a - G.h / 2) ** 2)
        num_perimeter = G.generate_geometry()[3]
        error = self.calc_error(analytical_perimeter, num_perimeter)
        assert error < 0.01

    def test_total_area_ver_model(self):
        ver_model = -0.0024705343591563712 # total area computed by verification model
        num_model = -G.A_total # total area computed by numerical model
        error = self.calc_error(ver_model, num_model)
        assert error < 0.01


    def test_centroid_ver_model(self):
        ver_model = -0.20728702965108006                 # z-coordinate of the centroid
        num_model = -G.x_centroid # Centroid location computed by numerical model
        error = self.calc_error(ver_model, num_model)
        assert error < 0.01

    def test_Iyy_ver_model(self):
        ver_model = 5.377416790820396e-05
        num_model = G.I_yy
        error = self.calc_error(ver_model, num_model)
        assert error < 0.05

    def test_Izz_ver_model(self):
        ver_model = 1.4221538884296291e-05
        num_model = G.I_zz
        error = self.calc_error(ver_model, num_model)
        assert error < 0.05

    def test_qb_vy_spar_analytical(self):
        """taking all values of the shear flow eqaution to be 1 and setting h of spar to 1 m,
        the value at the spar ends can easily be calculated analytically
        and compared to the value computed numerically"""
        analytical_qb_vy_top = -1 #N/m
        analytical_qb_vy_bottom = -1 #N/m
        h = 1  # Define h to be 2
        t_spar = 1
        I_zz = 1
        vy_xi = 1
        dx = 0.0001
        def qb_2_vy(s):
            return -vy_xi*t_spar*s/I_zz
        def qb_5_vy(s):
            return -vy_xi*t_spar*s/I_zz

        """Spar shear flow distribution"""
        y_spar = np.arange(-h, h, dx)
        sz_spar = np.size(y_spar)
        qb_spar_val_list = np.zeros(sz_spar)  # runs from -h to h
        qb_lastval = np.zeros(6) # last values of qb1, qb2, qb3, qb4, qb5 and qb6
        for i in range(sz_spar):
            s_current = y_spar[i]
            if y_spar[i] < 0:
                s_0 = 0
                qb_current = def_integral(qb_5_vy, s_0, s_current, num_var_integrals=1)
            else:
                s_0 = 0
                qb_current = def_integral(qb_2_vy, s_0, s_current, num_var_integrals=1)

            qb_spar_val_list[i] = qb_current

        qb_lastval[1] = qb_spar_val_list[-1]
        qb_lastval[4] = qb_spar_val_list[0]

        error1 = self.calc_error(analytical_qb_vy_bottom, qb_lastval[4])
        error2 = self.calc_error(analytical_qb_vy_top, qb_lastval[1])
        assert error1 < 0.05
        assert error2 < 0.05

