#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
title: aircraft_load_cases
project: 
created: 27/02/2020
author: lmaio
"""
# Imports
import sys
import os
import numpy as np

from definitions import AERO_LOADING_DATA_Do228, AERO_LOADING_DATA_B737

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports


class DornierDo228():
    def __init__(self, filename=AERO_LOADING_DATA_Do228):
        self.aircraft = 'dornier_do228'
        self.stringers = True
        self.C_a = 0.515  # [m]
        self.l_a = 2.691  # [m]
        self.x1 = 0.174  # [m]
        self.x2 = 1.051  # [m]
        self.x3 = 2.512  # [m]
        self.x_a = 30.0 / 100  # [m]
        self.x_a_1 = self.x2 - self.x_a / 2  # [m]
        self.x_a_2 = self.x2 + self.x_a / 2  # [m]
        self.z_tilde = -0.35375  # [m] distance to shear center #FIXME: THIS IS FROM VERIFICATION MODEL
        self.h = 24.8 / 100  # [m]
        self.d_1 = 1.034 / 100  # [m]
        self.d_3 = 2.066 / 100  # [m]
        self.z_h = -1 * self.h / 2  # [m]
        self.y_p = self.h / 2  # [m]
        self.theta = np.deg2rad(25)  # [rad]

        self.t_sk = 1.1 / 1000 # [m]
        self.t_sp = 2.2 / 1000 # [m]
        self.t_st = 1.2 / 1000 # [m]
        self.h_st = 1.5 / 100 # [m]
        self.w_st = 3.0 / 100 # [m]
        self.n_st = 11  # [-]

        self.I_xx = None
        self.I_yy = 5.643650631210155e-05  # [m^4]
        self.I_zz = 1.4221372629975417e-05  # [m^4]
        self.J = 0.00020500801555445113
        self.J_polar = self.I_yy + self.I_zz  # [m^4]

        ### Constants:
        self.E = 73.1 * 10 ** 9  # [Pa]
        self.G = 28 * 10 ** 9  # [Pa]

        ### Loads:
        self.P = 20.6 * 10 ** 3  # [N]
        self.pressure = np.genfromtxt(filename, delimiter=',') * 10 ** 3


# Used for validation
class Boeing737():
    def __init__(self, filename=AERO_LOADING_DATA_B737):
        self.aircraft = 'boeing_737'
        self.stringers = False
        self.C_a = 0.605  # [m]
        self.l_a = 2.661  # [m]
        self.x1 = 0.172  # [m]
        self.x2 = 1.211  # [m]
        self.x3 = 2.591  # [m]
        self.x_a = 35.0 / 100  # [m]
        self.x_a_1 = self.x2 - self.x_a / 2  # [m]
        self.x_a_2 = self.x2 + self.x_a / 2  # [m]
        self.z_tilde = -0.35375  # [m] distance to shear center #FIXME: THIS IS FROM VERIFICATION MODEL
        self.h = 20.5 / 100  # [m]
        self.d_1 = 1.154 / 100  # [m]
        self.d_3 = 1.840 / 100  # [m]
        self.z_h = -1 * self.h / 2  # [m]
        self.y_p = self.h / 2  # [m]
        self.theta = np.deg2rad(28)  # [rad]

        self.t_sk = 1.1 /1000  # [mm]
        self.t_sp = 2.2 /1000 # [mm]
        self.t_st = 1.2 /1000 # [mm]
        self.h_st = 1.5 / 100  # [mm]
        self.w_st = 3.0 / 100  # [mm]
        self.n_st = 0  # [-]

        self.I_xx = None
        self.I_yy = 6.327309776449187e-05  # [m^4]
        self.I_zz = 7.391448868699719e-06  # [m^4]

        self.J = 1.1758542029324649e-05
        self.J_polar = self.I_yy + self.I_zz  # [m^4]

        ### Constants:
        self.E = 73.1 * 10 ** 9  # [Pa]
        self.G = 28 * 10 ** 9  # [Pa]

        ### Loads:
        self.P = 20.6 * 10 ** 3  # [N]
        self.pressure = np.genfromtxt(filename, delimiter=',') * 10 ** 3
