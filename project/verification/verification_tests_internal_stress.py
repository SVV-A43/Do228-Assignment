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

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports

class InternalStress( unittest.TestCase ):
    Vy = -36718.115  # verification model
    Vz = 61222.16  # verification model
    My = 19633.84 # verification model
    Mz = -12318.44 # verification model
    def test_qb_spar_vy(self):
        """The base shear flow of vy at both ends should be equal."""
    def test_qb_spar_vz(self):
        """The base shear flow of vz at both ends should be equal."""
    def test_qb_spar_vy(self):
        """The base shear flow of vy at both ends should be equal."""
    def test_qb_spar_vy(self):
        """The base shear flow of vy at both ends should be equal."""
    def test_qb_spar_vy(self):
        """The base shear flow of vy at both ends should be equal."""