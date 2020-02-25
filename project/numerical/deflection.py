#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
title: deflection
project: 
created: 25/02/2020
author: lmaio
"""
# Main Structural Analysis model we need to develop

# Imports
import numpy as np
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports
from project.numerical.Loading.aileron_geometry import AileronGeometry
from project.numerical.Loading.integration import def_integral
from project.numerical.reaction_forces import reaction_forces


# CODE...


# Equations of motion


def deflection_y(x):
    G = AileronGeometry()
    r, r_names = reaction_forces()

    R1y = r[0]
    R2y = r[1]
    R3y = r[2]
    F   = r[6]
    C1  = r[7]
    C2  = r[8]

    terms = np.zeros((6,))

    if (x-G.x1) > 0:
        terms[0] = R1y * (x-G.x1)**3 / 6

    terms[1] = -1* def_integral(G.q_tilde().interpolate, )




