#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
title: reaction_forces
project: Do228-Assignment
date: 2/23/2020
author: lmaio
"""

import os
import sys
import numpy as np
from tqdm import tqdm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports
from project.numerical.Loading.interpolation import InterpolateRBF
from project.numerical.Loading.integration import def_integral, variable_integral
from project.numerical.Loading.aileron_geometry import AileronGeometry


### Constants:
E = None
P = 20.6 * 10**3    # [N]

# Geometry
G = AileronGeometry()




# Known input variables:
def q_tilde(**kwargs):
    num_bins = kwargs.pop('num_bins', 100)
    aileron = AileronGeometry()
    q_x = []
    x_coords = []

    for station in range(aileron.num_span_stations):
        x, z, p = aileron.station_data(station)
        int_fn = InterpolateRBF(z, p)
        station_load = def_integral(int_fn.interpolate, min(z), max(z), num_bins=num_bins)
        q_x.append(station_load)
        x_coords.append(x[0])
    q_x = np.array(q_x)
    x_coords = np.array(x_coords)

    # Create interpolation function
    q_tilde_x = InterpolateRBF(x_coords, q_x)

    return q_tilde_x








if __name__ == '__main__':
    x, q = q_tilde()

    import matplotlib.pyplot as plt
    plt.scatter(x, q)
    plt.show()
