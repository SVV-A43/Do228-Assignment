#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
title: reaction_forces
project: Do228-Assignment
date: 2/20/2020
author: lmaio
"""

import os
import sys

import numpy as np


def def_integral(fn, start, stop, num_var_integrals=1, **kwargs):
    num_bins = kwargs.pop('num_bins', 100)

    if isinstance(start, (float, int)):
        start_ar = np.array([start])
    else:
        start_ar = start
    if isinstance(stop, (float, int)):
        stop_ar = np.array([stop])
    else:
        stop_ar = stop

    del start, stop

    # New axis is added to the front (axis0)
    steps, step_width = np.linspace(start_ar, stop_ar, num_bins, axis=0, retstep=True) # Need to find integral of each column (axis0)

    s_diff = np.diff(steps, axis=0)
    width = np.take(s_diff, 0, axis=0)  # Take slice along axis 0 of the widths
    if width.ndim == 1:
        width = np.array([width]) # Make 2D array

    if steps.ndim > 2:
        steps = np.squeeze(steps)  # Remove extra dimensions with shape 1

    if num_var_integrals > 1:
        data = def_integral(fn, start_ar, steps, num_var_integrals=num_var_integrals - 1, num_bins=num_bins)
    else:
        data = fn(steps)

    # areas along axis 0
    def trapz(vals):
        return (vals[:-1] + vals[1:]) / 2

    trap_h = np.apply_along_axis(trapz, 0, data)    # Calc trapezoid heights
    areas = trap_h * width.T                        # Calc trapezoid areas
    areas = np.squeeze(areas)                       # Reduce dimensions to minimum (needed for next line)
    return areas.sum(axis=0)




def main():
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))  # This must come before the next imports
    from project.numerical.Loading.aileron_geometry import AileronGeometry
    from project.numerical.Loading.interpolation import InterpolateRBF

    ail = AileronGeometry()
    # min_z = min(ail.station_z_coords())
    # max_z = max(ail.station_z_coords())
    #
    station_id = 5
    _, z, p = ail.station_data(station_id)
    # xi = [0 for i in range(len(x))]

    int_fn = InterpolateRBF(z, p)

    second_int = def_integral(int_fn.interpolate, 0, -0.5, num_bins=10)
    print(second_int)

if __name__ == '__main__':
    main()
