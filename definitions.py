#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
title: definitions
project: 
date: 24/02/2020
author: lmaio
"""

import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

AERO_LOADING_DATA_Do228 = os.path.join(ROOT_DIR, 'project', 'numerical', 'Loading',
                                 'aero_loading_data', 'aerodynamicloaddo228.dat')
AERO_LOADING_DATA_B737 = os.path.join(ROOT_DIR, 'project', 'numerical', 'Loading',
                                 'aero_loading_data', 'aero_loading_737.dat')
