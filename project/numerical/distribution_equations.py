#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
title: deformations
project: Do228-Assignment
created: 25/02/2020
author: lmaio
"""
import os
import sys

# Imports
import numpy as np

# This must come before the next imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))
from project.numerical.Loading.aileron_geometry import AileronGeometry
from project.numerical.Loading.integration import def_integral
from project.numerical.reaction_forces import reaction_forces
from project.numerical.transformations import CoordinateTransforms


class DistributionEquations():
    def __init__(self):
        self.G = AileronGeometry()
        self.q_x = self.G.q_tilde()
        self.t_x = self.G.tau_tilde()
        self.r = reaction_forces()[0]
        
    def mac(self, x): # Macaulay step
        return np.maximum(0, x)
    
    def zmac(self, x):
        if x > 0:
            return 1
        else:
            return 0
    
    def deflection_y(self, x):
        R1y = self.r[0]
        R2y = self.r[1]
        R3y = self.r[2]
        F   = self.r[6]
        C1  = self.r[7]
        C2  = self.r[8]
    
        terms = np.zeros((6,))
        terms[0] = R1y /6 * self.mac(x - self.G.x1)** 3
        terms[1] = -1* def_integral(self.q_x, 0, x, num_var_integrals=4, num_bins=20)
        terms[2] = F * np.sin(self.G.theta) / 6 * self.mac(x-self.G.x_a_1)**3
        terms[3] = R2y / 6 * self.mac(x-self.G.x2)**3
        terms[4] = -1* self.G.P * np.sin(self.G.theta) / 6 * self.mac(x-self.G.x_a_2)**3
        terms[5] = R3y / 6 * self.mac(x-self.G.x3)**3
    
        v = (-1/(self.G.E * self.G.I_zz) * terms).sum() + C1*x + C2
        return v

    def deflection_z(self, x):
        R1z = self.r[3]
        R2z = self.r[4]
        R3z = self.r[5]
        F = self.r[6]
        C3 = self.r[9]
        C4 = self.r[10]
    
        terms = np.zeros((5,))
        terms[0] = -1* R1z/ 6 * self.mac(x - self.G.x1)**3
        terms[1] = -1*F * np.cos(self.G.theta) / 6 * self.mac(x - self.G.x_a_1)**3
        terms[2] = -1*R2z / 6 * self.mac(x - self.G.x2)**3
        terms[3] = self.G.P*np.cos(self.G.theta) / 6 * self.mac(x - self.G.x_a_2)**3
        terms[4] = -1*R3z / 6 * self.mac(x - self.G.x3)**3
    
        w = (-1/(self.G.E * self.G.I_yy) * terms).sum() + C3*x + C4
        return w

    def _theta_eq_variable_terms(self, x):
        # Terms of the eq excluding the +C5 value
        R1y = self.r[0]
        R2y = self.r[1]
        R3y = self.r[2]
        F = self.r[6]

        terms = np.zeros((8,))
        terms[0] = -1 * def_integral(self.t_x, 0, x, num_var_integrals=2)
        terms[1] = R1y * (self.G.z_h - self.G.z_tilde) * self.mac(x - self.G.x1)
        terms[2] = R2y * (self.G.z_h - self.G.z_tilde) * self.mac(x - self.G.x2)
        terms[3] = R3y * (self.G.z_h - self.G.z_tilde) * self.mac(x - self.G.x3)
        terms[4] = -1 * F * np.cos(self.G.theta) * self.G.y_p * self.mac(x - self.G.x_a_1)
        terms[5] = F * np.sin(self.G.theta) * (0 - self.G.z_tilde) * self.mac(x - self.G.x_a_1)
        terms[6] = self.G.P * np.cos(self.G.theta) * self.G.y_p * self.mac(x - self.G.x_a_2)
        terms[7] = -1 * self.G.P * np.sin(self.G.theta) * (0 - self.G.z_tilde) * self.mac(x - self.G.x_a_2)
        return (1 / (self.G.G * self.G.J) * terms).sum()

    
    def angle_of_twist(self, x):
        C5 = self._theta_eq_variable_terms(self.G.x_a_1)
        theta_rad = self._theta_eq_variable_terms(x) + -C5
        return np.degrees(theta_rad)


    # ---- BENDING MOMENTS --------
    def moment_about_z(self, x):
        R1y = self.r[0]
        R2y = self.r[1]
        R3y = self.r[2]
        F = self.r[6]
    
        terms = np.zeros((6,))
        terms[0] = R1y * self.mac(x - self.G.x1)
        terms[1] = -1* def_integral(self.q_x, 0, x, num_var_integrals=2)
        terms[2] = F* np.sin(self.G.theta) * self.mac(x - self.G.x_a_1)
        terms[3] = R2y * self.mac(x-self.G.x2)
        terms[4] = -1*self.G.P*np.sin(self.G.theta)*self.mac(x-self.G.x_a_2)
        terms[5] = R3y*self.mac(x-self.G.x3)
    
        return terms.sum()
    
    
    def moment_about_y(self, x):
        R1z = self.r[3]
        R2z = self.r[4]
        R3z = self.r[5]
        F = self.r[6]
    
        terms = np.zeros((5,))
        terms[0] = -1*R1z * self.mac(x - self.G.x1)
        terms[1] = -1*F * np.cos(self.G.theta) * self.mac(x - self.G.x_a_1)
        terms[2] = -1*R2z * self.mac(x - self.G.x2)
        terms[3] = self.G.P * np.cos(self.G.theta) * self.mac(x - self.G.x_a_2)
        terms[4] = -1*R3z * self.mac(x - self.G.x3)
    
        return terms.sum()
    
    
    def torsion_x(self, x):
        R1y = self.r[0]
        R2y = self.r[1]
        R3y = self.r[2]
        F = self.r[6]
    
        terms = np.zeros((8,))
        terms[0] = -1*def_integral(self.t_x, 0, x)
        terms[1] = R1y * (self.G.z_h - self.G.z_tilde) * self.zmac(x-self.G.x1)
        terms[2] = R2y * (self.G.z_h - self.G.z_tilde) * self.zmac(x-self.G.x2)
        terms[3] = R3y * (self.G.z_h - self.G.z_tilde) * self.zmac(x-self.G.x3)
        terms[4] = -1*F * np.cos(self.G.theta) * self.G.y_p * self.zmac(x - self.G.x_a_1)
        terms[5] = F * np.sin(self.G.theta) * (0-self.G.z_tilde)*self.zmac(x-self.G.x_a_1)
        terms[6] = self.G.P*np.cos(self.G.theta) * self.G.y_p * self.zmac(x-self.G.x_a_2)
        terms[7] = -1*self.G.P*np.sin(self.G.theta)*(0-self.G.z_tilde)*self.zmac(x-self.G.x_a_2)
    
        return terms.sum()
    
    
    def shear_y(self, x):
        R1y = self.r[0]
        R2y = self.r[1]
        R3y = self.r[2]
        F = self.r[6]
    
        terms = np.zeros((6,))
        terms[0] = R1y * self.zmac(x - self.G.x1)
        terms[1] = -1* def_integral(self.q_x, 0, x)
        terms[2] = F* np.sin(self.G.theta) * self.zmac(x - self.G.x_a_1)
        terms[3] = R2y * self.zmac(x-self.G.x2)
        terms[4] = -1*self.G.P*np.sin(self.G.theta)*self.zmac(x-self.G.x_a_2)
        terms[5] = R3y*self.zmac(x-self.G.x3)
    
        return terms.sum()
    
    
    def shear_z(self, x):
        R1z = self.r[3]
        R2z = self.r[4]
        R3z = self.r[5]
        F = self.r[6]
    
        terms = np.zeros((5,))
        terms[0] = -1*R1z * self.zmac(x - self.G.x1)
        terms[1] = -1*F * np.cos(self.G.theta) * self.zmac(x - self.G.x_a_1)
        terms[2] = -1*R2z * self.zmac(x - self.G.x2)
        terms[3] = self.G.P * np.cos(self.G.theta) * self.zmac(x - self.G.x_a_2)
        terms[4] = -1*R3z * self.zmac(x - self.G.x3)
    
        return -1 * terms.sum()
    
    
    

    def deformations_in_local_coords(self, steps=50):
        min_x, max_x = min(self.G.station_x_coords()), max(self.G.station_x_coords())
    
        deflect_data = np.zeros((steps, 4))
        deflect_data[:, 0] = np.linspace(min_x, max_x, steps)
    
        for i, x in enumerate(deflect_data[:, 0]):
            deflect_data[i, 1] = self.deflection_y(x)
            deflect_data[i, 2] = self.deflection_z(x)
            deflect_data[i, 3] = self.angle_of_twist(x)
    
        return deflect_data
    
    
    def deformations_in_global_coords(self, steps):
        CT = CoordinateTransforms()
        local_deform = self.deformations_in_local_coords(steps=steps)
    
        # Create array of [[x'], [y'], [z']]
        local_lateral = local_deform[:, :-1].T
        global_lateral = (CT.local2global @ local_lateral).T
        # Reapply angle of twist, since that is not affected by transform
        out = np.zeros((global_lateral.shape[0], global_lateral.shape[1] + 1))
        out[:, :-1] = global_lateral
        out[:, -1] = local_deform[:, -1]
        return out
