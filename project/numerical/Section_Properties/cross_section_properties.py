#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
title: Cross section properties
project: Do228-Assignment
date: 2/26/2020
author: kbislip
"""


# Imports
import sys
import os
import math as m
import numpy as np
import matplotlib.pyplot as plt

# This must come before the next imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))
from project.numerical.Loading.aileron_geometry import AileronGeometry


class Cross_section_properties(AileronGeometry):
    def __init__(self):
        AileronGeometry.__init__(self)

        self.dx = 0.0001

        areas = self.areas()
        self.A_stiff, self.A_total, self.Am1, self.Am2 = areas
        self.x_centroid = self.centroid()

        perimeters = self.perimeters()
        self.l_sk_triangle, self.l_halfcircle, self.l_spar, self.total_perimeter = perimeters

        geo = self.generate_geometry()
        self.x = geo[0]
        self.y = geo[1]
        self.stiff_loc = geo[2]

        self.I_zz = self.moment_of_inertia()[0]
        self.I_yy = self.moment_of_inertia()[1]

        J1 = (4*self.Am1**2)/(((m.pi*self.h/2) / self.t_sk) + (self.h / (0.5 * self.t_sp)))
        J2 = (4*self.Am2**2)/((2 * self.l_sk_triangle / self.t_sk) + (self.h / (0.5 * self.t_sp)))
        self.J = J1+J2

    def generate_geometry(self):
        x = np.arange(0, self.C_a + self.dx, self.dx)
        x = np.append(x, x[-2::-1])
        sz = np.size(x)
        y = np.zeros(sz)
        stiff_loc = np.zeros((self.n_st, 2))
        stiff_counter = 1

        # Stiffener locations and general shape

        perimeter = m.pi*(self.h/2) + 2*m.sqrt((self.h/2) ** 2 + (self.C_a - self.h / 2) ** 2)

        if self.stringers:
            spacing = perimeter / self.n_st

        perimeter_addition = 0
        perimeter_check = 0
        seg_i = 0
        for i in range(sz):
            # Switch segments
            if seg_i != 0 and x[i] == 0.:
                break

            if x[i] == self.h/2 or x[i] == self.C_a:
                s_current = 0
                seg_i += 1

            if x[i] <= self.h/2 and seg_i == 0:
                y[i] = m.sqrt((self.h/2)**2-(x[i]-self.h/2)**2)

            elif seg_i == 1:
                y[i] = -((self.h/2) / (self.C_a - self.h / 2)) * x[i] + (((self.h / 2) / (self.C_a - self.h / 2)) * self.C_a)

            elif seg_i == 2:
                y[i] = -(-((self.h/2) / (self.C_a - self.h / 2)) * x[i] + (((self.h / 2) / (self.C_a - self.h / 2)) * self.C_a))

            elif seg_i == 3:
                y[i] = -m.sqrt((self.h/2)**2-(x[i]-self.h/2)**2)

            perimeter_old = perimeter_addition
            perimeter_addition += m.hypot(self.dx, y[i]-y[i-1])
            perimeter_check += m.hypot(self.dx, y[i]-y[i-1])

            if self.stringers:
                if perimeter_addition > spacing > perimeter_old:
                    if perimeter_addition - spacing <= perimeter_old - spacing:
                        stiff_loc[stiff_counter,0] = x[i]
                        stiff_loc[stiff_counter,1] = y[i]
                    else:
                        stiff_loc[stiff_counter,0] = x[i-1]
                        stiff_loc[stiff_counter,1] = y[i-1]
                    perimeter_addition -= spacing
                    stiff_counter += 1

        if self.stringers:
            for i in range (int((self.n_st - 1) / 2)):
                stiff_loc[self.n_st - (i + 1), 0] = stiff_loc[i + 1, 0]
                stiff_loc[self.n_st - (i + 1), 1] = -stiff_loc[i + 1, 1]

        return x,y,stiff_loc, perimeter_check

    def areas(self):
        a_stiff = self.w_st*self.t_st + self.h_st*self.t_st
        a_total = self.n_st * a_stiff + self.h * self.t_sp + 2 * m.sqrt((self.h / 2) ** 2 + (self.C_a - self.h / 2) ** 2) * self.t_sk + m.pi * self.h / 2 * self.t_sk
        am1 = 0.5* m.pi*(self.h/2)**2
        am2 = 0.5 * (self.C_a - self.h / 2) * self.h
        return a_stiff, a_total, am1, am2

    def perimeters(self):
        l_sk_triangle = m.hypot(self.C_a - self.h / 2, self.h / 2)
        l_halfcircle = (2*m.pi*self.h/2)/2
        l_spar = self.h
        total_perimeter = m.pi*(self.h/2) + 2*m.sqrt((self.h/2) ** 2 + (self.C_a - self.h / 2) ** 2)
        return l_sk_triangle, l_halfcircle, l_spar, total_perimeter

    def centroid(self):

        stiff_cont = 0
        for i in range(self.n_st):
            stiff_cont += self.areas()[0]*self.generate_geometry()[2][i,0]
        spar_cont = self.h * self.t_sp * self.h / 2
        semi_cont = m.pi * self.h / 2 * self.t_sk * (self.h / 2 - self.h / m.pi)
        straight_cont = 2 * (m.sqrt((self.h/2) ** 2 + (self.C_a - self.h / 2) ** 2) * self.t_sk) * (self.h / 2 + (self.C_a - self.h / 2) / 2)

        x_centroid = (stiff_cont+spar_cont+semi_cont+straight_cont)/self.areas()[1]
        return x_centroid

    def moment_of_inertia(self):

        # Ixx

        sinalpha = (self.h/2)/(m.sqrt((self.h/2) ** 2 + (self.C_a - self.h / 2) ** 2))
        stiff_cont_xx = 0
        for i in range (self.n_st):
            stiff_cont_xx += self.A_stiff*(self.stiff_loc[i,1])**2
        spar_cont_xx = (self.t_sp * self.h ** 3) / 12
        plates_cont_xx = self.t_sk * m.sqrt((self.h / 2) ** 2 + (self.C_a - self.h / 2) ** 2) / 12 * (
            m.sqrt((self.h / 2) ** 2 + (self.C_a - self.h / 2) ** 2)) ** 2 * sinalpha ** 2 + m.sqrt(
            (self.h / 2) ** 2 + (self.C_a - self.h / 2) ** 2 ) * self.t_sk * (self.h / 4) ** 2
        semi_cont_xx = m.pi * (self.h/2) ** 3 * self.t_sk / 2

        Izz = stiff_cont_xx + spar_cont_xx + 2*plates_cont_xx + semi_cont_xx

        # Iyy (Still gives wrong value compared to Verification Model)
        #       Mistake most probable in semicircle or plate
        cosalpha = (self.C_a - self.h / 2) / (m.sqrt((self.h / 2) ** 2 + (self.C_a - self.h / 2) ** 2))
        stiff_cont_yy = 0
        for i in range (self.n_st):
            stiff_cont_yy += self.A_stiff*(self.stiff_loc[i,0]-self.x_centroid)**2
        spar_cont_yy = self.t_sp * self.h * (self.x_centroid - self.h / 2) ** 2
        plates_cont_yy = self.t_sk * m.sqrt((self.h / 2) ** 2 + (self.C_a - self.h / 2) ** 2) / 12 * (m.sqrt((self.h / 2) ** 2 + (self.C_a - self.h / 2) ** 2)) ** 2 * cosalpha ** 2 + m.sqrt((self.h / 2) ** 2 + (self.C_a - self.h / 2) ** 2) * self.t_sk * ((self.h / 2 + (self.C_a - self.h / 2) / 2) - self.x_centroid) ** 2
        semi_cont_yy = m.pi * (self.h/2) ** 3 * self.t_sk / 2 + m.pi * self.h / 2 * self.t_sk * ((self.h / 2 - self.h / m.pi) - self.x_centroid) ** 2

        Iyy = stiff_cont_yy + spar_cont_yy + 2*plates_cont_yy + semi_cont_yy

        # Iyy = 5.377416790820396*10**-05 #Correct value, take out later
        return Izz, Iyy





if __name__ == '__main__':
    props = Cross_section_properties()

    print(f'Iyy: {props.I_yy}')
    print(f'Izz: {props.I_zz}')
    print(f'J: {props.J}')


    plot = True
    if plot:
        plt.plot(props.stiff_loc[:,0],props.stiff_loc[:,1], marker="o")
        plt.plot(props.x,props.y)
        plt.plot(props.x_centroid,0, marker="o")
        plt.axis('equal')
        plt.show()

