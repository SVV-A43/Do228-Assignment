import math as m
import numpy as np
import matplotlib.pyplot as plt

import os
import sys

# This must come before the next imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))
from project.numerical.Loading.integration import def_integral
from project.numerical.reaction_forces import AileronGeometry
from project.numerical.Section_Properties.cross_section_properties import Cross_section_properties
from project.numerical.distribution_equations import DistributionEquations


# Initialize parameters

D = Cross_section_properties()
G = AileronGeometry()
E = DistributionEquations()

#############  OUR MODEL ################
sample_steps = 50
min_x, max_x = min(G.station_x_coords()), max(G.station_x_coords())
x_steps = np.linspace(min_x, max_x, sample_steps)[:1]

def Vy(xi):
    return E.shear_y(xi)
def Vz(xi):
    return E.shear_z(xi)
def My(xi):
    return E.moment_about_y(xi)
def Mz(xi):
    return E.moment_about_z(xi)
def Tx(xi):
    return E.torsion_x(xi)

sigma_vm_all = [0]*sample_steps

spar_segment = 0
for xi in x_steps:
    #########################################

    # Vy(xi) = -36718.115  # verification model
    # Vz(xi) = 61222.16  # verification model
    # My(xi) = 19633.84 # verification model
    # Mz(xi) = -12318.44 # verification model
    eta = G.z_tilde

C = D.C_a
h = D.h
n = D.n_st
w_st = D.w_st
h_st = D.h_st
t_st = D.t_st
t_skin = D.t_sk
t_spar = D.t_sp
A_stiff = D.A_stiff
l_sk = D.l_sk_triangle
I_zz = D.I_zz
I_yy = D.I_yy
x_centroid = D.x_centroid
dx = D.dx

    def shear_flow_Vy():

        h = D.h/2  # Define h to be h/2
        def qb_1_vy(s):
            return -Vy(xi)*t_skin*h*h*np.sin(s)/I_zz
        def qb_6_vy(s):
            return -Vy(xi)*t_skin*h*h*np.sin(s)/I_zz
        def qb_2_vy(s):
            return -Vy(xi)*t_spar*s/I_zz
        def qb_5_vy(s):
            return -Vy(xi)*t_spar*s/I_zz
        def qb_3_vy(s):
            return -Vy(xi)*t_skin*(h-h/l_sk)*s/I_zz
        def qb_4_vy(s):
            return -Vy(xi)*t_skin*(-h/l_sk)*s/I_zz

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

        """Outer segment shear flow distribution
        Order: 1,3,4,6
        """
        # Initialize simulation parameters
        x = np.arange(0, C+dx,dx)
        x = np.append(x, x[-2::-1])
        sz = np.size(x)
        y = np.zeros(sz)
        stiff_counter = 1
        qb_current = 0  # start at cut
        qb_stiff_current = 0
        qb_val_list = np.zeros(sz)

        perimeter = D.total_perimeter
        s_list = [m.pi/2, l_sk, l_sk, m.pi/2]  #list containing lengths of each skin 1,3,4,6
        i_list = []
        spacing = perimeter/n
        perimeter_addition = 0
        s_current = 0
        seg_i = 0

        for i in range(sz):  # go through the outer section

            # End
            if seg_i != 0 and x[i] == 0.:
                qb_lastval[5] = qb_current
                qb_val_list[i] = qb_current
                # print("Finished at segment {}, i={} with x,y {},{} and qb of {}".format(seg_i,i, x[i],y[i],qb_current))
                max_val = np.amax(qb_val_list)
                min_val = np.amin(qb_val_list)
                if abs(min_val) > max_val:
                    max_val = min_val
                max_val_i = np.where(qb_val_list == max_val)
                # print("Max qb value of {} found at x,y = {},{}".format(qb_val_list[max_val_i], x[max_val_i], y[max_val_i]))
                break

            # Switch segments
            if x[i] == h/2 or x[i] == C:
                s_current = 0
                if seg_i == 0:
                    qb_lastval[0] = qb_current
                elif seg_i == 1:
                    qb_lastval[2] = qb_current
                elif seg_i == 2:
                    qb_lastval[3] = qb_current

                seg_i += 1
                i_list.append(i)
                # print("Finished at segment {}, i={} with x,y {},{} and qb of {}".format(seg_i,i, x[i],y[i],qb_current))

            if seg_i == 0:
                y[i] = m.sqrt((h/2)**2-(x[i]-h/2)**2)
                split_rad = (h/2/dx)
                qb_current = def_integral(qb_1_vy, 0, (m.pi/2)*i/split_rad, num_var_integrals=1)

            elif seg_i == 1:
                y[i] = -((h/2)/(C-h/2))*x[i] + (((h/2)/(C-h/2))*C)
                qb_current = def_integral(qb_3_vy, 0, s_current, num_var_integrals=1) + qb_lastval[0] + qb_lastval[1]

            elif seg_i == 2:
                y[i] = -(-((h/2)/(C-h/2))*x[i] + (((h/2)/(C-h/2))*C))
                qb_current = def_integral(qb_4_vy, 0, s_current, num_var_integrals=1) + qb_lastval[2]

            elif seg_i == 3:
                y[i] = -m.sqrt((h/2)**2-(x[i]-h/2)**2)
                s_0 = -m.pi/2
                qb_current = def_integral(qb_6_vy, s_0, (m.pi/2)*(i-i_list[-1])/split_rad+s_0, num_var_integrals=1) + qb_lastval[3] - qb_lastval[4]

            #Stiffener addition to qb
            perimeter_old = perimeter_addition
            perimeter_addition += m.hypot(dx, y[i]-y[i-1])
            if perimeter_addition > spacing > perimeter_old:
                if perimeter_addition - spacing <= perimeter_old - spacing:
                    qb_stiff_current += (A_stiff * (x[i] + y[i]) )/I_zz
                else:
                    qb_stiff_current += (A_stiff * (x[i-1] + y[i-1]) )/I_zz
                perimeter_addition -= spacing
                stiff_counter += 1
                # qb_current += qb_stiff_current
            # print("Stiffener {} causes qb_stiff_current {}".format(stiff_counter,qb_stiff_current))
            qb_current += qb_stiff_current
            qb_val_list[i] = qb_current
            s_current += m.hypot(dx, y[i]-y[i-1])
        return qb_lastval, qb_val_list, qb_spar_val_list

    qb_lastval_vy, qbvy_outer_val_list, qbvy_spar_val_list = shear_flow_Vy()

    """Calculate shear center and redundant shear flows"""

    def shear_center():
        def qb_1_vy(s):
            return -Vy(xi)*t_skin*h*h*np.sin(s)/I_zz
        def qb_6_vy(s):
            return -Vy(xi)*t_skin*h*h*np.sin(s)/I_zz
        def qb_2_vy(s):
            return -Vy(xi)*t_spar*s/I_zz
        def qb_5_vy(s):
            return -Vy(xi)*t_spar*s/I_zz
        def qb_3_vy(s):
            return -Vy(xi)*t_skin*(h-h/l_sk)*s/I_zz
        def qb_4_vy(s):
            return -Vy(xi)*t_skin*(-h/l_sk)*s/I_zz
        qb_1 = qb_lastval_vy[0]
        qb_2 = qb_lastval_vy[1]
        qb_3 = qb_lastval_vy[2]
        qb_4 = qb_lastval_vy[3]
        qb_5 = qb_lastval_vy[4]
        qb_6 = qb_lastval_vy[5]
        # semicircle
        h = D.h
        p1 = m.pi*h/2 + h
        qs_0_1 = -(def_integral(qb_1_vy, 0, m.pi/2, num_var_integrals=2)*m.pi*h/4 -
                   def_integral(qb_2_vy, 0, h/2, num_var_integrals=2)*h/2 - def_integral(qb_5_vy, 0, -h/2, num_var_integrals=2)*h/2 +
                   def_integral(qb_6_vy, -m.pi/2, 0, num_var_integrals=2)*m.pi*h/4)/p1

        #triangle
        p2 = h + 2*l_sk
        qs_0_2 = -(def_integral(qb_3_vy, 0, l_sk, num_var_integrals=2)*l_sk +
                   def_integral(qb_4_vy, 0, l_sk, num_var_integrals=2)*l_sk + def_integral(qb_2_vy, 0, m.pi/2, num_var_integrals=2)*h/2 +
                   def_integral(qb_5_vy, 0, m.pi/2, num_var_integrals=2)*h/2)/p2

        """Shear flows are then added to the base shear flows to calculate shear center"""
        qs_0_1, qs_0_2 = 0,0 # Neglect redundant shear flow as calculation is faulty


        qb_1 = qb_1 + qs_0_1
        qb_2 = qb_2 + qs_0_2 - qs_0_1
        qb_3 = qb_3 + qs_0_2
        qb_4 = qb_4 + qs_0_2
        qb_5 = qb_5 + qs_0_2 - qs_0_1
        qb_6 = qb_6 + qs_0_1

        theta = m.asin((h/2)/l_sk)

        qb_1_h = qb_1 / 2
        qb_1_v = qb_1 / 2
        qb_2_h = 0
        qb_2_v = qb_2
        qb_3_h = qb_3 * m.cos(theta)
        qb_3_v = qb_3 * m.sin(theta)
        qb_4_h = qb_4 * m.cos(theta)
        qb_4_v = qb_4 * m.sin(theta)
        qb_5_h = 0
        qb_5_v = qb_5
        qb_6_h = qb_6 / 2
        qb_6_v = qb_6 / 2

        #moment arm from trailing edge [vertical distance, horizontal distance] 1,3,4,6
        #for 2 and 5 = [0,eta],[0,eta]
        moment_arm_list = [[0,C+h/2],[h/2,C-h/2],[0,0],[h/2,C-h/2]]

        m_b_1 = (h/2) * qb_1_v
        m_b_2 = 0
        m_b_3 = h/2 * qb_3_h
        m_b_4 = (C-h/2) * qb_3_v
        m_b_5 = 0
        m_b_6 = h/2 * qb_6_h
        m_b = m_b_1 + m_b_2 + m_b_3 + m_b_4 + m_b_5 + m_b_6

        eta = -m_b
        eta = eta - h/2
        # print("Shear center location from leading edge",eta)
        return eta, qs_0_1, qs_0_2

    eta, qs_0_1, qs_0_2 = shear_center()
    # print(eta, qs_0_1, qs_0_2)
    def shear_flow_Vz():
        x_centroid = -D.x_centroid
        h = D.h/2  # Define h to be h/2
        def qb_1_vz(s):
            return -Vz(xi) * t_skin * h * (-(1 - np.cos(s)) * h - x_centroid) / I_yy
        def qb_6_vz(s):
            return -Vz(xi) * t_skin * h * (-(1 - np.cos(s)) * h - x_centroid) / I_yy
        def qb_2_vz(s):
            val = -t_spar * (-h - x_centroid) / I_yy
            return np.ones_like(s) * val
        def qb_5_vz(s):
            val = -t_spar * (-h - x_centroid) / I_yy
            return np.ones_like(s) * val
        def qb_3_vz(s):
            return -Vz(xi) * t_skin * ((-h - x_centroid) - (((C - h) * s) / l_sk)) / I_yy
        def qb_4_vz(s):
            return -Vz(xi) * t_skin * ((-C - x_centroid) - ((C - h) * s)) / I_yy

        """Spar shear flow distribution"""
        y_spar = np.arange(-h, h, dx)
        sz_spar = np.size(y_spar)
        qb_spar_val_list = np.zeros(sz_spar)  # runs from -h to h
        qb_lastval = np.zeros(6) # last values of qb1, qb2, qb3, qb4, qb5 and qb6
        for i in range(sz_spar):
            s_current = y_spar[i]
            if y_spar[i] < 0:
                s_0 = 0
                qb_current = def_integral(qb_5_vz, s_0, s_current, num_var_integrals=1)
            else:
                s_0 = 0
                qb_current = def_integral(qb_2_vz, s_0, s_current, num_var_integrals=1)

            qb_spar_val_list[i] = qb_current
        qb_lastval[1] = qb_spar_val_list[-1]
        qb_lastval[4] = qb_spar_val_list[0]

        """Outer segment shear flow distribution
        Order: 1,3,4,6
        """
        # Initialize simulation parameters
        x = np.arange(0, C+dx,dx)
        x = np.append(x, x[-2::-1])
        sz = np.size(x)
        y = np.zeros(sz)
        stiff_counter = 1
        qb_current = 0  # start at cut
        qb_stiff_current = 0
        qb_val_list = np.zeros(sz)

        perimeter = D.total_perimeter
        s_list = [m.pi/2, l_sk, l_sk, m.pi/2]  #list containing lengths of each skin 1,3,4,6
        i_list = []
        spacing = perimeter/n
        perimeter_addition = 0
        s_current = 0
        seg_i = 0

        for i in range(sz):  # go through the outer section

            # End
            if seg_i != 0 and x[i] == 0.:
                qb_lastval[5] = qb_current
                qb_val_list[i] = qb_current
                # print("Finished calculation at segment {}, i={} with x,y {},{} and qb of {}".format(seg_i,i, x[i],y[i],qb_current))
                # print("Max qb value of {} found at x,y = {},{}".format(np.amax(qb_val_list), x[np.where(qb_val_list == np.amax(qb_val_list))], y[np.where(qb_val_list == np.amax(qb_val_list))]))
                break

            # Switch segments
            if x[i] == h/2 or x[i] == C:
                s_current = 0
                if seg_i == 0:
                    qb_lastval[0] = qb_current
                elif seg_i == 1:
                    qb_lastval[2] = qb_current
                elif seg_i == 2:
                    qb_lastval[3] = qb_current

                seg_i += 1
                i_list.append(i)
                # print("Finished at segment {}, i={} with x,y {},{} and qb of {}".format(seg_i,i, x[i],y[i],qb_current))

            if seg_i == 0:
                y[i] = m.sqrt((h/2)**2-(x[i]-h/2)**2)
                split_rad = (h/2/dx)
                si = (m.pi/2)*i/split_rad
                qb_current = def_integral(qb_1_vz, 0, si, num_var_integrals=1)

            elif seg_i == 1:
                y[i] = -((h/2)/(C-h/2))*x[i] + (((h/2)/(C-h/2))*C)
                qb_current = def_integral(qb_3_vz, 0, s_current, num_var_integrals=1) + qb_lastval[0] + qb_lastval[1]

            elif seg_i == 2:
                y[i] = -(-((h/2)/(C-h/2))*x[i] + (((h/2)/(C-h/2))*C))
                qb_current = def_integral(qb_4_vz, 0, s_current, num_var_integrals=1) + qb_lastval[2]

            elif seg_i == 3:
                y[i] = -m.sqrt((h/2)**2-(x[i]-h/2)**2)
                s_0 = -m.pi/2
                qb_current = def_integral(qb_6_vz, s_0, (m.pi/2)*(i-i_list[-1])/split_rad+s_0, num_var_integrals=1) + qb_lastval[3] - qb_lastval[4]

            #Stiffener addition to qb
            perimeter_old = perimeter_addition
            perimeter_addition += m.hypot(dx, y[i]-y[i-1])
            if perimeter_addition > spacing > perimeter_old:
                if perimeter_addition - spacing <= perimeter_old - spacing:
                    qb_stiff_current += (A_stiff * (x[i] + y[i]) )/I_yy
                else:
                    qb_stiff_current += (A_stiff * (x[i-1] + y[i-1]) )/I_yy
                perimeter_addition -= spacing
                stiff_counter += 1

            qb_current += qb_stiff_current
            qb_val_list[i] = qb_current
            s_current += m.hypot(dx, y[i]-y[i-1])
        return qb_lastval, qb_val_list, qb_spar_val_list


    qb_lastval_vz, qbvz_outer_val_list, qbvz_spar_val_list = shear_flow_Vz()

    """Final shear distribution is a linear combination of Vy(x) and Vz(x)"""

    def final_shear_flow():
        """Add qso to qb of Vy(x) and add qb of Vz(x) and the q from torque"""
        #Shear due to Vy(x) given by qbvy_outer_val_list in the outer section and qbvy_spar_val_list in the spar
        #Shear due to Vz(x) given by qbvz_outer_val_list in the outer section and qbvz_spar_val_list in the spar

        # Outer section
        x = np.arange(0, C+dx,dx)
        x = np.append(x, x[-2::-1])
        sz = np.size(x)
        y = np.zeros(sz)
        q_val_outer_list = np.array(qbvy_outer_val_list) + np.array(qbvz_outer_val_list)
        shear_stress_outer = np.zeros(sz)
        y_spar = np.arange(-h/2, h/2, dx)
        sz_spar = np.size(y_spar)
        q_val_spar_list = np.array(qbvy_spar_val_list) + np.array(qbvz_spar_val_list)
        shear_stress_spar = np.zeros(sz_spar)
        seg_i = 0
        sz = np.size(D.x)
        for i in range(sz):  # go through the outer section

            # Switch segments
            if seg_i != 0 and x[i] == 0.:
                q_val_outer_list[i] += qs_0_1
                shear_stress_outer[i] = q_val_outer_list[i]/t_skin
                max_val = np.amax(q_val_outer_list)
                min_val = np.amin(q_val_outer_list)
                if abs(min_val) > max_val:
                    max_val = min_val
                max_val_i = np.where(q_val_outer_list == max_val)
                # print("Max q value of {} found at outer section x,y = {},{}".format(q_val_outer_list[max_val_i], x[max_val_i], y[max_val_i]))
                break

            if x[i] == h/2 or x[i] == C:
                seg_i += 1

            if seg_i == 0:
                q_val_outer_list[i] += qs_0_1 - q_t_1

            elif seg_i == 1:
                q_val_outer_list[i] += qs_0_2  - q_t_2


            elif seg_i == 2:
                q_val_outer_list[i] += qs_0_2  - q_t_2

            elif seg_i == 3:
                q_val_outer_list[i] += qs_0_1  - q_t_1

            shear_stress_outer[i] = q_val_outer_list[i]/t_skin

        # Spar section
        for i in range(sz_spar):
            q_val_spar_list[i] += qs_0_2 - qs_0_1 + q_t_3
            shear_stress_spar[i] = q_val_spar_list[i]/t_spar

        max_val = np.amax(q_val_spar_list)
        min_val = np.amin(q_val_spar_list)
        if abs(min_val) > max_val:
            max_val = min_val
        max_val_i = np.where(q_val_spar_list == max_val)
        # print("Max q value of {} found at spar x,y = {},{}".format(q_val_spar_list[max_val_i], h/2, y_spar[max_val_i]))

        final_q_values = np.array([q_val_outer_list,q_val_spar_list])
        shear_stress_values = np.array([shear_stress_outer,shear_stress_spar])
        return final_q_values, shear_stress_values

    final_q_values, shear_stress_values = final_shear_flow()

    def direct_stress():
        """Direct stress distribution results from a linear combination of My(x) and Mz(x)"""
        x = D.x
        y = D.y
        x_centroid = D.x_centroid
        sz = len(x)
        sigma_xx_z = np.zeros(sz)
        sigma_xx_y = np.zeros(sz)
        # outer section
        for i in range(len(x)):
            sigma_xx_z[i] = My(xi)*(-x[i]+x_centroid)/I_yy
            sigma_xx_y[i] = Mz(xi)*y[i]/I_zz

        sigma_xx_outer = sigma_xx_z + sigma_xx_y

        # Spar section
        h = D.h
        y_spar = np.arange(-h/2, h/2, dx)
        sz_spar = np.size(y_spar)
        x_spar = np.ones(sz_spar)*h/2
        sigma_xx_z = np.zeros(sz_spar)
        sigma_xx_y = np.zeros(sz_spar)
        # outer section
        for i in range(sz_spar):
            sigma_xx_z[i] = My(xi)*(-x_spar[i]+x_centroid)/I_yy
            sigma_xx_y[i] = Mz(xi)*y_spar[i]/I_zz

        sigma_xx_spar = sigma_xx_z + sigma_xx_y
        sigma_xx_values = np.array([sigma_xx_outer,sigma_xx_spar])
        return sigma_xx_outer, sigma_xx_spar, sigma_xx_values

    sigma_xx_outer, sigma_xx_spar, sigma_xx_values = direct_stress()

    def von_mises_stress():
        sigma_vm_outer = np.sqrt( 0.5*(sigma_xx_outer**2) + 3*(shear_stress_values[0]**2))
        sigma_vm_spar = np.sqrt( 0.5*(sigma_xx_spar**2) + 3*(shear_stress_values[1]**2))
        sigma_vm_values = np.array([sigma_vm_outer,sigma_vm_spar])
        return sigma_vm_outer, sigma_vm_spar, sigma_vm_values

    sigma_vm_outer, sigma_vm_spar, sigma_vm_values = von_mises_stress()

    sigma_vm_all[spar_segment] = sigma_vm_values
    spar_segment += 1

sigma_vm_max_all = []
for i in range(len(sigma_vm_all)):
    # Find maximum of this cross section
    sigma_vm_max_all.append(np.maximum(sigma_vm_all[i]))

max_vm_value = max(sigma_vm_max_all)
max_vm_index = [i for i, j in enumerate(sigma_vm_max_all) if j == max_vm_value]


if __name__ == '__main__':
    final_q_values, shear_stress_values = final_shear_flow()
    q_vy_values = np.array([qbvy_outer_val_list,qbvy_spar_val_list])
    q_vz_values = np.array([qbvz_outer_val_list,qbvz_spar_val_list])
    plot_shear_flow = False
    plot_shear_stress = False
    plot_direct_stress = False
    plot_vonmises = False
    if plot_shear_flow:
        y_spar = np.arange(-h/2, h/2, dx)
        fig, axs = plt.subplots(3)
        fig.suptitle('qb due to Vy(xi) Final shear flow values')
        axs[0].plot(D.x[0:len(q_vy_values[0])//2], q_vy_values[0][0:len(q_vy_values[0])//2])
        axs[0].set_title('Upper outer section')

        axs[1].plot(D.x[len(q_vy_values[0])//2:], q_vy_values[0][len(q_vy_values[0])//2:])
        axs[1].set_title('Lower outer section')

        axs[2].plot(y_spar, q_vy_values[1])
        axs[2].set_title('Spar section')

        plt.show()
        y_spar = np.arange(-h/2, h/2, dx)
        fig, axs = plt.subplots(3)
        fig.suptitle('qb due to Vz(x) shear flow values')
        axs[0].plot(D.x[0:len(q_vz_values[0])//2], q_vz_values[0][0:len(q_vz_values[0])//2])
        axs[0].set_title('Upper outer section')

        axs[1].plot(D.x[len(q_vz_values[0])//2:], q_vz_values[0][len(q_vz_values[0])//2:])
        axs[1].set_title('Lower outer section')

        axs[2].plot(y_spar, q_vz_values[1])
        axs[2].set_title('Spar section')

        plt.show()

        y_spar = np.arange(-h/2, h/2, dx)
        fig, axs = plt.subplots(3)
        fig.suptitle('Final shear flow values')
        axs[0].plot(D.x[0:len(final_q_values[0])//2], final_q_values[0][0:len(final_q_values[0])//2])
        axs[0].set_title('Upper outer section')

        axs[1].plot(D.x[len(final_q_values[0])//2:], final_q_values[0][len(final_q_values[0])//2:])
        axs[1].set_title('Lower outer section')

        axs[2].plot(y_spar, final_q_values[1])
        axs[2].set_title('Spar section')

        plt.show()

    if plot_shear_stress:
        y_spar = np.arange(-h/2, h/2, dx)
        fig, axs = plt.subplots(3)
        shear_stress_values = shear_stress_values/(10**6)
        fig.suptitle('Final shear stress values in MPA')
        axs[0].plot(D.x[0:len(shear_stress_values[0])//2], shear_stress_values[0][0:len(shear_stress_values[0])//2])
        axs[0].set_title('Upper outer section')

        axs[1].plot(D.x[len(shear_stress_values[0])//2:], shear_stress_values[0][len(shear_stress_values[0])//2:])
        axs[1].set_title('Lower outer section')

        axs[2].plot(y_spar, shear_stress_values[1])
        axs[2].set_title('Spar section')

        plt.show()
    if plot_direct_stress:
        y_spar = np.arange(-h/2, h/2, dx)
        fig, axs = plt.subplots(3)
        shear_stress_values = shear_stress_values/(10**6)
        fig.suptitle('Direct stress distribution')
        axs[0].plot(D.x[0:len(sigma_xx_values[0])//2], sigma_xx_values[0][0:len(sigma_xx_values[0])//2])
        axs[0].set_title('Upper outer section')

        axs[1].plot(D.x[len(sigma_xx_values[0])//2:], sigma_xx_values[0][len(sigma_xx_values[0])//2:])
        axs[1].set_title('Lower outer section')

        axs[2].plot(y_spar, sigma_xx_values[1])
        axs[2].set_title('Spar section')

    if plot_vonmises:
        y_spar = np.arange(-h/2, h/2, dx)
        fig, axs = plt.subplots(3)
        fig.suptitle('Von Mises stress distribution')
        axs[0].plot(D.x[0:len(sigma_vm_values[0])//2], sigma_vm_values[0][0:len(sigma_vm_values[0])//2])
        axs[0].set_title('Upper outer section')

        axs[1].plot(D.x[len(sigma_vm_values[0])//2:], sigma_vm_values[0][len(sigma_vm_values[0])//2:])
        axs[1].set_title('Lower outer section')

        axs[2].plot(y_spar, sigma_vm_values[1])
        axs[2].set_title('Spar section')
