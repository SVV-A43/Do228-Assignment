import math as m
import numpy as np
import matplotlib.pyplot as plt
from project.numerical.Loading.integration import def_integral
from project.numerical.reaction_forces import AileronGeometry
from project.numerical.Section_Properties.cross_section_properties import Cross_section_properties


# Initialize parameters

D = Cross_section_properties()
G = AileronGeometry()

Vy = -36718.115  # verification model
Vz = 61222.16  # verification model
My = 19633.84 # verification model
Mz = -12318.44 # verification model
eta = G.z_tilde

C = D.C
h = D.h
n = D.n
w_st = D.w_st
h_st = D.h_st
t_st = D.t_st
t_skin = D.t_skin
t_spar = D.t_spar
A_stiff = D.A_stiff
l_sk = D.l_sk_triangle
I_zz = D.I_zz
I_yy = D.I_yy
x_centroid = D.x_centroid
dx = D.dx


# Test values

# h = 1
# t_skin = 1
# t_spar = 1
# Vy = 1
# Vz = 1
# I_zz = 1
# I_yy = 1
# eta = 1
# C = 3


def shear_flow_Vy():

    h = D.h/2  # Define h to be h/2
    def qb_1_vy(s):
        return -Vy*t_skin*h*h*np.sin(s)/I_zz
    def qb_6_vy(s):
        return -Vy*t_skin*h*h*np.sin(s)/I_zz
    def qb_2_vy(s):
        return -Vy*t_spar*s/I_zz
    def qb_5_vy(s):
        return -Vy*t_spar*s/I_zz
    def qb_3_vy(s):
        return -Vy*t_skin*(h-h/l_sk)*s/I_zz
    def qb_4_vy(s):
        return -Vy*t_skin*(-h/l_sk)*s/I_zz

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
            print("Finished at segment {}, i={} with x,y {},{} and qb of {}".format(seg_i,i, x[i],y[i],qb_current))
            print("Max qb value of {} found at x,y = {},{}".format(np.amax(qb_val_list), x[np.where(qb_val_list == np.amax(qb_val_list))], y[np.where(qb_val_list == np.amax(qb_val_list))]))
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
            print("Finished at segment {}, i={} with x,y {},{} and qb of {}".format(seg_i,i, x[i],y[i],qb_current))

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

        qb_current += qb_stiff_current
        qb_val_list[i] = qb_current
        s_current += m.hypot(dx, y[i]-y[i-1])
    return qb_lastval, qb_val_list, qb_spar_val_list

qb_lastval_vy, qbvy_outer_val_list, qbvy_spar_val_list = shear_flow_Vy()

"""Calculate shear center and redundant shear flows"""

def shear_center():
    qb_1 = qb_lastval_vy[0]
    qb_2 = qb_lastval_vy[1]
    qb_3 = qb_lastval_vy[2]
    qb_4 = qb_lastval_vy[3]
    qb_5 = qb_lastval_vy[4]
    qb_6 = qb_lastval_vy[5]
    # semicircle
    h = D.h
    p1 = m.pi*h/2 + h
    qs_0_1 = -(qb_1*m.pi*h/4 - qb_2*h/2 - qb_5*h/2 + qb_6*m.pi*h/4)/p1

    #triangle
    p2 = h + 2*l_sk
    qs_0_2 = -(qb_3*l_sk + qb_4*l_sk + qb_2*h/2 + qb_5*h/2)/p2

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
    print("Shear center location from leading edge",eta)
    return eta, qs_0_1, qs_0_2

eta, qs_0_1, qs_0_2 = shear_center()

def shear_flow_Vz():

    h = D.h/2  # Define h to be h/2
    def qb_1_vz(s):
        return -Vz * t_skin * h * (-1 * (1 - np.cos(s) * h - x_centroid)) / I_yy
    def qb_6_vz(s):
        return -Vz * t_skin * h * (-1 * (1 - np.cos(s) * h - x_centroid)) / I_yy
    def qb_2_vz(s):
        val = -t_spar * (-h - x_centroid) / I_yy
        return np.ones_like(s) * val
    def qb_5_vz(s):
        val = -t_spar * (-h - x_centroid) / I_yy
        return np.ones_like(s) * val
    def qb_3_vz(s):
        return -Vz * t_skin * ((-h - x_centroid) - ((C - h) * s / l_sk)) / I_yy
    def qb_4_vz(s):
        return -Vz * t_skin * ((-C - x_centroid) - ((C - h) * s)) / I_yy

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
            print("Finished at segment {}, i={} with x,y {},{} and qb of {}".format(seg_i,i, x[i],y[i],qb_current))
            print("Max qb value of {} found at x,y = {},{}".format(np.amax(qb_val_list), x[np.where(qb_val_list == np.amax(qb_val_list))], y[np.where(qb_val_list == np.amax(qb_val_list))]))
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
            print("Finished at segment {}, i={} with x,y {},{} and qb of {}".format(seg_i,i, x[i],y[i],qb_current))

        if seg_i == 0:
            y[i] = m.sqrt((h/2)**2-(x[i]-h/2)**2)
            split_rad = (h/2/dx)
            qb_current = def_integral(qb_1_vz, 0, (m.pi/2)*i/split_rad, num_var_integrals=1)

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
                qb_stiff_current += (A_stiff * (x[i] + y[i]) )/I_zz
            else:
                qb_stiff_current += (A_stiff * (x[i-1] + y[i-1]) )/I_zz
            perimeter_addition -= spacing
            stiff_counter += 1

        qb_current += qb_stiff_current
        qb_val_list[i] = qb_current
        s_current += m.hypot(dx, y[i]-y[i-1])
    return qb_lastval, qb_val_list, qb_spar_val_list


qb_lastval_vz, qbvz_outer_val_list, qbvz_spar_val_list = shear_flow_Vz()

"""Final shear distribution is a linear combination of Vy and Vz"""

def final_shear_flow():
    """Add qso to qb of vy and add qb of vz"""
    #Shear due to vy given by qbvy_outer_val_list in the outer section and qbvy_spar_val_list in the spar
    #Shear due to vz given by qbvz_outer_val_list in the outer section and qbvz_spar_val_list in the spar

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
            print("Max q value of {} found at outer section x,y = {},{}".format(q_val_outer_list[max_val_i], x[max_val_i], y[max_val_i]))
            break

        if x[i] == h/2 or x[i] == C:
            seg_i += 1

        if seg_i == 0:
            q_val_outer_list[i] += qs_0_1

        elif seg_i == 1:
            q_val_outer_list[i] += qs_0_2

        elif seg_i == 2:
            q_val_outer_list[i] += qs_0_2

        elif seg_i == 3:
            q_val_outer_list[i] += qs_0_1

        shear_stress_outer[i] = q_val_outer_list[i]/t_skin

    # Spar section
    for i in range(sz_spar):
        q_val_spar_list[i] += qs_0_2 - qs_0_1
        shear_stress_spar[i] = q_val_spar_list[i]/t_spar

    max_val = np.amax(q_val_spar_list)
    min_val = np.amin(q_val_spar_list)
    if abs(min_val) > max_val:
        max_val = min_val
    max_val_i = np.where(q_val_spar_list == max_val)
    print("Max q value of {} found at spar x,y = {},{}".format(q_val_spar_list[max_val_i], h/2, y_spar[max_val_i]))

    final_q_values = np.array([q_val_outer_list,q_val_spar_list])
    shear_stress_values = np.array([shear_stress_outer,shear_stress_spar])
    return final_q_values, shear_stress_values

final_q_values, shear_stress_values = final_shear_flow()

def direct_stress():
    """Direct stress distribution results from a linear combination of My and Mz"""
    x = D.x
    y = D.y
    x_centroid = D.x_centroid
    sz = len(x)
    sigma_xx_z = np.zeros(sz)
    sigma_xx_y = np.zeros(sz)
    # outer section
    for i in range(len(x)):
        sigma_xx_z[i] = My*(x[i]-x_centroid)/I_yy
        sigma_xx_y[i] = Mz*y[i]/I_zz

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
        sigma_xx_z[i] = My*(x_spar[i]-x_centroid)/I_yy
        sigma_xx_y[i] = Mz*y_spar[i]/I_zz

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

if __name__ == '__main__':
    final_q_values, shear_stress_values = final_shear_flow()
    plot_shear_flow = True
    plot_shear_stress = False
    plot_direct_stress = True
    plot_vonmises = True
    if plot_shear_flow:
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




# def shear_flow():
#     """
#     Calculate the shear flow distribution along the cross-section
#     and return an array with the distribution, with the rows being the segment: 1-2-3-4-5-6
#     """
#
#     # Initialize simulation parameters
#     x = np.arange(0, C+dx,dx)
#     x = np.append(x, x[-2::-1])
#     sz = np.size(x)
#     y = np.zeros(sz)
#     stiff_loc = np.zeros((n,2))
#     stiff_counter = 1
#     qb_current = 0  # start at cut
#     qb_stiff_current = 0
#     qb_prev = 0
#     qb_val_list = np.zeros(sz)
#     qb_lastval = [] # last values of qb1, qb3, qb4 and qb6
#
#     perimeter = D.total_perimeter
#     s_list = [m.pi/2, l_sk, l_sk, m.pi/2]  #list containing lengths of each skin 1,3,4,6
#     spacing = perimeter/n
#     perimeter_addition = 0
#     s_current = 0
#     seg_i = 0
#
#     # Split base shear flow calculation in segments
#     # First segment 2 is calculated
#     # Then the outer section segments 1,3,4,6 are calculated
#     # Finally segment 5 is calculated
#
#     def qb_2_1(s):
#         return -s*Vy*t_spar/I_zz
#     def qb_2_2(s):
#         val = -eta*Vz*t_spar/I_yy
#         return np.ones_like(s) * val
#
#     qb_2_val = def_integral(qb_2_1, 0, h/2, num_var_integrals=2) + def_integral(qb_2_2, 0, h/2, num_var_integrals=1)
#
#     for i in range(sz):  # go through the outer section
#
#         # Switch segments
#         if seg_i != 0 and x[i] == 0.:
#             qb_lastval.append(qb_current)
#             qb_val_list[i] = qb_current
#             print("Finished at segment {}, i={} with x,y {},{} and qb of {}".format(seg_i,i, x[i],y[i],qb_current))
#             print("Max qb value of {} found at x,y = {},{}".format(np.amax(qb_val_list), x[np.where(qb_val_list == np.amax(qb_val_list))], y[np.where(qb_val_list == np.amax(qb_val_list))]))
#             break
#
#         if x[i] == h/2 or x[i] == C:
#             s_current = 0
#             # print("qb_current: ", qb_current)
#             qb_lastval.append(qb_current)
#             seg_i += 1
#             if seg_i == 1:
#                 qb_current += qb_2_val
#                 qb_prev += qb_current - qb_stiff_current # to not add stiffeners twice.
#             else:
#                 qb_prev = qb_current
#             print("Finished at segment {}, i={} with x,y {},{} and qb of {}".format(seg_i,i, x[i],y[i],qb_current))
#
#         if x[i] <= h/2 and seg_i == 0:
#             y[i] = m.sqrt((h/2)**2-(x[i]-h/2)**2)
#             n_int1 = 2
#             n_int2 = 1
#
#         elif seg_i == 1:
#             y[i] = -((h/2)/(C-h/2))*x[i] + (((h/2)/(C-h/2))*C)
#             n_int1 = 2
#             n_int2 = 2
#
#         elif seg_i == 2:
#             y[i] = -(-((h/2)/(C-h/2))*x[i] + (((h/2)/(C-h/2))*C))
#             n_int1 = 2
#             n_int2 = 2
#
#         elif seg_i == 3:
#             y[i] = -m.sqrt((h/2)**2-(x[i]-h/2)**2)
#             n_int1 = 1
#             n_int2 = 1
#
#         qb_i_1 = -Vy*t_skin*y[i]/I_zz
#         qb_i_2 = -Vz*t_skin*(eta - x[i])/I_yy
#         qb_i = [qb_i_1,qb_i_2]
#
#         def integrate_func1(s):
#             return np.ones_like(s) * qb_i[0]
#
#         def integrate_func2(s):
#             return np.ones_like(s) * qb_i[1]
#
#         # Numerically integrate qb for each point
#         qb_current = def_integral(integrate_func1, 0, s_current, num_var_integrals=n_int1) + def_integral(integrate_func2, 0, s_current, num_var_integrals=n_int2)
#
#         #Stiffener addition to qb
#         perimeter_old = perimeter_addition
#         perimeter_addition += m.hypot(dx, y[i]-y[i-1])
#         if perimeter_addition > spacing and perimeter_old < spacing:
#             if perimeter_addition - spacing <= perimeter_old - spacing:
#                 qb_stiff_current += A_stiff * (x[i] + y[i])
#             else:
#                 qb_stiff_current += A_stiff * (x[i-1] + y[i-1])
#             print("added stiff{} at i = {}".format(stiff_counter,i))
#             perimeter_addition -= spacing
#             stiff_counter += 1
#
#         qb_current += qb_prev + qb_stiff_current
#         qb_val_list[i] = qb_current
#         s_current += m.hypot(dx, y[i]-y[i-1])
#
#
#     # base shear flow at 5
#     def qb_5_1(s):
#         return -Vy*t_spar*s/I_zz
#     def qb_5_2(s):
#         val = -eta*Vz*t_spar/I_yy
#         return np.ones_like(s) * val
#
#     qb_5_val = def_integral(qb_5_1, -h/2, 0, num_var_integrals=1) + def_integral(qb_5_2, -h/2, 0, num_var_integrals=1) + qb_lastval[2]
#
#     """Spar shear flow distribution"""
#     y_spar = np.arange(-h/2, h/2, dx)
#     sz_spar = np.size(y_spar)
#     x_spar = np.zeros(sz_spar)
#     qb_spar_val_list = np.zeros(sz_spar)
#     qb_last = qb_lastval[2] # qb value before entering 5
#     for i in range(sz_spar):
#         if dx/10 >= y_spar[i] >= -dx/10:
#             qb_last = 0  # cut at mid of spar
#         if y_spar[i] <= 0:
#             s_0 = -h/2
#             n_int1 = 2
#             n_int2 = 1
#         else:
#             n_int1 = 2
#             n_int2 = 1
#             s_0 = 0
#
#         s_current = y_spar[i]
#
#         qb_i_1 = -Vy*t_spar*y_spar[i]/I_zz
#         qb_i_2 = -Vz*t_spar*(eta - x_spar[i])/I_yy
#         qb_i = [qb_i_1,qb_i_2]
#
#         def integrate_func1(s):
#             return np.ones_like(s) * qb_i[0]
#
#         def integrate_func2(s):
#             return np.ones_like(s) * qb_i[1]
#
#         qb_current = def_integral(integrate_func1, s_0, s_current, num_var_integrals=n_int1) + def_integral(integrate_func2, s_0, s_current, num_var_integrals=n_int2)  #integrating qb_2_1 twice gives correct value...
#         qb_current += qb_last
#         qb_spar_val_list[i] = qb_current
#
#     # plt.subplot(2, 1, 1)
#     # plt.plot(D.x[0:len(qb_val_list)//2], qb_val_list[0:len(qb_val_list)//2])
#     # plt.xlabel('Upper outer section')
#     # plt.title('Shear flow values for outer section')
#     # plt.subplot(2, 1, 2)
#     # plt.plot(D.x[len(qb_val_list)//2:], qb_val_list[len(qb_val_list)//2:])
#     # plt.xlabel('Lower outer section')
#     # plt.ylabel('qb_values')
#
#     # y_spar = np.arange(-h/2, h/2, dx)
#     # fig, axs = plt.subplots(3)
#     # fig.suptitle('Base shear flow values for outer section')
#     # axs[0].plot(D.x[0:len(qb_val_list)//2], qb_val_list[0:len(qb_val_list)//2])
#     # axs[0].set_title('Upper outer section')
#     #
#     # axs[1].plot(D.x[len(qb_val_list)//2:], qb_val_list[len(qb_val_list)//2:])
#     # axs[1].set_title('Lower outer section')
#     #
#     # axs[2].plot(y_spar, qb_spar_val_list)
#     # axs[2].set_title('Spar section')
#
#
#     plt.show()
#
#     qb_1 = qb_lastval[0]
#     qb_2 = qb_2_val
#     qb_3 = qb_lastval[1]
#     qb_4 = qb_lastval[2]
#     qb_5 = qb_5_val
#     qb_6 = qb_lastval[3]
#
#     theta = m.asin((h/2)/l_sk)
#
#     qb_1_h = qb_1 / 2
#     qb_1_v = qb_1 / 2
#
#     qb_2_h = 0
#     qb_2_v = qb_2
#
#     qb_3_h = qb_3 * m.cos(theta)
#     qb_3_v = qb_3 * m.sin(theta)
#
#     qb_4_h = qb_4 * m.cos(theta)
#     qb_4_v = qb_4 * m.sin(theta)
#
#     qb_5_h = 0
#     qb_5_v = qb_5
#
#     qb_6_h = qb_6 / 2
#     qb_6_v = qb_6 / 2
#
#     #moment arm from trailing edge [vertical distance, horizontal distance] 1,3,4,6
#     #for 2 and 5 = [0,eta],[0,eta]
#     moment_arm_list = [[0,C+h/2],[h/2,C-h/2],[0,0],[h/2,C-h/2]]
#
#     m_b_1 = (C) * qb_1_v
#     m_b_2 = (C - h/2) * qb_2_v
#     m_b_3 = h/2 * qb_3_h + (C-h/2) * -qb_3_v
#     m_b_4 = 0
#     m_b_5 = (C - h/2) * qb_5_v
#     m_b_6 = h/2 * qb_6_h + (C-h/2) * qb_6_v
#
#     m_vy = (C-eta) * Vy
#     m_vz = 0
#
#     #semicircle
#     p1 = m.pi*h/2 + h
#     qs_0_1 = -(qb_1*m.pi*h/4 + qb_2*h/2 + qb_5*h/2 + qb_6*m.pi*h/4)/p1
#
#     #triangle
#     p2 = h + 2*l_sk
#     qs_0_2 = -(qb_3*l_sk + qb_4*l_sk + qb_2*h/2 + qb_5*h/2)/p2
#     print(qs_0_1,qs_0_2)
#
#     """Add qso to qb for final shear flow distribution"""
#
#     # Outer section
#     q_val_list = qb_val_list.copy()
#     seg_i = 0
#     for i in range(sz):  # go through the outer section
#
#         # Switch segments
#         if seg_i != 0 and x[i] == 0.:
#             q_val_list[i] += qs_0_1
#             print("Finished at segment {} with x,y {},{} and final q of {}".format(seg_i,x[i],y[i],qb_current))
#             max_val = np.amax(q_val_list)
#             min_val = np.amin(q_val_list)
#             if abs(min_val) > max_val:
#                 max_val = min_val
#             max_val_i = np.where(q_val_list == max_val)
#             print("Max q value of {} found at outer section x,y = {},{}".format(q_val_list[max_val_i], x[max_val_i], y[max_val_i]))
#             break
#
#         if x[i] == h/2 or x[i] == C:
#             seg_i += 1
#
#         if x[i] <= h/2 and seg_i == 0:
#             q_val_list[i] += qs_0_1
#
#         elif seg_i == 1:
#             q_val_list[i] += qs_0_2
#
#         elif seg_i == 2:
#             q_val_list[i] += qs_0_2
#
#         elif seg_i == 3:
#             q_val_list[i] += qs_0_1
#
#     # Spar section
#     q_spar_val_list = qb_spar_val_list.copy()
#     for i in range(sz_spar):
#         q_spar_val_list[i] += qs_0_2 - qs_0_1
#
#     max_val = np.amax(q_spar_val_list)
#     min_val = np.amin(q_spar_val_list)
#     if abs(min_val) > max_val:
#         max_val = min_val
#     max_val_i = np.where(q_spar_val_list == max_val)
#     print("Max q value of {} found at spar x,y = {},{}".format(q_spar_val_list[max_val_i], h/2, y_spar[max_val_i]))
#
#     # q_1 = qb_1 + qs_0_1
#     # q_2 = qb_2 + qs_0_2 - qs_0_1
#     # q_3 = qb_3 + qs_0_2
#     # q_4 = qb_4 + qs_0_2
#     # q_5 = qb_5 + qs_0_2 - qs_0_1
#     # q_6 = qb_6 + qs_0_1
#
#     final_q_values = np.array([q_val_list,q_spar_val_list])
#     return final_q_values
#