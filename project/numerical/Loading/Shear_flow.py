import os
import sys
import math as m
import numpy as np
import matplotlib.pyplot as plt

# This must come before the next imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))
from project.numerical.distribution_equations import DistributionEquations
from project.numerical.Loading.integration import def_integral
from project.numerical.reaction_forces import AileronGeometry

G = AileronGeometry()
E = DistributionEquations()


# Real values
C = 0.515
h = 0.248
n = 11
w_st = 0.03
h_st = 0.015
t_st = 0.0012
t_skin = 0.0011
t_spar = 0.0022
A_stiff = w_st*t_st + h_st*t_st
l_sk = m.hypot(C-h/2,h/2)

I_zz = 1.4221372629975417e-5
I_yy = 5.643650631210155e-5
Vy = 1
Vz = 1
eta = G.z_tilde
# Test values

# h = 2
# t_sk = 1
# t_sp = 1
# Vy = 1
# Vz = 1
# I_zz = 1
# I_yy = 1
# eta = 1
# C_a = 3

# Segment 2, then 1,3,4,6,5
# define qb_segment#_1 as Vy and  qbsegment#_2 as Vz
# Shear flow 2

def qb_2_1(s):
    return -s*Vy*t_spar/I_zz
def qb_2_2(s):
    val = -eta*Vz*t_spar/I_yy
    return np.ones_like(s) * val

qb_2_val = def_integral(qb_2_1, 0, h/2, num_var_integrals=1) + def_integral(qb_2_2, 0, h/2, num_var_integrals=1)

dx = 0.0001
x = np.arange(0, C+dx,dx)
x = np.append(x, x[-2::-1])
sz = np.size(x)
y = np.zeros(sz)
stiff_loc = np.zeros((n,2))
stiff_counter = 1
qb_val_list = np.zeros(sz)

s_list = [m.pi/2, l_sk, l_sk, m.pi/2]  #list containing lengths of each skin 1,3,4,6
s_current = 0
seg_i = 0

"""
Go through outer section, using x, y find if the current location is equal to a stiffener position, if true add to the correct segment.
Change segments by keeping track of perimeter, if perimeter exceeds the length a segment, then switch to next segment.
"""


qb_lastval = [] # last values of qb1, qb3, qb4 and qb6

# Stiffener locations and general shape

perimeter = m.pi*(h/2) + 2*m.sqrt((h/2)**2+(C-h/2)**2)
spacing = perimeter/n
perimeter_addition = 0

qb_current = 0  # start at cut
qb_stiff_current = 0
for i in range(sz):  # go through the outer section

    # Switch segments
    if seg_i != 0 and x[i] == 0.:
        qb_lastval.append(qb_current)
        qb_val_list[i] = qb_current
        print("Finished at segment {} with x,y {},{} and qb of {}".format(seg_i,x[i],y[i],qb_current))
        print("Max qb value of {} found at x,y = {},{}".format(np.amax(qb_val_list), x[np.where(qb_val_list == np.amax(qb_val_list))], y[np.where(qb_val_list == np.amax(qb_val_list))]))
        break

    if x[i] == h/2 or x[i] == C:
        s_current = 0
        # print("qb_current: ", qb_current)
        qb_lastval.append(qb_current)
        seg_i += 1
        if seg_i == 1:
            qb_current += qb_2_val
        else:
            qb_current += qb_lastval[-1]
        # print("segment ", seg_i, "with x,y", x[i], y[i-1])


    if x[i] <= h/2 and seg_i == 0:
        y[i] = m.sqrt((h/2)**2-(x[i]-h/2)**2)

    elif seg_i == 1:
        y[i] = -((h/2)/(C-h/2))*x[i] + (((h/2)/(C-h/2))*C)

    elif seg_i == 2:
        y[i] = -(-((h/2)/(C-h/2))*x[i] + (((h/2)/(C-h/2))*C))

    elif seg_i == 3:
        y[i] = -m.sqrt((h/2)**2-(x[i]-h/2)**2)

    qb_i_1 = -Vy*t_skin*y[i]/I_zz
    qb_i_2 = -Vz*t_skin*(eta - x[i])/I_yy
    qb_i = [qb_i_1,qb_i_2]

    def integrate_func1(s):
        return np.ones_like(s) * qb_i[0]

    def integrate_func2(s):
        return np.ones_like(s) * qb_i[1]

    # Numerically integrate qb for each point

    qb_current = def_integral(integrate_func1, 0, s_current, num_var_integrals=1) + def_integral(integrate_func2, 0, s_current, num_var_integrals=1)

    #Stiffener addition to qb
    perimeter_old = perimeter_addition
    perimeter_addition += m.hypot(dx, y[i]-y[i-1])
    if perimeter_addition > spacing and perimeter_old < spacing and stiff_counter <= 10:
        if perimeter_addition - spacing <= perimeter_old - spacing:
            stiff_loc[stiff_counter,0] = x[i]
            stiff_loc[stiff_counter,1] = y[i]
            qb_stiff_current += A_stiff * (x[i] + y[i])
        else:
            stiff_loc[stiff_counter,0] = x[i-1]
            stiff_loc[stiff_counter,1] = y[i-1]
            qb_stiff_current += A_stiff * (x[i - 1] + y[i - 1])
        perimeter_addition -= spacing
        stiff_counter += 1

    qb_current = qb_current + qb_stiff_current
    qb_val_list[i] = qb_current
    s_current += m.hypot(dx, y[i]-y[i-1])

for i in range (int((n-1)/2)):
    stiff_loc[n-(i+1),0] = stiff_loc[i+1,0]
    stiff_loc[n-(i+1),1] = -stiff_loc[i+1,1]


# # shear flow at 5
def qb_5_1(s):
    return -Vy*t_spar*s/I_zz
def qb_5_2(s):
    val = -eta*Vz*t_spar/I_yy
    return np.ones_like(s) * val

qb_2_val = def_integral(qb_2_1, 0, h/2, num_var_integrals=1) + def_integral(qb_2_2, 0, h/2, num_var_integrals=1)
qb_5_val = def_integral(qb_5_1, -h/2, 0, num_var_integrals=1) + def_integral(qb_5_2, -h/2, 0, num_var_integrals=1)

"""Spar shear flow distribution"""
y_spar = np.arange(-h/2, h/2, dx)
sz_spar = np.size(y_spar)
x_spar = np.zeros(sz_spar)
qb_spar_val_list = np.zeros(sz_spar)
qb_last = qb_lastval[2] # qb value before entering 5
for i in range(sz_spar):
    if dx/10 >= y_spar[i] >= -dx/10:
        qb_last = 0  # cut at mid of spar
    if y_spar[i] <= 0:
        s_0 = -h/2
    else:
        s_0 = 0

    s_current = y_spar[i]

    qb_i_1 = -Vy*t_spar*y_spar[i]/I_zz
    qb_i_2 = -Vz*t_spar*(eta - x_spar[i])/I_yy
    qb_i = [qb_i_1,qb_i_2]

    def integrate_func1(s):
        return np.ones_like(s) * qb_i[0]

    def integrate_func2(s):
        return np.ones_like(s) * qb_i[1]

    qb_current = def_integral(integrate_func1, s_0, s_current, num_var_integrals=2) + def_integral(integrate_func2, s_0, s_current, num_var_integrals=1)  #integrating qb_2_1 twice gives correct value...
    qb_current += qb_last
    qb_spar_val_list[i] = qb_current


print("Value of qb2 at spar",qb_2_val)

qb_1 = qb_lastval[0]
qb_2 = qb_2_val
qb_3 = qb_lastval[1]
qb_4 = qb_lastval[2]
qb_5 = qb_5_val
qb_6 = qb_lastval[3]


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

#moment arm from traling edge [vertical distance, horizontal distance] 1,3,4,6
#for 2 and 5 = [0,eta],[0,eta]
moment_arm_list = [[0,C+h/2],[h/2,C-h/2],[0,0],[h/2,C-h/2]]

M_b_1 = (C) * qb_1_v
M_b_2 = (C - h/2) * qb_2_v
M_b_3 = h/2 * qb_3_h + (C-h/2) * -qb_3_v
M_b_4 = 0
M_b_5 = (C - h/2) * qb_5_v
M_b_6 = h/2 * qb_6_h + (C-h/2) * qb_6_v

M_vy = (C-eta) * Vy
M_vz = 0

#semicircle
p1 = m.pi*h/2 + h
qs_0_1 = -(qb_1*m.pi*h/4 + qb_2*h/2 + qb_5*h/2 + qb_6*m.pi*h/4)/p1

#triangle
p2 = h + 2*l_sk
qs_0_2 = -(qb_3*l_sk + qb_4*l_sk + qb_2*h/2 + qb_5*h/2)/p2


"""Add qso to qb for final shear flow distribution"""

# Outer section
q_val_list = qb_val_list.copy()
seg_i = 0
for i in range(sz):  # go through the outer section

    # Switch segments
    if seg_i != 0 and x[i] == 0.:
        q_val_list[i] += qs_0_1
        print("Finished at segment {} with x,y {},{} and final q of {}".format(seg_i,x[i],y[i],qb_current))
        print("Max q value of {} found at x,y = {},{}".format(np.amax(q_val_list), x[np.where(q_val_list == np.amax(q_val_list))], y[np.where(q_val_list == np.amax(q_val_list))]))
        break

    if x[i] == h/2 or x[i] == C:
        seg_i += 1

    if x[i] <= h/2 and seg_i == 0:
        q_val_list[i] += qs_0_1

    elif seg_i == 1:
        q_val_list[i] += qs_0_2

    elif seg_i == 2:
        q_val_list[i] += qs_0_2

    elif seg_i == 3:
        q_val_list[i] += qs_0_1

# Spar section
q_spar_val_list = qb_spar_val_list.copy()
for i in range(sz_spar):
    if y_spar[i] == 0.:
        qb_current = 0  # cut at mid of spar
    q_spar_val_list[i] += qs_0_2 - qs_0_1
    s_current = y_spar[i]

q_1 = qb_1 + qs_0_1
q_2 = qb_2 + qs_0_2 - qs_0_1
q_3 = qb_3 + qs_0_2
q_4 = qb_4 + qs_0_2
q_5 = qb_5 + qs_0_2 - qs_0_1
q_6 = qb_6 + qs_0_1

# Total torque
Am1 = 0.5* m.pi*(h/2)**2
Am2 = 0.5 * (C-h/2)*h
T = 2*Am1*qs_0_1 + 2*Am2*qs_0_2

comp_eq_1 = 1/(2*Am1)*((q_1*m.pi*h/4)/t_skin + (qb_2*h/2)/t_spar + (qb_5*h/2)/t_spar + (qb_6*m.pi*h/4)/t_skin)
comp_eq_2 = 1/(2*Am2)*((qb_3*l_sk)/t_skin + (qb_4*l_sk)/t_skin + (qb_2*h/2)/t_spar + (qb_5*h/2)/t_spar)
T = 1
J1 = T/(comp_eq_1)
J2 = T/(comp_eq_2)
print(J1, J2)
# print(qb_val_list)
# print(q_val_list)

# print(stiff_loc)
# print(qb_lastval)

# plt.plot(stiff_loc[:,0],stiff_loc[:,1], marker="o")
# plt.plot(x,y)
plt.plot(x, q_val_list)
# plt.plot(x,-y)
# plt.axis('equal')
plt.show()
