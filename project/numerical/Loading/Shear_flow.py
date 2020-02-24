import math as m
import numpy as np
import matplotlib.pyplot as plt
from project.numerical.Loading.integration import def_integral, variable_integral


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

I_zz = 1.4221372629975417e-5
I_yy = 5.643650631210155e-5

# Test values

h = 2
t_skin = 1
Vy = 1
Vz = 1
I_zz = 1
I_yy = 1
eta = 1
C = 3
l_sk = m.hypot(C-h/2,h/2)

# Segment 2, then 1,3,4,6,5
# define qb_segment#_1 as Vy and  qbsegment#_2 as Vz
# Shear flow 2


def qb_2_1(s):
    return -s*Vy*t_skin/I_zz
def qb_2_2(s):
    val = -eta*Vz*t_skin/I_yy
    return np.ones_like(s) * val


qb_2_val = def_integral(qb_2_1, 0, 1, num_bins=100) + def_integral(qb_2_2, 0, 1, num_bins=100)



# Shear flow 1,3,4,6
def qb_1_1(s):
    return -Vy*t_skin*(h/2)*np.sin(s)*(h/2)/I_zz
def qb_1_2(s):
    return -Vz*t_skin*(eta+(h/2)*(1-np.cos(s)))*(h/2)/I_yy

def qb_3_1(s):
    return -Vy*t_skin*((h/2)-h*s/2*l_sk)/I_zz
def qb_3_2(s):
    return -Vz*t_skin*(eta-s)/I_yy

def qb_4_1(s):
    return -Vy*t_skin*-(h*s/2*l_sk)/I_zz
def qb_4_2(s):
    return -Vz*t_skin*(eta-s)/I_yy

def qb_6_1(s):
    return -Vy*t_skin*(h/2)*-m.sin(s)*(h/2)/I_zz
def qb_6_2(s):
    return -Vz*t_skin*(eta+(h/2)*(1-m.cos(s)))*(h/2)/I_yy

dx = 0.0001
x = np.arange(0, C+dx,dx)
x = np.append(x, x[-2::-1])
sz = np.size(x)
y = np.zeros(sz)
stiff_loc = np.array([[ 0.        ,  0.        ],
                      [ 0.0455    ,  0.09598828],
                      [ 0.1479    ,  0.11642046],
                      [ 0.2527    ,  0.08318465],
                      [ 0.3576    ,  0.04991714],
                      [ 0.4624    ,  0.01668133],
                      [ 0.4624    , -0.01668133],
                      [ 0.3576    , -0.04991714],
                      [ 0.2527    , -0.08318465],
                      [ 0.1479    , -0.11642046],
                      [ 0.0455    , -0.09598828]])

s_list = [m.pi/2, l_sk, l_sk, m.pi/2]  #list containing lengths of each skin 1,3,4,6
s_current = 0
seg_i = 0  #0: 1, 1: 3, 2: 4, 3: 6
s_last = -10**4

"""
Go through outer section, using x, y find if the current location is equal to a stiffener position, if true add to the correct segment.
Change segments by keeping track of perimeter, if perimeter exceeds the length a segment, then switch to next segment.
"""


func_list = [[qb_1_1,qb_1_2],[], [], []]
qb_val_list = [qb_2_val]
#moment arm from traling edge [vertical distance, horizontal distance] 1,3,4,6
#for 2 and 5 = [0,eta],[0,eta]
moment_arm_list = [[0,C+h/2],[h/2,C-h/2],[0,0],[h/2,C-h/2]]

for i in range(sz):  # go through the outer section

    # Switch segments
    if s_current > s_list[seg_i]:
        seg_i += 1
        s_current = 0

    if seg_i == 0:
        y[i] = m.sqrt((h/2)**2-(x[i]-h/2)**2)
    elif seg_i == 1:
        y[i] = -((h/2)/(C-h/2))*x[i] + (((h/2)/(C-h/2))*C)
    elif seg_i == 2:
        y[i] = -(-((h/2)/(C-h/2))*x[i] + (((h/2)/(C-h/2))*C))
    elif seg_i == 3:
        y[i] = -m.sqrt((h/2)**2-(x[i]-h/2)**2)

    # Numerically integrate qb
    qb_current = def_integral(func_list[seg_i][0], 0, s_list[seg_i], num_bins=100) + def_integral(func_list[seg_i][1], 0, s_list[seg_i], num_bins=100)

    current_loc = np.array(x[i], y[i])
    if current_loc in stiff_loc:
        qb_current += A_stiff*y[i] + A_stiff * x[i]
        print("stiffener added at", current_loc)


    qb_val_list.append(qb_current)
    s_current += m.hypot(dx, y[i]-y[i-1])

# # shear flow at 5
def qb_5_1(s):
    return -Vy*t_skin*(s-h/2)/I_zz
def qb_5_2(s):
    val = -eta*Vz*t_skin/I_yy
    return np.ones_like(s) * val
