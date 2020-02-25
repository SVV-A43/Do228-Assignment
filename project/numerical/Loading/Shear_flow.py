import math as m
import numpy as np
import matplotlib.pyplot as plt
from project.numerical.Loading.integration import def_integral


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

# Test values

# h = 2
# t_skin = 1
# t_spar = 1
Vy = 1
Vz = 1
# I_zz = 1
# I_yy = 1
eta = 0.2
# C = 3

# Segment 2, then 1,3,4,6,5
# define qb_segment#_1 as Vy and  qbsegment#_2 as Vz
# Shear flow 2

def qb_2_1(s):
    return -s*Vy*t_spar/I_zz
def qb_2_2(s):
    val = -eta*Vz*t_spar/I_yy
    return np.ones_like(s) * val


qb_2_val = def_integral(qb_2_1, 0, 1, num_var_integrals=1) + def_integral(qb_2_2, 0, 1, num_var_integrals=1)

print(qb_2_val)

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
#moment arm from traling edge [vertical distance, horizontal distance] 1,3,4,6
#for 2 and 5 = [0,eta],[0,eta]

# Stiffener locations and general shape

perimeter = m.pi*(h/2) + 2*m.sqrt((h/2)**2+(C-h/2)**2)
spacing = perimeter/n
perimeter_addition = 0

qb_current = 0  # start at cut

for i in range(sz):  # go through the outer section

    # Switch segments
    if seg_i != 0 and x[i] == 0.:
        qb_lastval.append(qb_current)
        print("Finished at segment {} with x,y {},{} and qb of {}".format(seg_i,x[i],y[i],qb_current))
        break

    if x[i] == h/2 or x[i] == C:
        s_current = 0
        print("qb_current: ", qb_current)
        qb_lastval.append(qb_current)
        seg_i += 1
        if seg_i == 1:
            qb_current += qb_2_val
        else:
            qb_current += qb_lastval[-1]
        print("segment ", seg_i, "with x,y", x[i], y[i-1])

    # Shear flow 1,3,4,6
    # def qb_3_1(s):
    #     return -Vy*t_skin*((h/2)-h*s/2*l_sk)/I_zz
    # def qb_3_2(s):
    #     return -Vz*t_skin*(eta-s)/I_yy
    #
    # def qb_4_1(s):
    #     return -Vy*t_skin*-(h*s/2*l_sk)/I_zz
    # def qb_4_2(s):
    #     return -Vz*t_skin*(eta-s)/I_yy
    #
    # def qb_6_1(s):
    #     return -Vy*t_skin*(h/2)*-np.sin(s)*(h/2)/I_zz
    # def qb_6_2(s):
    #     return -Vz*t_skin*(eta+(h/2)*(1-np.cos(s)))*(h/2)/I_yy

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

    perimeter_old = perimeter_addition
    perimeter_addition += m.hypot(dx, y[i]-y[i-1])
    if perimeter_addition > spacing and perimeter_old < spacing and stiff_counter <= 10:
        if perimeter_addition - spacing <= perimeter_old - spacing:
            stiff_loc[stiff_counter,0] = x[i]
            stiff_loc[stiff_counter,1] = y[i]
            qb_current += A_stiff*(x[i] + y[i])
            print("stiffener added at", x[i], y[i])
        else:
            stiff_loc[stiff_counter,0] = x[i-1]
            stiff_loc[stiff_counter,1] = y[i-1]
            qb_current += A_stiff*(x[i-1] + y[i-1])
            print("stiffener added at", x[i-1], y[i-1])
        perimeter_addition -= spacing
        stiff_counter += 1

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


moment_arm_list = [[0,C+h/2],[h/2,C-h/2],[0,0],[h/2,C-h/2]]

M_b_1 = (C) * qb_1_v
M_b_2 = (C - h/2) * qb_2_v
M_b_3 = h/2 * qb_3_h + (C-h/2) * -qb_3_v
M_b_4 = 0
M_b_5 = (C - h/2) * qb_5_v
M_b_6 = h/2 * qb_6_h + (C-h/2) * qb_6_v

M_vy = (C-eta) * Vy
M_vz = 0


print(stiff_loc)
print(qb_lastval)

plt.plot(stiff_loc[:,0],stiff_loc[:,1], marker="o")
plt.plot(x,y)
plt.plot(x,-y)
plt.axis('equal')
plt.show()
