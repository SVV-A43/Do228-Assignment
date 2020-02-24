import math as m
import numpy as np
import matplotlib.pyplot as plt
from project.numerical.Loading.integration import def_integral, variable_integral




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
    return s

def_integral(fn_test1, start, end, num_bins=100)

qb_2_Vy = value_integrate(0,h/2,[[-Vy*t_skin/I_zz],[1]])
qb_2_Vz = value_integrate(0,h/2,[[-Vz*t_skin*eta/I_zz],[0]])
qb_2 = Func(qb_2_Vy[0],qb_2_Vy[1]).__add__(qb_2_Vz)
qb_2_val = qb_2.get_const()
qb_2_var = qb_2.get_var()
print(qb_2)
print(qb_2_val,qb_2_var)

# Shear flow 1,3,4,6
def qb_2_1(s):
    return s
def qb_2_2(s):
    return s

def qb_2_1(s):
    return s
def qb_2_2(s):
    return s

def qb_2_1(s):
    return s
def qb_2_2(s):
    return s

def qb_2_1(s):
    return s
def qb_2_2(s):
    return s

# shear flow at 5
def qb_2_1(s):
    return s
def qb_2_2(s):
    return s


# # lists are ordered as 1,3,4,6,5, which the numbers are given in cross section from verification model
s_list = [h/2, m.pi/2]  #list containing lengths of each skin
func_list = [func1,func2  ]
final_qb_values_list = []
final_qb_variables_list = []
perimeter_addition = 0
new_qb = Func([],[])
s_current = 0
seg_i = 0
s_last = -10**4
for i in range(sz):  # go through the outer section
    if x[i] <= h/2:
        y[i] = m.sqrt((h/2)**2-(x[i]-h/2)**2)
    else:
        y[i] = -((h/2)/(C-h/2))*x[i] + (((h/2)/(C-h/2))*C)

    perimeter_old = perimeter_addition
    perimeter_addition += m.hypot(dx, y[i]-y[i-1])
    s_current += m.hypot(dx, y[i]-y[i-1])
    if perimeter_addition > spacing > perimeter_old:
        if perimeter_addition - spacing <= perimeter_old - spacing:
            new_qb.__add__([[A_stiff*y[i]],[0]])
        else:
            new_qb.__add__([[A_stiff*y[i-1]],[0]])

        perimeter_addition -= spacing
        new_qb.__add__(qb_i(s_current, s_current+dx, func_list[i]))
    else:
        new_qb.__add__(qb_i(s_current, s_current+dx, func_list[i]))

        perimeter_addition -= spacing
    if s_current - s_list[seg_i] < dx*2 and s_current - s_last > dx*4:
        s_last = s_current
        seg_i += 1
        final_qb_values_list.append(new_qb.get_func()[0])
        final_qb_variables_list.append(new_qb.get_func()[1])


print(final_qb_variables_list, final_qb_values_list)