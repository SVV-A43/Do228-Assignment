import math as m
import numpy as np
import matplotlib.pyplot as plt

C = 0.515
h = 0.248
n = 11
w_st = 0.03
h_st = 0.015
t_st = 0.0012
t_skin = 0.0011
t_spar = 0.0022

dx = 0.0001

x = np.arange(0,C+dx,dx)
sz = np.size(x)
y = np.zeros(sz)
stiff_loc = np.zeros((n,2))
stiff_counter = 1

# Stiffener locations and general shape

perimeter = m.pi*(h/2) + 2*m.sqrt((h/2)**2+(C-h/2)**2)
spacing = perimeter/n
perimeter_addition = 0
for i in range(sz):
    if x[i] <= h/2:
        y[i] = m.sqrt((h/2)**2-(x[i]-h/2)**2)
    else:
        y[i] = -((h/2)/(C-h/2))*x[i] + (((h/2)/(C-h/2))*C)

    perimeter_old = perimeter_addition
    perimeter_addition += m.hypot(dx, y[i]-y[i-1])
    if perimeter_addition > spacing and perimeter_old < spacing:
        if perimeter_addition - spacing <= perimeter_old - spacing:
            stiff_loc[stiff_counter,0] = x[i]
            stiff_loc[stiff_counter,1] = y[i]
        else:
            stiff_loc[stiff_counter,0] = x[i-1]
            stiff_loc[stiff_counter,1] = y[i-1]
        perimeter_addition -= spacing
        stiff_counter += 1

for i in range (int((n-1)/2)):
    stiff_loc[n-(i+1),0] = stiff_loc[i+1,0]
    stiff_loc[n-(i+1),1] = -stiff_loc[i+1,1]

# Centroid

A_stiff = w_st*t_st + h_st*t_st
A_total = n*A_stiff + h*t_spar + 2*m.sqrt((h/2)**2+(C-h/2)**2)*t_skin + m.pi*h/2*t_skin

stiff_cont = 0
for i in range (n):
    stiff_cont += A_stiff*stiff_loc[i,0]
spar_cont = h*t_spar * h/2
semi_cont = m.pi*h/2*t_skin * (h/2-h/m.pi)
straight_cont = 2* (m.sqrt((h/2)**2+(C-h/2)**2)*t_skin) * (h/2+(C-h/2)/2)

x_centroid = (stiff_cont+spar_cont+semi_cont+straight_cont)/A_total

# Moments of inertia


# Shear centre
"""
Due to symmetry Izy = 0. Shear centre must be located on the z-axis. To find the shear center location eta in the z-axis
A unit load Sy and Sz = 0 is applied. One cut is made in the spar and one in the circular section and qb = 0 at the cuts


"""
I_zy = 0
I_zz = 1
I_yy = 1
Vy = 1
Vx = 0

def variable_integrator(func):
    """func is a list with lists containing the constants in the first list, and coefficients in the second list
    if coefficient is not a number, i.e. sin it should be written as s, if it cosine it is c"""
    constants = func[0]
    coefficients = func[1]
    for i in range(len(coefficients)):
        if type(coefficients[i]) == int:
            constants[i] *= 1/(coefficients[i]+1)
            coefficients[i] += 1

        elif coefficients[i] == 's':
            constants[i] *= -1
            coefficients[i] = 'c'

        elif coefficients[i] == 'c':
            constants[i] *= 1
            coefficients[i] = 's'
    return [constants,coefficients]

def value_integrate(a,b,func):
    """Integrates a function going from a to b"""
    constants = func[0]
    coefficients = func[1]
    result = 0
    for i in range(len(coefficients)):
        if type(coefficients[i]) == int:
            constants[i] *= 1/(coefficients[i]+1)
            coefficients[i] += 1
            result += (constants[i]*b**coefficients[i] - constants[i]*a**coefficients[i])

        elif coefficients[i] == 's':
            constants[i] *= -1
            coefficients[i] = 'c'
            result += constants[i]*(m.cos(b) - m.cos(a))


        elif coefficients[i] == 'c':
            constants[i] *= 1
            coefficients[i] = 's'
            result += constants[i]*(m.sin(b) - m.sin(a))
    return result


a = variable_integrator([[2,4],[1,2]]) #2s + 4s^2
b = value_integrate(0,m.pi/2,[[2],['c']]) #2costheta from 0 to pi/2
# qb of section 1:
# def qb_in_variable_form():


# def res_stringer(s, stiff_locs):
#     sum = 0
#     for stringer in stiff_locs:
#         if s > stringer[si]:
#             sum += stringer[area] * stringer[yi]
#
#
# qb = np.zeros(sz)  # contains all qbs for each
# for i in range(1, sz):
#     qb[i] = -(1/Izz)*( t*num_int(y[i] + res_stringer(s, stiff_locs)) + qb[i-1] )

# Moments of inertia


# plt.plot(stiff_loc[:,0],stiff_loc[:,1], marker="o")
# plt.plot(x,y)
# plt.plot(x,-y)
# plt.plot(x_centroid,0, marker="o")
# plt.axis('equal')
# plt.show()


print("perimeter =", perimeter)
print(stiff_loc)
print(A_total)
print(x_centroid)

