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

def res_stringer(s, stiff_locs):
    sum = 0
    for stringer in stiff_locs:
        if s > stringer[si]:
            sum += stringer[area] * stringer[yi]


qb = np.zeros(sz)  # contains all qbs for each
for i in range(1, sz):
    qb[i] = -(1/Izz)*( t*num_int(y[i] + res_stringer(s, stiff_locs)) + qb[i-1] )

# Moments of inertia

'''    
plt.plot(stiff_loc[:,0],stiff_loc[:,1], marker="o")    
plt.plot(x,y)
plt.plot(x_centroid,0, marker="o")
plt.axis('equal')
plt.show()
'''

print("perimeter =", perimeter)
print(stiff_loc)
print(A_total)
print(x_centroid)

