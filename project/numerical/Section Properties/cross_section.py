import math as m
import numpy as np
import matplotlib.pyplot as plt

C = 0.515
h = 0.248
dx = 0.001

x = np.arange(0,C+dx,dx)
sz = np.size(x)
y = np.zeros(sz)

total_perimeter = m.pi*(h/2) + 2*m.sqrt((h/2)**2+(C-h/2)**2)
spacing = total_perimeter/11
perimeter = 0
for i in range(sz):
    if x[i] <= h/2:
        y[i] = m.sqrt((h/2)**2-(x[i]-h/2)**2)
    else:
        y[i] = -((h/2)/(C-h/2))*x[i] + (((h/2)/(C-h/2))*C)
    perimeter += m.hypot(dx, y[i]-y[i-1])

    if perimeter % spacing == 0:
        stiffener_loc = (x[i],y[i])

#Moment of inertia



plt.plot(x,y)
plt.axis('equal')
plt.show()

print("perimeter =", perimeter)
