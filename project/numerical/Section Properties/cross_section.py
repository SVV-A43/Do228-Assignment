from math import sqrt,hypot
import numpy as np
import matplotlib.pyplot as plt

C = 0.515
h = 0.248
dx = 0.00001

x = np.arange(0,C+dx,dx)
sz = np.size(x)
y = np.zeros(sz)

perimeter = 0
for i in range(sz):
    if x[i] <= h/2:
        y[i] = sqrt((h/2)**2-(x[i]-h/2)**2)
    else:
        y[i] = -((h/2)/(C-h/2))*x[i] + (((h/2)/(C-h/2))*0.515)
    perimeter += hypot(dx, y[i]-y[i-1])

print(perimeter)

plt.plot(x,y)
plt.axis('equal')
plt.show()

