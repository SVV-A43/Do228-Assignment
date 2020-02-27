import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cmx


#open and read the table with x, y and z location of all the nodes used for the input data
nodes = pd.read_csv("B737.inp", header=None, skiprows=9, nrows=6588, error_bad_lines=False)
elements = pd.read_csv("B737.inp", header=None, skiprows=6598, nrows=6588, error_bad_lines=False)

#open and read the table with displacements in x, y and z direction of all nodes
data = np.genfromtxt("B737.rpt", skip_header=6705, max_rows=5778)
stress = pd.DataFrame(data)

#rename the colums of the input and output tables
nodes.columns = ['node', 'x', 'y', 'z']
stress.columns = ['element', 'int', 'mises1', 'mises2', 's12_1', 's12_2']
elements.columns = ['element_number', 'node 1', 'node 2', 'node 3', 'node 4']

avgstress = pd.concat([stress['element'], (stress['mises1']+stress['mises2'])/2], axis=1)
avgstress.columns = ['element', 'stress']

check = elements.element_number.isin(avgstress['element']) 
new = elements[check]

test = elements.values.tolist()
ncoor = nodes.values.tolist()
x = []
y = []
z = []

for i in range(len(elements)):
    j = 1
    
    xlocal = []
    ylocal = []
    zlocal = []
    
    while j < 5:
        loc = ncoor[test[i][j]-1]
        xlocal.append(loc[1])
        ylocal.append(loc[2])
        zlocal.append(loc[3])
        j = j+1
    
    x.append(sum(xlocal)/4)
    y.append(sum(ylocal)/4)
    z.append(sum(zlocal)/4)




fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, z, y, s=8, cmap='hot', marker = 'o', c=avgstress['stress'])
ax.view_init(25,150)
ax.set_ylim3d(-502.5,102.5)
ax.set_zlim3d(-200,200)

