import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

#open and read the table with x, y and z location of all the nodes used for the input data
nodes = pd.read_csv("B737.inp", header=None, skiprows=9, nrows=6588, error_bad_lines=False)
elements = pd.read_csv("B737.inp", header=None, skiprows=6598, nrows=6588, error_bad_lines=False)

#open and read the table with displacements in x, y and z direction of all nodes
data = np.genfromtxt("B737.rpt", skip_header=26724, max_rows=6588)
displacement = pd.DataFrame(data)

#rename the colums of the input and output tables
nodes.columns = ['node', 'x', 'y', 'z']
displacement.columns = ['node', 'd_tot', 'dx', 'dy', 'dz']

#create dataframe for which contain only the nodes on the leading edge
nodes_le = nodes[nodes.z.eq(102.5)]
nodes_hl = nodes[nodes.z.eq(0)& nodes.y.eq(0)]
nodes_te = nodes[nodes.z.eq(-502.5)]

#make a Series of True and False values in order to filter the dataframe 
#containing the displacements of the nodes, so that we are left with a dataframe
#only containing the displacents of the nodes on the leading edge
checknodes = displacement.node.isin(nodes_le['node']) 
checknodes2 = displacement.node.isin(nodes_hl['node'])
checknodes3 = displacement.node.isin(nodes_te['node'])

#making dataframe which only contains the displacements of the nodes on the leading edge
displacement_le = displacement[checknodes]
displacement_hl = displacement[checknodes2]
displacement_te = displacement[checknodes3]

#plotting the nodes to check if all nodes create the aileron
x = nodes['x']
y = nodes['y']
z = nodes['z']

x_le = nodes_le['x']
y_le = nodes_le['y']
z_le = nodes_le['z']

x_hl = nodes_hl['x']
y_hl = nodes_hl['y']
z_hl = nodes_hl['z']

x_te = nodes_te['x']
y_te = nodes_te['y']
z_te = nodes_te['z']

dy_le = displacement_le['dy']
dy_hl = displacement_hl['dy']
dy_te = displacement_te['dy']

plt.subplot(311)
plt.scatter(x_le, dy_le)
plt.title('leading edge')

plt.subplot(312)
plt.scatter(x_hl, dy_hl)
plt.title('hinge line')

plt.subplot(313)
plt.scatter(x_te, dy_te)
plt.title('trailing edge')

plt.plot()


'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_le, y_le, z_le)
'''