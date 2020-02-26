import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

#open and read the table with x, y and z location of all the nodes used for the input data
nodes = pd.read_csv("B737.inp", header=None, skiprows=9, nrows=6588, error_bad_lines=False)
elements = pd.read_csv("B737.inp", header=None, skiprows=6598, nrows=6588, error_bad_lines=False)

#open and read the table with displacements in x, y and z direction of all nodes
displacement = pd.read_fwf("B737.rpt", header=None, skiprows=20074, nrows=6588, error_bad_lines=False)

data = np.genfromtxt("B737.rpt", skip_header=20074, max_rows=6588) #Gebruik dit maar

#rename the colums of the input and output tables
nodes.columns = ['node', 'x', 'y', 'z']
displacement.columns = ['node', 'd_tot', 'dx', 'dy', 'dz']

#create dataframe for which contain only the nodes on the leading edge
nodes_le = nodes[nodes.z.eq(102.5)]

#make a Series of True and False values in order to filter the dataframe 
#containing the displacements of the nodes, so that we are left with a dataframe
#only containing the displacents of the nodes on the leading edge
checknodes = displacement.node.isin(nodes_le['node']) 

#making dataframe which only contains the displacements of the nodes on the leading edge
displacement_le = displacement[checknodes]


#plotting the nodes to check if all nodes create the aileron
x = nodes['x']
y = nodes['y']
z = nodes['z']

x_le = nodes_le['x']
y_le = nodes_le['y']
z_le = nodes_le['z']


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_le, y_le, z_le)
