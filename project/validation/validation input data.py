import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

nodes = pd.read_csv("B737.inp", header=None, skiprows=9, nrows=6588, error_bad_lines=False)
elements = pd.read_csv("B737.inp", header=None, skiprows=6598, nrows=6588, error_bad_lines=False)

nodes.columns = ['node', 'x', 'y', 'z']

nodes_le = nodes[nodes.z.eq(0)]

x = nodes['x']
y = nodes['y']
z = nodes['z']

x_le = nodes_le['x']
y_le = nodes_le['y']
z_le = nodes_le['z']
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z)
