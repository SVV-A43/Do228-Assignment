# Main Structural Analysis model we need to develop

# Imports
import numpy as np
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next imports
from project.numerical.Loading.interpolation import InterpolateRBF
from project.numerical.reaction_forces import reaction_forces
from project.numerical.Loading.aileron_geometry import AileronGeometry


# CODE...

def vec_local2global(local_vec):
    '''
    :param local_vec: array [x',y',z']
    :type local_vec: np.ndarray
    :return: [x,y,z]
    '''
    in_vec = np.asarray(local_vec)
    if in_vec.ndim == 1:
        in_vec = np.expand_dims(in_vec, axis=1)
    G = AileronGeometry()

    R = np.array([[1, 0, 0],
                  [0, np.cos(G.theta), -1*np.sin(G.theta)],
                  [0, np.sin(G.theta), np.cos(G.theta)]])
    out_vec = R@in_vec

    return np.squeeze(out_vec)



if __name__ == '__main__':
    t_vec = [1,2,3]
    rotated = vec_local2global(t_vec)
    print(rotated)

