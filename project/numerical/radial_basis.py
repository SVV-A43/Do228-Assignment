import numpy as np
import sys
import os
from matplotlib import pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next two imports
from project.numerical.aileron_geometry import AileronGeometry


def dist_r(p1, p2):
    '''
    :param x1: point 1 [x1, z1]
    :param x2: point 2 [x2, z2]
    :return:
    '''

    return np.sqrt(((p2 - p1)**2).sum())

def _euclidean_norm(x1, x2):
    return np.sqrt(((x1 - x2) ** 2).sum(axis=0))

def phi(r, basis):
    '''
    Radial basis of euclidean distance
    :param r: euclidean distance
    :param basis: 'linear'
    :return:
    '''
    if basis=='linear':
        return r
    else:
        raise NotImplementedError

def coeffs(X, Z, F):
    X = X.flatten()
    Z = Z.flatten()
    coords = np.array([X, Z])
    n = len(X)

    A = np.zeros((n,n))

    for i in range(n):
        for j in range(n):
            r = dist_r(coords[:,i], coords[:,j])
            A[i][j] = phi(r, 'linear')

    coefficients = np.linalg.solve(A, F)

    return coefficients, coords, A

def interpolate(X_in, Z_in, known_coords, coeffs):
    X_in = np.asarray(X_in).flatten()
    Z_in = np.asarray(Z_in).flatten()
    # Array of points to interpolate, row 0 = X values, row 1 = Y values
    pts_in = np.array([X_in, Z_in])

    n = len(coeffs) # Number of terms in the final interpolant

    phi_x = np.zeros((len(X_in), n)) # Basis terms

    for i in range(len(X_in)):
        for j in range(n):
            r = dist_r(pts_in[:,i], known_coords[:,j])
            phi_x[i,j] = phi(r, 'linear')



    F = np.dot(phi_x, coeffs)

    return F




if __name__ == '__main__':
    from scipy.interpolate import Rbf
    x, z, d = np.random.rand(3, 50)

    interp_coef, pt_coords, A = coeffs(x, z, d)

    rbfi = Rbf(x, z, d, function='linear')  # radial basis function interpolator instance
    xi = zi = np.linspace(0, 1, 20)
    # xi = zi = 0.5
    di = rbfi(xi, zi)  # interpolated values
    fi = interpolate(xi, zi, pt_coords, interp_coef)

    print(x)