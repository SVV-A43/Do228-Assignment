# Imports
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next two imports
from project.numerical.aileron_geometry import AileronGeometry


# Attempt 1 Radial Basis Function - Gaussian Kernel

class InterpolationRBF():
    def __init__(self, data_xyf, basis='linear', epsilon=None, coeffs_path=None):
        '''
        :type data_xyf: list
        :param data_xyf: List of tuples of 3D data [(x1, y1, f1), (x2, y2, f2), etc...]
        :type epsilon: float
        :param epsilon: Constant to control shape of interp fn, closer to 0 = flatter radial fn
        :param basis: Select from: ('gaussian', 'inv_quad', 'inv_multi_quad', 'multi_quad', 'linear')
        '''
        self.__data = data_xyf
        self.__basis = basis
        if epsilon is None:
            self.__epsilon = 0.1  # TODO find an appropriate epsilon
        else:
            self.__epsilon = epsilon

        # Calculate coefficients of interpolating function
        if coeffs_path is None:
            self.__coeffs_calc()
        else:
            self.__coeffs = np.genfromtxt(coeffs_path, delimiter=',')



    def __phi(self, x, xj):
        '''
            RBF Basis function
            :type x: tuple
            :param x: Vector (x,y) you are interested in
            :type xj: tuple
            :param xj: data point at vector (x_j, y_j)
            :return: Basis funciton
            '''
        r = np.sqrt((x[0] - xj[0]) ** 2 + (x[1] - xj[1]) ** 2)

        if self.__basis == 'gaussian':
            return np.e ** ((-1 * self.__epsilon * r) ** 2)
        elif self.__basis == 'inv_quad':
            return 1 / (1 + (self.__epsilon * r)**2 )
        elif self.__basis == 'inv_multi_quad':
            return 1 / np.sqrt(1 + (self.__epsilon * r)**2 )
        elif self.__basis == 'multi_quad':
            return np.sqrt(1 + (self.__epsilon * r)**2 )
        else: # Linear
            return r

    def __coeffs_calc(self):
        '''
        Calculates coefficients of interpolating function using Gaussian RBF
        :return: array of coefficients
        '''
        n = len(self.__data)
        A = np.zeros((n, n))
        f = np.zeros((n, 1))

        for i in range(n):
            f[i] = self.__data[i][2]
            for j in range(n):
                A[i, j] = self.__phi(self.__data[i], self.__data[j])

        self.__coeffs = np.linalg.solve(A, f)

        coef_fname = 'rbf_coefficients_' + self.__basis + '.csv'
        np.savetxt(coef_fname, self.__coeffs, delimiter=',')


    def function_coefficients(self):
        return self.__coeffs

    def interpolate(self, point_coords):
        '''
        Interpolating Polynomial for the RBF
        :type point_coords: tuple
        :param point_coords: Vector (x,y) you are interested in
        :return: Interpolated value
        '''
        s = 0
        for i in range(len(self.__coeffs)):
            s += self.__coeffs[i][0] * self.__phi(point_coords, self.__data[i])

        return s


def load_press_data():
    filename = './aero_loading_data/aerodynamicloaddo228.dat'
    loading_data = np.genfromtxt(filename, delimiter=',')

    geometry = AileronGeometry()
    loading_coords = geometry.loading_data_coords()

    # Flatten array to list using method 'F', which is column-major order. This results in the first 81 entries being from the first cross section
    data_fl = list(loading_data.flatten('F'))
    coords_fl = list(loading_coords.flatten('F'))

    combined = list(zip(coords_fl, data_fl))

    # Unpack nested tuples, and re-order coordinates to more logical order
    coords_data_combined = [(x, z, f) for (z, x), f in combined]

    return coords_data_combined




if __name__ == '__main__':
    # data = data_to_xyf()
    data = load_press_data()

    # Initialize Interpolation Object
    inter_fn = InterpolationRBF(data, basis='linear')

    print(inter_fn.function_coefficients())

    # Test an arbitrary point:
    pt  = (data[1][0], data[1][1])
    pt2 = (data[0][0], data[0][1] - (data[0][1] - data[1][1])/2)

    val = inter_fn.interpolate(pt)
    val2 = inter_fn.interpolate(pt2)
    print(val)


