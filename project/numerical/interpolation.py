# Imports
import numpy as np
import sys
import os
from tqdm import tqdm
from scipy import linalg
from matplotlib import pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))  # This must come before the next two imports
from project.numerical.aileron_geometry import AileronGeometry


# Attempt 1 Radial Basis Function - Gaussian Kernel

class InterpolationRBF():
    def __init__(self, data_xyf, basis='linear', coeffs_path=None, epsilon=None,):
        '''
        :type data_xyf: list
        :param data_xyf: List of tuples of 3D loading_data_prepped [(x1, y1, f1), (x2, y2, f2), etc...]
        :type epsilon: float
        :param epsilon: Constant to control shape of interp fn, closer to 0 = flatter radial fn
        :param basis: Select from: ('gaussian', 'inv_quad', 'inv_multi_quad', 'multi_quad', 'linear')
        '''
        self.__full_data = data_xyf
        self.__data = data_xyf
        self.__basis = basis

        if epsilon is None:
            self.__epsilon = 0.1  # TODO find an appropriate epsilon
        else:
            self.__epsilon = epsilon

        # Load coefficients if given
        if coeffs_path is None:
            self._coeffs = None
        else:
            loaded_coeff = np.genfromtxt(coeffs_path, delimiter=',')
            self._coeffs = np.array([loaded_coeff]).T # Convert to a 2d array with 1 column



    def single_station_data(self, i=1):
        first_idx = (i-1)*81
        stop_idx = (i)*81
        return self.__full_data[first_idx:stop_idx]


    def interp_function_coefficients(self):
        return self._coeffs


    def data_used(self):
        return self.__data

    def _dist_r(self, x1, x2):
        '''
        :param x1: array vector of columns x, y
        :param x2: array vector of columns x, y
        :return:
        '''
        return np.sqrt((x1[0] - x2[0]) ** 2 + (x1[1] - x2[1]) ** 2)

    def _phi(self, x, xj):
        '''
            RBF Basis function
            :type x: tuple
            :param x: Vector (x,y) you are interested in
            :type xj: tuple
            :param xj: loading_data_prepped point at vector (x_j, y_j)
            :return: Basis funciton
            '''
        r = self._dist_r(x, xj)

        if self.__basis == 'gaussian':
            return np.e ** ((-1 * self.__epsilon * r) ** 2)
        elif self.__basis == 'inv_quad':
            return 1 / (1 + (self.__epsilon * r)**2 )
        elif self.__basis == 'inv_multi_quad':
            return 1 / np.sqrt(1 + (self.__epsilon * r)**2 )
        elif self.__basis == 'multiquadric':
            return np.sqrt(1 + (self.__epsilon * r)**2 )
        else: # Linear
            return r

    def compute_interp_coeffs(self, dimensions='2d', station_idx=1):
        '''
        Calculates coefficients of interpolating function
        :param station_idx: if dimentsions not 2d, integer from 1 to 41
        :return: array of coefficients
        '''
        if dimensions is not '2d':
            self.__data = self.single_station_data(station_idx)

        n = len(self.__data)
        A = np.zeros((n, n))
        f = np.zeros((n, 1))

        for i in range(n):
            f[i] = self.__data[i][2]
            for j in range(n):
                A[i, j] = self._phi(self.__data[i], self.__data[j])

        self._coeffs = np.linalg.solve(A, f).T
        # self._coeffs = linalg.solve(A, f)

        coef_fname = 'rbf_coefficients_' + self.__basis + '.csv'
        np.savetxt(coef_fname, self._coeffs, delimiter=',')

        return self._coeffs


    def interpolate(self, point_coords):
        '''
        Interpolating Polynomial for the RBF
        :type point_coords: tuple
        :param point_coords: Vector (x,y) you are interested in
        :return: Interpolated value
        '''
        s = 0
        for i in range(len(self._coeffs)): # FIXME: should be self._coeffs[i] ???
            s += self._coeffs[0][i] * self._phi(point_coords, self.__data[i])
            #TODO: check calculation of phi

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

def plt_function_1d(fn, min, max):
    pts = np.linspace(min, max, 100)


def plot_single_station_interpolation(load_data, station_id=30, epsilon=1):
    from sklearn.metrics import mean_squared_error
    from math import sqrt

    # Initialize Interpolation Object
    inter_fn = InterpolationRBF(load_data, basis='gaussian', epsilon=epsilon)

    coeffs = inter_fn.compute_interp_coeffs(dimensions='1d', station_idx=station_id)

    # ----- Plot 1D real and interpolated data ---
    # Real data
    data_used = inter_fn.data_used()
    data_used = np.array(data_used)
    # interpolated data
    # d_min = min(data_used[:, 1])
    # d_max = max(data_used[:, 1])
    x_coord = load_data[station_id * 81][0]
    # Z = np.linspace(d_min, d_max, 100)
    int_z = list(data_used[:, 1])
    int_f = [inter_fn.interpolate((x_coord, z)) for z in int_z]

    rms = sqrt(mean_squared_error(data_used[:, 2], int_f))


    plt.scatter(data_used[:, 1], data_used[:, 2], c='green')
    plt.plot(int_z, int_f)
    plt.title(f'Epsilon: {epsilon}')
    plt.show()




    return rms




def meshplot_RBF():
    loading_data_prepped = load_press_data()
    inter = InterpolationRBF(loading_data_prepped, coeffs_path='rbf_coefficients_linear.csv')
    data_known = np.array(loading_data_prepped)

    X, Z = np.meshgrid(data_known[:,0], data_known[:,1])
    F = inter.interpolate((X,Z))

    np.savetxt('linear_interp_mesh.csv', F, delimiter=',')

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.contour3D(X, Z, F, 50, cmap='binary')
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_zlabel('f')
    fig.show()


if __name__ == '__main__':
    # loading_data_prepped = data_to_xyf()
    loading_data_prepped = load_press_data()


    station_id = 1
    ep = 0.01
    base = 'linear'
    inter_fn = InterpolationRBF(loading_data_prepped, basis=base)

    coeffs = inter_fn.compute_interp_coeffs(dimensions='1d', station_idx=station_id)
    data_used = inter_fn.data_used()
    data_used = np.array(data_used)
    x_coord = loading_data_prepped[station_id * 81][0]
    Z = data_used[:, 1]
    F = np.array([inter_fn.interpolate((x_coord, z)) for z in Z])

    from sklearn.metrics import mean_squared_error
    from math import sqrt
    from scipy.interpolate import Rbf
    rbfi = Rbf(Z, data_used[:, 2], function=base, epsilon=ep)
    R = rbfi(Z)

    rms = sqrt(mean_squared_error(R, F))

    print((coeffs==rbfi.nodes).all())
    print(rms)

    plt.scatter(data_used[:, 1], data_used[:, 2], c='green')
    plt.plot(Z, F, c='blue')
    plt.plot(Z, R, c='red')
    plt.title(f'builtin comparison')
    plt.show()


    # # Test for correct epsilon
    # epsilon_range = np.linspace(2.2,2.9,20)
    # rmse_record = []
    #
    # for epsilon in tqdm(epsilon_range):
    #     rms = 0
    #     count = 0
    #     for i in range(1, 41):
    #         count += 1
    #         rms += plot_single_station_interpolation(loading_data_prepped, i, epsilon=epsilon)
    #     rmse_record.append(rms/count)
    #
    # plt.plot(epsilon_range, rmse_record)
    # plt.xlabel('epsilon')
    # plt.ylabel('RMSE')
    # plt.show()



    # meshplot_RBF()

    # print(inter_fn.interp_function_coefficients())

    # # Test an arbitrary point:
    # inter_fn = InterpolationRBF()
    # pt  = (loading_data_prepped[1][0], loading_data_prepped[1][1])
    # pt2 = (loading_data_prepped[0][0], loading_data_prepped[0][1] - (loading_data_prepped[0][1] - loading_data_prepped[1][1])/2)
    #
    # val = inter_fn.interpolate(pt)
    # val2 = inter_fn.interpolate(pt2)
    # print(f'Point1 on-point of known loading_data_prepped. Actual = {loading_data_prepped[1][2]}; Interpolated = {val}')
    #
    # estimate2 = (loading_data_prepped[0][2] + loading_data_prepped[1][2]) /2
    # print(f'Point2 mid-point between 1st and 2nd point. '
    #       f'Manually estimated (average) = {estimate2};  Interpolated = {val2}')









