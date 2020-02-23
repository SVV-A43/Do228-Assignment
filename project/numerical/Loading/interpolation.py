import os
import sys
import numpy as np
from matplotlib import pyplot as plt

# This must come before the next imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))
from project.numerical.Loading.aileron_geometry import AileronGeometry


class InterpolateRBF():
    def __init__(self, *data_arrays, basis='linear', coeff_path=None,
                 save_path=None):
        '''
        :param X: 1D array of x-coordinates
        :param Z: 1D array of z-coordinates
        :param F: 1D array of data at known data points
        :param basis: 'linear'
        '''
        self._basis = basis
        self.__save_path = save_path
        self.phi_x = None

        # Save initial data
        coords = []
        for ar in data_arrays[:-1]:
            coords.append(np.asarray(ar))
        self._known_coords = np.array(coords)
        self._known_data = np.asarray(data_arrays[-1])

        # Compute coefficients:
        if coeff_path is None:
            self._compute_coeffs()
        else:
            self.coefficients = np.genfromtxt(coeff_path)


    def _dist_r(self, p1, p2):
        '''
        :param x1: point 1 [x1, z1]
        :param x2: point 2 [x2, z2]
        :return:
        '''
        return np.sqrt(((p2 - p1)**2).sum())


    def _phi(self, r):
        '''
        Radial basis of euclidean distance
        :param r: euclidean distance
        :param basis: 'linear'
        :return:
        '''
        if self._basis == 'linear':
            return r

        raise NotImplementedError

    def _compute_coeffs(self):
        # self._known_coords = np.array([np.asarray(X), np.asarray(Z)])
        dim_n = len(self._known_coords[0, :])

        mat_A = np.zeros((dim_n, dim_n))

        for i in range(dim_n):
            for j in range(dim_n):
                r = self._dist_r(self._known_coords[:, i],
                                 self._known_coords[:, j])
                mat_A[i][j] = self._phi(r)

        self.coefficients = np.linalg.solve(mat_A, self._known_data)
        if self.__save_path:
            np.savetxt(self.__save_path, self.coefficients, delimiter=',')


    def interpolate(self, *coord_arrays):
        coord_ls = []
        for ar in coord_arrays:
            if isinstance(ar, (float, int)):
                ar = [ar]
            coord_ls.append(np.asarray(ar))
        pts_in = np.array(coord_ls)

        n = len(self.coefficients) # Number of terms in the final interpolant

        self.phi_x = np.zeros((pts_in.shape[1], n)) # Basis terms

        for i in range(len(pts_in[0, :])):
            for j in range(n):
                r = self._dist_r(pts_in[:, i], self._known_coords[:, j])
                self.phi_x[i, j] = self._phi(r)

        F = np.dot(self.phi_x, self.coefficients)

        return F


def select_station(station_id):
    '''
    :param station_id: station number
    :return:
    '''
    aileron = AileronGeometry()
    x, z, p = aileron.data_x_z_p()

    start = station_id * 81
    end = (station_id + 1) * 81
    if end > p.shape[0]:
        end = -1

    return x[start:end], z[start:end], p[start:end]




def main():
    # x_coords, z_coords, p_vals = aileron.data_x_z_p()

    # Use single station
    x_coords, z_coords, p_vals = select_station(6)

    # Create instance of our interpolation class
    my_i = InterpolateRBF(x_coords, z_coords, p_vals)  # , coeff_path='rbf_coefficients_linear.csv'

    # Create reference radial basis function
    # rbfi = Rbf(x_coords, z_coords, p_vals, function='linear')

    # Generate points to test
    z_i = np.linspace(0.001, -0.5, len(z_coords))
    x_i = np.array([x_coords[0] for i in range(len(x_coords))])
    # xi = [1.5]
    # zi = [-0.1]

    # di = rbfi(xi, zi)  # interpolated values
    # print(di)

    f_i = my_i.interpolate(x_i, z_i)
    print(f_i)

    # print(f'Check for equal: {(fi==di).all()}')
    #
    plt.scatter(z_coords, p_vals, c='green')
    # plt.plot(zi, di, c='red')
    plt.plot(z_i, f_i, c='blue')
    plt.show()
    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.scatter3D(x_coords, z_coords, p_vals, c=p_vals, cmap='Greens')
    # ax.plot3D(xi, zi, di, 'red')
    # ax.plot3D(xi, zi, fi, 'blue')
    # # ax.plot_trisurf(x_coords, z_coords, P_arr, cmap='viridis', edgecolor='none')
    # ax.set_xlabel('Spanwise coordinate x_coords [m]')
    # ax.set_ylabel('Chordwise coordinate z_coords [m] ')
    # ax.set_zlabel('Aerodynamic Load p_vals [kN/m^2]')
    # ax.view_init(10, 45)
    # ax.set_title('Aerodynamic loading')
    # plt.show()


if __name__ == '__main__':

    main()
