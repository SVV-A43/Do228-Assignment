import numpy as np
import sys
import os
from tqdm import tqdm
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))  # This must come before the next imports
from project.numerical.Loading.aileron_geometry import AileronGeometry


class InterpolateRBF():
    def __init__(self, *data_arrays, basis='linear', coeff_path=None):
        '''
        :param X: 1D array of x-coordinates
        :param Z: 1D array of z-coordinates
        :param F: 1D array of data at known data points
        :param basis: 'linear'
        '''
        self._basis = basis

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
        if self._basis=='linear':
            return r
        else:
            raise NotImplementedError

    def _compute_coeffs(self):
        # self._known_coords = np.array([np.asarray(X), np.asarray(Z)])
        n = len(self._known_coords[0,:])

        A = np.zeros((n,n))

        for i in range(n):
            for j in range(n):
                r = self._dist_r(self._known_coords[:, i], self._known_coords[:, j])
                A[i][j] = self._phi(r)

        self.coefficients = np.linalg.solve(A, self._known_data)
        # np.savetxt('rbf_coefficients_linear.csv', self.coefficients, delimiter=',')


    def interpolate(self, *coord_arrays):
        coord_ls = []
        for ar in coord_arrays:
            coord_ls.append(np.asarray(ar))
        pts_in = np.array(coord_ls)
        #
        # if isinstance(X_in, float) or isinstance(X_in, int):
        #     X_in = [X_in]
        # if isinstance(Z_in, float) or isinstance(Z_in, int):
        #     Z_in = [Z_in]
        # X_in = np.asarray(X_in)
        # Z_in = np.asarray(Z_in)
        # # Array of points to interpolate, row 0 = X values, row 1 = Y values
        # pts_in = np.array([X_in, Z_in])

        n = len(self.coefficients) # Number of terms in the final interpolant

        self.phi_x = np.zeros((len(pts_in[0,:]), n)) # Basis terms

        for i in range(len(pts_in[0, :])):
            for j in range(n):
                r = self._dist_r(pts_in[:, i], self._known_coords[:, j])
                self.phi_x[i,j] = self._phi(r)

        F = np.dot(self.phi_x, self.coefficients)

        return F


def select_station(id):
    aileron = AileronGeometry()
    x, z, p = aileron.data_x_z_p()

    start = id*81
    end = (id+1)*81
    if end > p.shape[0]:
        end = -1

    return x[start:end], z[start:end], p[start:end]



if __name__ == '__main__':
    from scipy.interpolate import Rbf

    aileron = AileronGeometry()
    # x, z, p = aileron.data_x_z_p()

    # Use single station
    x, z, p = select_station(6)




    # Create instance of our interpolation class
    my_i = InterpolateRBF(x, z, p) #, coeff_path='rbf_coefficients_linear.csv'

    # Create reference radial basis function
    # rbfi = Rbf(x, z, p, function='linear')

    # Generate points to test
    zi = np.linspace(0.001, -0.5, len(z))
    xi = np.array([x[0] for i in range(len(x))])
    # xi = [1.5]
    # zi = [-0.1]

    # di = rbfi(xi, zi)  # interpolated values
    # print(di)

    fi = my_i.interpolate(xi, zi)
    print(fi)

    # print(f'Check for equal: {(fi==di).all()}')
    #
    plt.scatter(z, p, c='green')
    # plt.plot(zi, di, c='red')
    plt.plot(zi, fi, c='blue')
    plt.show()
    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.scatter3D(x, z, p, c=p, cmap='Greens')
    # ax.plot3D(xi, zi, di, 'red')
    # ax.plot3D(xi, zi, fi, 'blue')
    # # ax.plot_trisurf(x, z, P_arr, cmap='viridis', edgecolor='none')
    # ax.set_xlabel('Spanwise coordinate x [m]')
    # ax.set_ylabel('Chordwise coordinate z [m] ')
    # ax.set_zlabel('Aerodynamic Load p [kN/m^2]')
    # ax.view_init(10, 45)
    # ax.set_title('Aerodynamic loading')
    # plt.show()
