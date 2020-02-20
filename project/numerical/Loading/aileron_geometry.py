# Imports
import numpy as np


# CODE...


class AileronGeometry():
    # Functions to calculate geometry
    def __th_xi(self, i):
        return (i - 1) / self.__num_pts_span * np.pi

    def __xi(self, i):
        return .5 * (self.__la / 2 * (1 - np.cos(self.__th_xi(i))) + self.__la / 2 * (1 - np.cos(self.__th_xi(i + 1))))

    def __th_zi(self, i):
        return (i - 1) / self.__num_pts_chord * np.pi

    def __zi(self, i):
        return -.5 * (self.__Ca / 2 * (1 - np.cos(self.__th_zi(i))) + self.__Ca / 2 * (1 - np.cos(self.__th_zi(i + 1))))





    def __init__(self, filename='./aero_loading_data/aerodynamicloaddo228.dat'):
        self.__Ca = 0.515
        self.__la = 2.961

        self.__pressure = np.genfromtxt(filename, delimiter=',')
        self.__num_pts_span = len(self.__pressure[0, :])
        self.__num_pts_chord = len(self.__pressure[:, 0])

        self.__x_coords = []
        self.__z_coords = []


        # Calculate x-point_coords of stations along span
        for i in range(self.__num_pts_span):
            self.__x_coords.append(self.__xi(i+1))

        # Calculate z-point coords of station along chord
        for j in range(self.__num_pts_chord):
            self.__z_coords.append(self.__zi(j+1))

        num_data_pts = self.__num_pts_span * self.__num_pts_chord
        self.load_data = np.zeros((num_data_pts, 3))

        n = 0
        for i, x in enumerate(self.__x_coords):
            for j, z in enumerate(self.__z_coords):
                self.load_data[n, 0] = x
                self.load_data[n, 1] = z
                self.load_data[n, 2] = self.__pressure[j,i]
                n += 1


        # Create Array of the (z,x coordinates where the load data is found) LOCAL COORD SYSTEM AS DEFINED IN READER fig3
        # The order [z, x] here is to match that of the loading data, which is also z-coord per row
        self.__load_coords = np.empty((self.__num_pts_chord, self.__num_pts_span), dtype=object)
        for z in range(self.__num_pts_chord):
            for x in range(self.__num_pts_span):
                self.__load_coords[z, x] = (self.__z_coords[z], self.__x_coords[x])

    def station_x_coords(self):
        return self.__x_coords

    def station_z_coords(self):
        return self.__z_coords

    def loading_data_coords(self):
        return self.__load_coords

    def data_x_z_p(self):
        return self.load_data[:,0], self.load_data[:,1], self.load_data[:,2]



if __name__ == '__main__':
    aileron = AileronGeometry()
    print(aileron.station_x_coords())