# Imports
import numpy as np


# CODE...


class geometry():
    def __init__(self, Nx=41, Nz=81, la=2.961):
        self.__num_span_bins = Nx
        self.__num_chord_bins = Nz
        self.__la = la

        self.__span_stations = []
        self.__chord_stations = []

        # Initialize x-coords of stations along span
        for i in range(self.__num_span_bins):
            self.__span_stations.append(self.xi(i))

    def span_stations(self):
        return self.__span_stations

    def chord_stations(self):
        return self.__chord_stations

    # Functions to calculate geometry
    def th_xi(self, i):
        return (i - 1) / self.__num_span_bins * np.pi

    def xi(self, i):
        return .5 * (self.__la / 2 * (1 - np.cos(self.th_xi(i))) + self.__la / 2 * (1 - np.cos(self.th_xi(i + 1))))


if __name__ == '__main__':
    aileron = geometry()
    print(aileron.span_stations())