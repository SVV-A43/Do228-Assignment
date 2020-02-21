import numpy as np
import os
import warnings
import sys
from tqdm import tqdm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))  # This must come before the next imports
from project.numerical.Loading.aileron_geometry import AileronGeometry
from project.numerical.Loading.interpolation import InterpolateRBF, select_station



def def_integral(fn, start, stop, num_bins=1000):
    '''
    Definite numerical integration of a 1 variable function
    :param start: start coordinate
    :param stop: end coordinate
    :param num_bins: number of bins to use
    :param fn: function of a single variable
    :return:
    '''
    steps = np.linspace(start, stop, num_bins)
    width = steps[1] - steps[0] # Uniform bin width

    # Change this to general function
    fi = fn(steps)

    areas = (fi[:-1] + fi[1:]) / 2 * width
    return areas.sum()

def indef_integral(fn, start, stop, num_bins=1000):
    '''
    Definite numerical integration of a 1 variable function
    :param start: start coordinate
    :param stop: end coordinate
    :param num_bins: number of bins to use
    :param fn: function of a single variable
    :return:
    '''
    steps = np.linspace(start, stop, num_bins)
    width = steps[1] - steps[0] # Uniform bin width

    fi = []
    for x in steps:
        fi.append(def_integral(fn, start, x, num_bins=num_bins))

    fi = np.array(fi)

    areas = (fi[:-1] + fi[1:]) / 2 * width
    return areas.sum()

def station_loads(**kwargs):
    num_bins = kwargs.pop('num_bins', 1000)
    q_x = []
    x_coord = []
    ail = AileronGeometry()
    for station in tqdm(range(ail.span_stations)):
        x, z, p = ail.station_data(station)
        int_fn = InterpolateRBF(z, p)
        station_load = def_integral(int_fn.interpolate, min(z), max(z), num_bins=num_bins)
        q_x.append(station_load)
        x_coord.append(x[0])
    q_x = np.array(q_x)
    x_coord = np.array(x_coord)

    return x_coord, q_x

def check_fn_simplification():
    warnings.warn('Using LOW BIN RESOLUTION for def_integration, update for final model')
    x_coords, q = station_loads(num_bins=10)

    int_fn = InterpolateRBF(x_coords, q)

    a = int_fn.coefficients

    def lin_fn(x):
        f = (x*a).sum() - np.dot(x_coords, a.T)
        return f

    pt = 0.5
    fi = lin_fn(pt)
    qi = int_fn.interpolate(pt)


    print( qi == fi)




if __name__ == '__main__':



    ail = AileronGeometry()
    # min_z = min(ail.station_z_coords())
    # max_z = max(ail.station_z_coords())
    #
    station_id = 5
    x, z, p = ail.station_data(station_id)
    # xi = [0 for i in range(len(x))]

    int_fn = InterpolateRBF(z, p)

    second_int = indef_integral(int_fn.interpolate, 0, -0.5, num_bins=10)

    #
    #
    # area = def_integral(min_z, max_z, 100, x_coord=x[0], fn=interpolant, fn2=rbfi)

    # station_loads(100)
    check_fn_simplification()
