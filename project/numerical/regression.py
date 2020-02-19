import numpy as np


def _euclidean_norm(x1, x2):
    return np.sqrt(((x1 - x2) ** 2).sum(axis=0))

def _dist_r(x1, x2):
    '''
    :param x1: array vector of columns x, y
    :param x2: array vector of columns x, y
    :return:
    '''
    return np.sqrt((x1[0] - x2[0]) ** 2 + (x1[1] - x2[1]) ** 2)

if __name__ == '__main__':
    p1 = np.array([2,3])
    p2 = np.array([1,1])

    mine = _dist_r(p1, p2)
    theirs = _euclidean_norm(p1, p2)

    print(mine == theirs)