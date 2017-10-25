# coding: utf-8
import sys
from numpy import array as a, matrix as m, arange, sqrt, isnan, pi, cos, sin, mean
from itertools import permutations
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt

from _hypervolume import hv


def show_pareto_front(filename='front.txt'):
    f = open(filename)
    pareto_front = []
    nadir = None
    for line in f:
        if 'nadir' not in line:
            xs, objs = line.split('(')
            pareto_front.append(xs.split() + [{'obj': [float(e.strip().strip(')')) for e in objs.split()]}])
        else:
            nadir = [float(e) for e in line.strip('nadir:').split()]

    f1 = pareto_front[0][-1]['obj'][0]
    f2 = pareto_front[0][-1]['obj'][1]
    min_f1, max_f1 = f1, f1
    min_f2, max_f2 = f2, f2

    for point in pareto_front:
        f1 = point[-1]['obj'][0]
        f2 = point[-1]['obj'][1]
        plt.plot([f1], [f2], 'bo')

        min_f1 = min([min_f1, f1])
        max_f1 = max([max_f1, f1])
        min_f2 = min([min_f2, f2])
        max_f2 = max([max_f2, f2])

    objs = [e[-1]['obj'] for e in pareto_front ]

    df1 = (max_f1 - min_f1) * 0.05
    df2 = (max_f2 - min_f2) * 0.05
    plt.xlabel('f1')
    plt.ylabel('f2')
    plt.axis([min_f1 - df1, max_f1 + df2, min_f2 - df2, max_f2 + df2])
    plt.show()


def uniformity(front):
    '''Return the uniformity of a *front*. Uniformity (UD) definition taken
    from Yu. Evtushenko article "Nonuniform covering method as applied to
    multicriteria" article.

    :param front: The Pareto front estimate consisting of non-dominated
    solutions in objective space.
    '''
    min_dists = []
    for p in front:
        min_dist = float('inf')
        for p2 in front:
            if p != p2:
                dist = sqrt(sum([(e-e2)**2 for e, e2 in zip(p, p2)]))
                if dist < min_dist:
                    min_dist = dist
        if min_dist == float('inf'):
            raise ValueError('Distance cannot be infinite')
        min_dists.append(min_dist)

    avg_dist = sum(min_dists) / len(min_dists)

    dist_squared_sum = 0
    for dist in min_dists:
        dist_squared_sum += (dist - avg_dist)**2

    return sqrt(dist_squared_sum)


def print_hv(filename='front.txt'):
    f = open(filename)
    pareto_front = []
    nadir = None
    for line in f:
        if 'nadir' not in line:
            xs, objs = line.split('(')
            pareto_front.append(xs.split() + [{'obj': [float(e.strip().strip(')')) for e in objs.split()]}])
        else:
            nadir = [float(e) for e in line.strip('nadir:').split()]

    objs = [e[-1]['obj'] for e in pareto_front ]
    print(hv.hypervolume(a(objs), nadir))

if __name__ == '__main__':
    if '-hv' in sys.argv:
        if sys.argv[1] != '-hv':
            print_hv(sys.argv[1])
        elif len(sys.argv) > 2 and sys.argv[2] != '-hv':
            print_hv(sys.argv[2])
        else:
            print_hv()

    elif len(sys.argv) == 2:
        show_pareto_front(sys.argv[1])
    else:
        show_pareto_front()
