# coding: utf-8
# should construct simplexes and call show potential function
import sys
from numpy import array as a, matrix as m, arange, sqrt, isnan, pi, cos, sin, mean
from itertools import permutations
from mpl_toolkits.mplot3d import axes3d
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib
import random
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import LightSource
# matplotlib.rcParams.update({'font.size': 22})
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)
plt.rc('axes', labelsize=24)


def show_partition(filename='log/partition.txt'):
    from matplotlib import pyplot as plt
    draw_from_iteration = 1
    iteration = 0
    f = open(filename)
    simplexes = []
    selected_mode = False
    wanted_mode = False
    selected = []
    wanted = []
    title = ''
    # ok = False
    for line in f:
        if 'Iteration' in line or 'Partition' in line:
            # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,6))
            # ax1.axis([-0.05, 1.05, -0.05, 1.05])
            # ax2.axis([-0.05, 1.05, -0.05, 1.05])
            # title = line
            iteration = int(line.strip().strip('Iteration:').strip('Partition:').strip())
            continue
        if 'Selected' in line:
            selected_mode = True
            continue
            # Ignore till empty line found
        if 'Wanted' in line:
            selected_mode = False
            wanted_mode = True
            continue
        parts = line.split(';')
        simplex = []

        if line == '\n':
            selected_mode = False
            if simplexes and iteration >= draw_from_iteration:
                show_potential(simplexes, selected, wanted, title=title)
            simplexes = []
            selected = []
            continue
                    # [[0.75, 0.5], [0.875, 0.375], [1.0, 0.5]]
                    # plt.plot([s[j-1][0], s[j][0]], [s[j-1][1], s[j][1]], 'b-')
        else:
            if not selected_mode and not wanted_mode:
                add_to = simplexes
            elif not wanted_mode:
                add_to = selected
            else:
                add_to = wanted
            for part in parts:
                if ',' not in part:
                    simplex.append([float(e) for e in part.split() if not ('(') in e])
                    if '(' in part:
                        simplex[-1].append({'obj': (float(part.split()[-1].strip().strip('()')),)})
                if '(' in part and ',' in part:
                    size, value = part.strip().strip('()').split(',')
                    simplex.append({'size': float(size), 'value': float(value)})
            if simplex:
                add_to.append(simplex)

    # title = title + 'Dalinimui pasirinkta simpleksu: ' + (str(len(selected)) + " is " + (str(len(simplexes))))
    show_potential(simplexes, selected, wanted, title=title)


def l2norm(a1, a2):
    '''Euclidean norm, which converts arguments to arrays automatically.'''
    if isinstance(a1, (int, float)):  # len(X) < 2:
        return abs(a(a1)-a(a2))
    return sqrt(sum([e**2 for e in (a(a1)-a(a2))]))


def show_potential(simplexes, selected=[], wanted=[], show=True, title=''):
    from matplotlib import pyplot as plt
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import random
    from matplotlib import pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import LightSource


    separate_images = True

    if separate_images:
        fig1 = plt.figure(figsize=(10,10))
        ax1 = fig1.add_subplot(111)
        fig2 = plt.figure(figsize=(10,10))
        ax2 = fig2.add_subplot(111)
        fig3 = plt.figure(figsize=(10,10))
        ax3 = fig3.add_subplot(111)
        fig4 = plt.figure(figsize=(10,10))
        ax4 = fig4.add_subplot(111)
    else:
        fig = plt.figure(figsize=(14,6))
        ax3 = fig.add_subplot(221)
        ax1 = fig.add_subplot(222)
        ax2 = fig.add_subplot(223)
        # ax4 = fig.add_subplot(224, projection='3d')
        ax4 = fig.add_subplot(224)

    ax1.set_xlim([-0.01, 1.+ 0.01])
    ax1.set_ylim([-0.01, 1.+ 0.01])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    ax4.set_xlim([-0.01, 1.+ 0.01])
    ax4.set_ylim([-0.01, 1.+ 0.01])


    # print('got title: ', title)
    # fig.suptitle(title)
    # plt.title(title)

    # ax2 = fig.add_subplot(111)


    ## Draw one plot
    # fig = plt.figure(figsize=(10,10))
    # fig.suptitle(title)
    # ax2 = fig.add_subplot(111)

    ## Convex-hull
    for i in range(len(selected[:-1])):
        ax2.plot([selected[i][-1]['size'], selected[i+1][-1]['size']],
                 [selected[i][-1]['value'], selected[i+1][-1]['value']], 'r-', linewidth=2)

    for simplex in simplexes:
        ax2.plot([simplex[-1]['size']], [simplex[-1]['value']], 'bo', markersize=6)
    # f = open('tmp.output', 'w')
    for simplex in selected:
        ax2.plot([simplex[-1]['size']], [simplex[-1]['value']], 'ro', markersize=10, markeredgewidth=1)
    #     f.write('G1: %f, G2: %f, vertices: %s' % (simplex[-1]['value'], simplex[-1]['size'], str(simplex[:-1])+'\n'))
    # f.close()



    ## Wanted simplices
    # for simplex in wanted:
    #     ax2.plot([simplex[-1]['size']], [simplex[-1]['value']], 'yo')

    ## Stairs rule
    # for i in range(len(selected[:-1])):
    #     mid_size = (selected[i][-1]['size'] + selected[i+1][-1]['size']) / 2.
    #     # mid_value = (selected[i][-1]['value'] + selected[i+1][-1]['value']) / 2.
    #     ax2.plot([selected[i][-1]['size'], mid_size, mid_size, selected[i+1][-1]['size']],
    #             [selected[i][-1]['value'], selected[i][-1]['value'], selected[i+1][-1]['value'], selected[i+1][-1]['value']], 'r-')

    ax2.set_ylabel('G1')
    ax2.set_xlabel('G2')
    ax1.set_ylabel('X2')
    ax1.set_xlabel('X1')
    ax4.set_ylabel('X2')
    ax4.set_xlabel('X1')
    ax3.set_ylabel('X2')
    ax3.set_xlabel('X1')

    for simplex in simplexes:
        s = simplex[:-1]
        for j in range(len(s)):
            if len(s) == 3:
                ax1.plot([s[j-1][0], s[j][0]], [s[j-1][1], s[j][1]], 'b-', linewidth=2)
            else:
                # should use permutations here
                ax1.plot([s[j-1][0], s[j][0]], [s[j-1][1], s[j][1]], [s[j-1][2], s[j][2]], 'b-', linewidth=2)

    for simplex in selected:
        s = sort_vertexes_longest_edge_first(simplex)[:-1]
        for i, j in permutations(range(len(s)), 2):
            if len(s) == 3:
                ax1.plot([s[i][0], s[j][0]], [s[i][1], s[j][1]], 'r-', linewidth=4)
            else:
                # should use permutations here
                ax1.plot([s[i][0], s[j][0]], [s[i][1], s[j][1]], [s[i][2], s[j][2]], 'r-', linewidth=4)

        edge_lengths = []   # [(vertex_index, vertex_index, edge_length),]
        for i, j in permutations(range(len(s)), 2):
            if j > i:
                edge_lengths.append((i, j, l2norm(s[i][:-1], s[j][:-1])))
        le_i, le_j, le_length = max(edge_lengths, key=lambda x: x[-1])

        if len(s) == 3:
            division_point = [(s[le_i][0]+s[le_j][0])/2., (s[le_i][1]+s[le_j][1])/2.]
            ax1.plot([division_point[0], s[2][0]], [division_point[1], s[2][1]], 'r--', linewidth=4)
            ax1.plot([division_point[0]], [division_point[1]], 'wo', markersize=10, markeredgewidth=2)
        else:
            division_point = [(s[le_i][0] + s[le_j][0])/2., (s[le_i][1] + s[le_j][1])/2., (s[le_i][2] + s[le_j][2])/2.]
            ax1.plot([division_point[0]], [division_point[1]], [division_point[2]], 'ro')
            for i in range(len(s)):
                if i != le_j and i != le_i:
                    ax1.plot([division_point[0], s[i][0]], [division_point[1], s[i][1]], [division_point[2], s[i][2]], 'r--')
                # ax1.plot([division_point[0], s[3][0]], [division_point[1], s[3][1]], [division_point[2], s[3][2]], 'r--')

    for simplex in simplexes:
        for j in range(len(s)):
            if len(s) == 3:
                ax1.plot([simplex[j][0]], [simplex[j][1]], 'bo', markersize=7)
            else:
                ax1.plot([simplex[j][0]], [simplex[j][1]], [simplex[j][2]], 'bo', markersize=7)

    # ax2.axis([min([simplexes]) -0.05, 1.05, -0.05, 1.05])
    # max_size = max([s[-1]['size'] for s in simplexes])
    # ax2.set_xlim([-0.05, max_size + 0.05])
    # ax1.axis([-0.05, 1.05, -0.05, 1.05])


    f = open('log/surface.txt')

    zs = []
    x = []
    y = []
    points = []
    for line in f:
        parts = line.split()
        parts = [float(p.strip().strip(':').strip(',')) for p in parts]
        points.append(parts)
        x.append(parts[0])
        y.append(parts[1])
        zs.append(parts[2])
    zs = np.array(zs)
    x = np.array(x)
    y = np.array(y)

    X = np.reshape(x, (-1, int(np.sqrt(x.shape[0]))))
    Y = np.reshape(y, (-1, int(np.sqrt(y.shape[0]))))

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    Z = zs.reshape(X.shape)

    ls = LightSource(270, 45)
    rgb = ls.shade(Z, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
    # ax4.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=rgb, linewidth=0, antialiased=False, shade=False)

    points = []
    def draw_points_from_partition():
        pf = open('log/final_partition.txt')
        for line in pf:
            if 'art' not in line and 'ele' not in line:
                verts = line.split(';')
                for vert in verts:
                    point = []
                    for part in vert.split():
                        if ',' not in part:
                            if '(' in part:
                                part = part.strip('()')
                            point.append(float(part))
                    points.append(point)
        for point in points:
            if len(point) == 3:
                # ax4.plot([point[0]], [point[1]], [point[2]], 'or', markersize=2)
                ax4.plot([point[0]], [point[1]], 'or', markersize=6)

    draw_points_from_partition()



    # Contour plot
    import matplotlib
    import numpy as np
    import matplotlib.cm as cm
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt

    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'

    delta = 0.025
    # contour_levels = [-1, -0.5, -0.2, 0.0, 0.2,
    #         0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2., 2.5, 3., 3.5, 4.]

    # contour_levels = np.arange(-1,4.1,0.3)
    contour_levels = [-1, -0.5,  0.0,  0.2,  0.5,  0.8,  1.1,  1.4,  1.7,  2. , 2.3,  2.6,  2.9,  3.2,  3.5,  3.8]
    #  -1.00000000e+00 , 
    #   -2.22044605e-16,   2.00000000e-01,   4.00000000e-01,
    #    6.00000000e-01,   8.00000000e-01,   1.00000000e+00,   1.20000000e+00,
    #    1.40000000e+00,   1.60000000e+00,   1.80000000e+00,   2.00000000e+00,
    #    2.20000000e+00,   2.40000000e+00,   2.60000000e+00,   2.80000000e+00,
    #    3.00000000e+00,   3.20000000e+00,   3.40000000e+00,   3.60000000e+00,
    #    3.80000000e+00,   4.00000000e+00]

    CS = ax3.contour(X, Y, Z, levels=contour_levels, linewidths=2)
    # CS = plt.contour(X, Y, Z)
    # ax3.xlabel('X1')
    # ax3.ylabel('X2')
    ax3.clabel(CS, inline=1, fontsize=10)
    ax3.plot([0.54198], [0.951363], '*r', markersize=12)


    # ax3.plot([0.1], [0.9], 'ro')

    if show:
        # fig1.savefig('fig1.eps', format='eps', dpi=1200)
        # fig2.savefig('fig2.eps', format='eps', dpi=1200)
        # fig3.savefig('fig3.eps', format='eps', dpi=1200)
        # fig4.savefig('fig4.eps', format='eps', dpi=1200)
        fig1.savefig('fig1.eps', format='eps', dpi=600)
        fig2.savefig('fig2.eps', format='eps', dpi=600)
        fig3.savefig('fig3.eps', format='eps', dpi=600)
        fig4.savefig('fig4.eps', format='eps', dpi=600)
        plt.show()


def sort_vertexes_longest_edge_first(simplex):
    '''nD->nD Moves longest edge vertexes to the simplex vertex list beginning.'''
    # Find simplex edges lengths
    edge_lengths = []   # [(vertex_index, vertex_index, edge_length),]
    for i, j in permutations(range(len(simplex[:-1])), 2):
        if j > i:
            edge_lengths.append((i, j, l2norm(simplex[i][:-1], simplex[j][:-1])))


    # Get longest edge vertexes ids
    le_i, le_j, le_length = max(edge_lengths, key=lambda x: x[-1])

    # Move longest edge vertexes to simplex vertex list beginning
    vi = simplex[le_i]
    vj = simplex[le_j]
    simplex.remove(vi)
    simplex.remove(vj)
    simplex.insert(0, vj)
    simplex.insert(0, vi)
    return simplex


def show_diff_indexes(f1, f2):
    f1_content = open(f1).read().split('\n')
    indexes = []
    not_matching_lines = []
    for i, l2 in enumerate(open(f2)):
        found_match = False
        for j, l1 in enumerate(f1_content):
            line_match = True
            for i2 in l2.split(';'):
                l1_elems = [e.strip() for e in l1.split(';')]
                if i2.strip() not in l1_elems:
                    line_match = False
            if line_match:
                found_match = True
                break
        if not found_match:
            indexes.append(i)
            not_matching_lines.append(l2)
    return indexes, not_matching_lines


if __name__ == '__main__':
    if len(sys.argv) == 3:
        indexes, lines = show_diff_indexes(sys.argv[1], sys.argv[2])
        print(indexes)
        print('No match for these lines in ' + sys.argv[1])
        for line in lines:
            print(line.strip())
    elif len(sys.argv) == 2:
        show_partition(sys.argv[1])
    else:
        show_partition()
