# coding: utf-8
# should construct simplexes and call show potential function
import sys
from numpy import array as a, matrix as m, arange, sqrt, isnan, pi, cos, sin, mean
from itertools import permutations
from mpl_toolkits.mplot3d import axes3d


def show_partition(filename='partition.txt'):
    from matplotlib import pyplot as plt
    draw_from_iteration = 1
    iteration = 0
    f = open(filename)
    simplexes = []
    selected_mode = False
    global_min_mode = False
    global_min = None
    D = None
    selected = []
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
        if 'Global' in line:
            global_min_mode = True
            selected_mode = False
            continue
            # Ignore till empty line found
        parts = line.split(';')
        simplex = []

        if line == '\n':
            selected_mode = False
            if simplexes and iteration >= draw_from_iteration:
                show_potential(simplexes, selected, title=title)
            simplexes = []
            selected = []
            continue
                    # [[0.75, 0.5], [0.875, 0.375], [1.0, 0.5]]
                    # plt.plot([s[j-1][0], s[j][0]], [s[j-1][1], s[j][1]], 'b-')
        else:
            if not selected_mode:
                add_to = simplexes
            else:
                add_to = selected
            if global_min_mode:
                global_min = [float(e) for e in parts[0].split() if not ('(') in e]
                D = len(simplexes[0]) -2
                continue
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

    title = title + 'Dalinimui pasirinkta simpleksu: ' + (str(len(selected)) + " is " + (str(len(simplexes))))
    show_potential(simplexes, selected, title=title, global_min=global_min, D=D)


def l2norm(a1, a2):
    '''Euclidean norm, which converts arguments to arrays automatically.'''
    if isinstance(a1, (int, float)):  # len(X) < 2:
        return abs(a(a1)-a(a2))
    return sqrt(sum([e**2 for e in (a(a1)-a(a2))]))


def show_potential(simplexes, selected=[], show=True, title='', global_min=None, D=None):
    from matplotlib import pyplot as plt
    ## Draw two plots
    fig = plt.figure(figsize=(14,6))
    if len(simplexes[0]) > 4:
        ax1 = fig.add_subplot(121, projection='3d')
    else:
        ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    # print('got title: ', title)
    fig.suptitle(title)
    plt.title(title)

    # ax2 = fig.add_subplot(111)


    ## Draw one plot
    # fig = plt.figure(figsize=(10,10))
    # fig.suptitle(title)
    # ax2 = fig.add_subplot(111)


    for simplex in simplexes:
        ax2.plot([simplex[-1]['size']], [simplex[-1]['value']], 'bo')
    for simplex in selected:
        ax2.plot([simplex[-1]['size']], [simplex[-1]['value']], 'ro')


    ## Convex-hull
    # for i in range(len(selected[:-1])):
    #     ax2.plot([selected[i][-1]['size'], selected[i+1][-1]['size']],
    #              [selected[i][-1]['value'], selected[i+1][-1]['value']], 'r-')

    ## Stairs rule
    for i in range(len(selected[:-1])):
        mid_size = (selected[i][-1]['size'] + selected[i+1][-1]['size']) / 2.
        # mid_value = (selected[i][-1]['value'] + selected[i+1][-1]['value']) / 2.
        ax2.plot([selected[i][-1]['size'], mid_size, mid_size, selected[i+1][-1]['size']],
                [selected[i][-1]['value'], selected[i][-1]['value'], selected[i+1][-1]['value'], selected[i+1][-1]['value']], 'r-')

    ax2.set_ylabel(u'Mažiausios funkcijos reikšm$\.{e}$s simplekse $\k{i}$vertis')
    ax2.set_xlabel('Simplekso diametras')
    ax1.set_ylabel('X2')
    ax1.set_xlabel('X1')

    for simplex in simplexes:
        s = simplex[:-1]
        for j in range(len(s)):
            if len(s) == 3:
                ax1.plot([s[j-1][0], s[j][0]], [s[j-1][1], s[j][1]], 'b-')
            else:
                # should use permutations here
                ax1.plot([s[j-1][0], s[j][0]], [s[j-1][1], s[j][1]], [s[j-1][2], s[j][2]], 'b-')

    for simplex in selected:
        s = sort_vertexes_longest_edge_first(simplex)[:-1]
        for i, j in permutations(range(len(s)), 2):
            if len(s) == 3:
                ax1.plot([s[i][0], s[j][0]], [s[i][1], s[j][1]], 'r-', linewidth=2)
            else:
                # should use permutations here
                ax1.plot([s[i][0], s[j][0]], [s[i][1], s[j][1]], [s[i][2], s[j][2]], 'r-', linewidth=2)

        edge_lengths = []   # [(vertex_index, vertex_index, edge_length),]
        for i, j in permutations(range(len(s)), 2):
            if j > i:
                edge_lengths.append((i, j, l2norm(s[i][:-1], s[j][:-1])))
        le_i, le_j, le_length = max(edge_lengths, key=lambda x: x[-1])

        if len(s) == 3:
            division_point = [(s[le_i][0]+s[le_j][0])/2., (s[le_i][1]+s[le_j][1])/2.]
            ax1.plot([division_point[0]], [division_point[1]], 'ro')
            ax1.plot([division_point[0], s[2][0]], [division_point[1], s[2][1]], 'r--')
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
                ax1.plot([simplex[j][0]], [simplex[j][1]], 'bo')
            else:
                ax1.plot([simplex[j][0]], [simplex[j][1]], [simplex[j][2]], 'bo')

    # ax2.axis([min([simplexes]) -0.05, 1.05, -0.05, 1.05])
    max_size = max([s[-1]['size'] for s in simplexes])
    ax2.set_xlim([-0.05, max_size + 0.05])
    ax1.axis([-0.05, 1.05, -0.05, 1.05])

    # Draw global minimum region
    deltas = {2: 1e-4, 3: 1e-6, 4: 1e-6, 5: 1e-7}
    delta = deltas[D]**(1./D)
    x = [global_min[0] - delta, global_min[0] + delta]
    y = [global_min[1] - delta, global_min[1] + delta]
    ax1.plot([x[0], x[1], x[1], x[0], x[0]], [y[0], y[0], y[1], y[1], y[0]], 'g', linewidth=2)
    # print([x[0], x[1], x[1], x[0], x[0]]+ [y[0], y[0], y[1], y[1], y[0]])

    if show:
        plt.show()
    # return ax1, ax2





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
