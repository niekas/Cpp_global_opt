import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import LightSource

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

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
Z = zs.reshape(X.shape)

ls = LightSource(270, 45)
rgb = ls.shade(Z, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=rgb, linewidth=0, antialiased=False, shade=False)

points = []
def draw_points_from_partition():
    pf = open('log/partition.txt')
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
            ax.plot([point[0]], [point[1]], [point[2]], 'or', markersize=2)

draw_points_from_partition()

plt.show()




import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

delta = 0.025
# x = np.arange(-3.0, 3.0, delta)
# y = np.arange(-2.0, 2.0, delta)
# X, Y = np.meshgrid(x, y)
# Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
# Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# # difference of Gaussians
# Z = 10.0 * (Z2 - Z1)


# Create a simple contour plot with labels using default colors.  The
# inline argument to clabel will control whether the labels are draw
# over the line segments of the contour, removing the lines beneath
# the label
plt.figure()
contour_levels =  np.arange(-1,4.1,0.2)
#  -1.00000000e+00 , 
#   -2.22044605e-16,   2.00000000e-01,   4.00000000e-01,
#    6.00000000e-01,   8.00000000e-01,   1.00000000e+00,   1.20000000e+00,
#    1.40000000e+00,   1.60000000e+00,   1.80000000e+00,   2.00000000e+00,
#    2.20000000e+00,   2.40000000e+00,   2.60000000e+00,   2.80000000e+00,
#    3.00000000e+00,   3.20000000e+00,   3.40000000e+00,   3.60000000e+00,
#    3.80000000e+00,   4.00000000e+00]

CS = plt.contour(X, Y, Z, levels=contour_levels)
# CS = plt.contour(X, Y, Z)
plt.xlabel('X1')
plt.ylabel('X2')
plt.clabel(CS, inline=1, fontsize=10)
# plt.title('Simplest default with labels')

plt.show()
