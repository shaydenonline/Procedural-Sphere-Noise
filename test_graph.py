import matplotlib.pyplot as plt
import numpy as np
import csv

from matplotlib import cm
from matplotlib.ticker import LinearLocator

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
X = []
Y = []
Z = []
X, Y, Z = np.loadtxt('perlin_2d_plane.txt', delimiter=',', unpack=True)

'''
with open('test_points.txt', 'r') as datafile:
    plotting = csv.reader(datafile, delimiter=',')
    
    for ROWS in plotting:
        X.append(int(ROWS[0]))

        Y.append(int(ROWS[1]))
        Z.append(int(ROWS[2]))
# Plot the surface.
    '''
surf = ax.plot_trisurf(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
