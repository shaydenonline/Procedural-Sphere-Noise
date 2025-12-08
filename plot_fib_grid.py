import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
points = np.loadtxt("polar_coords_107.txt", delimiter=" ")
x, y, z = points[:,0], points[:,1], points[:,2]

fig = plt.figure()

ax = fig.add_subplot(111, projection="3d")

ax.scatter(x, y, z, s=15, cmap='viridis')

ax.set_box_aspect([1, 1, 1])
plt.show()
