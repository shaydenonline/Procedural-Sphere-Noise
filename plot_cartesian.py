import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Create a dummy data file for demonstration
data = np.loadtxt('polar_coords_1007.txt')
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

# Set up the 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the coordinates
ax.scatter(x, y, z, c='r', marker='o') # 'c' for color, 'marker' for point style

# Add labels
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')

# Display the plot
plt.title('3D Coordinates from File')
plt.show()
