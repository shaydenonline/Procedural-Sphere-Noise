
import numpy as np
import matplotlib.pyplot as plt

# Skew factor for 2D simplex noise
F2 = 0.5 * (np.sqrt(3.0) - 1.0)

# Range of integer grid points (include negatives for clarity)
rng = np.arange(-2, 4)

# Generate square lattice points
grid_points = np.array([[i, j] for i in rng for j in rng])

# Skew transformation function
def skew_2d(x, y):
    s = (x + y) * F2
    return x + s, y + s

# Apply skew to all grid points
skewed_points = np.array([skew_2d(x, y) for x, y in grid_points])

# Pick a test point in Euclidean space (try changing it)
test_point = np.array([0.3, 0.7])
test_skewed = np.array(skew_2d(*test_point))

# Plot setup
fig, axes = plt.subplots(1, 2, figsize=(12, 6))
titles = ["Original (square grid)", "After skew (simplex lattice)"]

for ax, points, title in zip(axes, [grid_points, skewed_points], titles):
    ax.set_title(title)
    ax.set_aspect("equal", "box")

    # Plot grid points
    ax.scatter(points[:,0], points[:,1], color="black", s=10, zorder=2)

    # Draw grid lines (connect horizontal & vertical neighbors)
    for i in rng:
        for j in rng:
            # Horizontal neighbor
            if i + 1 in rng:
                p1 = np.array([i, j])
                p2 = np.array([i+1, j])
                if title.startswith("After"):
                    p1 = np.array(skew_2d(*p1))
                    p2 = np.array(skew_2d(*p2))
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], "gray", lw=0.8, zorder=1)

            # Vertical neighbor
            if j + 1 in rng:
                p1 = np.array([i, j])
                p2 = np.array([i, j+1])
                if title.startswith("After"):
                    p1 = np.array(skew_2d(*p1))
                    p2 = np.array(skew_2d(*p2))
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], "gray", lw=0.8, zorder=1)

    # Plot test point
    if title.startswith("After"):
        ax.scatter(*test_skewed, color="red", s=50, label="Skewed point", zorder=3)
    else:
        ax.scatter(*test_point, color="red", s=50, label="Original point", zorder=3)

    ax.legend()

plt.tight_layout()
plt.show()
