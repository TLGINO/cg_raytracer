import numpy as np

# Define vectors x and y (replace these with your actual vectors)
x = np.array([2**0.5, 1, 0])
y = np.array([1, 1, 1])

# Calculate the cross product of x and y
cross_product = np.cross(x, y)
print(cross_product)

# Normalize the cross product vector to obtain z
magnitude = np.linalg.norm(cross_product)
z = cross_product / magnitude

# Ensure that |z| = 1
z_normalized = z / np.linalg.norm(z)

print("Vector z perpendicular to x and y with magnitude 1:")
print(z_normalized)
