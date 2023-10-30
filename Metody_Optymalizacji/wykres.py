import numpy as np
import matplotlib.pyplot as plt

# Define the function
def f(x, y):
    return x**2 + y**2 - np.cos(2.5 * np.pi * x) - np.cos(2.5 * np.pi * y) + 2

# Generate x and y values
x = np.linspace(-1, 1, 400)
y = np.linspace(-1, 1, 400)
X, Y = np.meshgrid(x, y)

# Compute the function values for the grid points
Z = f(X, Y)

# Create a contour plot
plt.figure(figsize=(8, 6))
contours = plt.contour(X, Y, Z, levels=20, colors='k')
plt.clabel(contours, inline=True, fontsize=8)

# Add a color map to visualize the function surface
plt.imshow(Z, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower',
           cmap='viridis', alpha=0.5)
plt.colorbar()

# Add labels and title
plt.xlabel('x')
plt.ylabel('y')
plt.title('$x^2 + y^2 - \cos(2.5\pi x) - \cos(2.5\pi y) + 2$\nCzerwony - HJ, Niebieski - Rosen')

with open('Metody_Optymalizacji\\hj.txt', 'r') as file:
    for line in file:
        x, y = map(float, line.strip().split(';'))
        plt.scatter(x, y, color='red', marker='o', s=15)
        
with open('Metody_Optymalizacji\\rosen.txt', 'r') as file:
    for line in file:
        x, y = map(float, line.strip().split(';'))
        plt.scatter(x, y, color='blue', marker='o', s=15)


# Show the plot
plt.show()