import numpy as np
import matplotlib.pyplot as plt

# Define the function
def f(x, y):
    return (x + 2*y -7) * (x + 2*y -7) + (2*x + y -5) * (2*x + y -5)

# Generate x and y values
x = np.linspace(-10, 10, 500)
y = np.linspace(-10, 10, 500)
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
plt.xlabel('x1')
plt.ylabel('x2')
plt.title('Wykres 6\n$(x1 + 2*x2 -7)^2 + (2*x1 + x2 -5)^2$\nCzerwony - NEWTON 0,12, Niebieski - NEWTON, Zielony - NEWTON 0,05')

with open('NEWTONh5.000000.txt', 'r') as file:
    for line in file:
        x, y = map(float, line.strip().split(';'))
        plt.scatter(x, y, color='green', marker='.', s=15, alpha=0.75)

with open('NEWTONh12.000000.txt', 'r') as file:
    for line in file:
        x, y = map(float, line.strip().split(';'))
        plt.scatter(x, y, color='red', marker='.', s=15, alpha=0.75)
        
with open('NEWTONh-100.000000.txt', 'r') as file:
    for line in file:
        x, y = map(float, line.strip().split(';'))
        plt.scatter(x, y, color='blue', marker='.', s=15, alpha=0.75)

# Show the plot
plt.show()