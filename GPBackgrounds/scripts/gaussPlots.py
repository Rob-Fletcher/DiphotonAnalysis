'''
======================
3D surface (color map)
======================

Demonstrates plotting a 3D surface colored with the coolwarm color map.
The surface is made opaque by using antialiased=False.

Also demonstrates using the LinearLocator and custom formatting for the
z axis tick labels.
'''

#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.mlab as mlab
import numpy as np


fig = plt.figure()
#ax = fig.gca(projection='3d')

# Make data.
X = np.arange(-3, 3, 0.25)
Y = np.arange(-3, 3, 0.25)
X, Y = np.meshgrid(X, Y)
Z = 5/(2*np.pi) * np.exp(0.5*(-(X+Y)**2 - (Y)**2))


# Plot the surface.
"""
ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

cont = ax.contour(X,Y,Z,[0,1,2],zdir='x',offset=-3)
# Customize the z axis.
ax.set_zlim(0, 1.0)
ax.grid(False)
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.savefig('Gauss_Surface.pdf')

# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)
"""
x = np.linspace(-3,-3,100)
plt.plot(x, mlab.normpdf(x, 0, 1), lw=2)
plt

plt.show()
