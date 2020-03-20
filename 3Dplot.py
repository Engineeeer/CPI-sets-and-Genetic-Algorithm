import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x = np.arange(-2.0,2.1,0.1)
y = np.arange(-2.0,2.1,0.1)
z = np.arange(-3.0,3.1,0.1)
z = 2 * x + 2 * y
x_zero = np.zeros(len(x))
x_one = np.ones(len(x))
y_zero = np.zeros(len(x))
y_one = np.ones(len(x))
z_zero = np.zeros(len(x))

X,Y = np.meshgrid(x,y)

#Z_1 = 2 * X + 2 * Y - 5
#Z_2 = 2 * X + 2 * Y + 5

Z_1 = (25 * X + 8 * Y -5)/1
Z_2 = (25 * X + 8 * Y +5)/1

fig = plt.figure()

ax = Axes3D(fig)
ax.set_xlabel("x")
ax.set_ylabel("dx")
ax.set_zlabel("e")
ax.plot_wireframe(X,Y,Z_1)
ax.plot_wireframe(X,Y,Z_2)
#ax.plot(x,y,z)
ax.plot(x_zero,y,color = "k",linestyle = "--")
ax.plot(x,y_zero,color = "k",linestyle = "--")
ax.plot(x_zero,y_zero,z,color = "k",linestyle = "--")
ax.plot(x_one,y_zero,z,color = "k")

plt.show()
