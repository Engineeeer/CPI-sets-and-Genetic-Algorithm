import numpy as np
import matplotlib.pyplot as plt

x_1 = np.linspace(-15,15,301)

#if|u|<= 10

x_2 = (10 -4 * x_1)/5
x_3 = (10 -(-4) * x_1)/(-5)

#if|w|<= 2

x_4 = (8 -4 * x_1)/5
x_5 = (8 -(-4) * x_1)/(-5)

#if|u| <= 20

x_6 = (20 -4* x_1)/5
x_7 = (20 -(-4) * x_1)/(-5)

plt.figure()
plt.plot(x_1,x_2)
plt.plot(x_1,x_3)
plt.plot(x_1,x_4,color = "k")
plt.plot(x_1,x_5,color = "k")
plt.plot(0,0,"o")

plt.figure()
plt.plot(x_1,x_2)
plt.plot(x_1,x_3)
plt.plot(x_1,x_6,color = "k")
plt.plot(x_1,x_7,color = "k")
plt.plot(0,0,"o")
plt.show()
