import numpy as np
import matplotlib.pyplot as plt


A = np.array([[0.1,2.0],[-0.8,2.0]])
A_c = A - np.array([[2],[1]]) * np.array([[0.2,0.5]])
M_z = np.array([[-1],[1]])
C = np.array([[-0.2,-0.5]])

K_0 = M_z @ C
K_1 = K_0 @ A_c
K_2 = K_1 @ A_c
K_3 = K_2 @ A_c

x_1 = np.linspace(-15,15,301)

x_2 = (5 - K_0[0][0] * x_1)/K_0[0][1]
x_3 = (5 - K_0[1][0] * x_1)/K_0[1][1]
x_4 = (5 - K_1[0][0] * x_1)/K_1[0][1]
x_5 = (5 - K_1[1][0] * x_1)/K_1[1][1]
x_6 = (5 - K_2[0][0] * x_1)/K_2[0][1]
x_7 = (5 - K_2[1][0] * x_1)/K_2[1][1]
x_8 = (5 - K_3[0][0] * x_1)/K_2[0][1]
x_9 = (5 - K_3[1][0] * x_1)/K_2[1][1]
plt.figure()
plt.plot(x_1,x_2,color = "k")
plt.plot(x_1,x_3,color = "k")
plt.plot(x_1,x_4,color = "k")
plt.plot(x_1,x_5,color = "k")
plt.plot(x_1,x_6,color = "k")
plt.plot(x_1,x_7,color = "k")
plt.plot(x_1,x_8,color = "k")
plt.plot(x_1,x_9,color = "k")

plt.show()
