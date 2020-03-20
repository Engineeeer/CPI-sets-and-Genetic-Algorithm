import numpy as np
import matplotlib.pyplot as plt
from scipy import signal as sg
from scipy.linalg import expm

num = np.array([1])
den = np.array([1,1])
A = np.array([[0,1],[0,-1]])
B = np.array([[0],[1]])
C = np.array([1,0])
D = np.array([0])

A_t = np.array([[0,1,0],[-3,-3,1],[-1,0,0]])
B_t = np.array([[0],[0],[1]])
C_t = np.array([1,0,0])

b = sg.cont2discrete((num,den),dt = 0.01)
c = sg.cont2discrete((A,B,C,D),dt = 1)
d = sg.cont2discrete((A_t,B_t,C_t,D),dt = 0.05)

x_series = np.zeros((500,3),dtype = np.float64)
X_series = np.zeros((100,3),dtype = np.float64)

x = np.array([5,5,0],dtype = np.float64)
X = np.array([5,5,0],dtype = np.float64)
X_series[0] = X

dt = 0.01

for i in range(500):
    dx = (A_t @ x) * dt
    x += dx
    x_series[i] = x

for i in range(100-1):
    X = (d[0] @ X)
    X_series[i+1] = X

plt.figure()

plt.plot(x_series[::5,0])
plt.plot(X_series[:,0])

plt.show()


