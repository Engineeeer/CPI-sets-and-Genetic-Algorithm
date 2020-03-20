import numpy as np
import matplotlib.pyplot as plt
#plt.rcParams["font.size"] = 12

def simulate(x_0, r_series):
    x = np.array([0, 0], dtype = np.float64) # データ更新用
    A = np.array([[-0.3, 1.0],[-1.0, 1.5]])
    #A = np.array([[0,1.0],[-1.0,-0.5]])
    B = np.array([2.4,1.7])
    #B = np.array([0,1.0])
    #B = np.array([1, 1])
    Num = len(r_series)
    x_series = np.zeros((Num,2), dtype = np.float64)
    x_series[0] = x_0
    for i in range(Num-1):
        x = (A @ x_series[i] + B * r_series[i])
        x_series[i+1] = x

    return x_series


r_series = np.ones((100,1))
x_1 = np.array([0,0])
#x_1 = np.array([-4, 2])
x_2 = np.array([-3.5, -6.5])
x_3 = np.array([4,-2])
x_4 = np.array([3.5,6.5])

x_1_series = simulate(x_1, r_series)
x_2_series = simulate(x_2, r_series)
x_3_series = simulate(x_3, r_series)
x_4_series = simulate(x_4, r_series)

#CPI集合をプロット
A = np.array([[0.1, 2.0],[-0.8, 2.0]])
B = np.array([[1],[1]])
M_z = np.array([[-1],[1]])
C = np.array([[-0.2,-0.5]])
D = np.array([[0.7]])
A_c = A - np.array([[2],[1]]) * np.array([[0.2,0.5]])
B_c = B + np.array([[2],[1]]) @ D

#手計算で導出 |w|<= 1.0
#g_0 = 0
g_0 = 0.7
g_1 = 1.33
g_2 = 0.271
g_3 = 0.4063
g_4 = 0.63661

K_0 = M_z @ C
K_1 = K_0 @ A_c
K_2 = K_1 @ A_c
K_3 = K_2 @ A_c
K_4 = K_3 @ A_c
K_5 = K_4 @ A_c

x_series = np.linspace(-15,15,301)

x_0 = ((5 - g_0) - K_0[0][0] * x_series)/K_0[0][1]
x_1 = ((5 - g_0) - K_0[1][0] * x_series)/K_0[1][1]

x_2 = (5 - sum([g_0,g_1]) - K_1[0][0] * x_series)/K_1[0][1]
x_3 = (5 - sum([g_0,g_1]) - K_1[1][0] * x_series)/K_1[1][1]

x_4 = (5 - sum([g_0,g_1,g_2]) - K_2[0][0] * x_series)/K_2[0][1]
x_5 = (5 - sum([g_0,g_1,g_2]) - K_2[1][0] * x_series)/K_2[1][1]

x_6 = (5 - sum([g_0,g_1,g_2,g_3]) - K_3[0][0] * x_series)/K_3[0][1]
x_7 = (5 - sum([g_0,g_1,g_2,g_3]) - K_3[1][0] * x_series)/K_3[1][1]

x_8 = (5 - sum([g_0,g_1,g_2,g_3,g_4]) - K_4[0][0] * x_series)/K_4[0][1]
x_9 = (5 - sum([g_0,g_1,g_2,g_3,g_4]) - K_4[1][0] * x_series)/K_4[1][1]

plt.figure()
plt.plot(-4,2,"o",color = "grey")
plt.plot(-3.5,-6.5,"o",color = "k")
plt.plot(4,-2,"o",color = "grey")
plt.plot(3.5,6.5,"o",color = "k")
plt.plot(x_1_series[:,0],x_1_series[:,1],color = "grey",linewidth = "2")
plt.plot(x_2_series[:,0],x_2_series[:,1],color = "k",linewidth = "2")
plt.plot(x_3_series[:,0],x_3_series[:,1],color = "grey",linewidth = "2")
plt.plot(x_4_series[:,0],x_4_series[:,1],color = "k",linewidth = "2")
plt.plot(x_series,x_0,color = "k",linewidth = "4")
plt.plot(x_series,x_1,color = "k",linewidth = "4")
plt.plot(x_series,x_2,color = "k",linestyle = "--")
plt.plot(x_series,x_3,color = "k",linestyle = "--")
plt.plot(x_series,x_4,color = "k",linestyle = "--")
plt.plot(x_series,x_5,color = "k",linestyle = "--")
plt.plot(x_series,x_6,color = "k",linestyle = "--")
plt.plot(x_series,x_7,color = "k",linestyle = "--")
plt.plot(x_series,x_8,color = "k",linestyle = "--")
plt.plot(x_series,x_9,color = "k",linestyle = "--")
plt.xlabel("x_1")
plt.ylabel("x_2")
plt.xlim(-10,10)
plt.ylim(-10,10)
plt.show()
