import numpy as np
import matplotlib.pyplot as plt
from scipy import signal as sg
import sympy as sym
plt.rcParams["font.size"] = 24

def simulate(r_series,K_p,K_d,K_i):
    dt = 0.01
    A = np.array([[0,1,0],[-K_p - 1,K_d - 1,K_i],[-1,0,0]])
    B = np.array([0,K_p,1])

    x = np.array([0, 0, 0], dtype = np.float64)

    N = len(r_series)

    x_series = np.zeros((N, 3), dtype = np.float64)
    u_series = np.zeros((N, 1), dtype = np.float64)
    
    for i in range(N):
        dx = (A @ x + B * r_series[i]) * dt
        x += dx
        u = -K_p * x[0] + K_d * x[1] + K_i * x[2] + K_p * r_series[i]
        x_series[i] = x
        u_series[i] = u
    return x_series,u_series

def dis_simulate(r_series_d,K_p,K_d,K_i):
    Ts = 0.10
    A = np.array([[0,1,0],[-K_p - 1,K_d - 1,K_i],[-1,0,0]])
    B_w = np.array([[0],[K_p],[1]])
    C = np.array([1,0,0])
    D = np.array([0])

    c2d = sg.cont2discrete((A,B_w,C,D),dt = Ts) #C,Dは使わないフェイク

    A_d = c2d[0]
    B_w_d = c2d[1].reshape(3)

    x = np.array([0, 0, 0], dtype = np.float64)

    N = len(r_series_d)

    x_series_d = np.zeros((N, 3), dtype = np.float64)
    u_series_d = np.zeros((N, 1), dtype = np.float64)

    u_series_d[0] = -K_p * x[0] + K_d * x[1] + K_i * x[2] + K_p * r_series[0]
    
    for i in range(N-1):
        x = (A_d @ x_series_d[i] + B_w_d * r_series[i])
        x_series_d[i+1] = x 
        u = -K_p * x[0] + K_d * x[1] + K_i * x[2] + K_p * r_series[i+1]
        u_series_d[i+1] = u
    return x_series_d,u_series_d

def filter_simulate(r_series_f,K_p,K_d,K_i,u_max,u_min):
    dt = 0.01
    A = np.array([[0,1,0],[-1,-1,0],[-1,0,0]])
    B = np.array([0,1,0])
    C = np.array([0,0,1])
    F = np.array([-K_p, K_d, K_i])

    x = np.array([0, 0, 0], dtype = np.float64)

    N = len(r_series_f)

    x_series_f = np.zeros((N, 3), dtype = np.float64)
    u_series_f = np.zeros((N, 1), dtype = np.float64)
    
    for i in range(N):
        u = F @ x + K_p * r_series[i]
        #u = u_filter(u,u_max,u_min)
        dx = (A @ x + B * u + C * r_series[i] ) * dt
        x += dx
        x_series_f[i] = x
        u_series_f[i] = u
    return x_series_f,u_series_f

def dis_filter_simulate(r_series_f,K_p,K_d,K_i,u_max,u_min):
    dt = 0.01
    A = np.array([[0,1,0],[-1,-1,0],[-1,0,0]])
    B = np.array([0,1,0])
    C = np.array([0,0,1])
    F = np.array([-K_p, K_d, K_i])

    x = np.array([0, 0, 0], dtype = np.float64)

    N = len(r_series_f)

    x_series_f = np.zeros((N, 3), dtype = np.float64)
    u_series_f = np.zeros((N, 1), dtype = np.float64)
    
    for i in range(N):
        u = F @ x + K_p * r_series[i]
        #u = u_filter(u,u_max,u_min)
        dx = (A @ x + B * u + C * r_series[i] ) * dt
        x += dx
        x_series_f[i] = x
        u_series_f[i] = u
    return x_series_f,u_series_f

def u_filter(u,u_max,u_min):
    if u <= u_min:
        return u_min

    elif u >= u_max:
        return u_max

    else:
        return u

if __name__ == "__main__":
    #u,U = test_u_filter(1.2,-1.2)
    Ts_hi = 10 #Ts_hi = Ts/dt
    K_p = 2.0#0.5
    K_d = -1.5#-1.5
    K_i = 1.0#1.0
    u_max = 1.2
    u_min = -1.2
    time_lenth = 500
    dis_time_lenth = time_lenth / Ts_hi
     
    r_series = np.ones(time_lenth)
    r_series_d = np.ones(int(dis_time_lenth))
    
    t_series = np.arange(0,time_lenth)
    t_series_d = np.arange(0,time_lenth,Ts_hi)

    x_series,u_series = simulate(r_series,K_p,K_d,K_i)
    x_seriesf,u_seriesf = filter_simulate(r_series,K_p,K_d,K_i,u_max,u_min)
    x_series_d,u_series_d = dis_simulate(r_series_d,K_p,K_d,K_i)
    
    x_series_d_long,u_series_d_long = dis_simulate(r_series,K_p,K_d,K_i)
    
    plt.figure()
    plt.plot(np.arange(len(r_series)),x_seriesf[:,0])
    #plt.plot(np.arange(len(r_series)),x_seriesf[:,1])
    #plt.plot(np.arange(len(r_series)),x_seriesf[:,2])
    plt.plot(np.arange(len(r_series)),u_seriesf)
    plt.plot(np.arange(len(r_series)),r_series,color = "k",linestyle = "--")
    #plt.plot(np.arange(len(r_series)),np.ones(len(r_series)) * u_max) #サチリプロット
    #plt.plot(np.arange(len(r_series)),np.ones(len(r_series)) * u_min)
    #plt.ylim(-0.3,2.6)

    plt.figure()
    plt.plot(t_series,x_series[:,0],linewidth = 4)
    #plt.plot(t_series_d,x_series[::Ts_hi,0]) #間違えたけどカクカクになる
    #plt.plot(np.arange(len(r_series)),x_series[:,1])
    #plt.plot(np.arange(len(r_series)),x_series[:,2])
    plt.plot(t_series,u_series,linewidth = 4)
    #plt.plot(t_series_d,u_series[::Ts_hi]) #間違えたけどカクカクになる
    plt.plot(np.arange(len(r_series)),r_series,color = "k",linestyle = "--",linewidth = 4)
    

    #plt.plot(t_series_d,x_series_d[:,0])
    #plt.plot(t_series_d,u_series_d)
    plt.scatter(t_series_d,x_series_d[:,0],linewidth = 6)
    plt.scatter(t_series_d,u_series_d,linewidth = 6)
    #plt.scatter(t_series_d,x_series_d[:,2],linewidth = 6)
    plt.xlabel("step")

    #plt.plot(t_series,x_series_d_long[:,0])
    #plt.plot(t_series,u_series_d_long)
    plt.show()

    #x_series_d[:,0] / x_series[::5,0] で調べる
