import numpy as np
import matplotlib.pyplot as plt

def simulate(r_series,K_p,K_d,K_i):
    dt = 0.01
    A = np.array([[0,1,0],[K_p - 0.5,K_d - 0.2,K_i],[-1,0,0]])
    #B = np.array([[0],[0],[1]])
    B = np.array([0,0,1])

    x = np.array([0, 0, 0], dtype = np.float64)

    N = len(r_series)

    x_series = np.zeros((N, 3), dtype = np.float64)
    u_series = np.zeros((N, 1), dtype = np.float64)
    
    for i in range(N):
        dx = (A @ x + B * r_series[i]) * dt
        x += dx
        u = K_p * x[0] + K_d * x[1] + K_i * x[2]
        x_series[i] = x
        u_series[i] = u
    return x_series,u_series

def filter_simulate(r_series,K_p,K_d,K_i):
    dt = 0.01
    u_max = 1.1
    u_min = -1.1
    A = np.array([[0,1,0],[-0.5,-0.2,0],[-1,0,0]])
    B = np.array([0,1,0])
    C = np.array([0,0,1])

    x = np.array([0, 0, 0], dtype = np.float64)

    N = len(r_series)

    x_series = np.zeros((N, 3), dtype = np.float64)
    u_series = np.zeros((N, 1), dtype = np.float64)
    
    for i in range(N):
        u = K_p * x[0] + K_d * x[1] + K_i * x[2]
        u = u_filter(u,u_max,u_min)
        dx = (A @ x + B * u + C * r_series[i]) * dt
        x += dx
        x_series[i] = x
        u_series[i] = u
    return x_series,u_series

def u_filter(u,u_max,u_min):
    if u <= u_min:
        return u_min

    elif u >= u_max:
        return u_max

    else:
        return u

def test_u_filter(u_max,u_min):
    u = np.asarray([2,-2,1,-1,1.5,-1.5])
    U = []
    for i in u:
        a = u_filter(i,u_max,u_min)
        U.append(a)
    return np.asarray(u),np.asarray(U)

if __name__ == "__main__":
    #u,U = test_u_filter(1.2,-1.2)
    K_p = -5
    K_d = -5
    K_i = 2
    
    r_series = np.ones(2000)

    x_series,u_series = simulate(r_series,K_p,K_d,K_i)
    x_seriesf,u_seriesf = filter_simulate(r_series,K_p,K_d,K_i)

    
    plt.figure()
    plt.plot(np.arange(len(r_series)),x_seriesf[:,0])
    plt.plot(np.arange(len(r_series)),x_seriesf[:,1])
    #plt.plot(np.arange(len(r_series)),x_seriesf[:,2])
    plt.plot(np.arange(len(r_series)),u_seriesf)
    plt.plot(np.arange(len(r_series)),r_series,color = "k",linestyle = "--")
    #plt.ylim(-0.3,2.6)

    plt.figure()
    plt.plot(np.arange(len(r_series)),x_series[:,0])
    plt.plot(np.arange(len(r_series)),x_series[:,1])
    #plt.plot(np.arange(len(r_series)),x_series[:,2])
    plt.plot(np.arange(len(r_series)),u_series)
    plt.plot(np.arange(len(r_series)),r_series,color = "k",linestyle = "--")
    #plt.ylim(-0.3,2.6)

    plt.show()
