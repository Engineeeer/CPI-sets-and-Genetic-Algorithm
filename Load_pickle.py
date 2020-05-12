import numpy as np
import numpy.random
import random
import matplotlib.pyplot as plt
import pickle
from pathlib import Path
from scipy import signal as sg
#plt.rcParams["font.size"] = 36

#CPI用パラメータ部開始
Mz = np.array([[-1/5],[1/5]])
K_p = 1.0 
K_d = 0.5
A = np.array([[0,1],[-K_p,-0.5 - K_d]])
B = np.array([[0],[1.0]])
C = np.array([[-K_p,-K_d]])
D = np.array([K_p])

Ts = 0.05
c2d = sg.cont2discrete((A,B,C,D),dt = Ts)

A_d,B_d,C_d,D_d = c2d[0],c2d[1].reshape(2),c2d[2],c2d[3]

sys = np.array([A_d,B_d,C_d,D_d])
#CPI用パラメータ部終了

def simulate(r_series):
    Num = len(r_series)
    #GA用パラメータ部
    K_p = 1.0
    K_d = 0.5
    A = np.array([[0,1],[-K_p,-0.5 - K_d]])
    B = np.array([[0],[1.0]])
    C = np.array([[-K_p,-K_d]])
    D = np.array([K_p])
    
    Ts = 0.05
    c2d = sg.cont2discrete((A,B,C,D),dt = Ts)
    A_d,B_d,C_d,D_d = c2d[0],c2d[1].reshape(2),c2d[2],c2d[3]
    #GA用パラメータ部終了

    #保存用ログ確保
    x_series = np.zeros((Num, 2), dtype = np.float64)
    u_series = np.zeros((Num, 1), dtype = np.float64)

    for i in range(Num-1):
        x = A_d @ x_series[i] + B_d * r_series[i]
        u = C_d @ x_series[i] + D_d * r_series[i] 
        
        x_series[i+1] = x
        u_series[i+1] = u

    return x_series, u_series

def mkreference_series(gene):
    r_series = np.array(np.packbits(gene), dtype = np.float64) * 4.0 / 255.0 - 2.0
    return r_series

def long_mkreference_series(gene):
    n = 5
    r_series = mkreference_series(gene)
    lis = []
    for r in r_series:
        for _ in range(n):
            lis.append(r)
    return np.hstack([np.array(lis),np.ones( n * len(r_series))*1.5])

def evaluate_gene(gene):
    r_series = long_mkreference_series(gene)
    x_series, u_series = simulate(r_series)
    msae = np.abs(1.5 - x_series[:,0]).mean() + np.abs(1.5 - r_series).mean()
    return 1 / msae

def plot_result(gene):
    r_series = long_mkreference_series(gene)
    x_series, u_series = simulate(r_series)
    R_series = np.ones(len(r_series)) * 1.5
    X_series, U_series = simulate(R_series)
    f, a = plt.subplots()
    a.plot(np.arange(len(x_series[:,0]))/100,x_series[:,0],linewidth = 4)
    a.plot(np.arange(len(X_series[:,0]))/100,X_series[:,0],linewidth = 4)
    a.plot(np.arange(len(r_series))/100,R_series,linestyle = "--",color = "k",linewidth = 3)

    a.set_xlabel("time[s]")
    a.set_ylabel("Displacement x")

    f_1, a_1 = plt.subplots(2,1)
    a_1[0].plot(np.arange(len(r_series))/100,r_series,linewidth = "4")
    a_1[0].plot(np.arange(len(r_series))/100,np.ones(500) * 1.5,linewidth = "4")
    #a_1[0].set_xlabel("times[s]")
    a_1[0].set_ylabel("Reference")
    a_1[0].set_xlim(0,3)

    #f_2, a_2 = plt.subplots()
    a_1[1].plot(np.arange(len(r_series))/100,u_series[:,0],linewidth = "4")
    a_1[1].plot(np.arange(len(r_series))/100,U_series[:,0],linewidth = "4")
    a_1[1].plot(np.arange(len(r_series))/100,np.ones(len(r_series))*5,linestyle = "--",color = "k",linewidth = 3)
    a_1[1].plot(np.arange(len(r_series))/100,np.ones(len(r_series))*(-5),linestyle = "--",color = "k",linewidth = 3)
    #a_1[1].set_ylim(-3.5,3.5)
    a_1[1].set_xlabel("times[s]")
    a_1[1].set_ylabel("input")
    a_1[1].set_xlim(0,3)

def plot_CPI(CPI_set,gene):
    #CPI_set部分
    plt.figure()
    x1 = np.linspace(-15,15,301)
    for i in CPI_set:
        x2 = (1.0 - i[0] * (x1-1.5))/i[1]
        plt.plot(x1,x2,color = "k")

    #gene部分
    r_series = long_mkreference_series(gene)
    x_series, u_series = simulate(r_series)
    plt.plot(x_series[:,0],x_series[:,1],linewidth="7")

    plt.xlim(-2.0,5.0)
    plt.ylim(-3.5,3.5)
    plt.xlabel("x")
    plt.ylabel("dx")
    plt.plot(1.5,0,"o",color = "orange")

def plot_onlyCPI(CPI_set):
    #CPI_set部分
    plt.figure()
    x1 = np.linspace(-15,15,301)
    for i in CPI_set:
        x2 = (1.0 - i[0] * (x1-1.5))/i[1]
        plt.plot(x1,x2,color = "k")
    plt.xlim(-3.5,6.5)
    plt.ylim(-5,5)
    plt.xlabel("x")
    plt.ylabel("dx")
    plt.plot(1.5,0,"o")
    
def EoM_simulate(x_0, r_series, sys):
    #離散状態方程式の更新式
    A = sys[0]
    B = sys[1]
    Num = len(r_series)
    x_series = np.zeros((Num,2), dtype = np.float64)
    x_series[0] = x_0
    for i in range(Num-1):
        x = (A @ x_series[i] + B * r_series[i])
        x_series[i+1] = x
    return x_series    

if __name__ == "__main__":
    with open(Path(__file__).absolute().parent / "pickle_yard" / "genes100.pickle", "rb") as fp:
        genes = pickle.load(fp)
    with open(Path(__file__).absolute().parent / "pickle_yard" / "CPI_set.pickle", "rb") as ffp:
        CPI_set = pickle.load(ffp)
    evaluations = [evaluate_gene(g) for g in genes]
    elite_gene,elite_evaluation = sorted(zip(genes, evaluations), key = lambda item: item[1])[-1]
    plot_CPI(CPI_set,elite_gene)
    #plot_onlyCPI(CPI_set)
    X_0 = np.array([0,0])
    R_series = np.ones(1000) * 1.5
    X_series = EoM_simulate(X_0,R_series,sys)
    plt.plot(X_series[:,0],X_series[:,1],linewidth = "6",color = "orange")

    plot_result(elite_gene)
    
    plt.show()
