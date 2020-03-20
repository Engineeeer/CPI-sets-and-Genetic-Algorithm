import matplotlib.pyplot as plt
from scipy import signal as sg
from scipy.linalg import expm
import numpy as np
import sympy as sym

sym.init_printing()
x = sym.Symbol("x")
Ts = 0.1
K_p = 2.0
K_d = -1.5
K_i = 1.0
A_ch = np.array([[0,1,0],[-K_p - 1,K_d - 1,K_i],[-1,0,0]])
B_ch = np.array([[0],[K_p],[1]])
C = np.array([1,0,0])
D = np.array([0])

A = np.array([[0,1,0],[-1,-1,0],[-1,0,0]])
B = np.array([[0],[0],[1]])

B_b = np.array([[0],[1],[0]])
F = np.array([[-K_p, K_d, K_i]])
A_bf = B_b @ F
B_bk = np.array([[0],[K_p],[0]])

a = sg.cont2discrete((A_ch,B_ch,C,D),dt = Ts)
b = sg.cont2discrete((A,B,C,D),dt = Ts)
c = sg.cont2discrete((A_bf,B_bk,C,D),dt = Ts)

A_ch_d = expm(A_ch * Ts)
A_d = expm(A * Ts)
A_bf_d = expm(A_bf * Ts)
#print("A_d = expm(A*Ts) = ", expm(A_t*Ts))
#print(d[0] - expm(A_t * Ts))
A_s = sym.Matrix(A_ch * x).exp()
#A_s_int = sym.integrate(A_s,(x,0,Ts))
#B_s = A_s_int @ B_ch


