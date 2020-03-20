import numpy as np
import numpy.linalg as nl
import pulp

Mz = np.array([[-1],[1]])
A = np.array([[-0.3,1.0],[-1.0,1.5]])
B = np.array([2.4,1.7])
C = np.array([[-0.2,-0.5]])
D = np.array([0.7])
sys = np.array([A,B,C,D])

def make_right(Mz, sys,idx, w_max, w_min):
    problem = pulp.LpProblem("g_0", pulp.LpMaximize)
    w = pulp.LpVariable("w", w_min, w_max, "Continuous")
    
    g_stur = np.zeros(len(Mz))
    for j in range(len(g_stur)):
        problem += (Mz @ sys[2] @ nl.matrix_power(sys[0],idx) @ sys[1])[j] * w
        problem.solve()
        #print(problem)
        #print(problem.solve())
        g_stur[j] += (Mz @ sys[2] @ nl.matrix_power(sys[0],idx) @ sys[1])[j] * w.value()
        #print(w.value())
    #print(g_stur)
    return g_stur

def make_right_2(Mz, sys, idx, w_max, w_min):
    problem_0 = pulp.LpProblem("g_0", pulp.LpMaximize)
    problem = pulp.LpProblem("g_*", pulp.LpMaximize)
    w = pulp.LpVariable("w", w_min, w_max, "Continuous")
    
    g_stur = np.zeros(len(Mz))

    for a in range(len(g_stur)):
            problem_0 += (Mz @ sys[3])[a] * w
            problem_0.solve()
            g_stur[a] += (Mz @ sys[3])[a] * w.value()

    for i in range(idx):
        for j in range(len(g_stur)):
            problem += (Mz @ sys[2] @ nl.matrix_power(sys[0],i) @ sys[1])[j] * w
            problem.solve()
            g_stur[j] += (Mz @ sys[2] @ nl.matrix_power(sys[0],i) @ sys[1])[j] * w.value()
    return g_stur

if __name__ == "__main__":
    a = [make_right(Mz, sys,i,1,-1) for i in range(5)]
    b = [make_right_2(Mz,sys,i,1,-1) for i in range(5)]
