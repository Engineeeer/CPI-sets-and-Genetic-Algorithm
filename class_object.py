import numpy as np
import numpy.linalg as nl
import pulp
import matplotlib.pyplot as plt

"""
Objectクラスで逐次のK_iを管理
その後CPIクラスを作って[K_0,K_1,...,K_n]を管理し,ELIMを試してゆく
"""
Mz = np.array([[-1/5],[1/5]])
A = np.array([[-0.3,1.0],[-1.0,1.5]])
B = np.array([2.4,1.7])
C = np.array([[-0.2,-0.5]])
D = np.array([0.7])

sys = np.array([A,B,C,D])

def make_left_matrix(Mz, sys, idx):
    """
    sys = np.array([A,B,C,D])
    Mz@C@A^(idx-1)
    """
    return Mz @ sys[2] @ nl.matrix_power(sys[0],idx)

def make_right_raw(Mz, sys, idx, w_max, w_min):
    """
    Mz: Sz x p の行列
    lenは行の長さを返すので,len(Mz) = Szである
    
    pulp:
    変数定義には(名前,定義域(最小,最大),型)が必要
    floatはContinuousで表現される
    """

    problem_0 = pulp.LpProblem("g_0", pulp.LpMaximize)
    problem = pulp.LpProblem("g_*", pulp.LpMaximize)
    w = pulp.LpVariable("w", w_min, w_max, "Continuous")
    
    #メモリ確保
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
    return 1-g_stur

class Object:
    def __init__(self, left_matrix, right_raw):
        """
            left_matrix: A = np.array([[a11, a12, ...],[a21,...]...])
            right_raw: B = np.array([[b1],[b2],...[bn]])

        """
        self.left = left_matrix
        self.right = right_raw

    def check(self):
    #1-g* > 0 のcheck
        cnt = 0
        for i in range(len(self.right)):
            if self.right[i].all() > 0:
                cnt += 1
            #print(cnt)
        if cnt == len(self.right):
            return True
        else:
            return False

    def plot_to_figure(self):
        x_1 = np.linspace(-15,15,301)
        for i in range(len(self.left)):
            x_2 = (self.right[i] - self.left[i][0]*x_1)/self.left[i][1]
            plt.plot(x_1,x_2,color = "k")
        plt.xlim(-8.0,8.0)
        plt.ylim(-10,10)
        #plt.show()
            

if __name__ == "__main__":
    plt.figure()
    c = [Object(make_left_matrix(Mz,sys,i), make_right_raw(Mz,sys,i,1,-1)) for i in range(5)]
    for i in range(len(c)):
        c[i].plot_to_figure()
    plt.show()
    #a = [make_left_matrix(Mz, sys, i) for i in range(5)]
    #b = [make_right_raw(Mz,sys,i,1,-1) for i in range(5)]
