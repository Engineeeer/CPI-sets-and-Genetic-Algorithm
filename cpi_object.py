import numpy as np
import numpy.linalg as nl
import pulp
import matplotlib.pyplot as plt
import itertools
import scipy.signal as sg
import pickle
from pathlib import Path

#システムパラメータ部開始

#入力制限は左右(上下)対称を仮定
Mz = np.array([[-1/5.0],[1/5.0]])
K_p = 1.0 
K_d = 0.5
A = np.array([[0,1],[-K_p,-0.5 - K_d]])
B = np.array([[0],[1.0]])
C = np.array([[-K_p,-K_d]])
D = np.array([K_p])

Ts = 0.05
c2d = sg.cont2discrete([A,B,C,D],dt = Ts)
A_d,B_d,C_d,D_d = c2d[0],c2d[1].reshape(2),c2d[2],c2d[3]

sys = np.array([A_d,B_d,C_d,D_d])

#システムパラメータ部終了
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

def make_left_matrix(Mz, sys, idx):
    return Mz @ sys[2] @ nl.matrix_power(sys[0],idx)

def make_right_raw(Mz, sys, idx, w_max, w_min):
    problem_0 = pulp.LpProblem("g_0", pulp.LpMaximize)
    problem = pulp.LpProblem("g_*", pulp.LpMaximize)
    w = pulp.LpVariable("w", w_min, w_max, "Continuous")
    #メモリ確保
    g_stur = np.zeros(len(Mz))
    #D != 0 の時の処理
    if sys[3] != 0:
        for a in range(len(g_stur)):
                problem_0 += (Mz @ sys[3])[a] * w
                problem_0.solve()
                g_stur[a] += (Mz @ sys[3])[a] * w.value()

    for i in range(idx):
        for j in range(len(g_stur)):
            problem += (Mz @ sys[2] @ nl.matrix_power(sys[0],i) 
                    @ sys[1])[j] * w
            problem.solve()
            g_stur[j] += (Mz @ sys[2] @ nl.matrix_power(sys[0],i) @ sys[1])[j] * w.value()
    return 1-g_stur

class Object:
    #Ax = B の形であり,A = left_matrix,B = right_rawのイメージ
    def __init__(self, left_matrix, right_raw):
        self.left = left_matrix
        self.right = right_raw

    def check(self):
    #1-g* > 0 のcheck
        cnt = 0
        for i in range(len(self.right)):
            if self.right[i] > 0:
                cnt += 1
        if cnt != len(self.right):
            return True

    def plot_to_figure(self):
        x_1 = np.linspace(-15,15,301)
        for i in range(len(self.left)):
            x_2 = (self.right[i] - self.left[i][0]*x_1)/self.left[i][1]
            plt.plot(x_1,x_2,color = "k")
        plt.xlim(-10,10)
        plt.ylim(-10,10)

    def normalization(self):
        #len(left_matrix)はleft_matrixの行の長さを返す
        a = [self.left[i]/self.right[i] for i in range(len(self.left))]
        return np.vstack(a)

class CPI:
    def __init__(self, Object):
        self.obj = Object
        self.unnormal_list = []
        self.M_i_set = []
        self.CPI_set = []
        #self.CPI_setの中にM_0, M_1=[M_1/2 M_0]^T, ... , M_inf=[M_(inf-1)/2 M_0] が入る.
        
    def set_M_i_set(self, idx):
        a = [self.unnormal_list[i].normalization() for i in range(idx)]
        self.M_i_set.append(np.vstack(a))
        #self.M_i_setの要素(即ちM_i)を受け取って非冗長なM_iを返す

    def ELIM(self,M_r,th = 1E-4):
        i_list = []
        for i in range(len(M_r)):
            problem = pulp.LpProblem("h", pulp.LpMaximize)
            x1 = pulp.LpVariable("x1",-100,100,"Continuous")
            x2 = pulp.LpVariable("x2",-100,100,"Continuous")
            #目的関数
            problem += M_r[i][0] * x1 + M_r[i][1] * x2
            #制約条件x<=set(M_r)を示す
            for j in range(len(M_r)):
                    problem += M_r[j][0] * x1 + M_r[j][1] * x2 <= 1
            problem.solve()
            #print(problem)
            print("x1",x1.value())
            print("x2",x2.value())
            print()
            h = M_r[i][0] * x1.value() + M_r[i][1] * x2.value()
            print(h)
            if th < 1 - h:
               i_list.append(i)
            print(i_list)
        self.CPI_set.append(np.delete(M_r,i_list,0))
        
    #staticmethodはクラス内でのみselfを参照できない関数
    #この様に定義することでグローバルで使わない事が自明となる
    #while,if文等で条件分岐を行う際,staticmethodを作ると良い
    @staticmethod
    def check_CPI_identity(cpi1, cpi2):
        return cpi1.shape[0] == cpi2.shape[0]

    def simulate(self):
        #目標値幅は左右(上下)対象を仮定する
        w_min = -2.0
        w_max = 2.0
        idx = 0
        self.unnormal_list.append(self.obj(make_left_matrix(Mz,sys,idx),make_right_raw(Mz,sys,idx,w_max,w_min)))
        self.set_M_i_set(1)
        self.ELIM(self.M_i_set[0])
        for i in itertools.count(1):
            idx += 1
            self.unnormal_list.append(self.obj(make_left_matrix(Mz,sys,i),make_right_raw(Mz,sys,i,w_max,w_min)))
            if self.unnormal_list[i].check == True:
                print("check is False")
                break
            self.set_M_i_set(i+1)
            self.ELIM(self.M_i_set[i])
            if CPI.check_CPI_identity(self.CPI_set[i], self.CPI_set[i-1]):
                break
        return print("Please search self.CPI_set[-1] & idx == ",idx)
            
    def plot_CPI_set(self):
        plt.figure()
        x1 = np.linspace(-15,15,301)
        for i in self.CPI_set[-1]:
            x2 = (1.0 - i[0] * (x1-1.5))/i[1]
            plt.plot(x1,x2,color = "k")
        plt.xlim(-10,10)
        plt.ylim(-10,10)   
        
        plt.plot(1.5,0,"o")
        #plt.show()

    def Intersection(self,M_i):
        Intersection_list = []
        b = np.ones(2)#1ベクトルが2行
        #A_listは全ての組み合わせを含む
        A_list = list(itertools.combinations(M_i,2))
        for i in range(len(A_list)):
            a = np.vstack([A_list[i][0],A_list[i][1]])
            #print("a",a) #printデバック用
            try:
                Ans = nl.solve(a,b)
                #print("Ans",Ans) #printデバック用
                Intersection_list.append(Ans)
            except:
                pass
        return Intersection_list

if __name__ == "__main__":
    cpi = CPI(Object)
    cpi.simulate()
    cpi.plot_CPI_set()
    
    r_series = np.ones(1000) * 1.5
    x_0 = np.array([0,0])

    x_series = EoM_simulate(x_0,r_series,sys)
    plt.plot(x_series[:,0],x_series[:,1])
    X_0 = np.array([-1.5,0])
    R_series = np.zeros(1000)
    X_series = EoM_simulate(X_0,R_series,sys)
    plt.plot(X_series[:,0],X_series[:,1])
    plt.show()
    
    with open(Path(__file__).absolute().parent / "pickle_yard" / "CPI_set.pickle", "wb") as f:
        pickle.dump(cpi.CPI_set[-1],f)
    #計算された最大CPI集合を"CPI_set.pickle"として保存，Load_pickle.pyで読み取りが可能
