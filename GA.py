import numpy as np
import numpy.random
import random
import matplotlib.pyplot as plt
import pickle
from pathlib import Path
from scipy import signal as sg

def mkgene(n):
    return np.random.randint(0, 2, n, dtype = np.uint8)

def mutate_one_index(gene):
    #geneは1次元array
    gene = gene.copy()
    idx = np.random.randint(0,len(gene))
    gene[idx] = 1 - gene[idx]
    return gene

def crossover(gene1, gene2):
    a = gene1.copy()
    idx = np.random.randint(1, len(gene1))
    cross_gene1 = np.hstack((a[:idx], gene2[idx:]))
    cross_gene2 = np.hstack((gene2[:idx], a[idx:]))
    return cross_gene1, cross_gene2

def test_crossover():
    n = 10
    g1 = np.zeros(n,dtype = np.uint8)
    g2 = np.ones(n,dtype = np.uint8)
    h1, h2 = crossover(g1,g2)
    print(h1)
    print(h2)
    return 0

def mkreference_series(gene):
    #np.packbitsは8bit2進数=>10進数を生成，この計算によって |-1.0~-1/128~1/128~1.0|(r_series)の積分値(np.cumsum)が目標値となる
    r_series = np.array(np.packbits(gene), dtype = np.float64) * 4.0 / 255.0 - 2.0
    #return np.cumsum(r_series)
    return r_series

def long_mkreference_series(gene):
    n = 5
    r_series = mkreference_series(gene)
    lis = []
    for r in r_series:
        for _ in range(n):
            lis.append(r)
    return np.hstack([np.array(lis),np.ones(n * len(r_series))*1.5])

def simulate(r_series):
    
    Num = len(r_series)   
    Ts = 0.05
    K_p = 1.0
    K_d = 0.5
    #保存用ログ確保
    x_series = np.zeros((Num, 2), dtype = np.float64)
    u_series = np.zeros((Num, 1), dtype = np.float64)
    
    A_c = np.array([[0,1],[-K_p,-0.5 - K_d]])
    B_c = np.array([[0],[1.0]])
    C_c = np.array([[-K_p,-K_d]])
    D_c = np.array([K_p])

    c2d = sg.cont2discrete((A_c,B_c,C_c,D_c),dt = Ts)
    A_d,B_d,C_d,D_d = c2d[0],c2d[1].reshape(2),c2d[2],c2d[3]
    #sys = np.array([A_d,B_d,C_d,D_d])

    for i in range(Num-1):
        x = A_d @ x_series[i] + B_d * r_series[i]
        u = C_d @ x_series[i] + D_d * r_series[i] 
        
        x_series[i+1] = x
        u_series[i+1] = u

    return x_series, u_series
#水曜日はこのへんまでやってます
def absfunc(u):
    k = 1.0
    return k * np.abs(u)

#u_barrier(list) -> return score(int)
def u_barrier(u):
    L = -5.0
    R = 5.0
    return 0.0 if L <= u <= R else absfunc(u)

def test_u_barrier():
    u_series = (np.random.random_sample(10) * 2 - 1) * 10
    print(u_series)
    #vectorizeは関数をベクトル化,即ちu_listの各要素に全てに対して関数を当てる.
    ub = np.vectorize(u_barrier)
    
    return print(ub(u_series))

def evaluate_gene(gene):
    r_series = long_mkreference_series(gene)
    x_series, u_series = simulate(r_series)
    msae = 10 * np.abs(1.5 - x_series[:,0]).mean() + np.abs(1.5 - r_series).mean()
    #msae = np.abs(np.pi/2 - x_series[:,0]).mean() + np.abs(np.pi/2 - r_series).mean() #+ np.abs(r_series[1:] - r_series[:-1]).mean() + 0.1 * np.sum(evaluate_gene.ub(u_series))
    return 1 / msae
evaluate_gene.ub = np.vectorize(u_barrier)

if __name__ == "__main__":
    if  "genes" not in globals():
        individual = 100
        mutation_prob = 0.1
        genes = [mkgene(8*50) for _ in range(individual)]
        cnt = 0

    while True:
        cnt += 1
        evaluations = [evaluate_gene(g) for g in genes]
        elite_gene, elite_evaluation = sorted(zip(genes, evaluations), key=lambda item: item[1])[-1]
        if cnt % 100 == 0:
            with open(Path(__file__).absolute().parent / "pickle_yard" / "genes{}.pickle".format(cnt), "wb") as f:
                pickle.dump(genes, f)

        print("\r" + str(cnt) + ": " + str(elite_evaluation), end="")
        new_genes = [elite_gene]
        while len(new_genes) < len(genes):
            a, = random.choices(genes, weights=evaluations, k=1)
            b, = random.choices(genes, weights=evaluations, k=1)
            while id(a) == id(b):
                b, = random.choices(genes, weights=evaluations, k=1)

            a,b = crossover(a,b)
            new_genes.append(a)
            new_genes.append(b)
        genes = new_genes[:len(genes)]

        for i, _ in enumerate(genes[1:], start = 1):
            if random.random() < mutation_prob:
                genes[i] = mutate_one_index(genes[i])

