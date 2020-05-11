import numpy as np
import numpy.random
import random
import matplotlib.pyplot as plt
import pickle
from pathlib import Path
from scipy import signal as sg

class GA:
    def __init__(self,genes_lenth):
        #indivi_numで個体数を決定
        self.indivi_num = 100
        #mutate_probで突然変異確率を決定
        self.mutate_prob = 0.1
        #初期個体を生成
        self.individual = [self.mkgene(genes_lenth) for _ in range(self.indivi_num)]
    
    def mkgene(self,genes_lenth):
        #遺伝子{x \in R^{len(genes_lenth)} | x = 0 or 1} を生成
        return np.random.randint(0, 2, genes_lenth, dtype = np.uint8)
    
    def crossover(self,gene1,gene2):
        #一点交叉を実装
        a = gene1.copy()
        idx = np.random.randint(1, len(gene1))
        cross_gene1 = np.hstack((a[:idx], gene2[idx:]))
        cross_gene2 = np.hstack((gene2[:idx], a[idx:]))
        return cross_gene1, cross_gene2

    def mutate_one_index(self,gene):
        #突然変異を実装
        #遺伝子の一点を0->1 or 1->0に変換
        gene = gene.copy()
        idx = np.random.randint(0,len(gene))
        gene[idx] = 1 - gene[idx]
        return gene

    def mkreference_series(self,gene):
        #目標値を256で量子化及び正規化
        #今回は{x|-2<=x<=2}を生成する為に
        #4/255.0で{x|0<=x<=4}を作り,平行移動(-2.0)で{x|-2<=x<=2}を生成
        #np.packbitsは8bitのバイナリを0-255の10進数に対応させる
        r_series = np.array(np.packbits(gene), dtype = np.float64) * 4.0 / 255.0 - 2.0
        return r_series

    def long_mkreference_series(self,gene):
        n = 5
        r_series = self.mkreference_series(gene)
        lis = []
        for r in r_series:
            for _ in range(n):
                lis.append(r)
        return np.hstack([np.array(lis), np.ones(n * len(r_series))*1.5])

    def simulate(self,r_series):
        #実際の状態方程式におけるシミュレーション
        #印加される目標値に対して状態と入力を出力
        Num = len(r_series)
        #パラメータ
        Ts = 0.05
        K_p = 1.0
        K_d = 0.5

        x_series = np.zeros((Num, 2), dtype = np.float64)
        u_series = np.zeros((Num, 1), dtype = np.float64)
    
        A_c = np.array([[0,1],[-K_p,-0.5 - K_d]])
        B_c = np.array([[0],[1.0]])
        C_c = np.array([[-K_p,-K_d]])
        D_c = np.array([K_p])

        c2d = sg.cont2discrete((A_c,B_c,C_c,D_c),dt = Ts)
        A_d,B_d,C_d,D_d = c2d[0],c2d[1].reshape(2),c2d[2],c2d[3]

        for i in range(Num-1):
            x = A_d @ x_series[i] + B_d * r_series[i]
            u = C_d @ x_series[i] + D_d * r_series[i] 
        
            x_series[i+1] = x
            u_series[i+1] = u

        return x_series, u_series

    def evaluate_individual(self,ind,reference):
        #各個体の評価を実装
        r_series = self.long_mkreference_series(ind)
        x_series, u_series = self.simulate(r_series)
        #評価関数
        msae = 1 * np.abs(reference - x_series[:,0]).mean() + 1 * np.abs(reference - r_series).mean()
        return 1 / msae

    def GA_LOOP(self,Reference):
        #遺伝的アルゴリズムのメインループ
        cnt = 0
        while True:
            cnt += 1
            #評価点がついたリストを生成
            evaluations = [self.evaluate_individual(ind,Reference) for ind in self.individual]
            #評価点に基づいたソート
            elite_individual, elite_evaluation = sorted(zip(self.individual, evaluations), key=lambda item: item[1])[-1]
            if cnt % 100 == 0:
            #100ループ毎のelite_individualを"pickle_yard"に保存
                with open(Path(__file__).absolute().parent / "pickle_yard" / "genes{}.pickle".format(cnt), "wb") as f:
                    pickle.dump(self.individual, f)
            print("\r" + str(cnt) + ": " + str(elite_evaluation), end="")
            #次世代を生成,elite_individualを最初の次世代へ
            new_individual = [elite_individual]
            #次世代が現世代以上になるまで交叉を繰り返し,それらを次世代へ
            #random.choicesはweightsに依存した要素をself.individualよりk個取得
            while len(new_individual) <= len(self.individual):
                a, = random.choices(self.individual, weights=evaluations, k=1)
                b, = random.choices(self.individual, weights=evaluations, k=1)
                #同じ個体が選択された場合,片方を選択し直す
                while id(a) == id(b):
                    b, = random.choices(self.individual, weights=evaluations, k=1)
                a,b = self.crossover(a,b)
                new_individual.append(a)
                new_individual.append(b)
            self.individual = new_individual[:len(self.individual)]
            
            for i, _ in enumerate(self.individual[1:], start = 1):
                #elite_individual以外に対して突然変異を実施
                #self.individual[1:]はelite_individual以外を選択
                if random.random() <= self.mutate_prob:
                    self.individual[i] = self.mutate_one_index(self.individual[i])

if __name__ == "__main__":
    #8bitで1ステップなので,50ステップ生成
    genes_lenth = 8 * 50
    #ステップ状の目標値を1.5とする
    Step_Reference = 1.5
    #インスタンス生成
    ga = GA(genes_lenth)
    #GA_LOOPを走らせる
    ga.GA_LOOP(Step_Reference)
