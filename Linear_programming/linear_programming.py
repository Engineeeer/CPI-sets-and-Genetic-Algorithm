import pulp
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    #線形計画問題を定義
    problem = pulp.LpProblem("x-y", pulp.LpMaximize)

    #変数定義
    x = pulp.LpVariable("x",0,10,"Integer")
    y = pulp.LpVariable("y",0,10,"Integer")

    #評価関数
    problem += 2*x + 3*y

    #制約条件
    problem += x + 2*y <= 10
    problem += 2*x + y <= 8

    #最適化問題確認
    print(problem)
    
    #解く
    status = problem.solve()
    print(pulp.LpStatus[status])

    #結果表示
    print("Result")
    print("x:",x.value())
    print("y:",y.value())

    #グラフ生成
    #時系列
    x_1 = np.linspace(0,499,500) * 0.01
    
    y_1 = 5 - 0.5 * x_1
    y_2 = 8 - 2 * x_1
    #要素が全て0のlen(x_1)の配列を生成
    y_3 = np.zeros_like(x_1)
    y_4 = np.minimum(y_1,y_2)

    #y_4の領域を確認
    #plt.figure()
    #plt.plot(x_1,y_4)
    #plt.show()
    
    plt.figure()
    plt.plot(x_1,y_1, label = "x_1 + 2x_2 <= 10")
    plt.plot(x_1,y_2, label = "2x_1 + x_2 <= 8")
    plt.plot(x.value(),y.value(),"ro")
    plt.fill_between(x_1, y_3, y_4, where=y_4>y_3, facecolor = "yellow", alpha=0.3)
    plt.ylim(0,8.5)
    plt.xlim(0,4.5)
    plt.legend()
    plt.show()


