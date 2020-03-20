import pulp
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    #線形計画問題を定義
    problem = pulp.LpProblem("x-y", pulp.LpMinimize)

    #変数定義
    x = pulp.LpVariable("x",-1,1,"Continuous")
    y = pulp.LpVariable("y",-1,1,"Continuous")

    #評価関数
    problem += x * 2 - y

    #制約条件
    problem += x >= -1
    problem += x <= 1

    #最適化問題確認
    print(problem)
    
    #解く
    status = problem.solve()
    print(pulp.LpStatus[status])

    #結果表示
    print("Result")
    print("x:",x.value())
    print("y:",y.value())
