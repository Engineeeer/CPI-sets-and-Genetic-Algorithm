import numpy as np
import matplotlib.pyplot as plt
import pulp

a = 3
problem = pulp.LpProblem("x-y", pulp.LpMaximize)

x = pulp.LpVariable("x",-1,1,"Integer")
y = pulp.LpVariable("y",-1,1,"Integer")

problem += (-x) + y + a

print(problem)

status = problem.solve()
print(pulp.LpStatus[status])

print("Result")
print("x:",x.value())
print("y:",y.value())
#当たり前やけど,1次関数一本でも線形計画は可能!!
