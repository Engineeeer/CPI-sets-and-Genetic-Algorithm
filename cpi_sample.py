import numpy as np

class Condition:
    """ 条件`ax <= b`を表現するクラス。 """
    def __init__(self, left_coefficients, right_value):
        """
        パラメータである左辺係数と右辺値を設定する。

            left_coefficients: a = [a1, a2, a3, ...]
            right_value: b
        """
        self.left_coefficients = left_coefficients
        self.right_value = right_value

    def check(self, x):
        """
        引数xが条件を満たすか判定する。
        """
        left_value = self.left_coefficients @ x
        return left_value <= self.right_value

class CPI:
    """ CPI集合クラス。 """
    def __init__(self, initial_conditions=None):
        """
        集合を初期化する。

            initial_conditions: CPIの初期条件群。条件なし、条件未定の場合は引数
                                を与えずに呼び出す。空リストを与えても良い。
        """
        # 初期コンディションがない場合は空リストを設定する。
        if initial_conditions is None:
            self.conditions = []
        else:
            self.conditions = initial_conditions

    def check(self, x):
        # xが全条件を満たせば真、1つでも満たさなければ偽を返す。
        for c in self.conditions:
            if not c.check(x):
                return False
        return True
