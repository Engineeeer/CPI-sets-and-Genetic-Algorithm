# はじめに
ここでは最大CPI集合と遺伝的アルゴリズムのPythonによる実装例を示します．

本コードは以下の環境で動作確認を行いました．

- Ubuntu 18.04 LTS
- python 3.6.9
- numpy 1.13.3
- matplotlib 2.1.1
- PuLP 1.6.0

## 構成ファイル
- cpi_object.py
- GA_class.py
- Load_pickle.py
- pickle_yard

### cpi_object.py
2次元最大CPI集合を計算するためのコードです．

パラメータ部のパラメータを指定し，CPIクラスのsimulate関数にある，
目標値幅w_min，w_maxを指定することで，2次元のシステムに対して最大CPI集合を計算できます．

計算結果はpickle_yardにCPI_set.pickleとして保存されます．

### GA_class.py
遺伝的アルゴリズムのコードです．

使用法は以下の通りです

- 初期化メソッドにある個体数(indivi_num)と突然変異確率(mutate_prob)を指定
- mkreference_series関数の目標値整形幅を指定
- long_mkreference_series関数のn及び返り値を指定
- simulate関数のパラメータ部を指定
- evaluate_individual関数内の評価関数を指定
- 実行部(if __name__ == "__main__"以下)で遺伝子の長さ，定常の目標値をそれぞれ指定

世代交代数100step毎にpickle_yardにgenes(step数).pickleとして保存されます．

### Load_pickle.py
pickle_yardにあるpickleをロードするためのコードです．

CPI用，GA用パラメータをそれぞれcpi_object,GA_classと同じものを指定することで，同様の結果をプロットします．

細かな設定はplot部を調節してください．

### pickle_yard
cpi_objectとGA_classより得られる最大CPI集合と優良遺伝子を保存するためのファイルです．


