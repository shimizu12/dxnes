## 進化計算ライブラリ
本ライブラリは複数の連続値の最適化ができます．
最適化アルゴリズムにはDX-NESが採用されています．


## How to use
`eval/evaluator.cpp`内の`evalFunc()`が評価関数になります，適宜上書きしてください．

`main.cpp`内で，`dxnes`のインスタンス`n`と最適化させたい連続値の数と同じ長さを持つベクトル`v`を宣言し，`n.dxnes_start( v )`

`$ cd build`<br>
`$ cmake ..`<br>
`$ make`<br>
`$ ./a.out`<br>
