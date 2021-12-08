# 進化計算ライブラリ
本ライブラリは複数の連続値の最適化ができます．
最適化アルゴリズムにはDX-NESが採用されています．


## How to use
`$ brew install eigen`<br>

`$ mkdir build`<br>
`$ cd build`<br>
`$ cmake ..`<br>

`$ make`<br>
`$ ./a.out`<br>


## 構成
```
./
 +-- calVal/
 |     +-- calcurateVector : ベクトル・行列計算をまとめたもの
 |
 +-- dxnes/
 |     +-- dxnes : 進化計算アルゴリズムDXNESのアルゴリズムが実装してある
 |
 +-- eval/
       +-- evaluator : 評価関数といくつかのベンチマーク関数の実装
```
