# NMsolver.jl

Julia によるニュートン法ソルバーです。
スカラー方程式 $f(x) = 0$ および多変数連立方程式 $F(\mathbf{x}) = \mathbf{0}$ を解きます。

---

## 目次

1. [セットアップ](#セットアップ)
2. [クイックスタート](#クイックスタート)
3. [API リファレンス](#api-リファレンス)
   - [newton — 1 変数](#newton--1-変数)
   - [newton_system — 多変数](#newton_system--多変数)
4. [使用例](#使用例)
   - [1 変数：基本的な使い方](#1-変数基本的な使い方)
   - [1 変数：導関数の指定方法](#1-変数導関数の指定方法)
   - [1 変数：ダンピングと収束履歴](#1-変数ダンピングと収束履歴)
   - [多変数：基本的な使い方](#多変数基本的な使い方)
   - [多変数：ヤコビアンの指定方法](#多変数ヤコビアンの指定方法)
   - [多変数：ダンピングと収束履歴](#多変数ダンピングと収束履歴)
5. [収束の確認とエラー処理](#収束の確認とエラー処理)
6. [性能評価スクリプト](#性能評価スクリプト)
7. [テスト](#テスト)
8. [ヒントとトラブルシューティング](#ヒントとトラブルシューティング)

---

## セットアップ

### 動作要件

- Julia 1.6 以上
- 依存パッケージ：`ForwardDiff`, `LinearAlgebra`（標準ライブラリ）

### パッケージ環境の構築

```julia
# プロジェクトディレクトリへ移動してから Julia を起動
cd NMsolver/julia

# パッケージモードで依存関係をインストール
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### ロード

Julia のプロジェクト環境内で `using NMsolver` とするだけで使えます。

```julia
julia --project=.   # プロジェクト環境で Julia を起動

julia> using NMsolver
```

エクスポートされる関数は `newton` と `newton_system` の 2 つです。

---

## クイックスタート

```julia
using NMsolver

# 1 変数：x^2 - 2 = 0 の正の根 (= √2) を求める
f = x -> x^2 - 2
x, iter, converged = newton(f, 1.0)
# x ≈ 1.4142135623730951, iter = 5, converged = true

# 多変数：x^2 + y^2 = 1, x = y の解 (= 1/√2, 1/√2) を求める
F = v -> [v[1]^2 + v[2]^2 - 1, v[1] - v[2]]
x, iter, converged = newton_system(F, [1.0, 0.5])
# x ≈ [0.7071..., 0.7071...], iter = 4, converged = true
```

---

## API リファレンス

### `newton` — 1 変数

```julia
newton(f, x0; df=nothing, derivative=:finite, tol=1e-10,
       maxiter=100, damping=false, min_alpha=1e-8,
       return_history=false, verbose=false)
```

スカラー方程式 $f(x) = 0$ をニュートン法で解く。

#### 引数

| 引数 | 型 | デフォルト | 説明 |
|:-----|:---|:----------:|:-----|
| `f` | 関数 | （必須） | 解くべき関数 $f(x)$ |
| `x0` | `Real` | （必須） | 初期値 |
| `df` | 関数 または `nothing` | `nothing` | $f$ の導関数 $f'(x)$。省略時は `derivative` の方法で自動計算 |
| `derivative` | `:finite` / `:forwarddiff` | `:finite` | `df` 省略時の微分方法（後述） |
| `tol` | `Float64` | `1e-10` | 収束判定の許容誤差 $\varepsilon$ |
| `maxiter` | `Int` | `100` | 最大反復回数 |
| `damping` | `Bool` | `false` | `true` のとき，バックトラック直線探索でステップ幅を調整 |
| `min_alpha` | `Float64` | `1e-8` | ダンピング時のステップ幅の最小値 |
| `return_history` | `Bool` | `false` | `true` のとき，反復履歴を追加で返す |
| `verbose` | `Bool` | `false` | `true` のとき，各反復の情報を標準出力に表示 |

#### 返り値

`return_history=false`（デフォルト）：

```julia
x, iter, converged = newton(f, x0)
```

| 変数 | 型 | 説明 |
|:-----|:---|:-----|
| `x` | `Float64` | 収束解（または最終近似値） |
| `iter` | `Int` | 実際の反復回数 |
| `converged` | `Bool` | 収束したかどうか |

`return_history=true`：

```julia
x, iter, converged, history = newton(f, x0; return_history=true)
```

`history` は各反復の `NamedTuple` の配列。各要素は以下のフィールドを持つ。

| フィールド | 説明 |
|:-----------|:-----|
| `iter` | 反復番号 |
| `x` | そのステップの近似解 |
| `residual` | $|f(x)|$ |
| `step` | $|\alpha \cdot \Delta x|$ |
| `alpha` | ダンピング係数 $\alpha$（通常は `1.0`） |

#### `derivative` オプション

| 値 | 説明 | 精度 | 速度 |
|:---|:-----|:----:|:----:|
| `:finite`（デフォルト） | 中心差分 $[f(x+h)-f(x-h)]/(2h)$ | $O(h^2) \approx 10^{-8}$ | 速い |
| `:forwarddiff` | ForwardDiff.jl による自動微分 | 機械精度（$\approx 10^{-16}$） | やや遅い |

---

### `newton_system` — 多変数

```julia
newton_system(F, x0; J=nothing, jacobian=:finite, tol=1e-10,
              maxiter=100, damping=false, min_alpha=1e-8,
              return_history=false, verbose=false)
```

$n$ 次元連立方程式 $F(\mathbf{x}) = \mathbf{0}$ をニュートン法で解く。

#### 引数

| 引数 | 型 | デフォルト | 説明 |
|:-----|:---|:----------:|:-----|
| `F` | 関数 | （必須） | 解くべきベクトル値関数 $F(\mathbf{x})$ |
| `x0` | `AbstractVector` | （必須） | 初期ベクトル |
| `J` | 関数 または `nothing` | `nothing` | ヤコビアン行列を返す関数。省略時は `jacobian` の方法で自動計算 |
| `jacobian` | `:finite` / `:forwarddiff` | `:finite` | `J` 省略時のヤコビアン計算方法 |
| `tol` | `Float64` | `1e-10` | 収束判定の許容誤差（$\|F(\mathbf{x})\|$ で評価） |
| `maxiter` | `Int` | `100` | 最大反復回数 |
| `damping` | `Bool` | `false` | `true` のとき，バックトラック直線探索でステップ幅を調整 |
| `min_alpha` | `Float64` | `1e-8` | ダンピング時のステップ幅の最小値 |
| `return_history` | `Bool` | `false` | `true` のとき，反復履歴を追加で返す |
| `verbose` | `Bool` | `false` | `true` のとき，各反復の情報を標準出力に表示 |

#### 返り値

`return_history=false`（デフォルト）：

```julia
x, iter, converged = newton_system(F, x0)
```

| 変数 | 型 | 説明 |
|:-----|:---|:-----|
| `x` | `Vector{Float64}` | 収束解（または最終近似値） |
| `iter` | `Int` | 実際の反復回数 |
| `converged` | `Bool` | 収束したかどうか |

`return_history=true` の場合は `newton` と同様に 4 値で返す。
`history` の各要素の `x` フィールドは `Vector{Float64}` になる。

#### 収束判定

| 判定条件 | 説明 |
|:---------|:-----|
| $\lVert F(\mathbf{x}_k)\rVert < \varepsilon$ | 残差ノルムが閾値以下 |
| $\lVert\alpha\,\Delta\mathbf{x}_k\rVert < \varepsilon$ | 更新量ノルムが閾値以下 |

2 つの条件を **両方** 満たしたとき収束と判定する。

---

## 使用例

### 1 変数：基本的な使い方

#### 最もシンプルな呼び出し（数値微分・自動）

導関数を渡さなければ，中心差分で自動的に近似する。

```julia
using NMsolver

# 例1: x^2 - 2 = 0 → x = √2
f = x -> x^2 - 2
x, iter, converged = newton(f, 1.0)

println("解: $x")          # 1.4142135623730951
println("反復回数: $iter")  # 5
println("収束: $converged") # true
```

```julia
# 例2: cos(x) - x = 0（不動点）
x, iter, converged = newton(x -> cos(x) - x, 0.5)
println("解: $x")  # 0.7390851332151607
```

```julia
# 例3: 組み込み関数もそのまま使える
# sin(x) = 0 の解（x ≈ π）
x, iter, converged = newton(sin, 3.0)
println("解: $x")  # 3.141592653589793
```

---

### 1 変数：導関数の指定方法

#### パターン A：解析的導関数を渡す（最速・最高精度）

```julia
f  = x -> x^3 - x - 2
df = x -> 3x^2 - 1          # 解析的に求めた導関数

x, iter, converged = newton(f, 1.5; df=df)
println("解: $x")    # 1.5213797068045675
println("反復: $iter") # 5
```

#### パターン B：数値微分（デフォルト）

```julia
f = x -> x^3 - x - 2
x, iter, converged = newton(f, 1.5)          # derivative=:finite が省略値
# または明示的に
x, iter, converged = newton(f, 1.5; derivative=:finite)
```

#### パターン C：ForwardDiff による自動微分（高精度・手間なし）

```julia
f = x -> x^3 - x - 2
x, iter, converged = newton(f, 1.5; derivative=:forwarddiff)
println("解: $x")  # 1.5213797068045675（機械精度）
```

---

### 1 変数：ダンピングと収束履歴

#### ダンピング付きニュートン法

初期値が解から離れているときの発散を防ぐ。

```julia
f = x -> x^2 - 2

# damping=true でバックトラック直線探索を有効化
x, iter, converged = newton(f, 10.0; damping=true)
println("解: $x")           # √2
println("収束: $converged")  # true
```

#### 収束履歴の取得

```julia
f = x -> x^2 - 2
x, iter, converged, history = newton(f, 1.5; return_history=true)

# history は各反復の NamedTuple 配列
for h in history
    @printf("iter=%d  x=%.10f  residual=%.2e  step=%.2e  alpha=%.4f\n",
            h.iter, h.x, h.residual, h.step, h.alpha)
end
```

出力例：

```
iter=1  x=1.4166666667  residual=7.64e-03  step=8.33e-02  alpha=1.0000
iter=2  x=1.4142156863  residual=6.01e-06  step=2.51e-04  alpha=1.0000
iter=3  x=1.4142135624  residual=3.72e-12  step=2.13e-07  alpha=1.0000
iter=4  x=1.4142135624  residual=0.00e+00  step=1.57e-15  alpha=1.0000
iter=5  x=1.4142135624  residual=0.00e+00  step=0.00e+00  alpha=1.0000
```

#### 各反復の情報をリアルタイム表示

```julia
x, iter, converged = newton(x -> x^2 - 2, 1.5; verbose=true)
```

出力例：

```
iter=1, x=1.4166666666666667, residual=0.0076388..., step=0.0833..., alpha=1.0
iter=2, x=1.4142156862745099, residual=6.007e-06, step=0.0024..., alpha=1.0
...
```

---

### 多変数：基本的な使い方

`newton_system` に関数とベクトル型の初期値を渡す。

```julia
using NMsolver

# 例1: x^2 + y^2 = 1, x = y → (1/√2, 1/√2)
F = v -> [v[1]^2 + v[2]^2 - 1,
          v[1] - v[2]]

x, iter, converged = newton_system(F, [1.0, 0.5])
println("解: $x")           # [0.7071..., 0.7071...]
println("残差: $(norm(F(x)))")  # ≈ 0.0
```

```julia
# 例2: 3 変数 x+y+z=6, x^2+y^2+z^2=14, xyz=6 → (1, 2, 3)
F = v -> [
    v[1] + v[2] + v[3] - 6,
    v[1]^2 + v[2]^2 + v[3]^2 - 14,
    v[1] * v[2] * v[3] - 6
]
x, iter, converged = newton_system(F, [0.5, 1.5, 3.5])
println("解: $x")    # [1.0, 2.0, 3.0]
println("反復: $iter") # 6
```

```julia
# 例3: n 変数 x[i]^2 = i → x[i] = √i（n=5 の例）
n = 5
F = v -> [v[i]^2 - Float64(i) for i in 1:n]
x, iter, converged = newton_system(F, 2.0 * ones(n))
println("解: $x")  # [1.0, 1.414..., 1.732..., 2.0, 2.236...]
```

---

### 多変数：ヤコビアンの指定方法

#### パターン A：解析的ヤコビアンを渡す（最速）

```julia
F = v -> [v[1]^2 + v[2]^2 - 1,
          v[1] - v[2]]

# J(v) は n×n 行列を返す関数。J[i,j] = ∂F_i/∂x_j
Jfunc = v -> [2v[1]  2v[2]
              1.0    -1.0]

x, iter, converged = newton_system(F, [1.0, 0.5]; J=Jfunc)
```

3 変数の例：

```julia
F = v -> [
    v[1] + v[2] + v[3] - 6,
    v[1]^2 + v[2]^2 + v[3]^2 - 14,
    v[1] * v[2] * v[3] - 6
]

Jfunc = v -> [
    1.0        1.0        1.0
    2v[1]      2v[2]      2v[3]
    v[2]*v[3]  v[1]*v[3]  v[1]*v[2]
]

x, iter, converged = newton_system(F, [0.5, 1.5, 3.5]; J=Jfunc)
```

#### パターン B：数値ヤコビアン（デフォルト）

```julia
F = v -> [v[1]^2 + v[2]^2 - 1, v[1] - v[2]]
x, iter, converged = newton_system(F, [1.0, 0.5])          # jacobian=:finite が省略値
```

#### パターン C：ForwardDiff による自動ヤコビアン

```julia
F = v -> [v[1]^2 + v[2]^2 - 1, v[1] - v[2]]
x, iter, converged = newton_system(F, [1.0, 0.5]; jacobian=:forwarddiff)
```

---

### 多変数：ダンピングと収束履歴

#### ダンピング付き

```julia
F = v -> [v[1]^2 + v[2]^2 - 1, v[1] - v[2]]

x, iter, converged = newton_system(F, [5.0, 0.1]; damping=true)
println("解: $x")  # [0.7071..., 0.7071...]
```

#### 収束履歴の取得と可視化

```julia
using NMsolver, LinearAlgebra

F = v -> [v[1]^2 + v[2]^2 - 1, v[1] - v[2]]
x, iter, converged, history = newton_system(F, [1.0, 0.5]; return_history=true)

# 残差の変化を確認
for h in history
    println("iter=$(h.iter)  residual=$(h.residual)  step=$(h.step)")
end
```

出力例（2次収束の確認）：

```
iter=1  residual=0.08838...  step=0.20710...
iter=2  residual=0.00312...  step=0.06050...
iter=3  residual=3.90e-06    step=0.00249...
iter=4  residual=6.10e-12    step=3.90e-06
```

残差が反復ごとに約 2 乗されて小さくなっており，ニュートン法の **2 次収束** が確認できる。

---

## 収束の確認とエラー処理

### 返り値 `converged` の確認

```julia
x, iter, converged = newton(f, x0)

if converged
    println("収束しました: x = $x（$iter 回）")
else
    println("収束しませんでした（最大反復回数 $iter 回に達しました）")
    println("最終近似値: x = $x")
end
```

### よくある失敗パターンと対処

#### 1. 収束しない（`converged = false`）

```julia
# 初期値が悪い場合
x, iter, converged = newton(x -> x^2 - 2, 0.0)
# → 導関数 f'(0) = 0 でエラー

# 対処：初期値をずらすか，ダンピングを使う
x, iter, converged = newton(x -> x^2 - 2, 0.1; damping=true)
```

#### 2. 反復回数を増やしたい

```julia
x, iter, converged = newton(f, x0; maxiter=500)
```

#### 3. 許容誤差を緩くしたい（収束しやすくする）

```julia
x, iter, converged = newton(f, x0; tol=1e-6)
```

#### 4. 許容誤差を厳しくしたい（より高精度）

```julia
# tol を機械イプシロンに近づける
x, iter, converged = newton(f, x0; tol=1e-14)
```

#### 5. 複数の解が存在する場合

ニュートン法は**初期値に最も近い解**に収束する傾向がある。

```julia
f = x -> x^2 - 1    # 解は x = +1 と x = -1

x1, _, _ = newton(f,  0.5)   # x1 ≈  1.0
x2, _, _ = newton(f, -0.5)   # x2 ≈ -1.0
```

---

## 性能評価スクリプト

`julia/scripts/evaluate.jl` を実行すると，1〜10 変数の問題について
4 種類の手法（数値微分・解析的・自動微分・ダンピング）を比較し，
結果を CSV とグラフとして `julia/results/` に保存する。

### 実行方法

```bash
cd NMsolver/julia
julia --project=. scripts/evaluate.jl
```

### 評価内容

| 評価項目 | 詳細 |
|:--------|:-----|
| **テスト問題 1** | 分離型二乗系：$F_i(x) = x_i^2 - i = 0$，解 $x_i = \sqrt{i}$ |
| **テスト問題 2** | 三対角型：$x_{i-1} + x_i^2 + x_{i+1} = c_i$，解 $x_i = 1$ |
| **評価指標** | 解の誤差 L2 ノルム，残差 L2 ノルム，反復回数，計算時間（μs） |
| **計測方法** | 500 回試行の最小値（JIT ウォームアップ後） |

### 出力ファイル

| ファイル | 内容 |
|:--------|:-----|
| `results/comparison.csv` | 全結果の数値データ |
| `results/error_vs_n.png` | 誤差 vs 変数数（分離型） |
| `results/iters_vs_n.png` | 反復回数 vs 変数数 |
| `results/time_vs_n.png` | 計算時間 vs 変数数 |
| `results/tridiagonal_error_vs_n.png` | 誤差 vs 変数数（三対角型） |
| `results/summary.png` | 上記 4 グラフのまとめ |

---

## テスト

### テストの実行

```bash
cd NMsolver/julia
julia --project=. -e 'using Pkg; Pkg.test()'
```

### テスト内容

| テストセット | 内容 |
|:------------|:-----|
| スカラー（解析的導関数） | $x^2-2=0$，$x^3-x-2=0$ |
| スカラー（数値微分） | $x^2-2=0$，$\sin x=0$，$e^x-3=0$ |
| 収束しない場合 | `maxiter=2` で打ち切り確認 |
| 2 変数（数値微分） | 円と直線の交点，対称な 2 元連立方程式 |
| 2 変数（解析的ヤコビアン） | 同上 |
| 自動微分（ForwardDiff） | 1 変数・2 変数の両方 |
| ダンピングと収束履歴 | alpha の範囲，履歴の長さの確認 |
| 3 変数 | $x+y+z=6,\ x^2+y^2+z^2=14,\ xyz=6$ |
| 4〜10 変数 | 分離型二乗系 $x_i^2 = i$ |
| 10 変数（三対角型） | $x_{i-1}+x_i^2+x_{i+1}=c_i$ |

---

## ヒントとトラブルシューティング

### 導関数の指定方法の選び方

```
導関数を手で求められる？
  Yes → df=（解析的導関数）  ← 最速・最高精度
  No  →
      精度より手軽さ優先？
        Yes → derivative=:finite（デフォルト）  ← 簡単・十分な精度
        No  → derivative=:forwarddiff  ← 機械精度・コード変更不要
```

### ヤコビアンの指定方法の選び方

```
ヤコビアンを手で求められる？
  Yes → J=（解析的ヤコビアン）  ← 最速
  No  →
      変数が多い（n > 10）？
        Yes → jacobian=:forwarddiff  ← 高精度・スケーラブル
        No  → jacobian=:finite       ← デフォルトで十分
```

### `F` の定義で注意すること

`newton_system` に渡す関数 `F` は `Vector` を受け取り `Vector` を返す必要がある。

```julia
# 正しい例
F = v -> [v[1]^2 + v[2]^2 - 1,
          v[1] - v[2]]

# 誤りやすい例：スカラーを返してしまう
F_wrong = v -> v[1]^2 + v[2]^2 - 1   # これは newton_system には使えない
```

### ヤコビアン行列の書き方

$J[i, j] = \partial F_i / \partial x_j$ となるように定義する。

```julia
# F(v) = [F1, F2] の場合
# J = [∂F1/∂x1  ∂F1/∂x2
#      ∂F2/∂x1  ∂F2/∂x2]

Jfunc = v -> [∂F1/∂x1(v)  ∂F1/∂x2(v)
              ∂F2/∂x1(v)  ∂F2/∂x2(v)]

# Julia の行列リテラル：スペース区切りで列，改行または ; で行
Jfunc = v -> [2v[1]  2v[2];
              1.0    -1.0]
```

### ForwardDiff が使えない場合

`ForwardDiff.jl` は関数が **型安定** （`dual number` を透過できる）である必要がある。
条件分岐で異なる型を返すような関数では動作しないことがある。その場合は `:finite` を使う。

```julia
# ForwardDiff が苦手な例（条件分岐で型が変わる場合）
f_bad = x -> x > 0 ? x^2 - 2 : -(x^2 - 2)

# この場合は数値微分を使う
x, iter, converged = newton(f_bad, 1.0; derivative=:finite)
```

### 解の検証方法

収束後は必ず残差（$|f(x)|$ または $\|F(\mathbf{x})\|$）を確認することを推奨する。

```julia
using LinearAlgebra

# 1 変数
x, iter, converged = newton(f, x0)
println("残差 |f(x)| = $(abs(f(x)))")

# 多変数
x, iter, converged = newton_system(F, x0)
println("残差 ||F(x)|| = $(norm(F(x)))")
```
