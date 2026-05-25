# Fortran ニュートン法ソルバー

## 概要

ニュートン法を用いた方程式 $f(x) = 0$ の根を求めるソルバーを Fortran で実装したものです。
将来的には Hartree-Fock (HF) 方程式の数値解法への応用を目指しています。

---

## ニュートン法

### アルゴリズム

初期値 $x_0$ から出発し、以下の漸化式で根に収束させます。

$$x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$$

### 収束条件

以下の両方を同時に満たしたとき収束とみなします。

$$|f(x_n)| < \varepsilon \quad \text{かつ} \quad |\Delta x| < \varepsilon$$

### 注意点

- $f'(x_n) \approx 0$ の場合は発散するため、閾値（$10^{-15}$）以下で停止します。
- 収束は初期値の選択に強く依存します。

---

## ファイル構成

```
Fortran/
├── newton_solver.f90   # newton_module（ソルバー本体）
└── test_newton.f90     # テストプログラム
```

---

## モジュール・サブルーチン仕様

### `newton_module`（`newton_solver.f90`）

#### `newton_solve` サブルーチン

```fortran
subroutine newton_solve(f, df, x0, tol, max_iter, root, converged, iter)
```

| 引数 | 種別 | 説明 |
|---|---|---|
| `f` | `procedure(func)` | 対象関数 $f(x)$ |
| `df` | `procedure(func)` | 導関数 $f'(x)$ |
| `x0` | `real(dp), in` | 初期値 |
| `tol` | `real(dp), in` | 収束判定閾値 $\varepsilon$ |
| `max_iter` | `integer, in` | 最大反復回数 |
| `root` | `real(dp), out` | 求めた根 |
| `converged` | `logical, out` | 収束したか否か |
| `iter` | `integer, out` | 実際の反復回数 |

関数インターフェース `func` の形式：

```fortran
function func(x) result(y)
  real(dp), intent(in) :: x
  real(dp) :: y
end function func
```

---

## ビルド・実行方法

```bash
gfortran -o test_newton newton_solver.f90 test_newton.f90
./test_newton
```

---

## テスト内容

`test_newton.f90` には以下の3ケースが含まれています。

| # | 方程式 | 初期値 | 真の解 |
|---|---|---|---|
| 1 | $x^2 - 2 = 0$ | 1.5 | $\sqrt{2} \approx 1.41421356...$ |
| 2 | $x^3 - x - 2 = 0$ | 1.5 | $\approx 1.52137970680457$ |
| 3 | $\cos(x) - x = 0$ | 0.5 | $\approx 0.73908513321516$ |

---

## 今後の課題：HF方程式への応用

現在の実装はスカラー（1変数）の方程式専用です。
HF方程式は行列・固有値問題であるため、以下の拡張が必要になります。

- スカラー版 → **ベクトル版**への拡張（$\mathbf{F}(\mathbf{x}) = \mathbf{0}$）
- 導関数 → **ヤコビアン行列** $J_{ij} = \partial F_i / \partial x_j$ の導入
- 線形方程式 $J \Delta \mathbf{x} = -\mathbf{F}$ を各反復で解く処理（LAPACK 利用を検討）
- 行列版の収束判定（残差ノルム $\|\mathbf{F}\|$ など）

