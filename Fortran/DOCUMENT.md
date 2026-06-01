# Fortran ニュートン法ソルバー ドキュメント

---

## 1. 概要

`newton_module` は，Fortran で実装した**ニュートン法**による方程式ソルバーです。  
1次元スカラー方程式と N 次元連立方程式の両方に対応しています。

| ファイル | 役割 |
|---|---|
| `newton_solver.f90` | ソルバー本体（モジュール） |
| `test_newton.f90` | テストプログラム |

---

## 2. コンパイル・実行方法

```bash
# コンパイル
gfortran -o test_newton newton_solver.f90 test_newton.f90

# 実行
./test_newton
```

> **要件**: gfortran（Fortran 2003 以上対応）

---

## 3. アルゴリズム

### 3.1 ニュートン法（1次元）

スカラー方程式 $f(x) = 0$ の根を反復で求めます。

$$
x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)}
$$

**収束判定**（以下を両方満たしたとき収束）:

$$
|f(x_k)| < \varepsilon \quad \text{かつ} \quad |\Delta x_k| < \varepsilon
$$

---

### 3.2 ニュートン法（N次元）

連立方程式 $\mathbf{F}(\mathbf{x}) = \mathbf{0}$ の根を反復で求めます。

各反復で次の線形方程式を解き，解ベクトルを更新します。

$$
J(\mathbf{x}_k)\, \Delta\mathbf{x}_k = -\mathbf{F}(\mathbf{x}_k)
$$

$$
\mathbf{x}_{k+1} = \mathbf{x}_k + \Delta\mathbf{x}_k
$$

ここで $J(\mathbf{x})$ は**ヤコビアン行列**です。

$$
J_{ij} = \frac{\partial F_i}{\partial x_j}
$$

**収束判定**（以下を両方満たしたとき収束）:

$$
\|\mathbf{F}(\mathbf{x}_k)\| < \varepsilon \quad \text{かつ} \quad \|\Delta\mathbf{x}_k\| < \varepsilon
$$

ここで $\|\cdot\|$ はユークリッドノルム（`norm2`）。

---

### 3.3 線形方程式ソルバー（`lu_solve`）

N次元ニュートン法の内部で使用する $A\mathbf{x} = \mathbf{b}$ のソルバーです。  
**部分ピボット選択付きガウス消去法**を用います。

**手順**:
1. **前進消去**: 各ステップで絶対値最大の要素をピボットとして行を入れ替え，下三角部分をゼロにする。
2. **後退代入**: 上三角行列から解を逆順に求める。

ピボット要素が $10^{-15}$ 未満の場合，行列が特異に近いと判断してフラグを返します。

---

## 4. サブルーチン・関数の説明

### 4.1 `newton_solve`（1次元）

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
| `root` | `real(dp), out` | 計算した根 |
| `converged` | `logical, out` | 収束したか（`.true.` / `.false.`） |
| `iter` | `integer, out` | 実際の反復回数 |

**使用例**（`test_newton.f90` より）:

```fortran
call newton_solve(f1, df1, 1.5_dp, 1.0d-10, 100, root, converged, iter)
```

---

### 4.2 `newton_solve_nd`（N次元）

```fortran
subroutine newton_solve_nd(F_func, J_func, x0, n, tol, max_iter, root, converged, iter)
```

| 引数 | 種別 | 説明 |
|---|---|---|
| `F_func` | `procedure(func_nd)` | 連立方程式 $\mathbf{F}(\mathbf{x})$ を計算するサブルーチン |
| `J_func` | `procedure(jac_nd)` | ヤコビアン $J(\mathbf{x})$ を計算するサブルーチン |
| `x0(n)` | `real(dp), in` | 初期ベクトル |
| `n` | `integer, in` | 次元数 |
| `tol` | `real(dp), in` | 収束判定閾値 $\varepsilon$ |
| `max_iter` | `integer, in` | 最大反復回数 |
| `root(n)` | `real(dp), out` | 計算した解ベクトル |
| `converged` | `logical, out` | 収束したか（`.true.` / `.false.`） |
| `iter` | `integer, out` | 実際の反復回数 |

**使用例**（`test_newton.f90` より，2次元）:

```fortran
real(dp) :: x0(2), root(2)
x0 = [1.0_dp, 1.5_dp]
call newton_solve_nd(F2d, J2d, x0, 2, 1.0d-10, 100, root, converged, iter)
```

---

### 4.3 `lu_solve`（内部ルーチン）

```fortran
subroutine lu_solve(A_in, b_in, n, x, ok)
```

| 引数 | 種別 | 説明 |
|---|---|---|
| `A_in(n,n)` | `real(dp), in` | 係数行列 $A$ |
| `b_in(n)` | `real(dp), in` | 右辺ベクトル $\mathbf{b}$ |
| `n` | `integer, in` | 次元数 |
| `x(n)` | `real(dp), out` | 解ベクトル $\mathbf{x}$ |
| `ok` | `logical, out` | 正常に解けたか（行列が特異な場合 `.false.`） |

> このサブルーチンは `newton_solve_nd` の内部で自動的に呼ばれます。ユーザーが直接呼ぶ必要はありません。

---

## 5. インターフェースの説明

`newton_solve_nd` に渡す関数・ヤコビアンは，次のインターフェースに一致させる必要があります。

### `func_nd`（連立方程式）

```fortran
subroutine my_F(x, F)
  real(dp), intent(in)  :: x(:)   ! 入力ベクトル（長さ n）
  real(dp), intent(out) :: F(:)   ! 出力ベクトル（長さ n）
  F(1) = ...
  F(2) = ...
end subroutine my_F
```

### `jac_nd`（ヤコビアン）

```fortran
subroutine my_J(x, J)
  real(dp), intent(in)  :: x(:)    ! 入力ベクトル（長さ n）
  real(dp), intent(out) :: J(:,:)  ! ヤコビアン行列（n×n）
  J(1,1) = ...  ! ∂F1/∂x1
  J(1,2) = ...  ! ∂F1/∂x2
  ...
end subroutine my_J
```

> **注意**: `J(i,j)` は $\partial F_i / \partial x_j$ を格納します（行 = 方程式番号，列 = 変数番号）。

---

## 6. 変数の説明

### `newton_module` 内の共通変数

| 変数 | 型 | 説明 |
|---|---|---|
| `dp` | `integer, parameter` | 倍精度実数の種別値（`kind(1.0d0)`） |

### `newton_solve` の局所変数

| 変数 | 型 | 説明 |
|---|---|---|
| `x` | `real(dp)` | 現在の近似解 |
| `fx` | `real(dp)` | $f(x_k)$：現在点での関数値 |
| `dfx` | `real(dp)` | $f'(x_k)$：現在点での導関数値 |
| `dx` | `real(dp)` | 更新量 $\Delta x = -f(x)/f'(x)$ |
| `i` | `integer` | 反復カウンタ |

### `newton_solve_nd` の局所変数

| 変数 | 型 | 説明 |
|---|---|---|
| `x(n)` | `real(dp)` | 現在の近似解ベクトル |
| `Fx(n)` | `real(dp)` | $\mathbf{F}(\mathbf{x}_k)$：現在点での残差ベクトル |
| `Jx(n,n)` | `real(dp)` | $J(\mathbf{x}_k)$：現在点でのヤコビアン行列 |
| `dx(n)` | `real(dp)` | 更新量ベクトル $\Delta\mathbf{x}$（線形方程式の解） |
| `norm_Fx` | `real(dp)` | $\lVert\mathbf{F}(\mathbf{x}_k)\rVert$：残差のノルム |
| `norm_dx` | `real(dp)` | $\lVert\Delta\mathbf{x}_k\rVert$：更新量のノルム |
| `ok` | `logical` | `lu_solve` の成否フラグ |
| `i` | `integer` | 反復カウンタ |

### `lu_solve` の局所変数

| 変数 | 型 | 説明 |
|---|---|---|
| `A(n,n)` | `real(dp)` | 係数行列の作業コピー（消去過程で上書き） |
| `b(n)` | `real(dp)` | 右辺ベクトルの作業コピー |
| `factor` | `real(dp)` | 消去のための乗数 $A_{ik}/A_{kk}$ |
| `tmp` | `real(dp)` | スカラー行交換の一時変数 |
| `tmp_row(n)` | `real(dp)` | 行交換の一時ベクトル |
| `prow` | `integer` | ピボット行のインデックス |
| `pivot_idx(1)` | `integer` | `maxloc` の結果を受け取る一時配列 |
| `eps` | `real(dp), parameter` | 特異判定の閾値（$10^{-15}$） |
| `i, k` | `integer` | ループカウンタ |

---

## 7. 次元を増やす方法

3次元以上に拡張する場合，`test_newton.f90` に関数とヤコビアンを追加するだけで対応できます。  
ソルバー本体（`newton_solver.f90`）の変更は不要です。

**例: 4次元の場合**

```fortran
! 方程式の定義
subroutine F4d(x, F)
  real(dp), intent(in)  :: x(:)
  real(dp), intent(out) :: F(:)
  F(1) = ...  ! F1(x1, x2, x3, x4) = 0
  F(2) = ...  ! F2(x1, x2, x3, x4) = 0
  F(3) = ...  ! F3(x1, x2, x3, x4) = 0
  F(4) = ...  ! F4(x1, x2, x3, x4) = 0
end subroutine F4d

! ヤコビアンの定義
subroutine J4d(x, J)
  real(dp), intent(in)  :: x(:)
  real(dp), intent(out) :: J(:,:)
  J(1,1) = ...  ! ∂F1/∂x1
  ...           ! （4×4 = 16要素を埋める）
end subroutine J4d

! 呼び出し
real(dp) :: x0(4), root(4)
x0 = [...]
call newton_solve_nd(F4d, J4d, x0, 4, tol, max_iter, root, converged, iter)
```
