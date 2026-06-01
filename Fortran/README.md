# NMsolver — Fortran 版

Fortran によるニュートン法ソルバーです。
スカラー方程式 $f(x) = 0$ および多変数連立方程式 $F(\mathbf{x}) = \mathbf{0}$ を解きます。

---

## 目次

1. [セットアップ](#セットアップ)
2. [クイックスタート](#クイックスタート)
3. [API リファレンス](#api-リファレンス)
   - [newton_solve — 1 変数](#newton_solve--1-変数)
   - [newton_solve_nd — 多変数](#newton_solve_nd--多変数)
   - [lu_solve — 内部線形ソルバー](#lu_solve--内部線形ソルバー)
4. [インターフェース要件](#インターフェース要件)
5. [使用例](#使用例)
   - [1 変数：基本的な使い方](#1-変数基本的な使い方)
   - [2 変数の連立方程式](#2-変数の連立方程式)
   - [3 変数の連立方程式](#3-変数の連立方程式)
   - [次元を増やす方法](#次元を増やす方法)
6. [収束の確認](#収束の確認)
7. [テスト例題の実行](#テスト例題の実行)
8. [ヒントとトラブルシューティング](#ヒントとトラブルシューティング)

---

## セットアップ

### 動作要件

- gfortran（Fortran 2003 以上対応）
- 外部ライブラリ不要（標準 Fortran のみ）

### コンパイルと実行

```bash
# Fortran/ ディレクトリへ移動
cd NMsolver/Fortran

# コンパイル（モジュールファイルを先に，次にテストプログラムを渡す）
gfortran -o test_newton newton_solver.f90 test_newton.f90

# 実行
./test_newton
```

最適化オプションを付けてコンパイルする場合：

```bash
gfortran -O2 -o test_newton newton_solver.f90 test_newton.f90
```

### ファイル構成

| ファイル | 役割 |
|:--------|:-----|
| `newton_solver.f90` | ソルバー本体（`newton_module` モジュール） |
| `test_newton.f90` | テストプログラム（例題 5 問） |
| `DOCUMENT.md` | 詳細アルゴリズムドキュメント |

---

## クイックスタート

```fortran
program my_program
  use newton_module
  implicit none

  real(dp) :: root
  logical  :: converged
  integer  :: iter

  ! x^2 - 2 = 0 の根（= √2）を求める
  call newton_solve(my_f, my_df, 1.5_dp, 1.0d-10, 100, root, converged, iter)
  write(*,*) "root =", root        ! 1.41421356237...
  write(*,*) "iter =", iter        ! 5
  write(*,*) "converged =", converged  ! T

contains

  function my_f(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = x**2 - 2.0_dp
  end function my_f

  function my_df(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = 2.0_dp * x
  end function my_df

end program my_program
```

---

## API リファレンス

### `newton_solve` — 1 変数

```fortran
subroutine newton_solve(f, df, x0, tol, max_iter, root, converged, iter)
```

スカラー方程式 $f(x) = 0$ をニュートン法で解く。

#### 引数

| 引数 | 種別 | 入出力 | 説明 |
|:-----|:-----|:------:|:-----|
| `f` | `procedure(func)` | 入力 | 対象関数 $f(x)$ |
| `df` | `procedure(func)` | 入力 | 導関数 $f'(x)$（**必須**） |
| `x0` | `real(dp)` | 入力 | 初期値 |
| `tol` | `real(dp)` | 入力 | 収束判定閾値 $\varepsilon$ |
| `max_iter` | `integer` | 入力 | 最大反復回数 |
| `root` | `real(dp)` | 出力 | 計算した根（収束しなかった場合は最終近似値） |
| `converged` | `logical` | 出力 | 収束したか（`.true.` / `.false.`） |
| `iter` | `integer` | 出力 | 実際の反復回数 |

#### 収束判定

以下の **2 条件を両方** 満たしたとき収束と判定する。

$$|f(x_k)| < \varepsilon \quad \text{かつ} \quad |\Delta x_k| < \varepsilon$$

#### 出力例

各反復でステップ情報が標準出力に表示される。

```
  iter=   1  x=  1.416666666667E+00  f(x)=  6.3611E-02
  iter=   2  x=  1.414215686275E+00  f(x)=  6.0073E-06
  iter=   3  x=  1.414213562375E+00  f(x)=  5.2915E-12
  iter=   4  x=  1.414213562373E+00  f(x)= -4.4409E-16
  iter=   5  x=  1.414213562373E+00  f(x)= -4.4409E-16
```

---

### `newton_solve_nd` — 多変数

```fortran
subroutine newton_solve_nd(F_func, J_func, x0, n, tol, max_iter, root, converged, iter)
```

$n$ 次元連立方程式 $F(\mathbf{x}) = \mathbf{0}$ をニュートン法で解く。

#### 引数

| 引数 | 種別 | 入出力 | 説明 |
|:-----|:-----|:------:|:-----|
| `F_func` | `procedure(func_nd)` | 入力 | 残差ベクトル $F(\mathbf{x})$ を計算するサブルーチン |
| `J_func` | `procedure(jac_nd)` | 入力 | ヤコビアン $J(\mathbf{x})$ を計算するサブルーチン（**必須**） |
| `x0(n)` | `real(dp)` | 入力 | 初期ベクトル |
| `n` | `integer` | 入力 | 次元数 |
| `tol` | `real(dp)` | 入力 | 収束判定閾値 $\varepsilon$ |
| `max_iter` | `integer` | 入力 | 最大反復回数 |
| `root(n)` | `real(dp)` | 出力 | 計算した解ベクトル |
| `converged` | `logical` | 出力 | 収束したか（`.true.` / `.false.`） |
| `iter` | `integer` | 出力 | 実際の反復回数 |

#### 収束判定

$$\|F(\mathbf{x}_k)\| < \varepsilon \quad \text{かつ} \quad \|\Delta\mathbf{x}_k\| < \varepsilon$$

$\|\cdot\|$ はユークリッドノルム（組み込み関数 `norm2`）。

#### 出力例

```
  iter=   1  ||F(x)||=  2.5000E-01  ||dx||=  2.0711E-01
  iter=   2  ||F(x)||=  3.1250E-02  ||dx||=  6.0660E-02
  iter=   3  ||F(x)||=  4.8828E-04  ||dx||=  3.9052E-03
  iter=   4  ||F(x)||=  1.1921E-07  ||dx||=  9.5368E-07
```

---

### `lu_solve` — 内部線形ソルバー

```fortran
subroutine lu_solve(A_in, b_in, n, x, ok)
```

$A\mathbf{x} = \mathbf{b}$ を部分ピボット選択付きガウス消去法で解く。
**`newton_solve_nd` の内部で自動的に呼ばれるため，ユーザーが直接呼ぶ必要はない。**

#### 引数

| 引数 | 種別 | 入出力 | 説明 |
|:-----|:-----|:------:|:-----|
| `A_in(n,n)` | `real(dp)` | 入力 | 係数行列 $A$（内部でコピーされるため破壊されない） |
| `b_in(n)` | `real(dp)` | 入力 | 右辺ベクトル $\mathbf{b}$ |
| `n` | `integer` | 入力 | 次元数 |
| `x(n)` | `real(dp)` | 出力 | 解ベクトル $\mathbf{x}$ |
| `ok` | `logical` | 出力 | 正常に解けたか。ピボット要素 $< 10^{-15}$ のとき `.false.` |

---

## インターフェース要件

`newton_solve` および `newton_solve_nd` に渡す関数・サブルーチンは，
モジュール内で定義された **抽象インターフェース** に適合させる必要がある。

### `func`（1 変数：`newton_solve` 用）

```fortran
! f および df に渡す関数のテンプレート
function my_func(x) result(y)
  use newton_module, only: dp
  real(dp), intent(in) :: x   ! スカラー入力
  real(dp)             :: y   ! スカラー出力
  y = ...
end function my_func
```

### `func_nd`（多変数：`F_func` 用）

```fortran
! F_func に渡すサブルーチンのテンプレート
subroutine my_F(x, F)
  use newton_module, only: dp
  real(dp), intent(in)  :: x(:)   ! 入力ベクトル（長さ n）
  real(dp), intent(out) :: F(:)   ! 出力ベクトル（長さ n）
  F(1) = ...   ! 第 1 式の残差
  F(2) = ...   ! 第 2 式の残差
  ! ...
end subroutine my_F
```

### `jac_nd`（多変数：`J_func` 用）

```fortran
! J_func に渡すサブルーチンのテンプレート
subroutine my_J(x, J)
  use newton_module, only: dp
  real(dp), intent(in)  :: x(:)    ! 入力ベクトル（長さ n）
  real(dp), intent(out) :: J(:,:)  ! ヤコビアン行列（n×n）
  J(1,1) = ...   ! ∂F1/∂x1
  J(1,2) = ...   ! ∂F1/∂x2
  J(2,1) = ...   ! ∂F2/∂x1
  J(2,2) = ...   ! ∂F2/∂x2
  ! ...
end subroutine my_J
```

> **インデックスの規則**：`J(i,j)` は $\partial F_i / \partial x_j$（行 = 方程式番号，列 = 変数番号）。

---

## 使用例

### 1 変数：基本的な使い方

#### 例 1：$x^2 - 2 = 0$（解：$x = \sqrt{2}$）

```fortran
program example1
  use newton_module
  implicit none

  real(dp) :: root
  logical  :: converged
  integer  :: iter
  real(dp), parameter :: tol      = 1.0d-10
  integer,  parameter :: max_iter = 100

  call newton_solve(f, df, 1.5_dp, tol, max_iter, root, converged, iter)

  write(*,'(A,ES20.12)') "root      = ", root       ! 1.414213562373...
  write(*,'(A,ES12.4)')  "error     = ", abs(root - sqrt(2.0_dp))  ! ~1e-15
  write(*,'(A,I4)')      "iter      = ", iter        ! 5
  write(*,'(A,L1)')      "converged = ", converged   ! T

contains

  function f(x) result(y)
    real(dp), intent(in) :: x; real(dp) :: y
    y = x**2 - 2.0_dp
  end function f

  function df(x) result(y)
    real(dp), intent(in) :: x; real(dp) :: y
    y = 2.0_dp * x
  end function df

end program example1
```

#### 例 2：$x^3 - x - 2 = 0$（解：$x \approx 1.5214$）

```fortran
! 関数と導関数の定義
function f(x) result(y)
  real(dp), intent(in) :: x; real(dp) :: y
  y = x**3 - x - 2.0_dp
end function f

function df(x) result(y)
  real(dp), intent(in) :: x; real(dp) :: y
  y = 3.0_dp * x**2 - 1.0_dp
end function df

! 呼び出し
call newton_solve(f, df, 1.5_dp, 1.0d-10, 100, root, converged, iter)
```

#### 例 3：$\cos(x) - x = 0$（解：$x \approx 0.7391$）

```fortran
function f(x) result(y)
  real(dp), intent(in) :: x; real(dp) :: y
  y = cos(x) - x
end function f

function df(x) result(y)
  real(dp), intent(in) :: x; real(dp) :: y
  y = -sin(x) - 1.0_dp
end function df

call newton_solve(f, df, 0.5_dp, 1.0d-10, 100, root, converged, iter)
```

---

### 2 変数の連立方程式

#### 例：$x^2 + y^2 = 4$，$x = y$（解：$(\sqrt{2},\, \sqrt{2})$）

```fortran
program example2d
  use newton_module
  implicit none

  real(dp), parameter :: tol      = 1.0d-10
  integer,  parameter :: max_iter = 100
  integer,  parameter :: n        = 2

  real(dp) :: x0(n), root(n)
  logical  :: converged
  integer  :: iter

  x0 = [1.0_dp, 1.5_dp]   ! 初期ベクトル
  call newton_solve_nd(my_F, my_J, x0, n, tol, max_iter, root, converged, iter)

  write(*,'(A,2ES20.12)') "root = ", root     ! [1.4142..., 1.4142...]
  write(*,'(A,I4)')       "iter = ", iter      ! 4
  write(*,'(A,L1)')       "converged = ", converged  ! T

contains

  ! 残差ベクトル
  subroutine my_F(x, F)
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: F(:)
    F(1) = x(1)**2 + x(2)**2 - 4.0_dp   ! F1 = x^2 + y^2 - 4
    F(2) = x(1) - x(2)                   ! F2 = x - y
  end subroutine my_F

  ! ヤコビアン行列
  ! J = [∂F1/∂x  ∂F1/∂y]   [2x  2y]
  !     [∂F2/∂x  ∂F2/∂y] = [ 1  -1]
  subroutine my_J(x, J)
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: J(:,:)
    J(1,1) = 2.0_dp * x(1);  J(1,2) = 2.0_dp * x(2)
    J(2,1) = 1.0_dp;          J(2,2) = -1.0_dp
  end subroutine my_J

end program example2d
```

---

### 3 変数の連立方程式

#### 例：$x+y+z=6$，$x^2+y^2+z^2=14$，$xyz=6$（解：$(1,\,2,\,3)$）

```fortran
program example3d
  use newton_module
  implicit none

  real(dp), parameter :: tol      = 1.0d-10
  integer,  parameter :: max_iter = 100
  integer,  parameter :: n        = 3

  real(dp) :: x0(n), root(n)
  logical  :: converged
  integer  :: iter

  x0 = [1.5_dp, 2.0_dp, 2.5_dp]
  call newton_solve_nd(my_F, my_J, x0, n, tol, max_iter, root, converged, iter)

  write(*,'(A,3ES20.12)') "root = ", root   ! [1.0, 2.0, 3.0]

contains

  subroutine my_F(x, F)
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: F(:)
    F(1) = x(1) + x(2) + x(3) - 6.0_dp
    F(2) = x(1)**2 + x(2)**2 + x(3)**2 - 14.0_dp
    F(3) = x(1) * x(2) * x(3) - 6.0_dp
  end subroutine my_F

  ! ヤコビアン J(i,j) = ∂F_i/∂x_j
  subroutine my_J(x, J)
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: J(:,:)
    J(1,1) = 1.0_dp;          J(1,2) = 1.0_dp;          J(1,3) = 1.0_dp
    J(2,1) = 2.0_dp * x(1);  J(2,2) = 2.0_dp * x(2);  J(2,3) = 2.0_dp * x(3)
    J(3,1) = x(2) * x(3);    J(3,2) = x(1) * x(3);    J(3,3) = x(1) * x(2)
  end subroutine my_J

end program example3d
```

---

### 次元を増やす方法

ソルバー本体（`newton_solver.f90`）の変更は**不要**。
`my_F`・`my_J` を追加して `newton_solve_nd` に渡すだけでよい。

#### 4 変数への拡張例：$x_i^2 = i$（解：$x_i = \sqrt{i}$）

```fortran
integer, parameter :: n = 4

! 残差ベクトル: F(i) = x(i)^2 - i
subroutine my_F(x, F)
  real(dp), intent(in)  :: x(:)
  real(dp), intent(out) :: F(:)
  integer :: i
  do i = 1, 4
    F(i) = x(i)**2 - real(i, dp)
  end do
end subroutine my_F

! ヤコビアン（対角行列）: J(i,i) = 2*x(i)，その他は 0
subroutine my_J(x, J)
  real(dp), intent(in)  :: x(:)
  real(dp), intent(out) :: J(:,:)
  integer :: i
  J = 0.0_dp
  do i = 1, 4
    J(i,i) = 2.0_dp * x(i)
  end do
end subroutine my_J

! 呼び出し
real(dp) :: x0(n), root(n)
x0 = [2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp]
call newton_solve_nd(my_F, my_J, x0, n, 1.0d-10, 100, root, converged, iter)
! root ≈ [1.0, 1.4142, 1.7321, 2.0]
```

#### n 変数の一般化ルール

| ステップ | 作業内容 |
|:--------|:--------|
| 1 | 変数の個数 `n` を決める |
| 2 | `my_F(x, F)` で `F(1)` 〜 `F(n)` を定義する |
| 3 | `my_J(x, J)` で `J(i,j) = ∂F_i/∂x_j` を計算する（$n^2$ 要素） |
| 4 | `x0(n)` に初期ベクトルを設定する |
| 5 | `call newton_solve_nd(my_F, my_J, x0, n, tol, max_iter, root, converged, iter)` を呼ぶ |

---

## 収束の確認

### `converged` フラグの確認

```fortran
call newton_solve(f, df, x0, tol, max_iter, root, converged, iter)

if (converged) then
  write(*,'(A,ES20.12)') "収束しました: root = ", root
  write(*,'(A,I4)')      "反復回数: ", iter
else
  write(*,'(A)')         "収束しませんでした（最大反復回数に達しました）"
  write(*,'(A,ES20.12)') "最終近似値: root = ", root
end if
```

### 残差の確認

収束後は関数値を評価して残差を直接確認することを推奨する。

```fortran
! 1 変数の場合
write(*,'(A,ES12.4)') "残差 |f(root)|  = ", abs(f(root))

! 多変数の場合
real(dp) :: Fval(n)
call my_F(root, Fval)
write(*,'(A,ES12.4)') "残差 ||F(root)|| = ", norm2(Fval)
```

---

## テスト例題の実行

`test_newton.f90` には以下の 5 つの例題が含まれている。

```bash
gfortran -o test_newton newton_solver.f90 test_newton.f90
./test_newton
```

### 例題一覧

| テスト | 問題 | 初期値 | 真の解 |
|:------:|:-----|:------:|:------:|
| Test 1 (1D) | $x^2 - 2 = 0$ | 1.5 | $\sqrt{2} \approx 1.41421$ |
| Test 2 (1D) | $x^3 - x - 2 = 0$ | 1.5 | $\approx 1.52138$ |
| Test 3 (1D) | $\cos x - x = 0$ | 0.5 | $\approx 0.73909$ |
| Test 4 (2D) | $x^2+y^2=4,\ x=y$ | (1.0, 1.5) | $(\sqrt{2}, \sqrt{2})$ |
| Test 5 (3D) | $x+y+z=6,\ x^2+y^2+z^2=14,\ xyz=6$ | (1.5, 2.0, 2.5) | $(1, 2, 3)$ |

### 実行結果の例（Test 1）

```
=== Test 1 (1D): x^2 - 2 = 0 ===
  x0 =  1.50000
  iter=   1  x=  1.416666666667E+00  f(x)=  6.3611E-02
  iter=   2  x=  1.414215686275E+00  f(x)=  6.0073E-06
  iter=   3  x=  1.414213562375E+00  f(x)=  5.2915E-12
  iter=   4  x=  1.414213562373E+00  f(x)= -4.4409E-16
  iter=   5  x=  1.414213562373E+00  f(x)= -4.4409E-16
  root      =   1.414213562373E+00
  exact     =   1.414213562373E+00
  error     =   4.4409E-16
  iter      =    5
  [Result] Converged.
```

---

## ヒントとトラブルシューティング

### ヤコビアン行列の書き方

$F = (F_1, F_2)^T$，$\mathbf{x} = (x_1, x_2)^T$ の場合：

$$J = \begin{pmatrix} \partial F_1/\partial x_1 & \partial F_1/\partial x_2 \\ \partial F_2/\partial x_1 & \partial F_2/\partial x_2 \end{pmatrix}$$

```fortran
! 対応する Fortran コード
J(1,1) = 偏微分 ∂F1/∂x1  ! 第1方程式を x1 で微分
J(1,2) = 偏微分 ∂F1/∂x2  ! 第1方程式を x2 で微分
J(2,1) = 偏微分 ∂F2/∂x1  ! 第2方程式を x1 で微分
J(2,2) = 偏微分 ∂F2/∂x2  ! 第2方程式を x2 で微分
```

**よくある間違い**：`J(i,j)` と `J(j,i)` の転置ミス。
ヤコビアンが正しければ数回で収束するが，転置していると収束が遅くなるか発散する。

### 倍精度定数の書き方

精度を落とさないため，定数には必ず `_dp` または `d` 表記を付ける。

```fortran
! 良い例（倍精度）
y = x**2 - 2.0_dp
y = 3.14159265358979_dp * x

! 悪い例（単精度定数が混入する）
y = x**2 - 2.0    ! 2.0 は単精度リテラル
y = 3.14 * x      ! π が単精度精度になる
```

### 収束しない場合の対処

| 症状 | 原因の候補 | 対処 |
|:-----|:----------|:-----|
| `converged = .false.` で終了 | 初期値が解から遠い | `x0` を真の解に近い値に変更する |
| `[Warning] df/dx is nearly zero` | 初期値が極値付近 | `x0` をずらす |
| `[Warning] Jacobian is nearly singular` | 解のヤコビアンが特異（または初期値が悪い） | `x0` を変更する，問題を確認する |
| 発散・大きな値になる | 初期値と解が遠くオーバーシュート | `x0` を解に近い値にするか，ステップ幅を手動で小さくする |

### 反復回数・許容誤差の調整

```fortran
! 収束しにくい場合：反復回数を増やす
integer, parameter :: max_iter = 500

! より高精度が必要な場合：tol を小さくする
real(dp), parameter :: tol = 1.0d-14

! 収束しやすくしたい場合：tol を大きくする
real(dp), parameter :: tol = 1.0d-6
```

> デフォルト設定 `tol = 1.0d-10`，`max_iter = 100` は，倍精度の範囲で通常の問題に対し
> 十分な精度と反復余裕を持つ値である。

### `use newton_module` の記述場所

`newton_module` は `newton_solver.f90` で定義されており，
これを利用するすべてのプログラム単位の先頭で `use newton_module` を記述する。

```fortran
program my_program
  use newton_module          ! ← 必ずここで読み込む
  implicit none
  ...
end program my_program
```

`contains` 内のサブルーチン・関数からは，ホストの `use` 文を引き継ぐため，
個別に `use newton_module` を書く必要はない（ホスト結合）。
ただし独立したモジュールや外部サブルーチンでは個別に `use` が必要。
