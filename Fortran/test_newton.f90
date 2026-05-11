program test_newton
  use newton_module
  implicit none

  real(dp) :: root
  logical  :: converged
  integer  :: iter
  real(dp), parameter :: tol      = 1.0d-10
  integer,  parameter :: max_iter = 100

  ! ---- テスト1: x^2 - 2 = 0  (真の解: sqrt(2) ≈ 1.41421356...) ----
  write(*,'(/,A)') "=== テスト1: x^2 - 2 = 0 ==="
  write(*,'(A,F8.5)') "  初期値 x0 = ", 1.5_dp
  call newton_solve(f1, df1, 1.5_dp, tol, max_iter, root, converged, iter)
  call print_result(root, sqrt(2.0_dp), converged, iter)

  ! ---- テスト2: x^3 - x - 2 = 0  (真の解: ≈ 1.52137970680457) ----
  write(*,'(/,A)') "=== テスト2: x^3 - x - 2 = 0 ==="
  write(*,'(A,F8.5)') "  初期値 x0 = ", 1.5_dp
  call newton_solve(f2, df2, 1.5_dp, tol, max_iter, root, converged, iter)
  call print_result(root, 1.5213797068045675_dp, converged, iter)

  ! ---- テスト3: cos(x) - x = 0  (真の解: ≈ 0.73908513321516) ----
  write(*,'(/,A)') "=== テスト3: cos(x) - x = 0 ==="
  write(*,'(A,F8.5)') "  初期値 x0 = ", 0.5_dp
  call newton_solve(f3, df3, 0.5_dp, tol, max_iter, root, converged, iter)
  call print_result(root, 0.7390851332151607_dp, converged, iter)

contains

  ! ---------- テスト1: f(x) = x^2 - 2 ----------
  function f1(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = x**2 - 2.0_dp
  end function f1

  function df1(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = 2.0_dp * x
  end function df1

  ! ---------- テスト2: f(x) = x^3 - x - 2 ----------
  function f2(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = x**3 - x - 2.0_dp
  end function f2

  function df2(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = 3.0_dp * x**2 - 1.0_dp
  end function df2

  ! ---------- テスト3: f(x) = cos(x) - x ----------
  function f3(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = cos(x) - x
  end function f3

  function df3(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = -sin(x) - 1.0_dp
  end function df3

  ! ---------- 結果の表示 ----------
  subroutine print_result(root, exact, converged, iter)
    real(dp), intent(in) :: root, exact
    logical,  intent(in) :: converged
    integer,  intent(in) :: iter
    write(*,'(A,ES20.12)') "  計算した根  = ", root
    write(*,'(A,ES20.12)') "  真の解      = ", exact
    write(*,'(A,ES12.4)')  "  誤差        = ", abs(root - exact)
    write(*,'(A,I4)')      "  反復回数    = ", iter
    if (converged) then
      write(*,'(A)') "  [結果] 収束しました"
    else
      write(*,'(A)') "  [結果] 収束しませんでした"
    end if
  end subroutine print_result

end program test_newton
