program test_newton
  use newton_module
  implicit none

  real(dp) :: root
  logical  :: converged
  integer  :: iter
  real(dp), parameter :: tol      = 1.0d-10
  integer,  parameter :: max_iter = 100

  real(dp) :: x0_2d(2), root_2d(2)
  real(dp) :: x0_3d(3), root_3d(3)

  ! ---- Test 1: x^2 - 2 = 0  (exact: sqrt(2) = 1.41421356...) ----
  write(*,'(/,A)') "=== Test 1 (1D): x^2 - 2 = 0 ==="
  write(*,'(A,F8.5)') "  x0 = ", 1.5_dp
  call newton_solve(f1, df1, 1.5_dp, tol, max_iter, root, converged, iter)
  call print_result(root, sqrt(2.0_dp), converged, iter)

  ! ---- Test 2: x^3 - x - 2 = 0  (exact: 1.52137970680457) ----
  write(*,'(/,A)') "=== Test 2 (1D): x^3 - x - 2 = 0 ==="
  write(*,'(A,F8.5)') "  x0 = ", 1.5_dp
  call newton_solve(f2, df2, 1.5_dp, tol, max_iter, root, converged, iter)
  call print_result(root, 1.5213797068045675_dp, converged, iter)

  ! ---- Test 3: cos(x) - x = 0  (exact: 0.73908513321516) ----
  write(*,'(/,A)') "=== Test 3 (1D): cos(x) - x = 0 ==="
  write(*,'(A,F8.5)') "  x0 = ", 0.5_dp
  call newton_solve(f3, df3, 0.5_dp, tol, max_iter, root, converged, iter)
  call print_result(root, 0.7390851332151607_dp, converged, iter)

  ! ---- Test 4 (2D): x^2 + y^2 = 4, x = y  (exact: (sqrt(2), sqrt(2))) ----
  write(*,'(/,A)') "=== Test 4 (2D): x^2 + y^2 = 4,  x = y ==="
  x0_2d = [1.0_dp, 1.5_dp]
  write(*,'(A,2F8.5)') "  x0 = ", x0_2d
  call newton_solve_nd(F2d, J2d, x0_2d, 2, tol, max_iter, root_2d, converged, iter)
  call print_result_nd(root_2d, [sqrt(2.0_dp), sqrt(2.0_dp)], converged, iter)

  ! ---- Test 5 (3D): x+y+z=6, x^2+y^2+z^2=14, xyz=6  (exact: (1,2,3)) ----
  write(*,'(/,A)') "=== Test 5 (3D): x+y+z=6,  x^2+y^2+z^2=14,  xyz=6 ==="
  x0_3d = [1.5_dp, 2.0_dp, 2.5_dp]
  write(*,'(A,3F8.5)') "  x0 = ", x0_3d
  call newton_solve_nd(F3d, J3d, x0_3d, 3, tol, max_iter, root_3d, converged, iter)
  call print_result_nd(root_3d, [1.0_dp, 2.0_dp, 3.0_dp], converged, iter)

contains

  ! ---------- Test 1: f(x) = x^2 - 2 ----------
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

  ! ---------- Test 2: f(x) = x^3 - x - 2 ----------
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

  ! ---------- Test 3: f(x) = cos(x) - x ----------
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

  ! ---------- Test 4 (2D): F1 = x^2 + y^2 - 4,  F2 = x - y ----------
  subroutine F2d(x, F)
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: F(:)
    F(1) = x(1)**2 + x(2)**2 - 4.0_dp
    F(2) = x(1) - x(2)
  end subroutine F2d

  subroutine J2d(x, J)
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: J(:,:)
    J(1,1) = 2.0_dp * x(1);  J(1,2) = 2.0_dp * x(2)
    J(2,1) = 1.0_dp;         J(2,2) = -1.0_dp
  end subroutine J2d

  ! ---------- Test 5 (3D): F1=x+y+z-6, F2=x^2+y^2+z^2-14, F3=xyz-6 ----------
  subroutine F3d(x, F)
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: F(:)
    F(1) = x(1) + x(2) + x(3) - 6.0_dp
    F(2) = x(1)**2 + x(2)**2 + x(3)**2 - 14.0_dp
    F(3) = x(1) * x(2) * x(3) - 6.0_dp
  end subroutine F3d

  subroutine J3d(x, J)
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: J(:,:)
    J(1,1) = 1.0_dp;          J(1,2) = 1.0_dp;          J(1,3) = 1.0_dp
    J(2,1) = 2.0_dp * x(1);  J(2,2) = 2.0_dp * x(2);  J(2,3) = 2.0_dp * x(3)
    J(3,1) = x(2) * x(3);    J(3,2) = x(1) * x(3);    J(3,3) = x(1) * x(2)
  end subroutine J3d

  ! ---------- print result (1D) ----------
  subroutine print_result(root, exact, converged, iter)
    real(dp), intent(in) :: root, exact
    logical,  intent(in) :: converged
    integer,  intent(in) :: iter
    write(*,'(A,ES20.12)') "  root      = ", root
    write(*,'(A,ES20.12)') "  exact     = ", exact
    write(*,'(A,ES12.4)')  "  error     = ", abs(root - exact)
    write(*,'(A,I4)')      "  iter      = ", iter
    if (converged) then
      write(*,'(A)') "  [Result] Converged."
    else
      write(*,'(A)') "  [Result] Did not converge."
    end if
  end subroutine print_result

  ! ---------- print result (ND) ----------
  subroutine print_result_nd(root, exact, converged, iter)
    real(dp), intent(in) :: root(:), exact(:)
    logical,  intent(in) :: converged
    integer,  intent(in) :: iter
    integer :: k
    do k = 1, size(root)
      write(*,'(A,I1,A,ES20.12,A,ES20.12,A,ES12.4)') &
        "  x(", k, ")    = ", root(k), &
        "  exact = ", exact(k), &
        "  error = ", abs(root(k) - exact(k))
    end do
    write(*,'(A,ES12.4)') "  ||error|| = ", norm2(root - exact)
    write(*,'(A,I4)')     "  iter      = ", iter
    if (converged) then
      write(*,'(A)') "  [Result] Converged."
    else
      write(*,'(A)') "  [Result] Did not converge."
    end if
  end subroutine print_result_nd

end program test_newton
