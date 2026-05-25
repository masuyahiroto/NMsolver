module newton_module
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  abstract interface
    function func(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function func

    subroutine func_nd(x, F)
      import dp
      real(dp), intent(in)  :: x(:)
      real(dp), intent(out) :: F(:)
    end subroutine func_nd

    subroutine jac_nd(x, J)
      import dp
      real(dp), intent(in)  :: x(:)
      real(dp), intent(out) :: J(:,:)
    end subroutine jac_nd
  end interface

contains

  subroutine newton_solve(f, df, x0, tol, max_iter, root, converged, iter)
    procedure(func) :: f, df
    real(dp), intent(in)  :: x0, tol
    integer,  intent(in)  :: max_iter
    real(dp), intent(out) :: root
    logical,  intent(out) :: converged
    integer,  intent(out) :: iter

    real(dp) :: x, fx, dfx, dx
    integer  :: i

    x         = x0
    converged = .false.
    iter      = 0

    do i = 1, max_iter
      fx  = f(x)
      dfx = df(x)

      if (abs(dfx) < 1.0d-15) then
        write(*,'(A)') "  [Warning] df/dx is nearly zero. Stopping."
        exit
      end if

      dx = -fx / dfx
      x  = x + dx
      iter = i

      write(*,'(A,I4,A,ES20.12,A,ES12.4)') &
        "  iter=", i, "  x=", x, "  f(x)=", f(x)

      if (abs(fx) < tol .and. abs(dx) < tol) then
        converged = .true.
        exit
      end if
    end do

    root = x
  end subroutine newton_solve

  subroutine newton_solve_nd(F_func, J_func, x0, n, tol, max_iter, root, converged, iter)
    procedure(func_nd) :: F_func
    procedure(jac_nd)  :: J_func
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: x0(n), tol
    integer,  intent(in)  :: max_iter
    real(dp), intent(out) :: root(n)
    logical,  intent(out) :: converged
    integer,  intent(out) :: iter

    real(dp) :: x(n), Fx(n), Jx(n,n), dx(n)
    real(dp) :: norm_Fx, norm_dx
    integer  :: i
    logical  :: ok

    x         = x0
    converged = .false.
    iter      = 0

    do i = 1, max_iter
      call F_func(x, Fx)
      call J_func(x, Jx)

      call lu_solve(Jx, -Fx, n, dx, ok)
      if (.not. ok) then
        write(*,'(A)') "  [Warning] Jacobian is nearly singular. Stopping."
        exit
      end if

      norm_Fx = norm2(Fx)
      norm_dx = norm2(dx)

      x    = x + dx
      iter = i

      write(*,'(A,I4,A,ES12.4,A,ES12.4)') &
        "  iter=", i, "  ||F(x)||=", norm_Fx, "  ||dx||=", norm_dx

      if (norm_Fx < tol .and. norm_dx < tol) then
        converged = .true.
        exit
      end if
    end do

    root = x
  end subroutine newton_solve_nd

  subroutine lu_solve(A_in, b_in, n, x, ok)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: A_in(n,n), b_in(n)
    real(dp), intent(out) :: x(n)
    logical,  intent(out) :: ok

    real(dp) :: A(n,n), b(n), factor, tmp, tmp_row(n)
    integer  :: i, k, prow, pivot_idx(1)
    real(dp), parameter :: eps = 1.0d-15

    A  = A_in
    b  = b_in
    ok = .true.

    do k = 1, n-1
      pivot_idx = maxloc(abs(A(k:n, k)))
      prow = k + pivot_idx(1) - 1

      if (abs(A(prow, k)) < eps) then
        ok = .false.
        return
      end if

      if (prow /= k) then
        tmp_row  = A(k,:);  A(k,:)  = A(prow,:);  A(prow,:) = tmp_row
        tmp      = b(k);    b(k)    = b(prow);     b(prow)   = tmp
      end if

      do i = k+1, n
        factor      = A(i,k) / A(k,k)
        A(i,k+1:n) = A(i,k+1:n) - factor * A(k,k+1:n)
        b(i)        = b(i) - factor * b(k)
      end do
    end do

    if (abs(A(n,n)) < eps) then
      ok = .false.
      return
    end if

    x(n) = b(n) / A(n,n)
    do i = n-1, 1, -1
      x(i) = (b(i) - dot_product(A(i, i+1:n), x(i+1:n))) / A(i,i)
    end do
  end subroutine lu_solve

end module newton_module
