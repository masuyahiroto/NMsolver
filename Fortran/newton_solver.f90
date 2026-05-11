module newton_module
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  abstract interface
    function func(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function func
  end interface

contains

  ! ニュートン法による方程式 f(x)=0 の根を求める
  ! f       : 対象関数
  ! df      : f の導関数
  ! x0      : 初期値
  ! tol     : 収束判定閾値
  ! max_iter: 最大反復回数
  ! root    : 求めた根
  ! converged: 収束したかどうか
  ! iter    : 実際の反復回数
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
        write(*,'(A)') "  [警告] 微分値がゼロに近いため停止"
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

  

end module newton_module
