c
c  lusolv :  solve a square linear system, given an LU factorization.
c
      subroutine lusolv (n, lu, lurow, lcolst, ucolst, perm, x, rwork)
c
c  Solve for X in a square linear system Ax = b, given the factorization
c  PA=LU.
c
c  Input parameters:
c    n                          dimension of matrix.
c    lu, lurow, lcolst,
c    ucolst, perm               PA=LU factorization (see lufact for format).
c
c  Modified parameter:
c    x                          Real array of length n.
c                               On entry, holds B.  On exit, holds X.
c
c  Work parameter:
c    rwork                      Real array of length n; holds intermediate
c                               solution.
c
      integer   n, lurow(1), lcolst(1), ucolst(1), perm(1)
      real*8    lu(1), x(1), rwork(1)
c
c
      call lsolve (n, lu, lurow, lcolst, ucolst, perm, x, rwork)
      call usolve (n, lu, lurow, lcolst, ucolst, rwork, x)
      return
      end
