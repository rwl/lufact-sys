c  Added error parameter, colunm perm. Sivan Toledo 2/97.

c
c  lsolve :  solve lower triangular system.
c
      subroutine lsolve (n, lu, lurow, lcolst, ucolst, rperm, cperm,
     $                   b, x, error)
c
c  This routine takes an LU factorization from lufact (i.e. P, L, U with
c  PA = LU) and solves Lx = Pb for x.  There is nothing clever at all
c  about sparse right-hand sides here; we always look at every nonzero
c  of L.  We do make some checks for consistency of the LU data
c  structure.
c
c  Input parameters:
c    n    Dimension of the system.
c    lu, lurow, lcolst, ucolst, rperm, cperm  LU factorization
c    b    Right-hand side, as a dense n-vector.
c
c  Output parameter:
c    x    Solution, as a dense n-vector.
c    error 0 if successful, 1 otherwise
c
      integer   n, lurow(1), lcolst(1), ucolst(1), rperm(1), cperm(1)
      integer   i, j, nzst, nzend, nzptr
      integer   error
      real*8    lu(1), b(1), x(1)
c
      if (n .le. 0) goto 800
c
c  Check that rperm is really a permutation.
c
      do 10 i = 1, n
          x(i) = 0.0
10        continue
      do 20 i = 1, n
          if (rperm(i) .lt. 1  .or.  rperm(i) .gt. n) goto 803
          if (x(rperm(i)) .ne. 0.0) goto 803
          x(rperm(i)) = 1.0
20        continue
c
c  Check that cperm is really a permutation.
c
      do 110 i = 1, n
          x(i) = 0.0
110       continue
      do 120 i = 1, n
          if (cperm(i) .lt. 1  .or.  cperm(i) .gt. n) goto 804
          if (x(cperm(i)) .ne. 0.0) goto 804
          x(cperm(i)) = 1.0
120        continue
c
c  Solve the system.
c
      do 50 i = 1, n
          x(rperm(i)) = b(i)
50        continue
c
      do 200 j = 1, n
          nzst =  lcolst(j)
          nzend = ucolst(j+1) - 1
          if (nzst .lt. 1  .or.  nzst .gt. nzend+1) goto 801
          if (nzst .gt. nzend) goto 150
          do 100 nzptr = nzst, nzend
              i = lurow(nzptr)
              if (i .le. j .or. i .gt. n) goto 802
              x(i) = x(i) - lu(nzptr)*x(j)
100           continue
150       continue
200       continue

      error = 0
      return
c
c  Bug traps.
c
800   continue
      write (6, 900) n
900   format ('lsolve called with nonpositive n =', i6)
      error = 1
      return
c
801   continue
      write (6, 901) j, nzst, nzend
901   format ('In lsolve, inconsistent column of L: j, nzst, nzend =',
     1        3i6)
      error = 1
      return
c
802   continue
      write (6, 902) i, j, nzptr
902   format ('In lsolve, illegal row i in column j of L:'/
     1        '  i, j, nzptr =', 3i6)
      error = 1
      return
c
803   continue
      print *,rperm(i)
      write (6, 903) i
903   format ('In lsolve, rpermutation is illegal in position i =', i6)
      error = 1
      return
c
804   continue
      print *,cperm(i)
      write (6, 904) i
904   format ('In lsolve, cpermutation is illegal in position i =', i6)
      error = 1
      return
c
      end
