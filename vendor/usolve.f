c  Added error parameter, permutations. Sivan Toledo 2/97.

c
c  usolve :  solve upper triangular system.
c
      subroutine usolve (n, lu, lurow, lcolst, ucolst, rperm, cperm,
     $                   b, x, error)
c
c  This routine takes an LU factorization from lufact (i.e. L, U
c  with PA = LU) and solves Ux = b for x.  Note that P is not used
c  and is not a parameter.  There is nothing clever at all about
c  sparse right-hand sides here; we always look at every nonzero of U.
c  We do make some checks for consistency of the LU data structure.
c
c  Input parameters:
c    n    Dimension of the system.
c    lu, lurow, lcolst, ucolst  LU factorization; see lufact for format.
c    b    Right-hand side, as a dense n-vector.
c
c  Output parameter:
c    x    Solution, as a dense n-vector.
c    error 0 if successful, 1 otherwise
c
      integer   n, lurow(1), lcolst(1), ucolst(1), rperm(1), cperm(1)
      integer   i, j, nzst, nzend, nzptr, jj
      integer   error
      real*8    lu(1), b(1), x(1)
c
      if (n .le. 0) goto 800
      do 50 i = 1, n
          x(i) = b(i)
50        continue
c
      do 200 jj = 1, n
          j = n + 1 - jj
          nzst =  ucolst(j)
          nzend = lcolst(j) - 1
          if (nzst .lt. 1  .or.  nzst .gt. nzend) goto 801
          if (lurow(nzend) .ne. j) goto 803
          if (lu(nzend) .eq. 0.0) goto 804
          x(j) = x(j) / lu(nzend)
          nzend = nzend - 1
          if (nzst .gt. nzend) goto 150
          do 100 nzptr = nzst, nzend
              i = lurow(nzptr)
              if (i .le. 0 .or. i .ge. j) goto 802
              x(i) = x(i) - lu(nzptr)*x(j)
100           continue
150       continue
200       continue

      do 250 i = 1, n
          b(i) = x(i)
250       continue
      do 260 i = 1, n
          x(cperm(i)) = b(i)
260       continue

      error = 0
      return
c
c  Bug traps.
c
800   continue
      write (6, 900) n
900   format ('usolve called with nonpositive n =', i6)
      error = 1
      return
c
801   continue
      write (6, 901) j, nzst, nzend
901   format ('In usolve, inconsistent column of U: j, nzst, nzend =',
     1        3i6)
      error = 1
      return
c
802   continue
      write (6, 902) i, j, nzptr
902   format ('In usolve, illegal row i in column j of U:'/
     1        '  i, j, nzptr =', 3i6)
      error = 1
      return
c
803   continue
      write (6, 903) j, nzend, lurow(nzend)
903   format ('In usolve, diagonal elt of col j is not in last place.'/
     1       'j, nzend, lurow(nzend) =', 3i6)
      error = 1
      return
c
804   continue
      write (6, 904) j
904   format ('In usolve, zero diagonal element in column j =', i6)
      error = 1
      return
      end
