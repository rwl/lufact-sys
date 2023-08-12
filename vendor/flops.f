c
c  flops :  count flops in an LU factorization
c
      real*8 function flops (nrow, ncol, lurow, lcolst, ucolst, ucount)
c
c  This counts the number of nonzero arithmetic operations that would
c  be needed to multiply L by U, which is the number of operations the
c  LUFACT code performed to get the factorization in the first place
c  (assuming that no coincidental zeros arose).  Notice that we count
c  both multiplications and additions.
c
c
c  Input parameters:
c    ncol                       columns in matrix.
c    lurow, lcolst, ucolst      structure of factorization 
c                               (see lufact for format)
c
c  Working parameter:
c    ucount(nrow)               used to count nonzeros in each row of U.
c
      integer   nrow, ncol, lurow(1), lcolst(1), ucolst(1)
      integer   ucount(1)
      integer   update, lcount, i, j, k
      real*8    sum
c
c
c  Count nonzeros in each row of U
c
      do 70 i = 1, nrow
          ucount(i) = 0
70        continue
      do 90 j = 1, ncol
          do 80 k = ucolst(j), lcolst(j)-1
              ucount(lurow(k)) = ucount(lurow(k)) + 1
80            continue
90        continue
c
c  Use nonzero counts for rows of U and columns of L to count flops.
c
      sum = 0.0d0
      do 100 k = 1, ncol
          lcount = ucolst(k+1) - lcolst(k)
          update = 2 * lcount * ucount(k) - lcount
          sum = sum + update
100       continue
      flops = sum
      return
      end
