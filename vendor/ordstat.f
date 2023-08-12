      subroutine dordstat (n, k, A, kth, info)

      implicit none

      integer          n, k
      double precision A(n)
      double precision kth
      integer          info

      integer          p,r,i,j,q
      double precision x,tmp

      integer          rnd
      save             rnd
      data             rnd / 1 /

      if (k .lt. 0 .or. k .gt. n) then
         info = -1
         return
      end if

      p  = 1
      r  = n

 100  continue

      if (p .eq. r) goto 900

      if (r - p .ge. 8) then
         rnd = mod(1366 * rnd + 150889,714025)
         q = p + mod(rnd,r-p+1)

         tmp  = A(p)
         A(p) = A(q)
         A(q) = tmp
      end if

      x = A(p)
      i = p-1
      j = r+1

 200  continue
      
 210  continue
      j = j-1
      if (A(j) .gt. x) goto 210

 220  continue
      i = i+1
      if (A(i) .lt. x) goto 220

      if (i .lt. j) then
         tmp  = A(i)
         A(i) = A(j)
         A(j) = tmp
         goto 200
      end if

      if (j .lt. k) then 
         p = j+1
      else
         r = j
      end if

      goto 100

 900  continue

      kth = A(p)

      info = 0

      return

      end
