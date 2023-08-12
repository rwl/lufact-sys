      subroutine cntrow(arow, lasta, rowcnt)
c
c  CNTROW fills its last argument with the nonzero row counts of the
c  matrix specified in the first two arguments.  (See lufact for matrix
c  format.)
c
      integer   arow(1), lasta, rowcnt(1)
      integer   i, j, k, maxk
c
c  MAXK marks the highest numbered row that has been seen.
c
      maxk = 0
      do 300 i = 1, lasta
          k = arow(i)
          if (k .le. maxk) goto 200
              do 100 j = maxk+1, k
                  rowcnt(j) = 0
100           continue
              maxk = k
200       continue
          rowcnt(k) = rowcnt(k) + 1
300   continue
      return
      end
