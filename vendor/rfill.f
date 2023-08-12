      subroutine rfill(a, la, rval)
c
c  RFILL fills a real*8 array with a given value.
c
      integer   la, i
      real*8    a(1), rval
c
      do 100 i = 1, la
          a(i) = rval
100   continue
      return
      end
