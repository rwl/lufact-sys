      subroutine ifill(a, la, ival)
c
c  IFILL fills an integer array with a given value.
c
      integer a(1), la, ival, i
c
      do 100 i = 1, la
          a(i) = ival
100   continue
      return
      end
