c  In the following routine for copying whole arrays, the direction
c  of iteration (which makes a difference if the arrays overlap) is
c  controlled by MODE, which is set false for backward displacement and
c  true for forward displacement.
c
      subroutine icopy(a, b, la, mode)
c
c  ICOPY copies an integer array A to another array B.
c
      integer   a(1), b(1), la, i
      logical   mode
c
      if (mode) goto 200
          do 100 i = 1, la
              b(i) = a(i)
100       continue
          return
c
200   continue
          do 300 i = la, 1, -1
              b(i) = a(i)
300       continue
          return
      end
