      logical function requiv(a, b)
c
c  REQUIV tests if two real*8 arrays start at the same address.
c
      real*8 a(1), b(1), temp
c
      requiv = .false.
      temp = a(1)
      a(1) = 0.0
      if (b(1) .ne. 0.0) goto 100
          a(1) = 1.0
          if (b(1) .ne. 1.0) goto 100
              requiv = .true.
100   continue
      a(1) = temp
      return
      end
