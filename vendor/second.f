c  Timing subroutine for use with Linpack benchmarks distributed 
c  by netlib.  The call t = second () returns elapsed running time
c  since something-or-other, in seconds.
c
      real function second ()
      real etime
      real tarray(2)
c
      second = etime (tarray)
      return
      end
