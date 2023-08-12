c Changes zpivot from logical to integer. Sivan Toledo 2/97.
c Cut out rowcnt and the corresponding code.
c Added flops

c
c  lucopy :  Copy dense column to sparse structure, pivot, and divide
c
      subroutine lucopy (
     $     pivot, pthresh, dthresh, nzcount,
     $     jcol,  ncol, 
     $     lastlu, lu, lurow, lcolst, ucolst, 
     $     rperm, cperm,
     $     dense, 
     $     pattern, 
     $     twork,
     $     flops,
     $     zpivot)
c
c  Copy column jcol from the dense vector to the sparse data structure,
c  zeroing out the dense vector.  Then find the diagonal element (either
c  by partial or threshold pivoting or by looking for row number=col
c  number), move it from L into U, and divide the column of L by it.
c
c  Input variables:
c    pivot   = -1 for columns with no diagonal element
c            = 0 for no pivoting
c            = 1 for partial (row) pivoting
c            = 2 for threshold (row) pivoting
c    pthresh  fraction of max pivot candidate acceptable for pivoting
c    jcol    Current column number.
c    ncol    Total number of columns; upper bound on row counts.
c
c  Modified variables:
c    lastlu                 Index of last nonzero in lu, updated here.
c    lu                     On entry, cols 1 through jcol-1 of Pt(L-I+U).
c                           On exit, cols 1 through jcol.
c    lurow, lcolst, ucolst  Nonzero structure of columns 1 through jcol
c                           of PtL and PtU.
c                           No pivoting has been done, so on entry the
c                           element that will be U(jcol,jcol) is still
c                           somewhere in L; on exit, it is the last
c                           nonzero in column jcol of U.
c    rperm                  The row permutation P.
c                           rperm(r) = s > 0 means row r of A is row s of PA.
c                           rperm(r) = 0 means row r of A has not yet been used
c                           as a pivot.  On input, perm reflects rows 1 through
c                           jcol-1 of PA; on output, rows 1 through jcol.
c    cperm                  The column permutation.
c    dense                  On entry, column jcol of Pt(U(jcol,jcol)*(L-I)+U).
c                           On exit, zero.
c    flops                  flop count
c
c  Output variable:
c    zpivot                 > 0 for success (pivot row), -1 for zero pivot element.
c  (was:
c    zpivot                 .false. for success, .true. for zero pivot element.
c
c  Local variables:
c    nzptr       Index into lurow of current nonzero.
c    nzst, nzend Loop bounds for nzptr.
c    irow        Row number of current nonzero (according to A, not PA).
c    pivrow      Pivot row number (according to A, not PA).
c    maxpiv      Temporary to find maximum element in column for pivoting.
c    utemp       Temporary for computing maxpiv.
c    ujj         Diagonal element U(jcol,jcol) = PtU(pivrow,jcol).
c    ujjptr      Index into lu and lurow of diagonal element.
c    dptr        Temporary index into lu and lurow.
c    diagptr     Index to diagonal element of QAQt
c    diagpiv     Value of diagonal element
c
      real*8    lu(1), dense(1), twork(1)
      real*8    pthresh, ujj, utemp, maxpiv, flops
      real*8    maxpivglb, dthresh, dthreshabs, kth
      real*8    udthreshabs, ldthreshabs
      integer   pivot, jcol, ncol, lurow(1), lcolst(1), ucolst(1)
      integer   rperm(1), cperm(1), pattern(1)
      integer   lastlu, nzcount
      integer   nzptr, nzst, nzend, irow, ujjptr, dptr, pivrow
      integer   nzcpy
      integer   zpivot
      integer   lastu
      integer   i

      real*8    diagpiv
      integer   diagptr

c
c
c     Copy column jcol from dense to sparse, recording the position of
c     the diagonal element.
c
      ujjptr = 0

      if (pivot .le. 0) then

c
c     No pivoting, diagonal element has irow = jcol.  Copy the
c     column elements of U and L, throwing out zeros.
c     
         if (ucolst(jcol+1)-1 .lt. ucolst(jcol)) goto 9001
         
c     Start with U
         
         nzcpy = ucolst(jcol)
         
         do nzptr = ucolst(jcol),lcolst(jcol)-1
            irow = lurow(nzptr)
            
            if (pattern(irow) .ne. 0 .or.
     $          irow .eq. cperm(jcol)) then
               lurow(nzcpy) = irow
               lu(nzcpy) = dense(irow)
               dense(irow) = 0.0
               nzcpy = nzcpy + 1
            else
               dense(irow) = 0.0
            end if
         end do
         lastu = nzcpy - 1
         
c     Now do L. Same action as U, except that we search for diagonal
         
         do nzptr = lcolst(jcol),ucolst(jcol+1)-1
            irow = lurow(nzptr)
c            if (irow .eq. cperm(jcol)) then
            if (pattern(irow) .eq. 2) then
               ujjptr = nzcpy
            end if
            if (pattern(irow) .ne. 0 .or. 
c     $           irow .eq. cperm(jcol)) then
     $          pattern(irow) .eq. 2) then
               lurow(nzcpy) = irow
               lu(nzcpy) = dense(irow)
               dense(irow) = 0.0
               nzcpy        = nzcpy + 1
            else
               dense(irow) = 0.0
            end if
         end do
         
         lcolst(jcol)   = lastu + 1
         ucolst(jcol+1) = nzcpy
         lastlu         = nzcpy - 1
         
         if (pivot .eq. -1) goto 8000
         
      else
c
c     Partial and threshold pivoting.  

         if (ucolst(jcol+1) - 1 .lt. lcolst(jcol)) goto 9002
c     
c             Partial pivoting, diagonal elt. has max. magnitude in L.
c

c     Compute the drop threshold for the column

         if (nzcount .le. 0) then

            maxpivglb = -1.0D0
            do nzptr = ucolst(jcol), lcolst(jcol) - 1
               irow = lurow(nzptr)
               utemp = abs(dense(irow))
               if (utemp .gt. maxpivglb) then
                  maxpivglb = utemp
               end if
            end do
            udthreshabs = dthresh * maxpivglb

            maxpivglb = -1.0D0
            do nzptr = lcolst(jcol), ucolst(jcol+1) - 1
               irow = lurow(nzptr)
               utemp = abs(dense(irow))
               if (utemp .gt. maxpivglb) then
                  maxpivglb = utemp
               end if
            end do
            ldthreshabs = dthresh * maxpivglb

         else

            i = 0
            do nzptr = ucolst(jcol), lcolst(jcol) - 1
               i = i + 1
               irow = lurow(nzptr)
               utemp = abs(dense(irow))
               twork(i) = utemp
            end do
            if (nzcount .lt. i) then
               call dordstat(i,i-nzcount+1,twork,kth,i)
               udthreshabs = kth
            else
               udthreshabs = 0.0D0
            end if

            i = 0
            do nzptr = lcolst(jcol), ucolst(jcol+1) - 1
               i = i + 1
               irow = lurow(nzptr)
               utemp = abs(dense(irow))
               twork(i) = utemp
            end do
            if (nzcount .lt. i) then
               call dordstat(i,i-nzcount+1,twork,kth,i)
               ldthreshabs = kth
            else
               ldthreshabs = 0.0D0
            end if

         end if

c     Copy the column elements of
c     U, throwing out zeros.
c

         nzcpy = ucolst(jcol)

         if (lcolst(jcol) - 1 .ge. ucolst(jcol)) then
            do nzptr = ucolst(jcol), lcolst(jcol) - 1
               irow = lurow(nzptr)

c               if (pattern(irow) .ne. 0 .or. 
c     $          pattern(irow) .eq. 2) then
               if (pattern(irow) .ne. 0 .or. 
     $             abs(dense(irow)) .ge. udthreshabs) then
                  lurow(nzcpy) = irow
                  lu(nzcpy) = dense(irow)
                  dense(irow) = 0.0
                  nzcpy = nzcpy + 1
               else
                  dense(irow) = 0.0
               end if
            end do
         end if
         lastu = nzcpy - 1
c     
c         Copy the column elements of L, throwing out zeros.  Keep track
c         of maximum magnitude element for pivot.
c
         if (ucolst(jcol+1) - 1 .lt. lcolst(jcol)) goto 9002
c     
c             Partial pivoting, diagonal elt. has max. magnitude in L.
c
         diagptr = 0
         diagpiv = 0.0D0

         ujjptr = 0
         maxpiv = -1.0D0
         maxpivglb = -1.0D0
         
         do nzptr = lcolst(jcol), ucolst(jcol+1) - 1
            irow = lurow(nzptr)
            utemp = abs(dense(irow))
            
c            if (irow .eq. cperm(jcol)) then
            if (pattern(irow) .eq. 2) then
               diagptr = irow
               diagpiv = utemp
c               if (diagpiv .eq. 0) 
c     $           print*,
c     $            'WARNING: Numerically zero diagonal element at col',
c     $            jcol
            end if

c     original            
c     if (utemp .gt. maxpiv) then

c     do not pivot outside the pattern
c            if (utemp .gt. maxpiv .and.
c     $          pattern(irow) .ne. 0) then
c     pivot outside pattern
           if (utemp .gt. maxpiv) then
               ujjptr = irow
               maxpiv = utemp
            end if

c     global pivot outside pattern
            if (utemp .gt. maxpivglb) then
               maxpivglb = utemp
            end if

         end do

c     threshold pivoting

         if (diagptr .ne. 0 .and.
     $       diagpiv .ge. (pthresh * maxpiv)) then
            ujjptr = diagptr
         end if

c         print *, ujjptr,diagptr,pattern(ujjptr)

         if (diagptr .eq. 0 .and. ujjptr .eq. 0)
     $        print *,'ERROR!!!',ucolst(jcol+1) - lcolst(jcol)

c         if (diagptr .ne. ujjptr) 
c     $        print*,'pivoting',pthresh,maxpiv,diagpiv,diagptr

         
         diagptr = ujjptr
         ujjptr  = 0

cccccccccccccccccccccccccccccccccccccccccc

         do nzptr = lcolst(jcol), ucolst(jcol+1) - 1
            irow = lurow(nzptr)
            utemp = abs(dense(irow))


c            if (irow .eq. cperm(jcol)) then
c               diagptr = nzcpy
c               diagpiv = utemp
c            end if

c            if (utemp .gt. maxpiv) then
c               ujjptr = nzcpy
c               maxpiv = utemp
c            end if


c     pattern dropping
c            if (pattern(irow) .eq. 0 .and.
c     $          irow .ne. diagptr) then
c     pattern + threshold dropping

            if (pattern(irow) .eq. 0 .and.
     $           irow .ne. diagptr .and.
     $           utemp .lt. ldthreshabs) then
               dense(irow) = 0.0
            else
               if (irow .eq. diagptr) ujjptr = nzcpy

               lurow(nzcpy) = irow
               lu(nzcpy) = dense(irow)
               dense(irow) = 0.0
               nzcpy        = nzcpy + 1
            end if
         end do

         lcolst(jcol)   = lastu + 1
         ucolst(jcol+1) = nzcpy
         lastlu         = nzcpy - 1
         
      end if
c
c     Diagonal element has been found.  Swap U(jcol,jcol) from L into U.
c
      if (ujjptr .eq. 0) goto 9003

      pivrow = lurow(ujjptr)
      ujj    = lu(ujjptr)

c      if (jcol.eq.58 .or. jcol.eq.85) print*,'diag',jcol,ujj,pivrow

      if (ujj .eq. 0.0) goto 9004
      dptr = lcolst(jcol)
      lurow(ujjptr) = lurow(dptr)
      lu(ujjptr)    = lu(dptr)
      lurow(dptr) = pivrow
      lu(dptr)    = ujj
      lcolst(jcol) = dptr + 1
c
c     Record the pivot in P.
c
      rperm(pivrow) = jcol
c      if (pivrow.eq.38) print*,'exchanging',jcol,pivrow
c
c     Divide column jcol of L by U(jcol,jcol).
c
      nzst = lcolst(jcol)
      nzend = ucolst(jcol+1) - 1
      if (nzst .gt. nzend) goto 8000
      flops = flops + dble( nzend - nzst + 1 )
      do 7000 nzptr = nzst, nzend
          lu(nzptr) = lu(nzptr) / ujj
7000  continue
c
c  Normal return.
c
8000  continue
      zpivot = pivrow
      return
c
c  Non-fatal error return: zero pivot found.
c
9001  continue
      print *,'zero length (U-I+L) column'
      zpivot = -1
      return

9002  continue
      print *,'zero length L column'
      zpivot = -1
      return

9003  continue
      print *,'ujjptr not set (1)'
      print *,diagptr
      print *,lcolst(jcol),ucolst(jcol+1) - 1
      zpivot = -1
      return

9004  continue
      print *,'numerically zero diagonal element at column',jcol
      zpivot = -1
      return

      end


