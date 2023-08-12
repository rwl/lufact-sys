c  Added flops. Sivan Toledo 2/97.

c
c  lucomp :  Compute one column of L and U in the dense vector.
c
      subroutine lucomp (jcol, lastlu, lu, lurow, lcolst, ucolst, 
     1                   rperm, cperm,
     1                   dense, found, pattern, flops)
c
c  This routine computes column jcol of L and U (except for dividing
c  through by U(jcol,jcol)) in the dense vector, which is equal to
c  column jcol of A on entry.  It also allocates space in the sparse
c  data structure for the fill entries in column jcol of L.
c
c  Input parameters:
c    jcol    current column number.
c    rperm   row permutation P.
c            rperm(r) = s > 0 means row r of A is row s < jcol of PA.
c            rperm(r) = 0 means row r of A has not yet been used as a
c            pivot and is therefore still below the diagonal.
c    cperm   column permutation.
c
c  Modified parameters:
c    lastlu  number of positions used in lurow array.
c    lu, lurow, lcolst, ucolst  nonzeros in Pt(L-I+U); see lufact for format.
c            On entry, columns 1 through jcol-1 are complete,
c            ucolst(jcol) and lcolst(jcol) point to the storage for
c            column jcol, and lurow has entries for column jcol of U and
c            the non-fill entries in column jcol of L.  The diagonal
c            element is allocated in L, not U.  On exit, storage has
c            been allocated for all of column jcol of L, and
c            ucolst(jcol+1) has been set up.  No values of column jcol
c            are in lu yet.
c    dense   column jcol as a dense vector.  On entry, column jcol of A;
c            on exit, column jcol of Pt*(U(jcol,jcol)*(L-I) + U).
c    found   found(i)=jcol if storage for position (i,jcol) has been
c            allocated in the sparse data structure.
c    flops   flop count
c
c            Both dense and found are indexed according to the row
c            numbering of A, not PA.
c
c  Local variables:
c    nzuptr                pointer to current nonzero PtU(krow,jcol).
c    nzuend, nnzu, nzuind  used to compute nzuptr.
c    krow                  row index of current nonzero (according to A,
c                          not PA), which is PtU(krow,jcol) = U(kcol,jcol).
c    kcol                  rperm(krow), that is, row index of current
c                          nonzero according to PA.
c    ukj                   value of PtU(krow,jcol).
c    nzlptr                pointer to PtL(irow,kcol), being used for update.
c    nzlst, nzlend         used to compute nzlptr.
c    irow                  row index in which update is taking place,
c                          according to A (not PA).
c
      real*8    lu(1), dense(1), ukj, flops
      integer   jcol, lastlu, lurow(1), lcolst(1), ucolst(1)
      integer   rperm(1), cperm(1), found(1), pattern(1)
      integer   error, nzuptr, nzuend, nnzu, nzuind
      integer   krow, kcol, nzlptr, nzlst, nzlend, irow
c
c
c     For each krow with PtU(krow,jcol) <> 0, in reverse postorder, use
c     column kcol = rperm(krow) of L to update the current column.
c

c      if (jcol.eq.58) print*,'A([39 37 38],58)=',
c     $     dense(39),dense(37),dense(38)
c     $     ,pattern(39),pattern(37),pattern(38)
c      if (jcol.eq.85) print*,'A([39 37 38],85)=',
c     $     dense(39),dense(37),dense(38)
c     $     ,pattern(39),pattern(37),pattern(38)
c      if (jcol.eq.129) print*,'A([39 37 38],129)=',
c     $     dense(39),dense(37),dense(38)
c     $     ,pattern(39),pattern(37),pattern(38)
c      if (dense(38) .ne. 0) print*,'A(',38,',',jcol,')=',dense(38)

      nzuend = lcolst(jcol)
      nnzu   = nzuend - ucolst(jcol)
      if (nnzu .eq. 0) goto 300
      do 200 nzuind = 1, nnzu
          nzuptr = nzuend - nzuind
          krow = lurow(nzuptr)
          kcol = rperm(krow)
          ukj = dense(krow)
ccc          if (pattern(rperm(krow)) .eq. 0) ukj=0.0D0
c
c         For each irow with PtL(irow,kcol) <> 0, update PtL(irow,jcol)
c         or PtU(irow,jcol)
c
          nzlst = lcolst(kcol)
          nzlend = ucolst(kcol+1) - 1
          if (nzlend .lt. nzlst) goto 200
          flops = flops + dble( 2 * (nzlend - nzlst + 1) )
          do 100 nzlptr = nzlst, nzlend
              irow = lurow(nzlptr)
c              if (jcol .eq. 129 .and. irow .eq. 38)

c              if (irow .eq. 38 .or. irow .eq. 39 .or. irow .eq. 37)
c     $             print *,irow,jcol,kcol,
c     $             '---',dense(irow),ukj,lu(nzlptr)
c     $             ,pattern(irow)

c              if (dense(irow) .ge. 1.0D-8 .and.
c     $            dabs(dense(irow) - ukj*lu(nzlptr)) .le. 1.0D-14)
c     $             print *,'WARNING, EXACT CANCELATION',irow,jcol,kcol,
c     $               '---',dense(irow),ukj,lu(nzlptr),':',rperm(irow)


              dense(irow) = dense(irow) - ukj*lu(nzlptr) 

c              if (irow .eq. 38 .or. irow .eq. 39 .or. irow .eq. 37)
c     $             print *,
c     $             'new dense ---',dense(irow)

c
c             If this is a new nonzero in L, allocate storage for it.
c
              if (found(irow) .eq. jcol) goto 100
                  found(irow) = jcol
                  lastlu = lastlu + 1
                  lurow(lastlu) = irow
100       continue
c
200   continue
c
300   continue
      ucolst(jcol+1) = lastlu + 1
      error = 0
      return
c
      end
