c Changed stop to return and added error argument, cperm. 
c Sivan Toledo 2/97.

c
c  ludfs :  Depth-first search to allocate storage for U
c
      subroutine ludfs (jcol, a, arow, acolst, lastlu, lurow, lcolst,
     1                  ucolst, rperm, cperm, dense, 
     1                  found, parent, child,
     1                  error)
c
c  Input parameters:
c    jcol             current column number.
c    a, arow, acolst  the matrix A; see lufact for format.
c    rperm            row permutation P.
c                     perm(r) = s > 0 means row r of A is row s < jcol of PA.
c                     perm(r) = 0 means row r of A has not yet been used as a
c                     pivot and is therefore still below the diagonal.
c    cperm            column permutation.
c
c  Modified parameters (see below for exit values):
c    lastlu           last used position in lurow array.
c    lurow, lcolst, ucolst  nonzero structure of Pt(L-I+U);
c                           see lufact for format.
c    dense            current column as a dense vector.
c    found            integer array for marking nonzeros in this column of
c                     Pt(L-I+U) that have been allocated space in lurow.
c                     Also, marks reached columns in depth-first search.
c                     found(i)=jcol if i was found in this column.
c    parent           parent(i) is the parent of vertex i in the dfs,
c                     or 0 if i is a root of the search.
c    child            child(i) is the index in lurow of the next unexplored
c                     child of vertex i.
c                     Note that parent and child are also indexed according to
c                     the vertex numbering of A, not PA; thus child(i) is
c                     the position of a nonzero in column rperm(i), 
c                     not column i.
c
c  Output parameters:
c    error            0 if successful, 1 otherwise
c
c  On entry:
c    found(*)<jcol
c    dense(*)=0.0
c    ucolst(jcol)=lastlu+1 is the first free index in lurow.
c
c  On exit:
c    found(i)=jcol iff i is a nonzero of column jcol of PtU or
c      a non-fill nonzero of column jcol of PtL.
c    dense(*)=column jcol of A.
c      Note that found and dense are kept according to the row
c      numbering of A, not PA.
c    lurow has the rows of the above-diagonal nonzeros of col jcol of U in
c      reverse topological order, followed by the non-fill nonzeros of col
c      jcol of L and the diagonal elt of U, in no particular order.
c      These rows also are numbered according to A, not PA.
c    lcolst(jcol) is the index of the first nonzero in col j of L.
c    lastlu is the index of the last non-fill nonzero in col j of L.
c
c  Local variables:
c    nzast, nzaend   range of indices in arow for column jcol of A.
c    nzaptr          pointer to current position in arow.
c    krow            current vertex in depth-first search (numbered
c                    according to A, not PA).
c    nextk           possible next vertex in depth-first search.
c    chdend          next index after last child of current vertex.
c    chdptr          index of current child of current vertex
c
      real*8    a(1), dense(1)
      integer   jcol, arow(1), acolst(1), lastlu, lurow(1)
      integer   lcolst(1), ucolst(1), rperm(1), cperm(1)
      integer   error, krow, nextk
      integer   nzast, nzaend, nzaptr, parent(1), child(1), found(1)
      integer   chdend, chdptr
c
c
c  Depth-first search through columns of L from each nonzero of
c  column jcol of A that is above the diagonal in PA.
c
c     For each krow such that A(krow,jcol) is nonzero do ...
c

ccc Old code (no cperm):
ccc      nzast = acolst(jcol)
ccc      nzaend = acolst(jcol+1)

      nzast  = acolst(cperm(jcol))
      nzaend = acolst(cperm(jcol)+1)

      if (nzaend .lt. nzast) goto 800
      nzaend = nzaend - 1
      do 500 nzaptr = nzast, nzaend
          krow = arow(nzaptr)
c
c         Copy A(krow,jcol) into the dense vector.  If above diagonal in
c         PA, start a depth-first search in column rperm(krow) of L.
c
          dense(krow) = a(nzaptr)
          if (rperm(krow) .eq. 0)    goto 500
          if (found(krow) .eq. jcol) goto 500
          if (dense(krow) .eq. 0.0)  goto 500
          parent(krow) = 0
          found(krow) = jcol
          chdptr = lcolst(rperm(krow))
c
c         The main depth-first search loop starts here.
c         repeat
c             if krow has a child that is not yet found
c             then step forward
c             else step back
c         until a step back leads to 0
c
100       continue
c
c             Look for an unfound child of krow.
c
              chdend = ucolst(rperm(krow)+1)
200           continue
                  if (chdptr .ge. chdend) goto 400
                  nextk = lurow(chdptr)
                  chdptr = chdptr + 1
                  if (rperm(nextk) .eq. 0) goto 200
                  if (found(nextk) .eq. jcol) goto 200
c
c             Take a step forward.
c
300           continue
              child(krow) = chdptr
              parent(nextk) = krow
              krow = nextk
              found(krow) = jcol
              chdptr = lcolst(rperm(krow))
              goto 100
c
c             Take a step back.
c
c             Allocate space for U(rperm(k),jcol) = PtU(krow,jcol) in the
c             sparse data structure.
c
400           continue
              lastlu = lastlu + 1
              lurow(lastlu) = krow
              krow = parent(krow)
              if (krow .eq. 0) goto 500
              chdptr = child(krow)
              goto 100
c
c         The main depth-first search loop ends here.
c
500       continue
c
c     Close off column jcol of U and allocate space for the non-fill
c     entries of column jcol of L.  The diagonal element goes in L, not
c     U, until we do the column division at the end of the major step.
c
      lcolst(jcol) = lastlu + 1
      do 600 nzaptr = nzast, nzaend
          krow = arow(nzaptr)
          if (rperm(krow) .ne. 0) goto 600
              found(krow) = jcol
              lastlu = lastlu + 1
              lurow(lastlu) = krow
600   continue
c
      error = 0
      return
c
c     Bug traps and error returns.
c
800   continue
      write (6, 900) jcol, nzast, nzaend
900   format ('In ludfs, negative length for column', i6,
     1        ' of A.  nzast, nzend =', 2i6)
      error = 1
      return
      end
