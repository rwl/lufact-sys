c
c  lufact :  Sparse LU factorization with partial pivoting.
c
      subroutine lufact (pivot, thresh, nrow, ncol, a, arow, acolst,
     1                   maxlu, lastlu, lu, lurow, lcolst, ucolst, perm,
     1                   error, rwork, iwork)
c
c  Given a matrix A in sparse format by columns, perform an LU
c  factorization, with partial or threshold pivoting, if desired.  The
c  factorization is PA = LU, where L and U are triangular.  P, L, and U
c  are returned; see below for format.  This subroutine uses the
c  Coleman-Gilbert-Peierls algorithm, in which total time is O(nonzero
c  multiplications).
c
c  If A is the same array as LU, the solution overwrites the original
c  values of A.
c
c  Input parameters:
c    pivot   = 0 for no pivoting
c            = 1 for partial (row) pivoting
c            = 2 for threshold (row) pivoting
c    thresh  Used only in threshold pivoting.  For each major step of the
c            algorithm, the pivot is chosen to be a nonzero below the
c            diagonal in the current column of A with the most nonzeros to
c            the right in its row, with absolute value at least
c            thresh*maxpiv, where maxpiv is the largest absolute value below
c            the diagonal in the current column.  Note that if thresh .le.
c            0.0, then the pivot is chosen purely on the basis of row
c            sparsity.  Also, if thresh .ge. 1.0, then the pivoting is
c            effectively partial pivoting with ties broken on the basis of
c            sparsity.
c    nrow    number of rows in A.
c    ncol    number of columns in A.
c    a       nonzeros in A.  Nonzeros in each column are contiguous and columns
c            are in left-to-right order, but nonzeros are not necessarily in
c            order by row within each column.
c    arow    arow(i) is the row number of the nonzero a(i).
c    acolst  acolst(j) is the index in a of the first nonzero in column j.
c            acolst(ncol+1) is one more than the number of nonzeros in A.
c    maxlu   size of arrays lu and lurow.
c
c  Output parameters:
c    lastlu  index of last nonzero in lu; number of nonzeros in L-I+U.
c    lu      nonzeros in L and U.  Nonzeros in each column are contiguous.
c            Columns of U and L alternate: u1, l1, u2, l2, ..., un, ln.
c            Nonzeros are not necessarily in order by row within columns.
c            The diagonal elements of L, which are all 1, are not stored.
c    lurow   lurow(i) is the row number of the nonzero lu(i).
c            During the computation, these correspond to row numbers in A,
c            so we really store the non-triangular PtL and PtU.
c            At the end we transform them to row numbers in PA, so the
c            L and U we return are really triangular.
c    lcolst  lcolst(j) is the index in lu of the first nonzero in col j of L.
c    ucolst  ucolst(j) is the index in lu of the first nonzero in col j of U.
c            The last nonzero of col j of U is in position lcolst(j)-1, and
c            the last nonzero of col j of L is in position ucolst(j+1)-1.
c            ucolst(ncol+1) is one more than the number of nonzeros in L-I+U.
c            Notice that ucolst has dimension ncol+1 and lcolst has dimension ncol,
c            although ucolst(ncol+1)=lcolst(ncol) because the last column of L
c            contains only the diagonal one, which is not stored.
c    perm    perm(r) = s means that row r of A is in position s in PA.
c    error   Success or failure indication:
c              0  success
c              1  zero pivot encountered (factorization completed)
c              2  out of space
c
c  Working parameters:
c    rwork   real work vector of length nrow
c    iwork   integer work vector of length 3*nrow (or 4*nrow);
c            partitioned by this routine into the following three (or
c            four) integer vectors of length nrow:
c
c    iwork(found)       integer vector used to control depth-first search.
c    iwork(parent)      also used by ludfs.
c    iwork(child)       also used by ludfs.
c    iwork(rowcnt)      rowcnt(i) is the number of remaining nonzeros in
c                       row i of A (not of PA); for threshold pivoting.
c
c  Local variables:
c    jcol    column number being computed
c    lasta   number of nonzeros in A/AROW
c    locpiv  local pivot code: equals -1 for columns with no diagonal,
c            otherwise same as pivot
c    maxcol  index of current column of LU/LUROW
c    overwr  set to .true. if A and LU as given are equivalent arrays
c    xa      starting index of nonzeros in A/AROW (initially 1)
c    zpivot  set to .true. if a zero pivot is found by lucopy()
c
      real*8    thresh, a(1), lu(1), rwork(1)
      integer   pivot, nrow, ncol, arow(1), acolst(1), maxlu
      integer   lastlu, lurow(1), lcolst(1), ucolst(1), perm(1), error
      integer   iwork(1), found, parent, child, rowcnt
      integer   i, jcol, lasta, locpiv, maxcol, xa
      logical   overwr, requiv, zpivot
c
c
c  Partition IWORK into FOUND, CHILD, PARENT, ROWCNT.
c  (ROWCNT is used only with threshold pivoting.)
c
      found  = 1
      child  = nrow + 1
      parent = 2*nrow + 1
      rowcnt = 3*nrow + 1
c
c  Initialize useful values and zero out the dense vectors.  If we are
c  threshold pivoting, get row counts.
c
      locpiv = pivot
      error = 0
      lastlu = 0
      lasta = acolst(ncol+1) - 1
      ucolst(1) = 1
c
      call ifill(iwork(found), nrow, 0)
      call ifill(perm, nrow, 0)
      call rfill(rwork, nrow, 0.0)
c
      if (pivot .eq. 2) call cntrow(arow, lasta, iwork(rowcnt))
c
c  If A and LU are equivalent, copy A to end of LU.  XA is pointer to
c  start of copied array elements.  OVERWR is set true if A and LU are
c  equivalent.
c
      xa = 1
      overwr = requiv(a, lu)
      if (.not. overwr) goto 100
          xa = maxlu - lasta + 1
          call rcopy(a, a(xa), lasta, .true.)
          call icopy(arow, arow(xa), lasta, .true.)
100   continue
c
c  Compute one column at a time.
c
      do 200 jcol = 1, ncol
c
c         If overwriting, set MAXCOL to just before start of column JCOL
c         of A.  Otherwise, MAXCOL is just MAXLU.  Make sure there is
c         room for up to NROW nonzeros in this column.
c
          maxcol = maxlu
          if (overwr)
     1        maxcol = xa + acolst(jcol) - 2
          if (lastlu + nrow .ge. maxcol) goto 800
c
c         Depth-first search from each above-diagonal nonzero of column
c         jcol of A, allocating storage for column jcol of U in
c         topological order and also for the non-fill part of column
c         jcol of L.
c
          call ludfs (jcol, a(xa), arow(xa), acolst, lastlu, lurow,
     1                lcolst, ucolst, perm, rwork, iwork(found),
     1                iwork(parent), iwork(child))
c
c         Compute the values of column jcol of L and U in the dense
c         vector, allocating storage for fill in L as necessary.
c
          call lucomp (jcol, lastlu, lu, lurow, lcolst, ucolst, perm,
     1                 rwork, iwork(found))
c
c         Copy the dense vector into the sparse data structure, find the
c         diagonal element (pivoting if specified), and divide the
c         column of L by it.
c
          call lucopy (locpiv, thresh, jcol, ncol, iwork(rowcnt),
     1                 lastlu, lu, lurow, lcolst, ucolst, perm, rwork,
     1                 iwork(found), zpivot)
          if (zpivot)
     1        error = 1
c
c         If there are no diagonal elements after this column, change
c         the pivot mode.
c
          if (jcol .eq. nrow)
     1        locpiv = -1
c
200   continue
c
c  Fill in the zero entries of the permutation vector, and renumber the
c  rows so the data structure represents L and U, not PtL and PtU.
c
      jcol = ncol + 1
      do 300 i = 1, nrow
          if (perm(i) .ne. 0) goto 300
              perm(i) = jcol
              jcol = jcol + 1
300   continue
c
      do 400 i = 1, lastlu
          lurow(i) = perm(lurow(i))
400   continue
c
c  Normal return
c
      return
c
c  Error return
c
800   continue
      write (6, 900) maxcol, jcol
900   format ('Out of space:  The limit of maxcol=', i6, 
     1        ' was exceeded at column jcol =', 2i6)
      error = 2
      return
      end
