use libc::{c_double, c_int};

extern "C" {
    /// Sparse LU factorization with partial pivoting.
    pub fn lufact_(
        pivot: *const c_int,
        thresh: *const c_double,
        nrow: *const c_int,
        ncol: *const c_int,
        a: *const c_double,
        arow: *const c_int,
        acolst: *const c_int,
        maxlu: *const c_int,
        lastlu: *mut c_int,
        lu: *mut c_double,
        lurow: *mut c_int,
        lcolst: *mut c_int,
        ucolst: *mut c_int,
        perm: *mut c_int,
        error: *mut c_int,
        rwork: *mut c_double,
        iwork: *mut c_int,
    );

    /// Find maximum matching.
    pub fn maxmatch_(
        nrows: *const c_int,
        ncols: *const c_int,
        colstr: *const c_int,
        rowind: *const c_int,
        prevcl: *mut c_int,
        prevrw: *mut c_int,
        marker: *mut c_int,
        tryrow: *mut c_int,
        nxtchp: *mut c_int,
        rowset: *mut c_int,
        colset: *mut c_int,
    );

    /// Fills an integer array with a given value.
    pub fn ifill_(a: *mut c_int, la: *const c_int, ival: *const c_int);
    /// Fills a real*8 array with a given value.
    pub fn rfill_(a: *mut c_double, la: *const c_double, rval: *const c_double);

    /// Depth-first search to allocate storage for U.
    pub fn ludfs_(
        jcol: *const c_int,
        a: *const c_double,
        arow: *const c_int,
        acolst: *const c_int,
        lastlu: *mut c_int,
        lurow: *mut c_int,
        lcolst: *mut c_int,
        ucolst: *mut c_int,
        rperm: *mut c_int,
        cperm: *mut c_int,
        dense: *mut c_double,
        found: *mut c_int,
        parent: *mut c_int,
        child: *mut c_int,
        error: *mut c_int,
    );

    /// Compute one column of L and U in the dense vector.
    pub fn lucomp_(
        jcol: *const c_int,
        lastlu: *mut c_int,
        lu: *mut c_double,
        lurow: *mut c_int,
        lcolst: *mut c_int,
        ucolst: *mut c_int,
        rperm: *const c_int,
        cperm: *const c_int,
        dense: *mut c_double,
        found: *mut c_int,
        pattern: *mut c_int,
        flops: *mut c_double,
    );

    /// Copy dense column to sparse structure, pivot, and divide.
    pub fn lucopy_(
        pivot: *const c_int,
        pthresh: *const c_double,
        dthresh: *mut c_double,
        nzcount: *mut c_int,
        jcol: *const c_int,
        ncol: *const c_int,
        lastlu: *mut c_int,
        lu: *mut c_double,
        lurow: *mut c_int,
        lcolst: *mut c_int,
        ucolst: *mut c_int,
        rperm: *mut c_int,
        cperm: *mut c_int,
        dense: *mut c_double,
        pattern: *mut c_int,
        twork: *mut c_double,
        flops: *mut c_int,
        zpivot: *mut c_int,
    );

    /// Solve lower triangular system.
    pub fn lsolve_(
        n: *const c_int,
        lu: *const c_double,
        lurow: *const c_int,
        lcolst: *const c_int,
        ucolst: *const c_int,
        rperm: *const c_int,
        cperm: *const c_int,
        b: *const c_double,
        x: *mut c_double,
        error: *mut c_int,
    );

    /// Solve upper triangular system.
    pub fn usolve_(
        n: *const c_int,
        lu: *const c_double,
        lurow: *const c_int,
        lcolst: *const c_int,
        ucolst: *const c_int,
        rperm: *const c_int,
        cperm: *const c_int,
        b: *const c_double,
        x: *mut c_double,
        error: *mut c_int,
    );

    /// Solve transpose of lower triangular system.
    pub fn ltsolve_(
        n: *const c_int,
        lu: *const c_double,
        lurow: *const c_int,
        lcolst: *const c_int,
        ucolst: *const c_int,
        rperm: *const c_int,
        cperm: *const c_int,
        b: *const c_double,
        x: *mut c_double,
        error: *mut c_int,
    );

    /// Solve transpose of upper triangular system.
    pub fn utsolve_(
        n: *const c_int,
        lu: *const c_double,
        lurow: *const c_int,
        lcolst: *const c_int,
        ucolst: *const c_int,
        rperm: *const c_int,
        cperm: *const c_int,
        b: *const c_double,
        x: *mut c_double,
        error: *mut c_int,
    );
}
