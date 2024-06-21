C COPYRIGHT (c) 1999 CCLRC Council for the Central Laboratory
*            of the Research Councils
C
C Version 2.3.1
C For full version history see ChangeLog
C

    SUBROUTINE ME57I(CNTL, ICNTL)
C
C****************************************************************
C
C  Purpose
C  =======
C
C  The entries of the arrays CNTL and ICNTL control the action of
C  ME57. Default values for the entries are set in this routine.
C
    REAL    CNTL(5)
    INTEGER         ICNTL(20)
C
C  Parameters
C  ==========
C
C  The entries of the arrays CNTL and ICNTL control the action of
C   ME57. Default values are set by ME57I.
C
C   CNTL(1) has default value 0.01 and is used for threshold pivoting
C    Values less than 0.0 will be treated as 0.0 and values greater
C    than 0.5 as 0.5. Values near 0.0 may perhaps give faster
C    factorization times and less entries in the factors but may
C    result in a less stable factorization.  This parameter is only
C    accessed if ICNTL(7) is equal to 1.
C
C   CNTL(2) has default value $10^{-20}$. ME57B/BD will treat any
C    pivot whose modulus is less than CNTL(2) as zero.
C    If ICNTL(16) = 1, then blocks of entries less than CNTL(2) can
C    be discarded during the factorization and the corresponding
C    pivots are placed at the end of the ordering. In this case,
C    a normal value for CNTL(2) could be $10{-12}$.
C
C   CNTL(3) has default value 0.5.  It is used by ME57D/DD to
C    monitor the convergence of the iterative refinement.
C    If the norm of the scaled residuals does not decrease by a
C    factor of at least CNTL(3), convergence is deemed to be too
C    slow and ME57D/DD terminates with INFO(1) set to -8.
C
C   CNTL(4) has default value 0.0. It is used by ME57B/BD to control
C    the static pivoting. If CNTL(4) > 0.0, then
C    small pivots may be replaced by entries of value CNTL(4) so that
C    the factorization could be inaccurate.  If CNTL(5) is also
C    greater than zero, then CNTL(4) is treated as zero (that is
C    uneliminated variables are delayed) until CNTL(5)*N fully summed
C    variables have been delayed.  If static pivots are used, it is
C    recommended that iterative refinement is used when computing the
C    solution.
C
C   CNTL(5)  has default value zero. Static pivoting is invoked if
C    CNTL(4) is greater than zero and the accumulated number of
C    delayed pivots exceeds CNTL(5)*N.
C
C  ICNTL(1) has default value 6.
C   It is the output stream for error messages. If it is negative,
C   these messages are suppressed.
C
C  ICNTL(2) has default value 6.
C   It is the output stream for warning messages. If it is negative,
C   these messages are suppressed.
C
C  ICNTL(3) has default value 6. It is the output stream for monitoring
C   printing. If it is negative, these messages are suppressed.
C
C  ICNTL(4) has default value -1.
C   It is the output stream for the printing of statistics.
C   If it is negative, the statistics are not printed.
C
C  ICNTL(5) is used by ME57 to control printing of error,
C   warning, and monitoring messages. It has default value 2.
C   Possible values are:
C
C    <1     No messages output.
C   1     Only error messages printed.
C   2     Errors and warnings printed.
C   3     Errors and warnings and terse monitoring
C         (only first ten entries of arrays printed).
C    >3     Errors and warnings and all information
C         on input and output parameters printed.
C
C  ICNTL(6) has default value 5 and must be set by the user to 1
C   if the pivot order in KEEP is to be used by ME57A. For any
C   other value of ICNTL(6), a suitable pivot order will be chosen
C   automatically.  The choices available so far are:
C  ICNTL(6) = 0   AMD using MC47 (with dense row detection disabled)
C  ICNTL(6) = 2   AMD using MC47 (this was previously the MC50 option)
C  ICNTL(6) = 3   MA27 minimum degree ordering
C  ICNTL(6) = 4   METIS_NODEND_ME57 ordering from MeTiS package.
C  ICNTL(6) = 5   Ordering chosen depending on matrix characteristics.
C             At the moment choices are MC47 or METIS.
C             INFO(36) is set to ordering used.
C  ICNTL(6) > 5   At the moment this is treated as 5 (the default).
C
C  ICNTL(7) is used by ME57B to control numerical pivoting.  It
C   has default value 1.  Values out of range cause an error return
C   with INFO(1) equal to -10.
C   Possible values are:
C
C  1   Numerical pivoting is performed using the threshold value in
C    CNTL(1).
C
C  2   No pivoting will be performed and an error exit will occur
C    immediately a sign change or a zero is detected among the pivots.
C    This is suitable for cases when A is thought to be definite
C    and is likely to decrease the factorization time while still
C    providing a stable decomposition.
C
C  3   No pivoting will be performed and an error exit will occur if a
C    zero pivot is detected. This is likely to decrease the
C    factorization time, but may be unstable if there is a sign
C    change among the pivots.
C
C  4   No pivoting will be performed but the matrix will be altered
C    so that all pivots are of the same sign.
C
C  ICNTL(8) has default value 0. If ME57B is called with ICNTL(8) NE 0,
C    then the factorization will discard factors and try
C    to continue the factorization to determine the amount of space
C    needed for a successful factorization.  In this case, a
C    factorization will not have been produced.  If the default value
C    of 0 is used and the factorization stops because of lack of
C    space, the user should reallocate the real or integer space for
C    FACT or IFACT, respectively and reset LFACT or LIFACT
C    appropriately, using ME57E before recalling ME57B.
C
C  ICNTL(9) has default value 10.  It corresponds to the maximum number
C    of steps of iterative refinement.
C
C  ICNTL(10) has default value 0. A positive value will return the
c    infinity norm of the input matrix, the computed solution, and
C    the scaled residual in RINFO(5) to RINFO(7), respectively,
C    a backward error estimate in RINFO(8) and RINFO(9), and an
C    estimate of the forward error in RINFO(10).  If ICNTL(10) is
C    negative or zero no estimates are returned.
C
C  ICNTL(11) The block size to be used by the Level 3 BLAS.
C
C  ICNTL(12) Two nodes of the assembly tree are merged only if both
C    involve less than ICNTL(12) eliminations.
C
C  ICNTL(13) Threshold on number of rows in a block for using BLAS2 in
C    ME57C.
C
C  ICNTL(14) Threshold on number of entries in a row for declaring row
C    full. Set as percentage of N.  So 100 means row must be full.
C
C  ICNTL(15) should be set to 1 (the default) if MC64 scaling is
C   requested.
C
C  ICNTL(16) should be set to 1 if "small" entries are to be removed
C   removed from the frontal matrices.  The default is 0.
C
C  ICNTL(17) to ICNTL(20) are set to zero by ME57I but are not
C   currently used by ME57.
C
C Local variables
    INTEGER I
    REAL ZERO
    PARAMETER (ZERO=0.0E0)
C===============================================
C Default values for variables in control arrays
C===============================================
C   Threshold value for pivoting
    CNTL(1)   = 0.01E0
C Test for zero pivot
    CNTL(2)   = 1.0E-12
C Iterative refinement convergence
    CNTL(3)   = 0.5E0
C Static pivoting control
    CNTL(4) = ZERO
C Control to allow some delayed pivots
    CNTL(5) = ZERO
C   Printing control
    ICNTL(1)  = 6
    ICNTL(2)  = 6
    ICNTL(3)  = 6
    ICNTL(4)  = -1
    ICNTL(5)  = 2
C   Provide pivot order (1=NO)
C   Set to make automatic choice between METIS and MC47
    ICNTL(6)  = 5
C   Pivoting control
    ICNTL(7)  = 1
C   Restart facility
    ICNTL(8)  = 0
C   IR steps
    ICNTL(9)  = 10
C   Error estimates
    ICNTL(10) = 0
C   Blocking for Level 3 BLAS
    ICNTL(11) = 16
C   Node amalgamation parameter (NEMIN)
    ICNTL(12) = 16
C   Switch for use of Level 2 BLAS in solve
    ICNTL(13) = 10
C Flag to indicate threshold will be set to N
    ICNTL(14) = 100
C Flag to invoke MC64 scaling (0 off, 1 on)
    ICNTL(15) = 1
C Flag to invoke dropping small entries from front
C Default is not to drop (set to 1 to drop)
    ICNTL(16) = 0

C Set unused parameters
    DO 110 I=17,20
      ICNTL(I) = 0
  110 CONTINUE

    RETURN
    END


    SUBROUTINE ME57A(KIND,N,NE,IRN,JCN,LKEEP,KEEP,IWORK,
   *            ICNTL,INFO,RINFO)
C This subroutine is a user-callable driver for the analysis phase of
C   ME57. It performs an ordering, a symbolic factorization and
C   computes information for the numerical factorization by ME57B/BD.
    INTEGER KIND,N,NE,IRN(NE),JCN(NE),IWORK(5*N),LKEEP,KEEP(LKEEP),
   *      ICNTL(20),INFO(40)
    REAL RINFO(20)
C
C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C
C  N is an INTEGER variable which must be set by the user to the
C    order n of the matrix A.  It is not altered by the subroutine.
C    Restriction: N > 0.
C
C  NE is an INTEGER variable which must be set by the user to the
C    number of entries being input.  It is not altered by the
C    subroutine. Restriction: NE >= 0.
C
C  IRN and JCN are INTEGER  arrays of length NE. The user
C    must set them so that each off-diagonal nonzero $a_{ij}$ is
C    represented by IRN(k)=i and JCN(k)=j or by IRN(k)=j
C    and JCN(k)=i.  Multiple entries are allowed and any with
C    IRN(k) or JCN(k) out of range are ignored. These arrays will
C    be unaltered by the subroutine.
C
C  IWORK is an INTEGER  array of length 5*N. This need not be set
C    by the user and is used as workspace by ME57A.
C
C  LKEEP is an INTEGER that must be set to length of array KEEP.
C    Restriction: LKEEP >= 5*N+NE+MAX(N,NE)+42
C
C  KEEP is an INTEGER  array of length LKEEP. It need not be set
C    by the user and must be preserved between a call to ME57A
C    and subsequent calls to ME57B and ME57C.
C    If the user wishes to input
C    the pivot sequence, the position of variable {i} in the pivot
C    order should be placed in KEEP(i), i = 1, 2,..., n and ICNTL(6)
C    should be set to 1.  The subroutine may replace the given order
C    by another that gives the same fill-in pattern and
C    virtually identical numerical results.
C
C ICNTL is an INTEGER array of length 10
C    that contains control parameters and must be set by the user.
C    Default values for the components may be set by a call to
C    ME57I/ID. Details of the control parameters are given in ME57I/ID.
C
C INFO is an INTEGER array of length 40 that need not be set by the
C    user.  On return from ME57A, a value of zero for INFO(1)
C    indicates that the subroutine has performed successfully.

C RINFO is a REAL (DOUBLE_PRECISION in the D version) array of length 20
C    that need not be set by the
C    user.  This array supplies information on the execution of ME57A.
C
C**** Still to be updated
C    INFO(1):
C     0  Successful entry.
C    -1 N < 1
C    -2 NE < 0.
C    -9 Invalid permutation supplied in KEEP.
C   -15 LKEEP < 5*N+NE+MAX(N,NE)+42
C    +1 One or more indices out of range.
C    +2 One or more duplicate entries found.
C    +3. Combination of warnings +1 and +2.
C    INFO(2):
C      if INFO(1) = -1, the value input for N.
C      if INFO(1) = -2, the value input for NE.
C      if INFO(1) = -9, index at which error first detected.
C      if INFO(1) = -15, the value input for LKEEP.
C    INFO(3) Number of entries with out-of-range indices.
C    INFO(4) Number of off-diagonal duplicate entries.
C    INFO(5) Forecast number of reals to hold the factorization.
C    INFO(6) Forecast number of integers to hold the factorization.
C    INFO(7) Forecast maximum front size.
C    INFO(8) Number of nodes in the assembly tree.
C    INFO(9) Minimum size for LA of ME57B (without compress).
C    INFO(10) Minimum size for LIFACT of ME57B (without compress).
C    INFO(11) Minimum size for LA of ME57B (with compress).
C    INFO(12) Minimum size for LIFACT of ME57B (with compress).
C    INFO(13) Number of compresses.
C    INFO(14:40) Not used.

C Procedures
    INTRINSIC MIN
    EXTERNAL ME57G,MC47I,MC47B,ME57V,ME57H,ME57J,ME57K,
   *       ME57L,ME57M,ME57N
C ME57G Expand representation to whole matrix and sort.
C MC47B Approximate Minimum Degree avoiding problems with dense rows.
C ME57V Is same as MA27GD. Sort for using MA27HD/ME57H.
C ME57H Is same as MA27HD. Minimum degree ordering from MA27.
C ME57J Sort to upper triangular form (pivot sequence given).
C ME57K Construct tree pointers given output from ME57J.
C ME57L Depth-first search of tree.
C ME57M Construct map.
C ME57N Calculate storage and operation counts.

C Local variables
    INTEGER I,IL,IN,IPE,IRNPRM,COUNT,FILS,FRERE,HOLD,IFCT,INVP,IPS,
   +      IW,IWFR,K,LDIAG,LP,LW,LROW,MAP,EXPNE,
   +      MP,NCMPA,NEMIN,NODE,NST,NSTEPS,NV,PERM,
   +      IW1,IW2,IW3,IW4,IW5,NSTK,ND,NELIM,NZE,ALENB,
   +      J,JJ,J1,J2,SIZE22,OXO
C Local variables for MA27 ordering call
    INTEGER IF27H
C Local variables for MeTiS call
    INTEGER METOPT(8),METFTN,ICNTL6,INF47(10),ICNT47(10)
    REAL ZERO,THRESH,AVNUM,MC47FI,RINF47(10)
    PARAMETER (ZERO=0.0E0)

C I      Temporary DO index.
C IPE    Displacement in array KEEP for array IPE of ME57G/HD/JD/KD.
C COUNT    Displacement in array KEEP for array COUNT of
C      ME57G/HD/JD/KD/LD/ND.
C IWFR   First unused location in IFCT(1:LW).
C K      Temporary variable.
C LDIAG    Control for amount of information output.
C LP     Stream number for error messages.
C LW     Length of IFACT when four arrays of length N+1 are excluded.
C LROW   Subscript in array KEEP for array LAST of ME57K/LD
C      and LROW of ME57M.
C MAP    Subscript in array KEEP for MAP array.
C MP     Stream number for monitoring.
C NODE   Subscript in array KEEP for array FLAG of ME57K, and
C      NODE of ME57L.
C NSTEPS   Number of nodes in the assembly tree.
C NV     Displacement in array IFACT for array NV of ME57G/HD, PERM
C      of ME57J/JD, IPR of ME57K/KD, and NE of ME57L\LD.
C PERM   Subscript in array KEEP for array PERM of ME57L/ND.
C SP     Stream number for statistics.
C SIZES    subscript in array KEEP.  On exit,
C      KEEP(SIZES) is number of faulty entries.
C      KEEP(SIZES+1) is the number of nodes in the assembly tree.


C Set local print variables
    LP = ICNTL(1)
    MP = ICNTL(3)
    LDIAG = ICNTL(5)

C Initialize information array.
    DO 10 I = 1,40
      INFO(I) = 0
   10 CONTINUE
    DO 11 I = 1,20
      RINFO(I) = ZERO
   11 CONTINUE

C Check N, NE, and LKEEP for obvious errors
    IF (N.LT.1)  GO TO 20
    IF (NE.LT.0) GO TO 30
    IF (LKEEP.LT.5*N+NE+MAX(N,NE)+42) GO TO 40

    IF (ICNTL(6).EQ.1) THEN
C Check permutation array
      DO 12 I = 1,N
      IWORK(I) = 0
   12   CONTINUE
      DO 14 I=1,N
      K = KEEP(I)
      IF (K.LE.0 .OR. K.GT.N) GO TO 80
      IF (IWORK(K).NE.0) GO TO 80
      IWORK(K) = I
   14   CONTINUE
    ENDIF

C If requested, print input variables.
    IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
      WRITE(MP,'(/A)') 'Entering analysis phase (ME57A) with ...'
      IF (KIND.EQ.1) THEN
      WRITE(MP,'(A)') 'Matrix is complex symmetric'
      ELSE
      WRITE(MP,'(A)') 'Matrix is Hermitian'
      ENDIF
      WRITE(MP,99980) N,NE,(ICNTL(I),I=1,7),ICNTL(12),ICNTL(15)
99980 FORMAT('N       Order of matrix             =',I12/
   2    'NE      Number of entries             =',I12/
   6    'ICNTL(1)  Stream for errors             =',I12/
   7    ' --- (2)  Stream for warnings           =',I12/
   8    ' --- (3)  Stream for monitoring         =',I12/
   9    ' --- (4)  Stream for statistics         =',I12/
   1    ' --- (5)  Level of diagnostic printing      =',I12/
   2    ' --- (6)  Flag for input pivot order      =',I12/
   2    ' --- (7)  Numerical pivoting control (st est) =',I12/
   2    ' --- (12) Node amalgamation parameter       =',I12/
   2    ' --- (15) Scaling control (storage estimate)  =',I12)
      K = MIN(10,NE)
      IF (LDIAG.GE.4) K = NE
C      WRITE (MP,'(/A/(3(I6,A,2I8,A)))') ' Matrix entries:',
C   +      (I,': (',IRN(I),JCN(I),')',I=1,K)
      IF (K.LT.NE) WRITE (MP,'(A)') '   . . .'

      IF (ICNTL(6).EQ.1) THEN
C Print out permutation array.
      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(A,10I6:/(7X,10I6))') ' KEEP =', (KEEP(I),I=1,K)
      IF (K.LT.N) WRITE (MP,'(7X,A)') '   . . .'
      END IF

    END IF

C Partition IWORK
    IW1 = 1
    IW2 = IW1 + N
    IW3 = IW2 + N
    IW4 = IW3 + N
    IW5 = IW4 + N
    FILS  = IW1
    FRERE = IW2
    ND    = IW3
    NELIM = IW4
    NV    = IW5

C Partition KEEP
    PERM = 1
    NSTEPS = PERM + N
    EXPNE  = NSTEPS + 1
    HOLD   = EXPNE + 1
    LROW = HOLD + 40
    NODE = LROW + N
    NSTK = NODE + N
    MAP  = NSTK + N
    IRNPRM = MAP + MAX(N,NE)
    INVP  = NODE
    IW    = NODE
    IPE   = LROW
    IFCT  = MAP
C This is set for ME57VD/ME57HD (MA27 orderings) only to allow more
C   space for expanded factors
    IF27H = NODE
    IPS   = MAP
    COUNT = NSTK

C Set HOLD(1)
    KEEP(HOLD) = 0

C Sort and order ... generate tree

C Set local value for ICNTL(6)
    ICNTL6 = ICNTL(6)
    IF (ICNTL(6).GT.5) ICNTL6 = 5

    IF (ICNTL6.EQ.4 .OR. ICNTL6.EQ.5) THEN
C MeTiS ordering requested.  Use dummy call to see if it has been
C   installed.
C Set flag for Fortran-style numbering of arrays
      METFTN    = 1
C Set default values for parameters.
      METOPT(1) = 0
C Dummy call with 1 by 1 matrix
      KEEP(IPE)   = 1
      KEEP(IPE+1) = 2
      KEEP(IFCT)  = 1
      CALL METIS_NODEND_ME57(1,KEEP(IPE),KEEP(IFCT),METFTN,METOPT,
   *            KEEP(NSTK),KEEP(PERM))
C Flag set if dummy code for METIS_NODEND has been used
      IF (KEEP(PERM).EQ.-1) THEN
      IF (ICNTL6 .EQ. 4) GO TO 90
C Reset ICNTL6 to use MC47
      ICNTL6 = 2
      ENDIF
    ENDIF
    IF (ICNTL6.NE.1) THEN
C Ordering is to be calculated by program

      CALL MC47I(ICNT47)

      IF (ICNTL6 .NE. 3) THEN
C ELSE clause (ICNTL6.EQ.3) for MA27 ordering
C MC47 (with dense row detection disabled) used if ICNTL6 equal to 0.
C MC47 used if ICNTL6 equal to 2.
C METIS used if ICNTL6 equal to 4.
C Automatic choice of METIS or MC47 if ICNTL6 equal to 5.

C Sort matrix to obtain complete pattern (upper and lower triangle)
C   but omitting diagonals, duplicates, and out-of-range entries.
C   On exit, sorted matrix is given by IPE(pointers), IFCT (indices),
C   and COUNT (lengths).  IWFR is position after last index in IFCT.
      CALL ME57G(KIND,N,NE,IRN,JCN,KEEP(IFCT),KEEP(IPE),
   +          KEEP(COUNT),KEEP(IW),IWFR,ICNTL,INFO)
      IF (INFO(1).LT.0) RETURN

      IF (ICNTL6.EQ.5) THEN
C Calculate matrix statistics to determine ordering
        IF (ICNTL(7).EQ.2) THEN
C Action if positive definite option has been chosen.
          AVNUM = FLOAT(IWFR+N-1)/FLOAT(N)
          IF (N.GE.50000) THEN
          ICNTL6 = 4
C         IF (AVNUM.LE.6.0) ICNTL6 = 2
          GO TO 97
          ENDIF
          IF (N.LE.30000) THEN
          ICNTL6 = 2
          IF (AVNUM.GT.100.0) ICNTL6 = 4
          GO TO 97
          ENDIF
          IF (N.GT.30000 .AND. N.LT.50000) THEN
          IF (AVNUM.GT.46.0) THEN
            ICNTL6 = 4
          ELSE
            ICNTL6 = 2
          ENDIF
          GO TO 97
          ENDIF
        ELSE
C Matrix has not been declared positive definite.
          AVNUM = FLOAT(IWFR+N-1)/FLOAT(N)
C Set flag for detection of OXO matrix
          OXO = 0
C Check for KKT and OXO.
C Calculate size of possible  trailing block
          J2 = IWFR - 1
          SIZE22 = 0
          DO 100 J = N,1,-1
          J1 = KEEP(IPE+J-1)
C Note that ME57G does not sort within order in columns
          DO  99 JJ = J1,J2
            IF (KEEP(IFCT+JJ-1).GT.J) GO TO 101
   99       CONTINUE
          SIZE22 = SIZE22 + 1
          J2 = J1-1
  100       CONTINUE
  101       IF (SIZE22 .GT. 0) THEN
C Check to see if there are no entries in (1,1) block.
          DO 98 I = 1,NE
            IF (IRN(I) .LE. N-SIZE22
   *      .AND. JCN(I) .LE. N-SIZE22) THEN
              AVNUM = FLOAT(IWFR+N-SIZE22-1)/FLOAT(N)
              GO TO 96
            ENDIF
   98       CONTINUE
C The (1,1) block is zero.
          OXO = 1
          AVNUM = FLOAT(IWFR-1)/FLOAT(N)
          ENDIF
   96       IF (N .GE. 100000) THEN
          IF (AVNUM.GT.5.42) THEN
            ICNTL6 = 4
          ELSE
            ICNTL6 = 2
          ENDIF
          GO TO 97
          ENDIF
C Logic for OXO matrices
          IF (OXO.EQ.1) THEN
          IF (FLOAT(N-SIZE22)/FLOAT(SIZE22) .GT. 1.8D0) THEN
            ICNTL6 = 2
          ELSE
            ICNTL6 = 4
          ENDIF
          GO TO 97
          ENDIF
C We can try further simple logic here ... then ...
C Call MC47 to test whether fill-in projected is large
          LW = LKEEP-IFCT+1
          CALL MC47B(N,LW,KEEP(IPE),IWFR,KEEP(COUNT),
   +            KEEP(IFCT),IWORK(NV),
   +            KEEP(INVP),KEEP(PERM),IWORK(IW1),
   +            IWORK(IW2),IWORK(IW3),IWORK(IW4),
   +            ICNT47,INF47,RINF47)
          INFO(13) = INF47(2)
          ICNTL6 = 2
          NEMIN    = ICNTL(12)
    CALL ME57L(N,KEEP(IPE),IWORK(NV),KEEP(IPS),IWORK(NELIM),
   +        KEEP(NSTK),KEEP(NODE),KEEP(PERM),
   +        KEEP(NSTEPS),IWORK(FILS),IWORK(FRERE),IWORK(ND),
   +        NEMIN,KEEP(IRNPRM))
          NST = KEEP(NSTEPS)
    CALL ME57M(KIND,N,NE,IRN,JCN,KEEP(MAP),KEEP(IRNPRM),
   +        KEEP(LROW),KEEP(PERM),
   +        IWORK(IW2),IWORK(IW5))
          KEEP(EXPNE) = IWORK(IW5)
    CALL ME57N(N,KEEP(LROW),KEEP(NSTK),IWORK(NELIM),
   +        IWORK(ND),NST,IWORK(IW1),IWORK(IW2),
   +        INFO,RINFO)
C Check relative fill-in of MC47
          IF (FLOAT(INFO(5))/FLOAT(NE) .LT. 10.0) THEN
C We will run with the MC47 ordering
          GO TO 93
          ELSE
C Save value of relative fill-in for testing against METIS
          MC47FI = FLOAT(INFO(5))/FLOAT(NE)
C Must test METIS ordering now ... ugh
      CALL ME57G(KIND,N,NE,IRN,JCN,KEEP(IFCT),KEEP(IPE),
   +          KEEP(COUNT),KEEP(IW),IWFR,ICNTL,INFO)
      IF (INFO(1).LT.0) RETURN
          KEEP(IPE+N) = IWFR
          METFTN    = 1
          METOPT(1) = 0
          IF (N.LT.50) GO TO 92
          DO 91 I = 1,N
            IF ((KEEP(IPE+I)-KEEP(IPE+I-1)) .GT. N/10) THEN
            METOPT(1) = 1
            METOPT(2) = 3
            METOPT(3) = 1
            METOPT(4) = 2
            METOPT(5) = 0
            METOPT(6) = 1
            METOPT(7) = 200
            METOPT(8) = 1
            GO TO 92
            ENDIF
   91       CONTINUE
   92   CALL METIS_NODEND_ME57(N,KEEP(IPE),KEEP(IFCT),METFTN,METOPT,
   *              KEEP(NSTK),KEEP(PERM))
      CALL ME57J(KIND,N,NE,IRN,JCN,KEEP(PERM),KEEP(IFCT),KEEP(IPE),
   +          KEEP(COUNT),IWORK(IW1),IWFR,ICNTL,INFO)
      IF (INFO(1).LT.0) RETURN
          LW = LKEEP - IFCT + 1
      CALL ME57K(N,KEEP(IPE),KEEP(IFCT),LW,IWFR,KEEP(PERM),
   +          KEEP(INVP),IWORK(NV),IWORK(IW1),NCMPA)
          INFO(13) = NCMPA
          NEMIN = ICNTL(12)
    CALL ME57L(N,KEEP(IPE),IWORK(NV),KEEP(IPS),IWORK(NELIM),
   +        KEEP(NSTK),KEEP(NODE),KEEP(PERM),
   +        KEEP(NSTEPS),IWORK(FILS),IWORK(FRERE),IWORK(ND),
   +        NEMIN,KEEP(IRNPRM))
          NST = KEEP(NSTEPS)
    CALL ME57M(KIND,N,NE,IRN,JCN,KEEP(MAP),KEEP(IRNPRM),
   +        KEEP(LROW),KEEP(PERM),
   +        IWORK(IW2),IWORK(IW5))
          KEEP(EXPNE) = IWORK(IW5)
    CALL ME57N(N,KEEP(LROW),KEEP(NSTK),IWORK(NELIM),
   +        IWORK(ND),NST,IWORK(IW1),IWORK(IW2),
   +        INFO,RINFO)
          IF (FLOAT(INFO(5))/FLOAT(NE).LT.MC47FI) THEN
            ICNTL6 = 4
            GO TO 93
          ELSE
C Double groan  ... we will run with MC47 after all
            ICNTL6=2
C KEEP(IPE) has been corrupted must reset it.
      CALL ME57G(KIND,N,NE,IRN,JCN,KEEP(IFCT),KEEP(IPE),
   +          KEEP(COUNT),KEEP(IW),IWFR,ICNTL,INFO)
      IF (INFO(1).LT.0) RETURN
            GO TO 97
          ENDIF
C End of METIS check
          ENDIF
C End of indef case calculation
        ENDIF
C End of logic for ICNTL6 = 5
      ENDIF

   97   IF (ICNTL6.EQ.4) THEN
C Set last pointer in IPE
        KEEP(IPE+N) = IWFR
C Use MeTiS ordering
C Set flag for Fortran-style numbering of arrays
        METFTN    = 1
C This would use only defaults
        METOPT(1) = 0
C Set options for METIS, particularly one for dense columns
C First determine if there are any dense columns
        IF (N.LT.50) GO TO 103
        DO 102 I = 1,N
          IF ((KEEP(IPE+I)-KEEP(IPE+I-1)) .GT. N/10) THEN
C The rest are set to default values
          METOPT(1) = 1
          METOPT(2) = 3
          METOPT(3) = 1
          METOPT(4) = 2
          METOPT(5) = 0
          METOPT(6) = 1
          METOPT(7) = 200
          METOPT(8) = 1
          GO TO 103
          ENDIF
  102     CONTINUE
  103   CALL METIS_NODEND_ME57(N,KEEP(IPE),KEEP(IFCT),METFTN,METOPT,
   *              KEEP(NSTK),KEEP(PERM))
        GO TO 111
      ENDIF


C Obtain ordering using approximate minimum degree ordering.
C Input IPE,IWFR,COUNT,IFCT as from ME57G.
C Output
C   IPE (- father pointer .. if NV(I) > 0, - subordinate variable
C   pointer if NV(I) = 0)
C   NV(I) for subordinate variables of supervariable, otherwise
C   degree when eliminated.
C   IWFR is set to length required by MC47B if no compresses.
C   COUNT, IFCT undefined on exit.
C   INVP is inverse permutation and PERM is permutation
C Length of LW set to maximum to avoid compresses in MC47B/BD
      LW = LKEEP-IFCT+1
C ICNTL6 =  0.  MC47 uses code for dealing with dense rows disabled
C In HSL 2002 it was only used in the F90 version.
C ICNTL6 =  2.  MC47 implements code for dealing with dense rows
      IF (ICNTL6 .EQ. 0) ICNT47(4) = -1
      CALL MC47B(N,LW,KEEP(IPE),IWFR,KEEP(COUNT),
   +          KEEP(IFCT),IWORK(NV),
   +          KEEP(INVP),KEEP(PERM),IWORK(IW1),
   +          IWORK(IW2),IWORK(IW3),IWORK(IW4),
   +          ICNT47,INF47,RINF47)
      INFO(13) = INF47(2)

      ELSE
C End of ICNTL6 .NE. 3

C MA27 ordering being used.  Must insert row lengths in KEEP(IFCT)
C Length of LW set to maximum to avoid compresses.
      LW = LKEEP-IF27H+1
      CALL ME57V(KIND,N,NE,IRN,JCN,KEEP(IF27H),LW,KEEP(IPE),
   *          IWORK(IW1),IWORK(IW2),IWFR,ICNTL,INFO)
      IF (INFO(1).LT.0) RETURN
C Analyse using minimum degree ordering
      THRESH = FLOAT(ICNTL(14))/100.0
      CALL ME57H(N,KEEP(IPE),KEEP(IF27H),LW,IWFR,IWORK(NV),
   *          IWORK(IW1),IWORK(IW2),IWORK(IW3),IWORK(IW4),
   +          2139062143,INFO(13),THRESH)
C Set IPE correctly
      DO 110 I = 1,N
        IF (IWORK(NV+I-1).NE.0) GO TO 110
        IN = I
  105     IL = IN
        IN = - KEEP(IPE+IL-1)
        IF (IWORK(NV+IN-1).EQ.0) GO TO 105
C Make subordinate node point to principal node
        KEEP(IPE+I-1) = -IN
  110   CONTINUE
      ENDIF

C End of block for generating ordering
    ENDIF

  111 IF (ICNTL6.EQ.1 .OR. ICNTL6.EQ.4) THEN
C If we have generated ordering using MeTiS then we need to feed
C   permutation as if it were coming from the user as we do not
C   have a tight coupling to MeTiS as for other orderings.

C Sort using given order
      CALL ME57J(KIND,N,NE,IRN,JCN,KEEP(PERM),KEEP(IFCT),
   +       KEEP(IPE),KEEP(COUNT),IWORK(IW1),IWFR,ICNTL,INFO)
      IF (INFO(1).LT.0) RETURN

C Generating tree using given ordering
C Input:  N,KEEP(IPE),KEEP(IFCT),LW,IWFR,KEEP(PERM)
C Output:  KEEP(IPE),IWORK(NV)
CCC For some reason LW was set to 2*NE which failed when NE < N
C     LW = 2*NE
      LW = LKEEP - IFCT + 1
      CALL ME57K(N,KEEP(IPE),KEEP(IFCT),LW,IWFR,KEEP(PERM),
   +          KEEP(INVP),IWORK(NV),IWORK(IW1),NCMPA)
      INFO(13) = NCMPA

    END IF


C Perform depth-first search of assembly tree
C Set NEMIN
    NEMIN = ICNTL(12)
C Input  IPE,NV,NEMIN
C Output
C   IPE .. father and younger brother pointer
C   NV  .. unchanged
C   NE/NSTK/ND defined for nodes of tree
C   PERM
C   IPS(I) position of node I in order
C   LROW(I) is size of frontal matrix at node I
    CALL ME57L(N,KEEP(IPE),IWORK(NV),KEEP(IPS),IWORK(NELIM),
   +        KEEP(NSTK),KEEP(NODE),KEEP(PERM),
   +        KEEP(NSTEPS),IWORK(FILS),IWORK(FRERE),IWORK(ND),
   +        NEMIN,KEEP(IRNPRM))
    NST = KEEP(NSTEPS)

C Construct map for storing the permuted upper triangle by rows.
C Input N,NE,IRN,JCN,PERM
C Output MAP,LROW,IRNPRM
    CALL ME57M(KIND,N,NE,IRN,JCN,KEEP(MAP),KEEP(IRNPRM),
   +        KEEP(LROW),KEEP(PERM),
   +        IWORK(IW2),IWORK(IW5))

C Set number of entries in expanded input matrix
    KEEP(EXPNE) = IWORK(IW5)

C Evaluate storage and operation counts.
C Input  LROW,NSTK,NELIM,ND
C Output LROW,NSTK (unchanged)
    CALL ME57N(N,KEEP(LROW),KEEP(NSTK),IWORK(NELIM),
   +        IWORK(ND),NST,IWORK(IW1),IWORK(IW2),
   +        INFO,RINFO)

C Set INFO entry to record ordering used
   93 INFO(36) = ICNTL6
C Add for BIGA
    ALENB    = 1
C Add for Schnabel-Eskow
    IF (ICNTL(7).EQ.4) ALENB = ALENB + N + 5
C Add for scaling
    IF (ICNTL(15).EQ.1) ALENB = ALENB + N

C Allow enough to get started
    INFO(9)  = MAX(INFO(9)+ALENB,ALENB+KEEP(EXPNE)+1)
    INFO(11) = MAX(INFO(11)+ALENB,ALENB+KEEP(EXPNE)+1)
C This is N+5 for starting the factorization, N for first row (max).
    INFO(10) = MAX(INFO(10),KEEP(EXPNE)+N+5)
    INFO(12) = MAX(INFO(12),KEEP(EXPNE)+N+5)

C Needed by ME57B/BD
    IF (ICNTL(15).EQ.1) THEN
      INFO(9) = MAX(INFO(9),ALENB+3*KEEP(EXPNE)+3*N)
      INFO(11) = MAX(INFO(11),ALENB+3*KEEP(EXPNE)+3*N)
C Allow space for integers in computing scaling factors
      INFO(10) = MAX(INFO(10),3*KEEP(EXPNE)+5*N+1)
      INFO(12) = MAX(INFO(12),3*KEEP(EXPNE)+5*N+1)
    ENDIF

C If requested, print parameter values on exit.
    IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
      NZE = KEEP(EXPNE)
      WRITE (MP,99999) INFO(1),NZE,
   *            (INFO(I),I=3,13),INFO(36),(RINFO(I),I=1,2)
99999 FORMAT (/'Leaving analysis phase (ME57A) with ...'/
   1    'INFO(1)  Error indicator              =',I12/
   2    'Number of entries in matrix incl diagonal   =',I12/
   2    'INFO(3)  Number of out-of-range indices     =',I12/
   2    'INFO(4)  Number of off-diagonal duplicates    =',I12/
   2    'INFO(5)  Forecast real storage for factors    =',I12/
   3    '----(6)  Forecast integer storage for factors =',I12/
   3    '----(7)  Forecast maximum front size      =',I12/
   4    '----(8)  Number of nodes in assembly tree   =',I12/
   5    '----(9)  Size of FACT without compress      =',I12/
   6    '----(10) Size of IFACT without compress     =',I12/
   5    '----(11) Size of FACT with compress       =',I12/
   5    '----(12) Size of IFACT with compress      =',I12/
   5    '----(13) Number of compresses           =',I12/
   5    '----(36) Ordering strategy used by code     =',I12/
   9    'RINFO(1) Forecast additions for assembly    =',1P,E12.5/
   9    'RINFO(2) Forecast ops for elimination       =',1P,E12.5)

      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(/A/(5I12))')  'Permutation array:',
   +            (KEEP(I),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
      WRITE (MP,'(/A/(5I12))')
   +      'Number of entries in rows of permuted matrix:',
   +      (KEEP(LROW+I-1),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
      K = MIN(10,NZE)
      IF (LDIAG.GE.4) K = NZE
      WRITE (MP,'(/A/(5I12))')
   *      'Column indices of permuted matrix:',
   *                 (KEEP(IRNPRM+I-1),I=1,K)
      IF (K.LT.NZE) WRITE (MP,'(16X,A)') '   . . .'
      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(/A/(5I12))')
   +    'Tree nodes at which variables eliminated:',
   +    (KEEP(NODE+I-1),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
      K = MIN(10,NE)
      IF (LDIAG.GE.4) K = NE
      WRITE (MP,'(/A/(5I12))') 'Map array:',
   *                   (KEEP(I),I=MAP,MAP+K-1)
      IF (K.LT.NE) WRITE (MP,'(16X,A)') ' . . .'
    END IF

    RETURN

C Error conditions.
   20 INFO(1) = -1
    INFO(2) = N
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A,I10)')
   +    '**** Error return from ME57A ****  INFO(1) =',INFO(1),
   +    'N has value ',INFO(2)
    RETURN

   30 INFO(1) = -2
    INFO(2) = NE
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A,I10)')
   +    '**** Error return from ME57A ****  INFO(1) =',INFO(1),
   +    'NE has value',INFO(2)
     RETURN

   40 INFO(1) = -15
    INFO(2) = LKEEP
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A,I10/A,I10)')
   +    '**** Error return from ME57A ****  INFO(1) =',INFO(1),
   +    'LKEEP has value    ',INFO(2),
   +    'Should be at least ',5*N+NE+MAX(N,NE)+42
     RETURN

   80 INFO(1) = -9
    INFO(2) = I
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A/A,I10,A)')
   +    '**** Error return from ME57A ****  INFO(1) =',INFO(1),
   +    'Invalid permutation supplied in KEEP',
   +    'Component',INFO(2),' is faulty'
    RETURN

   90 INFO(1) = -18
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A)')
   +    '**** Error return from ME57A ****  INFO(1) =',INFO(1),
   +    'MeTiS ordering requested but MeTiS not linked'

    END


C--------------------------------------------------------------------
C        HSL 2000 (2000)
C      --
C-         Copyright Rutherford Appleton Laboratory
C      --
C--------------------------------------------------------------------
    SUBROUTINE ME57B(KIND, N, NE, A, FACT, LFACT, IFACT, LIFACT,
   * LKEEP, KEEP, PPOS, WORK, ICNTL, CNTL, INFO, RINFO)
C
C Purpose
C =======
C
C
C This subroutine computes the factorization of the matrix input in
C   A using information (in KEEP and IFACT) from ME57A.
C
    INTEGER KIND,N,NE,LFACT,LIFACT,LKEEP
    COMPLEX A(NE),FACT(LFACT)
    REAL WORK(2*NE+3*N)
    REAL RINFO(20)
C
C   Control parameters: see description in ME57I
    REAL CNTL(5)
    INTEGER ICNTL(20), IFACT(LIFACT)
    INTEGER   INFO(40), KEEP(LKEEP), PPOS(N)
C
C Parameters
C ==========
C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C
C N is an INTEGER variable which must be set by the user to the
C    order n of the matrix A. It must be unchanged since the
C    last call to ME57A and is not altered by the
C    subroutine.  Restriction: N > 0
C
C NE is an INTEGER variable which must be set by the user to the
C    number of    entries in the matrix A.  It is not altered by
C    the subroutine.  Restriction: NE >= 0.
C
C A is a COMPLEX (COMPLEX in the D version) array of length NE.
C   It is not altered by the subroutine.
C
C FACT is a COMPLEX (COMPLEX in the D version)
C    array of length LFACT. It need not
C    be set by the user. On exit, entries 1 to INFO(15) of FACT hold
C    the real part of the factors and should be passed unchanged to
C    ME57C.
C
C LFACT is an INTEGER variable that must be set by the user to
C    the size of array FACT.
C    It should be passed unchanged to ME57C.
C
C IFACT is an INTEGER array of length LIFACT. It need not
C    be set by the user. On exit, entries 1 to INFO(16) of IFACT hold
C    the integer part of the factors and should be passed unchanged to
C    ME57C.
C
C LIFACT is an INTEGER variable that must be set by the user to
C    the size of array IFACT.
C    It should be passed unchanged to ME57C.
C
C LKEEP is an INTEGER that must be set to the length of array KEEP.
C    Restriction: LKEEP >= 5*N+NE+MAX(N,NE)+42.
C
C KEEP is an INTEGER  array of length LKEEP which must be
C     passed unchanged since the last call to ME57A.  It is not
C     altered by ME57B.
C
C PPOS is an INTEGER array of length N that is used as workspace.
C
C ICNTL is an INTEGER array of length 20
C        that contains control parameters and must be set
C        by the user.
C        Default values for the components may be set by a call
C        to ME57I. Details
C        of the control parameters are given in ME57I
C
C CNTL is a REAL (DOUBLE_PRECISION in the D version) array of length 5
C        that contains control parameters and must be set
C        by the user.
C        Default values for the components may be set by a call
C        to ME57I. Details
C        of the control parameters are given in ME57I
C
C RINFO is a REAL (DOUBLE_PRECISION in the D version) array of
C       length 20 which need not be set by the user.
C       We describe below the components of this array modified
C       in the subroutine.
C
C    ______(3)  Number of floating point operations involved
C       during the assembly process
C
C    ______(4)  Number of floating point operations involved
C       during the elimination process
C
C    ______(5)  Number of extra floating point operations caused by
C       use of Level 3 BLAS
C
C
C    .. Error Return ..
C       ============
C
C  A successful return from ME57B
C  is indicated by a value of INFO(1) positive.
C  Negative values of INFO(1) correspond to
C  error message whereas positive values correspond to
C  warning messages. A negative value of INFO(1) is associated with
C  an error message which will be output on unit ICNTL(1).
C
C
C     .. Local variables ..
C      ===============
C
    INTEGER EXPNE,HOLD,I,IRNPRM,K,LDIAG,LLFACT,LP,LROW,MAP,MM1,MM2,MP
    INTEGER J,JJ,KK,ISCALE,NUM,NE64,IDUP,IMAT,IPT,JLOOP,JNEW,NN,ISING
    INTEGER J1,J2
    INTEGER NSTEPS,NODE,NSTK,PERM,INEW,ALENB

    REAL ONE,ZERO,RINF,FD15A,FCT,SMAX,SMIN
    REAL BIGA
    PARAMETER (ONE = 1.0E0, ZERO=0.0E0)

C?? To identify bug
C   LOGICAL LCASE
C   COMMON /CCASE/LCASE

    INTRINSIC MIN

C
C EXPNE is number of entries of original matrix plus any missing
C    diagonals.
C HOLD points to position in IFACT for first entry of array that holds
C    values for restart.
C MM1,MM2 are used to define start point of arrays for matrix
C    modification.  Needed if LFACT < N + 5.
C NSTEPS holds the number of nodes in the tree.

C
C External Subroutines
C ====================
C

    EXTERNAL ME57O,ME57U,FD15A,MC34A,MC64W

C Set RINF to largest positive real number (infinity)
    RINF = FD15A('H')
C Set REPS to smallest number st 1.0+REPS > 1.0
C   REPS = FD15A('E')

C Set INFO for -3 and -4 error return. Ok that is done on every call
    INFO(17) = 0
    INFO(18) = 0

C
C
C Initialisation of printing controls.
C
    LP   = ICNTL(1)
    MP   = ICNTL(3)
    LDIAG  = ICNTL(5)
C
C??
C   Check if analysis has been effectively performed
C
    IF (N.LE.0)  GO TO 25
    IF (NE.LT.0) GO TO 30
    IF (LKEEP.LT.5*N+NE+MAX(N,NE)+42) GO TO 40
    IF (ICNTL(7).LT.1 .OR. ICNTL(7).GT.4) GO TO 35

C   Partition of array KEEP
    NSTEPS = KEEP(N+1)
    EXPNE  = KEEP(N+2)
    PERM = 1
    HOLD = PERM + N + 2
    LROW = HOLD + 40
    NODE = LROW + N
    NSTK = NODE + N
    MAP  = NSTK + N
    IRNPRM = MAP + MAX(NE,N)

    LLFACT = LFACT - 1

    IF (ICNTL(15).EQ.1) THEN
C Matrix is being scaled using MC64SYM
      ISCALE = LLFACT - N + 1
      LLFACT = ISCALE - 1
    ENDIF

    IF (ICNTL(7).EQ.4) THEN
C Schnabel-Eskow modification being used
C Reserve space in FACT for diagonal entries and controls.
      LLFACT = LLFACT - N - 5
C Set MM1 and MM2 to point to first entries in arrays.
      MM1 = LLFACT+6
      MM2 = LLFACT+1
    ELSE
      MM1 = 1
      MM2 = 1
    ENDIF

C One for BIGA
    ALENB = 1
    IF (ICNTL(7).EQ.4)  ALENB = ALENB + N + 5
    IF (ICNTL(15).EQ.1) ALENB = ALENB + N
C +1 because MAP of o-o-r maps to entry 0
    IF (LLFACT.LT.EXPNE+1)   GO TO 85
C The first five positions and room for a whole row are needed
C at the beginning of the ME57O/OD factorization
    IF (LIFACT.LT.EXPNE+N+5)  GO TO 95
C Check that there is enough space for scaling within ME57B/BD.
    IF (ICNTL(15).EQ.1)  THEN
      IF (LFACT .LT. ALENB + 3*EXPNE  + 3*N) GO TO 85
      IF (LIFACT .LT. 3*EXPNE + 5*N + 1) GO TO 95
    ENDIF

C
C PRINTING OF INPUT PARAMETERS
C*****************************
    IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
      WRITE (MP,99999)
99999 FORMAT (/'Entering factorization phase (ME57B) with ...')
      IF (KEEP(HOLD).GT.0) WRITE (MP,99998)
99998 FORMAT ('Re-entry call after call to ME57E')
      IF (KIND.EQ.1) THEN
      WRITE(MP,'(A)') 'Matrix is complex symmetric'
      ELSE
      WRITE(MP,'(A)') 'Matrix is Hermitian'
      ENDIF
      WRITE (MP,99997) N,NE,EXPNE,(ICNTL(I),I=1,5),ICNTL(7),ICNTL(8),
   +       ICNTL(11),ICNTL(15),ICNTL(16),LFACT,LIFACT, NSTEPS,
   +       CNTL(1), CNTL(2), CNTL(4), CNTL(5)
99997 FORMAT ('N     Order of input matrix         =',I12/
   2      'NE    Entries in input matrix         =',I12/
   2      '      Entries in input matrix (inc diags) =',I12/
   6      'ICNTL(1)  Stream for errors           =',I12/
   7      ' --- (2)  Stream for warnings         =',I12/
   8      ' --- (3)  Stream for monitoring         =',I12/
   9      ' --- (4)  Stream for statistics         =',I12/
   1      ' --- (5)  Level of diagnostic printing    =',I12/
   1      ' --- (7)  Numerical pivoting control      =',I12/
   1      ' --- (8)  Restart or discard factors      =',I12/
   1      ' --- (11) Block size for Level 3 BLAS     =',I12/
   1      ' --- (15) Scaling control (1 on)        =',I12/
   1      ' --- (16) Dropping control (1 on)       =',I12/
   4      'LFACT   Size of real working space      =',I12/
   5      'LIFACT  Size of integer working space     =',I12/
   7      '      Number nodes in assembly tree     =',I12/
   9      'CNTL(1) Value of threshold parameter      =',E12.5/
   9      'CNTL(2) Threshold for zero pivot        =',E12.5/
   9      'CNTL(4) Control for value of static pivots  =',E12.5/
   9      'CNTL(5) Control for number delayed pivots   =',E12.5)
      K = MIN(10,NE)
      IF (LDIAG.GE.4) K = NE
      IF (NE.GT.0) THEN
C      WRITE (MP,'(/A/(2(I4,A,1P,2E16.8,A)))') 'Matrix entries:',
C   +   (I,': (',A(I),')',I=1,K)
      IF (K.LT.NE) WRITE (MP,'(A)') '   . . .'
      END IF
      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(/A/(5I12))')  'Permutation array:',
   +            (KEEP(I),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
      WRITE (MP,'(/A/(5I12))')
   +      'Number of entries in rows of permuted matrix:',
   +      (KEEP(LROW+I-1),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
      WRITE (MP,'(/A/(5I12))')
   +    'Tree nodes at which variables eliminated:',
   +    (KEEP(NODE+I-1),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
      K = MIN(10,NSTEPS)
      IF (LDIAG.GE.4) K = NSTEPS
      IF (K.GT.0) WRITE (MP,'(/A/(5I12))')
   +   'Number of assemblies at each tree node:',
   +   (KEEP(NSTK+I-1),I=1,K)
      IF (K.LT.NSTEPS) WRITE (MP,'(16X,A)') ' . . .'
      K = MIN(10,NE)
      IF (LDIAG.GE.4) K = NE
      WRITE (MP,'(/A/(5I12))') 'Map array:',
   *                   (KEEP(I),I=MAP,MAP+K-1)
      IF (K.LT.NE) WRITE (MP,'(16X,A)') ' . . .'
      K = MIN(10,EXPNE)
      IF (LDIAG.GE.4) K = EXPNE
      WRITE (MP,'(/A/(5I12))')
   *      'Column indices of permuted matrix:',
   *                   (KEEP(IRNPRM+I-1),I=1,K)
      IF (K.LT.EXPNE) WRITE (MP,'(16X,A)') '   . . .'
    ENDIF

C Jump if it is reentry
    IF (KEEP(HOLD) .GT. 0) THEN
      BIGA = ABS(FACT(LFACT))
      GO TO 22
    ENDIF

C
C***************************************************
C MAP input nonzeros to appropriate position in FACT
C***************************************************
C
C?? For the moment to handle missing diagonals
    DO 19 K = 1,EXPNE
      FACT(LLFACT-EXPNE+K) = ZERO
   19 CONTINUE

    BIGA = ZERO
    IF (KIND.EQ.1) THEN
      DO 20 K = 1,NE
      BIGA = MAX(BIGA,ABS(A(K)))
      FACT(KEEP(MAP+K-1)+LLFACT-EXPNE) = A(K)
   20   CONTINUE
    ENDIF
    IF (KIND.NE.1) THEN
C Matrix is Hermitian
      DO 201 K = 1,NE
      BIGA = MAX(BIGA,ABS(A(K)))
      KK = KEEP(MAP+K-1)
      IF (KK.GT.0) THEN
        FACT(KK+LLFACT-EXPNE) = A(K)
      ELSE
C Entry is switching triangles
        FACT(-KK+LLFACT-EXPNE) = CONJG(A(K))
      ENDIF
  201   CONTINUE
    ENDIF
    RINFO(18) = BIGA
    DO 21 K = 1,EXPNE
      IFACT(LIFACT-EXPNE+K) = KEEP(IRNPRM+K-1)
   21 CONTINUE
C Invert array PERM
    DO 23 I = 1,N
      PPOS(KEEP(PERM+I-1)) = I
   23 CONTINUE

    IF (ICNTL(15).EQ.1) THEN
C Scaling using MC64.  Matrix must be generated in correct format.

    IPT = 1
    IDUP = IPT+N+1
    IMAT = IDUP+N
    ISING = IMAT + MAX(NE,EXPNE)

C Copy matrix, remove duplicates, and initialize IP array.
    DO 4444 I = 1,N
      IFACT(IDUP+I-1) = 0
 4444 CONTINUE
C Must use new coordinates to keep matrix (half) symmetric
    IFACT(IPT) = 1
    KK = 1
    K = 1
    DO 3333 J = 1,N
      DO 2222 JJ = 1,KEEP(LROW+J-1)
      I = KEEP(PERM+IFACT(LIFACT-EXPNE+K)-1)
      IF (IFACT(IDUP+I-1).GE.IFACT(IPT+J-1)) THEN
C Duplicate
        FACT(IFACT(IDUP+I-1)) =
   &      FACT(IFACT(IDUP+I-1)) + FACT(LLFACT-EXPNE+K)
      ELSE
C Remove explicit zeros
        IF (FACT(LLFACT-EXPNE+K).NE.ZERO) THEN
          IFACT(IDUP+I-1) = KK
          FACT(KK) = FACT(LLFACT-EXPNE+K)
          IFACT(IMAT-1+KK) = I
          KK = KK+1
        ENDIF
      ENDIF
      K = K + 1
 2222   CONTINUE
      IFACT(IPT+J) = KK
 3333 CONTINUE

C Check that Hermitian matrix has reals on diagonal
C Note that this check only applies when scaling is switched on
    IF (KIND.NE.1) THEN
      DO 3111 J = 1,N
      J1 = IFACT(IPT+J-1)
      J2 = IFACT(IPT+J)-1
      DO 2111 JJ = J1,J2
        I = IFACT(IMAT-1+JJ)
        IF (I.EQ.J) THEN
          IF (FACT(JJ) .NE. REAL(FACT(JJ))) THEN
          INFO(1) = -20
          DO 2110 K = 1,N
            IF (KEEP(PERM+K-1).EQ.I) THEN
            INFO(2) = K
            GO TO 2109
            ENDIF
 2110       CONTINUE
 2109 IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I8,A)')
   +   '**** Error return from ME57B ****  INFO(1) =',INFO(1),
   *   'Hermitian matrix but diagonal entry ',INFO(2),' is complex'
          RETURN
          ENDIF
        ENDIF
 2111   CONTINUE
 3111   CONTINUE
    ENDIF

    DO 3334 I=1,KK-1
      WORK(I) = ABS(FACT(I))
 3334 CONTINUE
C Expand matrix
    CALL MC34A(N,IFACT(IMAT),IFACT(IPT),.TRUE.,WORK,KEEP(PERM))
    NE64 = IFACT(IPT+N)-1
    DO 75 J = 1,N
      FCT = ZERO
      DO 60 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
C       FACT(K) = ABS(FACT(K))
      IF (WORK(K).GT.FCT) FCT = WORK(K)
   60   CONTINUE
      WORK(NE64+2*N+J) = FCT
      IF (FCT.NE.ZERO) THEN
      FCT = LOG(FCT)
      ELSE
C This can happen if only if column is null so matrix is singular.
      FCT = RINF/N
      ENDIF
      DO 70 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
C Note that zeros have been screened out so that FACT(K) always > 0.
      WORK(K) = FCT - LOG(WORK(K))
   70   CONTINUE
   75 CONTINUE
C Scale matrix
C B = DW(3N+1:3N+NE); IW(1:5N) and DW(1:2N) are workspaces
    CALL MC64W(N,NE64,IFACT(IPT),IFACT(IMAT),WORK,KEEP(PERM),NUM,
   &   IFACT(IDUP),IFACT(IMAT+NE64),IFACT(IMAT+NE64+N),
   &   IFACT(IMAT+NE64+2*N),IFACT(IMAT+NE64+3*N),
   &   WORK(NE64+1),WORK(NE64+N+1))
    IF (NUM.EQ.N) THEN
      DO 80 J = 1,N
      WORK(NE64+N+J) = WORK(NE64+N+J) - LOG(WORK(NE64+2*N+J))
   80   CONTINUE
C Check size of scaling factors
C   FCT = 0.5*LOG(RINF)
C   DO 86 J = 1,N
C     IF (FACT(NE64+J).LT.FCT .AND. FACT(NE64+N+J).LT.FCT) GO TO 86
C     INF64(1) = 2
C  86 CONTINUE

C Scaling is permuted to scaling on original matrix
      DO 5555 I=1,N
        FACT(ISCALE+PPOS(I)-1) =
   &      SQRT(EXP(WORK(NE64+I)+WORK(NE64+N+I)))
 5555   CONTINUE
      ELSE
C Matrix is singular
C Regenerate PERM and set PPOS to indicate nonsingular block
      K = 0
      DO 3501 I = 1,N
      IF (KEEP(PERM+I-1).LT.0) THEN
        PPOS(I) = -PPOS(I)
        IFACT(ISING+I-1) = 0
      ELSE
        K = K + 1
        IFACT(ISING+I-1) = K
      ENDIF
 3501   CONTINUE
      DO 3502 I = 1,N
      KEEP(PERM+ABS(PPOS(I))-1) = I
 3502   CONTINUE
C Copy matrix, remove duplicates, and initialize IP array.
      DO 3503 I = 1,N
      IFACT(IDUP+I-1) = 0
 3503   CONTINUE
C Must use new coordinates to keep matrix (half) symmetric
      IFACT(IPT) = 1
      KK = 1
      K = 1
      JNEW = 0
      NN = N
      DO 3505 J = 1,N
      IF (PPOS(J).LT.0) THEN
        NN = NN - 1
        K = K + KEEP(LROW+J-1)
        GO TO 3505
      ENDIF
      JNEW = JNEW + 1
      DO 3504 JJ = 1,KEEP(LROW+J-1)
        I = KEEP(PERM+IFACT(LIFACT-EXPNE+K)-1)
        IF (PPOS(I).GT.0) THEN
          IF (IFACT(IDUP+I-1).GE.IFACT(IPT+J-1)) THEN
C Duplicate
          FACT(IFACT(IDUP+I-1)) =
   &        FACT(IFACT(IDUP+I-1)) + FACT(LLFACT-EXPNE+K)
          ELSE
C Remove explicit zeros
          IF (FACT(LLFACT-EXPNE+K).NE.ZERO) THEN
            IFACT(IDUP+I-1) = KK
            FACT(KK) = FACT(LLFACT-EXPNE+K)
            IFACT(IMAT-1+KK) = IFACT(ISING+I-1)
            KK = KK+1
          ENDIF
          ENDIF
        ENDIF
        K = K + 1
 3504   CONTINUE
      IFACT(IPT+JNEW) = KK
 3505   CONTINUE
      DO 3500 I=1,KK-1
      WORK(I) = ABS(FACT(I))
 3500   CONTINUE
    NE64 = IFACT(IPT+NN)-1
      CALL MC34A(NN,IFACT(IMAT),IFACT(IPT),.TRUE.,WORK,KEEP(PERM))
      NE64 = IFACT(IPT+NN)-1
      DO 3508 J = 1,NN
      FCT = ZERO
      DO 3506 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
C       FACT(K) = ABS(FACT(K))
        IF (WORK(K).GT.FCT) FCT = WORK(K)
 3506   CONTINUE
      WORK(NE64+2*N+J) = FCT
      FCT = LOG(FCT)
      DO 3507 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
        WORK(K) = FCT - LOG(WORK(K))
 3507   CONTINUE
 3508   CONTINUE
      CALL MC64W(NN,NE64,IFACT(IPT),IFACT(IMAT),WORK,KEEP(PERM),NUM,
   &    IFACT(IDUP),IFACT(IMAT+NE64),IFACT(IMAT+NE64+N),
   &    IFACT(IMAT+NE64+2*N),IFACT(IMAT+NE64+3*N),
   &    WORK(NE64+1),WORK(NE64+N+1))
      DO 3509 J = 1,NN
      WORK(NE64+N+J) = WORK(NE64+N+J) - LOG(WORK(NE64+2*N+J))
 3509   CONTINUE
C Check size of scaling factors
C   FCT = 0.5*LOG(RINF)
C   DO 86 J = 1,N
C     IF (FACT(NE64+J).LT.FCT .AND. FACT(NE64+N+J).LT.FCT) GO TO 86
C     INF64(1) = 2
C  86 CONTINUE

C Scaling is permuted to scaling on original matrix for scale factors
C for nonsingular block
      K=0
C Loop is on new indices
      DO 3510 I=1,N
        IF (PPOS(I).LT.0) THEN
          K = K + 1
          FACT(ISCALE-PPOS(I)-1) = ZERO
        ELSE
          FACT(ISCALE+PPOS(I)-1) =
   &      SQRT(EXP(WORK(NE64+I-K)+WORK(NE64+N+I-K)))
        ENDIF
 3510   CONTINUE
C Compute scaling on nonsingular part
C Remember that PPOS maps from new to original but is flag on new
      DO 3516 I = 1,N
        KEEP(PERM+ABS(PPOS(I))-1) = I
 3516   CONTINUE
C Looping on new indices
      K = 1
      DO 3514 JJ = 1,N
        J = PPOS(JJ)
        IF (J.GT.0) THEN
          DO 3511 JLOOP = 1,KEEP(LROW+JJ-1)
          I = IFACT(LIFACT-EXPNE+K)
C I is original index so have to map to new to do PPOS test
          INEW = KEEP(PERM+I-1)
          IF (PPOS(INEW).LT.0)
   &        FACT(ISCALE+I-1) = MAX(ABS(FACT(ISCALE+I-1)),
   &           ABS(ABS(FACT(LLFACT-EXPNE+K))*FACT(ISCALE+J-1)))
          K = K + 1
 3511       CONTINUE
        ELSE
          DO 3512 JLOOP = 1,KEEP(LROW+JJ-1)
          I = IFACT(LIFACT-EXPNE+K)
C I is original index so have to map to new to do PPOS test
          INEW = KEEP(PERM+I-1)
C Shouldn't happen otherwise nonsingular block not maximum
C Sorry can happen but entry is implicit zero on diagonal
C Note that J is negative
          IF (I .NE. -J)  THEN
          FACT(ISCALE-J-1) =
   &          MAX(ABS(FACT(ISCALE-J-1)),
   &          ABS(ABS(FACT(LLFACT-EXPNE+K))*FACT(ISCALE+I-1)))
          ENDIF
          K = K + 1
 3512       CONTINUE
        ENDIF
 3514   CONTINUE
C Set scaling factors for singular block and reset PPOS
      DO 3513 I = 1,N
        INEW = KEEP(PERM+I-1)
        IF (PPOS(INEW) .LT. 0) THEN
          PPOS(INEW) = - PPOS(INEW)
          IF (FACT(ISCALE+I-1) .EQ. ZERO) THEN
          FACT(ISCALE+I-1) = ONE
          ELSE
          FACT(ISCALE+I-1) = ONE/FACT(ISCALE+I-1)
          ENDIF
        ENDIF
 3513   CONTINUE
      ENDIF
C End of logic for singular matrix

C       DO 8888 I = 1, N
C       FACT(ISCALE+I-1) = ONE
C8888   CONTINUE
      SMAX = ABS(FACT(ISCALE))
      SMIN = ABS(FACT(ISCALE))
      DO 5566 I = 1,N
        SMAX = MAX(SMAX,ABS(FACT(ISCALE+I-1)))
        SMIN = MIN(SMIN,ABS(FACT(ISCALE+I-1)))
 5566   CONTINUE
      RINFO(16) = SMIN
      RINFO(17) = SMAX
C Scale matrix
      K = 1
      BIGA = ZERO
      DO 6666 JJ = 1,N
        J = PPOS(JJ)
        DO 7777 JLOOP = 1,KEEP(LROW+JJ-1)
          I = IFACT(LIFACT-EXPNE+K)
          FACT(LLFACT-EXPNE+K) =
   &      FACT(ISCALE+I-1)*FACT(LLFACT-EXPNE+K)*FACT(ISCALE+J-1)
          BIGA = MAX(BIGA, ABS(FACT(LLFACT-EXPNE+K)))
          K = K + 1
 7777     CONTINUE
 6666   CONTINUE
CPRINT
C Only for testing will remove later
C       DO 6661 I = 1,EXPNE
C       IF (ABS(FACT(LLFACT-EXPNE+I)).GT.ONE+100*REPS) THEN
C         write(LP,*) '**** Matrix not scaled properly'
C         write(LP,*) 'Entry ',I,' has value ',FACT(LLFACT-EXPNE+I)
C         GO TO 6663
C       ENDIF
C6661   CONTINUE
C6663   CONTINUE
C End of scaling
    ELSE
      RINFO(16) = ONE
      RINFO(17) = ONE
    ENDIF
C
C**********************************
C Numerical Factorization
C**********************************
C Work arrays FACT(MM1/MM2), KEEP(PERM), IFACT(1)
    FACT(LFACT) = BIGA
   22 CALL ME57O(KIND, N, EXPNE, FACT, LLFACT, IFACT, LIFACT,
   *        KEEP(LROW), PPOS,
   *        NSTEPS, KEEP(NSTK), KEEP(NODE), FACT(MM1),
   *        FACT(MM2),
   *        KEEP(PERM),
   *        CNTL, ICNTL,
   *        INFO, RINFO, KEEP(HOLD), BIGA)
    IF (INFO(1).EQ.10 .OR. INFO(1).EQ.11) THEN
      IF (LDIAG.GT.2 .AND. MP.GE.0)  THEN
      IF (INFO(1).EQ.10) WRITE (MP,99982) INFO(1)
99982 FORMAT (/'Leaving factorization phase (ME57B) with ...'/
   1  'Factorization suspended because of lack of real space'/
   1  'INFO (1) = ',I3)
      IF (INFO(1).EQ.11) WRITE (MP,99983) INFO(1)
99983 FORMAT (/'Leaving factorization phase (ME57B) with ...'/
   1  'Factorization suspended because of lack of integer space'/
   1  'INFO (1) = ',I3)
      ENDIF
      RETURN
    ENDIF
C Regenerate array PERM
    DO 24 I = 1,N
      KEEP(PERM+PPOS(I)-1) = I
   24 CONTINUE
C Compute INFO(17-20)
      INFO(17) = ALENB + INFO(17)
      INFO(19) = ALENB + INFO(19)
C Allow space for scaling in ME57B/BD
    IF (ICNTL(15).EQ.1) THEN
      INFO(17) = MAX(INFO(17),ALENB + 3*EXPNE+3*N)
      INFO(19) = MAX(INFO(19),ALENB + 3*EXPNE+3*N)
      INFO(18) = MAX(INFO(18),3*EXPNE+5*N+1)
      INFO(20) = MAX(INFO(20),3*EXPNE+5*N+1)
    ENDIF
    IF (INFO(1).EQ.-3) GO TO 85
    IF (INFO(1).EQ.-4) GO TO 95
    IF (INFO(1).LT.0) RETURN
    GO TO 100
C************************
C **** Error returns ****
C************************
   25 INFO(1) = -1
    INFO(2) =  N
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
   +    '**** Error return from ME57B ****  INFO(1) =',INFO(1),
   +    'N has value ',INFO(2)
    RETURN
   30 INFO(1) = -2
    INFO(2) = NE
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
   +    '**** Error return from ME57B ****  INFO(1) =',INFO(1),
   +    'NE has value',INFO(2)
    RETURN
   40 INFO(1) = -15
    INFO(2) = LKEEP
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10/A,I10)')
   +    '**** Error return from ME57B ****  INFO(1) =',INFO(1),
   +    'LKEEP has value    ',INFO(2),
   +    'Should be at least ',5*N+NE+MAX(N,NE)+42
    RETURN
   35 INFO(1) = -10
    INFO(2) = ICNTL(7)
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
   +    '**** Error return from ME57B ****  INFO(1) =',INFO(1),
   +    'ICNTL(7) has value',ICNTL(7)
    RETURN

   85 INFO(1) = -3
    INFO(2) = LFACT
    INFO(17) = MAX(INFO(17), ALENB + EXPNE + 1)
    IF (ICNTL(15).EQ.1)
   *    INFO(17) = MAX(INFO(17), ALENB + 3*EXPNE + 3*N)
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
   +    '**** Error return from ME57B ****  INFO(1) =',INFO(1),
   +    'Insufficient real space in FACT, LFACT = ',INFO(2)
    RETURN

   95 INFO(1) = -4
    INFO(2) = LIFACT
    INFO(18) = MAX(INFO(18), EXPNE+N+5)
    IF (ICNTL(15).EQ.1)
   *    INFO(18) = MAX(INFO(18), 3*EXPNE + 5*N + 1)
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
   +    '**** Error return from ME57B ****  INFO(1) =',INFO(1),
   +    'Insufficient integer space in IFACT, LIFACT = ',INFO(2)
    RETURN

C****************
C Printing section
C****************
 100  IF (LDIAG.LE.2 .OR. MP.LT.0) RETURN
    WRITE (MP,99980) INFO(1), INFO(2),
   *    (INFO(I),I=14,25),INFO(28),INFO(29)
    WRITE (MP,99984) (INFO(I),I=31,35),RINFO(3), RINFO(4),
   *           RINFO(5), RINFO(18)
99980 FORMAT (/'Leaving factorization phase (ME57B) with ...'/
   1  'INFO (1)                        =',I12/
   2  ' --- (2)                        =',I12/
   3  ' --- (14) Number of entries in factors      =',I12/
   4  ' --- (15) Real storage for factors        =',I12/
   5  ' --- (16) Integer storage for factors       =',I12/
   6  ' --- (17) Min LFACT with compresses       =',I12/
   7  ' --- (18) Min LIFACT with compresses      =',I12/
   8  ' --- (19) Min LFACT without compresses      =',I12/
   9  ' --- (20) Min LIFACT without compresses     =',I12/
   *  ' --- (21) Order of largest frontal matrix   =',I12/
   1  ' --- (22) Number of 2x2 pivots          =',I12/
   2  ' --- (23) Number of delayed pivots        =',I12/
   3  ' --- (24) Number of negative eigenvalues    =',I12/
   4  ' --- (25) Rank of factorization         =',I12/
   5  ' --- (28) Number compresses on real data    =',I12/
   6  ' --- (29) Number compresses on integer data   =',I12)
    IF (ICNTL(15).EQ.1) WRITE (MP,99985) RINFO(16),RINFO(17)
99985 FORMAT (
   1  'RINFO(16) Minimum value of scaling factor   =  ',1PE10.3/
   2  '-----(17) Maximum value of scaling factor   =  ',1PE10.3)
99984 FORMAT (
   7  ' --- (31) Number of block pivots in factors   =',I12/
   7  ' --- (32) Number of zeros factors triangle    =',I12/
   7  ' --- (33) Number of zeros factors rectangle   =',I12/
   7  ' --- (34) Number of zero cols factors rect    =',I12/
   7  ' --- (35) Number of static pivots         =',I12/
   1  'RINFO(3)  Operations during node assembly   =  ',1PE10.3/
   2  '-----(4)  Operations during node elimination  =  ',1PE10.3/
   3  '-----(5)  Extra operations because of BLAS    =  ',1PE10.3/
   3  '-----(18) Largest modulus of entry in matrix  =  ',1PE10.3)
    IF (INFO(27).GT.0) WRITE (MP,99981) INFO(27),RINFO(14),RINFO(15)
99981 FORMAT (/'Matrix modification performed'/
   1  'INFO (27) Step at which matrix first modified =',I12/
   2  'RINFO(14) Maximum value added to diagonal   =  ',1PE10.3/
   2  'RINFO(15) Smallest pivot in modified matrix   =  ',1PE10.3)
C Print out matrix factors from ME57B.
    CALL ME57U(FACT,LFACT,IFACT,LIFACT,ICNTL)
C Print scaling factors
    IF (ICNTL(15).NE.1) RETURN
    K = MIN(10,N)
    IF (LDIAG.GE.4) K = N
    WRITE (MP,'(/A/(5E12.5))')  'Scaling factors:',
   +            (REAL(FACT(ISCALE+I-1)),I=1,K)
    IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'

    END
    SUBROUTINE ME57C(KIND,JOB,N,FACT,LFACT,IFACT,LIFACT,
   *            NRHS,RHS,LRHS,W,
   *            LW,IW1,ICNTL,INFO)
C This subroutine uses the factorisation of the matrix in FACT,IFACT to
C   solve a system of equations.
    INTEGER KIND,JOB,N,LFACT
    COMPLEX FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
    COMPLEX W(LW),RHS(LRHS,NRHS)
    INTEGER IW1(N),ICNTL(20),INFO(40)
CCC
C   LOGICAL LCASE
C   COMMON /CCASE/LCASE
C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C JOB must be set by the user to determine the coefficient matrix
C   of the equations being solved.  If the factorization is
C                  T  T
C         A =  P  L  D  L  P
C   then coefficient matrix is:
C
C   JOB <= 1   A
C               T
C   JOB  = 2   P  L  P
C               T
C   JOB  = 3   P  D  P
C            T  T
C   JOB >= 4   P  L  P
C
C N must be set to the order of the matrix and must be unchanged since
C   the call to ME57B. It is not altered.
C FACT holds information on the factors and must be unchanged since
C   the call to ME57B. It is not altered by ME57C.
C LFACT must be set to the length of FACT. It is not altered.
C IFACT holds information on the factors and must be unchanged since
C   the call to ME57B. It is not altered by ME57C.
C LIFACT must be set to the length of IFACT. It is not altered.
C NRHS is the number of right-hand sides being solved for.
C RHS must be set to the right hand sides for the equations being
C   solved. On exit, this array will hold the solutions.
C LHS must be set to the leading dimension of array RHS.
C W is used as a work array.
C LW  must be set to the length of array W.  It must be at least
C   as large as N*NRHS.  (Actually only INFO(21)*NRHS but no way to
C   check this).
C IW1 is used as a work array.
C ICNTL must be set by the user as follows and is not altered.
C   ICNTL(1)  must be set to the stream number for error messages.
C     A value less than 0 suppresses output.
C   ICNTL(2) must be set to the stream number for warning output.
C     A value less than 0 suppresses output.
C   ICNTL(3) must be set to the stream number for monitor output.
C     A value less than 0 suppresses output.
C   ICNTL(4) must be set to the stream number for statistics output.
C     A value less than 0 suppresses output.
C   ICNTL(5) must be set to control the amount of output:
C     0 None.
C     1 Error messages only.
C     2 Error and warning messages.
C     3 As 2, plus scalar parameters and a few entries of array
C       parameters on entry and exit.
C   > 3  As 2, plus all parameters on entry and exit.
C   ICNTL(6:12) Not referenced.
C   ICNTL(13) Threshold on number of columns in a block for direct
C       addressing using Level 2 and Level 3 BLAS.

C Procedures
    INTRINSIC MIN
    EXTERNAL ME57Q,ME57R,ME57S,ME57T,ME57U,ME57X,ME57Y

C
C Local variables
    REAL SCALE,ONE
    PARAMETER (ONE = 1.0E0)
    INTEGER I,J,K,LDIAG,LLW,LP,MP,ISCALE
C I  Temporary variable.
C J  Temporary variable.
C K  Temporary variable.
C LDIAG Control for amount of information output.
C LP Stream number for error printing.
C MP Stream number for monitor printing.
C

C Set local print control variables
    LP = ICNTL(1)
    MP = ICNTL(3)
    LDIAG = ICNTL(5)

    INFO(1) = 0

C Check input data
    IF (N.LE.0) THEN
      INFO(1) = -1
      INFO(2) = N
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
   +    '**** Error return from ME57C ****  INFO(1) =',INFO(1),
   +    'N has value',N
      GOTO 500
    ENDIF

    IF (NRHS.LT.1) THEN
      INFO(1) = -16
      INFO(2) = NRHS
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I10,A)')
   +    '**** Error return from ME57C ****  INFO(1) =',INFO(1),
   +    'value of NRHS =',NRHS,' is less than 1'
      GOTO 500
    ENDIF

    IF (LRHS.LT.N) THEN
      INFO(1) = -11
      INFO(2) = LRHS
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I10,A,I10)')
   +    '**** Error return from ME57C ****  INFO(1) =',INFO(1),
   +    'value of LRHS =',LRHS,' is less than N=',N
      GOTO 500
    ENDIF

    IF (LW.LT.N*NRHS) THEN
      INFO(1) = -17
      INFO(2) = N*NRHS
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I10,A,I10)')
   +    '**** Error return from ME57C ****  INFO(1) =',INFO(1),
   +    'value of LW =',LW,' is less than', N*NRHS
      GOTO 500
    ENDIF

C If requested, print input parameters
    IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
      WRITE (MP,'(/A/)') 'Entering solution phase (ME57C) with ...'
      IF (KIND.EQ.1) THEN
      WRITE(MP,'(A)') 'Matrix is complex symmetric'
      ELSE
      WRITE(MP,'(A)') 'Matrix is Hermitian'
      ENDIF
      WRITE (MP,99999) JOB,N,(ICNTL(I),I=1,5),LFACT,LIFACT,NRHS,
   +       LRHS,LW,ICNTL(13)
99999 FORMAT(
   +    'JOB     Control on coefficient matrix     =',I12/
   +    'N       Order of matrix             =',I12/
   6    'ICNTL(1)  Stream for errors             =',I12/
   7    ' --- (2)  Stream for warnings           =',I12/
   8    ' --- (3)  Stream for monitoring         =',I12/
   9    ' --- (4)  Stream for statistics         =',I12/
   1    ' --- (5)  Level of diagnostic printing      =',I12/
   +    'LFACT   Length of array FACT          =',I12/
   +    'LIFACT    Length of array IFACT         =',I12/
   +    'NRHS    Number of right-hand sides      =',I12/
   +    'LRHS    Leading dimension of RHS array    =',I12/
   +    'LW      Leading dimension of work array   =',I12/
   +    'ICNTL(13) Threshold for Level 2 and 3 BLAS    =',I12)

C Print out matrix factors.
      CALL ME57U(FACT,LFACT,IFACT,LIFACT,ICNTL)
C Print scaling factors
      IF (ICNTL(15).EQ.1) THEN
      ISCALE = LFACT-N
      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(/A/(5E12.5))')  'Scaling factors:',
   +                (REAL(FACT(ISCALE+I-1)),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
      ENDIF
      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      DO 10 J = 1,NRHS
      WRITE(MP,'(/A,I10)') 'Right-hand side',J
      WRITE (MP,'((1P,6E13.3))') (RHS(I,J),I=1,K)
      IF (K.LT.N) WRITE (MP,'(A)') '   . . .'
   10   CONTINUE
    END IF

    LLW = LW/NRHS


C Scale right-hand side
    IF (ICNTL(15).EQ.1) THEN
      ISCALE = LFACT-N
C   write(7,*) 'N,NRHS,ISCALE', N,NRHS,ISCALE
      DO 5555 I = 1, N
      SCALE = REAL(FACT(ISCALE+I-1))
C   write(7,*) 'SCALE',SCALE
      IF (JOB.GE.4) SCALE = ONE/SCALE
      DO 4444 J = 1, NRHS
C   write(7,*) 'I,J,RHS', I,J,RHS(I,J)
        RHS(I,J) = SCALE*RHS(I,J)
C   write(7,*) 'I,J,RHS', I,J,RHS(I,J)
 4444   CONTINUE
 5555   CONTINUE
    ENDIF

C Forward substitution
    IF (JOB.LE.2) THEN
      IF (NRHS.EQ.1) THEN
      CALL ME57X(N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
   *          W,LLW,IW1,ICNTL)
CCC
C     IF (LCASE) THEN
C       DO 201 J = 1,NRHS
C       WRITE(7,'(/A,I10)') 'Solution     ',J
C       WRITE (7,'(1P,6E13.3)') (RHS(I,J),I=1,N)
C 201   CONTINUE
C     ENDIF
      ELSE
      CALL ME57Q(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
   *          W,LLW,IW1,ICNTL)
CCC
C     IF (LCASE) THEN
C       DO 202 J = 1,NRHS
C       WRITE(7,'(/A,I10)') 'Solution     ',J
C       WRITE (7,'(1P,6E13.3)') (RHS(I,J),I=1,N)
C 202   CONTINUE
C     ENDIF
      ENDIF
      IF (JOB.EQ.2) GO TO 15
C Back substitution.
      IF (NRHS.EQ.1) THEN
      CALL ME57Y(KIND,N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
   *          W,LLW,IW1,ICNTL)
      ELSE
      CALL ME57R(KIND,N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
   *          W,LLW,IW1,ICNTL)
      ENDIF
    ENDIF
    IF (JOB.EQ.3)
   *  CALL ME57S(KIND,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
   *          W,LLW,ICNTL)
    IF (JOB.GE.4)
   *  CALL ME57T(KIND,N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
   *          W,LLW,IW1,ICNTL)

C Scale solution
   15 IF (ICNTL(15).EQ.1) THEN
      ISCALE = LFACT-N
      DO 6666 I = 1, N
      SCALE = REAL(FACT(ISCALE+I-1))
      IF (JOB.EQ.2) SCALE = ONE/SCALE
      DO 7777 J = 1, NRHS
        RHS(I,J) = SCALE*RHS(I,J)
 7777   CONTINUE
 6666   CONTINUE
    ENDIF

C
C If requested, print output parameters.
    IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
      WRITE (MP,'(//A)')
   *     'Leaving solution phase (ME57C) with ...'
      DO 20 J = 1,NRHS
      WRITE(MP,'(/A,I10)') 'Solution     ',J
      WRITE (MP,'(1P,6E13.3)') (RHS(I,J),I=1,K)
      IF (K.LT.N) WRITE (MP,'(A)') '   . . .'
   20   CONTINUE
    ENDIF

  500 RETURN

    END


    SUBROUTINE ME57Q(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
   *            W,LW,IW1,ICNTL)
C This subroutine performs forward elimination using the factors
C   stored in FACT/IFACT by ME57B.
C   Because of how the factors are stored this subroutine is the
C   same for complex symmetric and Hermitian matrices.
    INTEGER N,LFACT
    COMPLEX FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
    COMPLEX W(LW,NRHS),RHS(LRHS,NRHS)
    INTEGER IW1(N),ICNTL(20)
C N   must be set to the order of the matrix. It is not altered.
C FACT   must be set to hold the real values corresponding to the
C   factors. This must be unchanged since the preceding call to
C   ME57B. It is not altered.
C LFACT  length of array FACT. It is not altered.
C IFACT  holds the integer indexing information for the matrix factors
C   in FACT. This must be unchanged since the preceding call to
C   ME57B. It is not altered.
C LIFACT length of array IFACT. It is not altered.
C NRHS must be set to number of right-hand sides.
C RHS on input, must be set to hold the right hand side vector.  On
C   return, it will hold the modified vector following forward
C   elimination.
C LHS must be set to the leading dimension of array RHS.
C W   used as workspace to hold the components of the right hand
C   sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C   larger than INFO(21) as returned from ME57B.
C IW1 need not be set on entry. On exit IW1(I) (I = 1,NBLK), where
C   NBLK = IFACT(3) is the number of block pivots, will
C   hold pointers to the beginning of each block pivot in array IFACT.
C ICNTL Not referenced except:
C   ICNTL(13) Threshold on number of columns in a block for using
C       addressing using Level 2 and Level 3 BLAS.
C

C Procedures
    INTRINSIC ABS
    EXTERNAL CGEMM,CTPSV

C Constant
    REAL ONE
    PARAMETER (ONE=1.0E0)
    COMPLEX CONE
    PARAMETER (CONE=(1.0E0,0.0E0))
C
C Local variables
    INTEGER APOS,I,IBLK,II,IPIV,IRHS,IWPOS,J,J1,J2,K,
   +      NCOLS,NROWS
    COMPLEX W1
C
C APOS  Current position in array FACT.
C I   Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IWPOS Position in IFACT of start of current index list.
C J   Temporary DO index
C K   Temporary pointer to position in real array.
C J1    Position in IFACT of index of leading entry of row.
C J2    Position in IFACT of index of trailing entry of row.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.


    APOS = 1
    IWPOS = 4
    DO 270 IBLK = 1,IFACT(3)

C Find the number of rows and columns in the block.
      IW1(IBLK) = IWPOS
      NCOLS = IFACT(IWPOS)
      NROWS = IFACT(IWPOS+1)
      IWPOS = IWPOS + 2

      IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN


C Perform operations using direct addressing.

C Load appropriate components of right-hand sides into W.
      DO 10 I = 1,NCOLS
        II = ABS(IFACT(IWPOS+I-1))
        DO 11 J = 1,NRHS
          W(I,J) = RHS(II,J)
   11     CONTINUE
   10   CONTINUE


C Treat diagonal block (direct addressing)
      DO 12 J = 1,NRHS
        CALL CTPSV('L','N','U',NROWS,FACT(APOS),W(1,J),1)
   12   CONTINUE
      APOS = APOS + (NROWS* (NROWS+1))/2

C Treat off-diagonal block (direct addressing)
      IF (NCOLS.GT.NROWS) CALL CGEMM('N','N',NCOLS-NROWS,NRHS,NROWS,
   +                      CONE,FACT(APOS),NCOLS-NROWS,
   +                      W,LW,CONE,W(NROWS+1,1),LW)
      APOS = APOS + NROWS* (NCOLS-NROWS)

C Reload W back into RHS.
      DO 35 I = 1,NCOLS
        II = ABS(IFACT(IWPOS+I-1))
        DO 36 J = 1,NRHS
          RHS(II,J) = W(I,J)
   36     CONTINUE
   35   CONTINUE

      ELSE

C Perform operations using indirect addressing.

      J1 = IWPOS
      J2 = IWPOS + NROWS - 1


C Treat diagonal block (indirect addressing)
      DO 130 IPIV = 1,NROWS
      APOS = APOS + 1
      DO 101 II = 1,NRHS
        W1 = RHS(ABS(IFACT(J1)),II)
        K = APOS
        DO 100 J = J1+1,J2
          IRHS = ABS(IFACT(J))
          RHS(IRHS,II) = RHS(IRHS,II) - FACT(K)*W1
          K = K + 1
  100     CONTINUE
  101   CONTINUE
      APOS = K
      J1 = J1 + 1
  130   CONTINUE

C Treat off-diagonal block (indirect addressing)
      J2 = IWPOS + NCOLS - 1
      DO 136 IPIV = 1,NROWS
      DO 135 II = 1,NRHS
        K = APOS
        W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)),II)
        DO 133 J = J1,J2
          IRHS = ABS(IFACT(J))
          RHS(IRHS,II) = RHS(IRHS,II) + W1*FACT(K)
          K = K + 1
  133     CONTINUE
  135   CONTINUE
      APOS = K
  136   CONTINUE

    END IF

    IWPOS = IWPOS + NCOLS
  270 CONTINUE

    END


    SUBROUTINE ME57R(KIND,N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
   *            W,LW,IW1,ICNTL)
C This subroutine performs backward elimination operations
C   using the factors stored in FACT/IFACT by ME57B.
    INTEGER KIND,N,LFACT
    COMPLEX FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
    COMPLEX W(LW,NRHS),RHS(LRHS,NRHS)
    INTEGER IW1(N),ICNTL(20)
C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C N    must be set to the order of the matrix. It is not altered.
C FACT    must be set to hold the real values corresponding to the
C    factors. This must be unchanged since the
C    preceding call to ME57B. It is not altered.
C LFACT   length of array FACT. It is not altered.
C IFACT   holds the integer indexing information for the matrix factors
C    in FACT. This must be unchanged since the preceding call to
C    ME57B. It is not altered.
C LIFACT  length of array IFACT. It is not altered.
C NRHS must be set to number of right-hand sides.
C RHS  on entry, must be set to hold the right hand side modified by
C    the forward substitution operations. On exit, holds the
C    solution vector.
C LHS must be set to the leading dimension of array RHS.
C W    used as workspace to hold the components of the right hand
C    sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C   larger than INFO(21) as returned from ME57B.
C IW1  on entry IW1(I) (I = 1,NBLK), where  NBLK = IFACT(3) is the
C    number of block pivots, must hold pointers to the beginning of
C    each block pivot in array IFACT, as set by ME57Q/QD.
C    It is not altered.
C ICNTL Not referenced except:
C   ICNTL(13) Threshold on number of columns in a block for using
C       addressing using Level 2 and Level 3 BLAS.

C Procedures
    INTRINSIC ABS
    EXTERNAL CGEMM,CTPSV

C Constants
    REAL ONE
    PARAMETER (ONE=1.0E0)
    COMPLEX CONE
    PARAMETER (CONE=(1.0E0,0.0E0))
C
C Local variables.
    INTEGER APOS,APOS2,I,IBLK,II,IPIV,IRHS,IRHS1,
   +      IRHS2,IWPOS,J,JPIV,J1,J2,K,KK,LROW,NCOLS,NROWS
    COMPLEX W1
C APOS  Current position in array FACT.
C APOS2 Current position in array FACT for off-diagonal entry of 2x2
C     pivot.
C I   Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IRHS1 RHS index.
C IRHS2 RHS index.
C IWPOS Position in IFACT of start of current index list.
C J   Temporary DO index
C JPIV  Has the value -1 for the first row of a 2 by 2 pivot and 1 for
C     the second.
C K   Temporary pointer to position in real array.
C J1    Position in IFACT of index of leading entry of row.
C J2    Position in IFACT of index of trailing entry of row.
C K   Temporary variable.
C LROW  Length of current row.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.
C

    APOS = IFACT(1)
    APOS2 = IFACT(2)
C Run through block pivot rows in the reverse order.
    DO 380 IBLK = IFACT(3),1,-1

C Find the number of rows and columns in the block.
      IWPOS = IW1(IBLK)
      NCOLS = ABS(IFACT(IWPOS))
      NROWS = ABS(IFACT(IWPOS+1))
      APOS = APOS - NROWS* (NCOLS-NROWS)
      IWPOS = IWPOS + 2

      IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN


C Perform operations using direct addressing.

C Load latter part of right-hand side into W.
      DO 5 I = NROWS + 1,NCOLS
        II = ABS(IFACT(IWPOS+I-1))
        DO 3 J = 1,NRHS
          W(I,J) = RHS(II,J)
    3     CONTINUE
    5   CONTINUE


C Multiply by the diagonal matrix (direct addressing)
      DO 10 IPIV = NROWS,1,-1
        IRHS = ABS(IFACT(IWPOS+IPIV-1))
        APOS = APOS - (NROWS+1-IPIV)
        DO 9 J = 1,NRHS
          W(IPIV,J) = RHS(IRHS,J)*FACT(APOS)
    9     CONTINUE
   10   CONTINUE
      JPIV = -1
      DO 20 IPIV = NROWS,1,-1
        IRHS = IFACT(IWPOS+IPIV-1)
        IF (IRHS.LT.0) THEN
          IRHS1 = -IFACT(IWPOS+IPIV-1+JPIV)
          IF (KIND.EQ.1 .OR. JPIV.EQ.1) THEN
          DO 19 J = 1,NRHS
            W(IPIV,J) = RHS(IRHS1,J)*FACT(APOS2) + W(IPIV,J)
   19       CONTINUE
          ELSE
          DO 21 J = 1,NRHS
            W(IPIV,J) = RHS(IRHS1,J)*CONJG(FACT(APOS2))
   *                + W(IPIV,J)
   21       CONTINUE
          ENDIF
          IF (JPIV.EQ.1) APOS2 = APOS2 - 1
          JPIV = -JPIV
        END IF

   20   CONTINUE

C Treat off-diagonal block (direct addressing)
      K = NCOLS - NROWS
      IF (K.GT.0) THEN
        IF (KIND.EQ.1) CALL CGEMM('T','N',NROWS,NRHS,K,CONE,
   +                 FACT(APOS+(NROWS*(NROWS+1))/2),K,
   +                 W(NROWS+1,1),LW,CONE,W,LW)
        IF (KIND.NE.1) CALL CGEMM('C','N',NROWS,NRHS,K,CONE,
   +                 FACT(APOS+(NROWS*(NROWS+1))/2),K,
   +                 W(NROWS+1,1),LW,CONE,W,LW)
      ENDIF

C Treat diagonal block (direct addressing)
      IF (KIND.EQ.1) THEN
        DO 22 J = 1,NRHS
          CALL CTPSV('L','T','U',NROWS,FACT(APOS),W(1,J),1)
   22     CONTINUE
      ELSE
        DO 23 J = 1,NRHS
          CALL CTPSV('L','C','U',NROWS,FACT(APOS),W(1,J),1)
   23     CONTINUE
      ENDIF
C Reload W back into RHS.
      DO 60 I = 1,NROWS
        II = ABS(IFACT(IWPOS+I-1))
        DO 59 J = 1,NRHS
          RHS(II,J) = W(I,J)
   59     CONTINUE
   60   CONTINUE

      ELSE
C
C Perform operations using indirect addressing.
      J1 = IWPOS
      J2 = IWPOS + NCOLS - 1


C Multiply by the diagonal matrix (indirect addressing)
      JPIV = -1
      IF (KIND.EQ.1) THEN
        DO 210 IPIV = NROWS,1,-1
          IRHS = IFACT(IWPOS+IPIV-1)
          LROW = NROWS + 1 - IPIV
          IF (IRHS.GT.0) THEN
C 1 by 1 pivot.
          APOS = APOS - LROW
          DO 65 J = 1,NRHS
            RHS(IRHS,J) = RHS(IRHS,J)*FACT(APOS)
   65       CONTINUE
          ELSE
C 2 by 2 pivot
          IF (JPIV.EQ.-1) THEN
            IRHS1 = -IFACT(IWPOS+IPIV-2)
            IRHS2 = -IRHS
            APOS = APOS - LROW - LROW - 1
            DO 68 J = 1,NRHS
            W1 = RHS(IRHS1,J)*FACT(APOS) +
   +             RHS(IRHS2,J)*FACT(APOS2)
            RHS(IRHS2,J) = RHS(IRHS1,J)*FACT(APOS2) +
   +                   RHS(IRHS2,J)*FACT(APOS+LROW+1)
            RHS(IRHS1,J) = W1
   68         CONTINUE
            APOS2 = APOS2 - 1
          END IF
          JPIV = -JPIV
          END IF
  210     CONTINUE
      ELSE
        DO 211 IPIV = NROWS,1,-1
          IRHS = IFACT(IWPOS+IPIV-1)
          LROW = NROWS + 1 - IPIV
          IF (IRHS.GT.0) THEN
C 1 by 1 pivot.
          APOS = APOS - LROW
          DO 66 J = 1,NRHS
            RHS(IRHS,J) = RHS(IRHS,J)*FACT(APOS)
   66       CONTINUE
          ELSE
C 2 by 2 pivot
          IF (JPIV.EQ.-1) THEN
            IRHS1 = -IFACT(IWPOS+IPIV-2)
            IRHS2 = -IRHS
            APOS = APOS - LROW - LROW - 1
            DO 69 J = 1,NRHS
            W1 = RHS(IRHS1,J)*FACT(APOS) +
   +             RHS(IRHS2,J)*FACT(APOS2)
            RHS(IRHS2,J) = RHS(IRHS1,J)*CONJG(FACT(APOS2)) +
   +                   RHS(IRHS2,J)*FACT(APOS+LROW+1)
            RHS(IRHS1,J) = W1
   69         CONTINUE
            APOS2 = APOS2 - 1
          END IF
          JPIV = -JPIV
          END IF
  211     CONTINUE
      ENDIF
      APOS = APOS + (NROWS* (NROWS+1))/2

C Treat off-diagonal block (indirect addressing)
      KK = APOS
      J1 = IWPOS + NROWS
      IF (KIND.EQ.1) THEN
        DO 220 IPIV = 1,NROWS
          IRHS = ABS(IFACT(IWPOS+IPIV-1))
          DO 218 II = 1,NRHS
          W1 = RHS(IRHS,II)
          K = KK
          DO 215 J = J1,J2
            W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)),II)
            K = K + 1
  215       CONTINUE
          RHS(IRHS,II) = W1
  218       CONTINUE
          KK = K
  220     CONTINUE
      ELSE
        DO 221 IPIV = 1,NROWS
          IRHS = ABS(IFACT(IWPOS+IPIV-1))
          DO 219 II = 1,NRHS
          W1 = RHS(IRHS,II)
          K = KK
          DO 216 J = J1,J2
            W1 = W1 + CONJG(FACT(K))*RHS(ABS(IFACT(J)),II)
            K = K + 1
  216       CONTINUE
          RHS(IRHS,II) = W1
  219       CONTINUE
          KK = K
  221     CONTINUE
      ENDIF

C Treat diagonal block (indirect addressing)
      J2 = IWPOS + NROWS - 1
      DO 260 IPIV = 1,NROWS
        IRHS = ABS(IFACT(J1-1))
        APOS = APOS - IPIV
        DO 240 II = 1,NRHS
          W1 = RHS(IRHS,II)
          K = APOS + 1
          DO 230 J = J1,J2
          IF (KIND.EQ.1) THEN
            W1 = W1 - FACT(K)*RHS(ABS(IFACT(J)),II)
          ELSE
            W1 = W1 - CONJG(FACT(K))*RHS(ABS(IFACT(J)),II)
          ENDIF
          K = K + 1
  230       CONTINUE
          RHS(IRHS,II) = W1
  240     CONTINUE
        J1 = J1 - 1
  260   CONTINUE

      END IF

  380 CONTINUE

    END

    SUBROUTINE ME57U(FACT,LFACT,IFACT,LIFACT,ICNTL)
C   Print out matrix factors from ME57B or a symbolic representation
C   of them.
    INTEGER LFACT
    COMPLEX FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT),ICNTL(20)
C FACT   array holding the reals of the factorization.
C      It is not altered.
C LFACT  length of array FACT. It is not altered.
C IFACT  array holding the integers of the factorization. It is not
C   altered.
C LIFACT length of array IFACT. It is not altered.
C ICNTL is not referenced except:
C   ICNTL(3) must be set to the stream number for diagnostic output.
C     A value less than 1 suppresses output.
C   ICNTL(5) must be set to control the amount of output:
C    <3 None.
C     3 First block only.
C     4 All blocks.
C     5 All blocks, but each entry represented by a single character:
C        + for a positive integer
C        - for a negative integer
C        * for a nonzero entry
C        . for a zero entry

C Procedures
    INTRINSIC MIN,SIGN

C Local variables
    COMPLEX CZERO
    PARAMETER (CZERO = (0.0E0,0.0E0))
    CHARACTER*80 LINE
    INTEGER APOS,APOS2,IBLK,ILINE,IROW,IWPOS,J,JPIV,J1,J2,K,
   +      LDIAG,LEN,MP,NBLK,NCOLS,NROWS
C APOS Current position in FACT.
C APOS2 Position in FACT of next off-diagonal entry of 2x2 pivot.
C ILINE Current position in the line.
C IBLK  Current block.
C IROW  Current row.
C IWPOS Current position in IFACT.
C JPIV  has value 1 only for the first row of a 2x2 pivot.
C J   Column index.
C K   Temporary pointer to position in real array.
C J1    Position of last zero in leading part of row.
C J2    Position of last nonzero in leading part of row.
C K   Temporary DO index.
C LDIAG Control for diagnostic printing.
C LEN   1 if pattern only to be printed and 12 if values to be printed.
C LINE  Character variable in which an output line is built.
C MP    Stream number for warning or diagnostic messages
C NBLK  Number of blocks to be printed.
C NCOLS Number of columns in the block.
C NROWS Number of rows in the block.

    CHARACTER*1 PM(-2:2)
    DATA PM/'*','-','.','+','.'/
    REAL ZERO
    PARAMETER (ZERO=0.0E0)

C Initialize MP and LDIAG
    MP = ICNTL(3)
    LDIAG = ICNTL(5)

    APOS2 = IFACT(1)
    NBLK = IFACT(3)
    IF (LDIAG.EQ.3) NBLK = MIN(1,NBLK)
    LEN = 20
    IF (LDIAG.EQ.5) LEN = 1
    IF (LEN.EQ.20) THEN
      IF (NBLK.EQ.IFACT(3)) THEN
      WRITE (MP,'(/A)')
   +    'For each block, the following information is provided:'

      ELSE
      WRITE (MP,'(/A,A)') 'For the first block only,',
   +    ' the following information is provided:'
      END IF

    END IF

    IF (LEN.EQ.20) WRITE (MP,'(A)')
   +    '   1. Block number, number of rows, number of columns',
   +    '   2. List of indices for the pivot, each negated if part of'
   +    ,'    a 2x2 pivot',
   +    '   3. The factorized block pivot',
   +    '    It has the form',
   +    '        -1  T',
   +    '      L  D   L ',
   +    '               -1    T',
   +    '    and is printed as D and L  packed together.',
   +    '   4. List of indices for the non-pivot columns',
   +    '   5. The non-pivot part as rectangular block by rows'

    IWPOS = 4
    APOS = 1

    DO 300 IBLK = 1,NBLK
      NCOLS = IFACT(IWPOS)
      NROWS = IFACT(IWPOS+1)
      IWPOS = IWPOS + 2

      WRITE (MP,'(/4(A,I6))') 'Block pivot',IBLK,' with',NROWS,
   +      ' rows and', NCOLS,' columns'
C      IF (LEN.EQ.20) WRITE (MP,'(6I12)')
C   +               (IFACT(K),K=IWPOS,IWPOS+NROWS-1)
C      IF (LEN.EQ.1) WRITE (MP,'(80A1)') (PM(SIGN(1,IFACT(K))),
C   +    K=IWPOS,IWPOS+NROWS-1)

      JPIV = 0
      DO 30 IROW = 1,NROWS
      IF (JPIV.EQ.1) THEN
        JPIV = 0
      ELSE
        IF (IFACT(IWPOS+IROW-1).LT.0) JPIV = 1
      END IF

      ILINE = 1
      DO 10 J = 1,IROW - 1
        WRITE (LINE(ILINE:ILINE+LEN-1),'(A)') ' '
        ILINE = ILINE + LEN
        IF (ILINE.GT.80) THEN
          WRITE (MP,'(A)') LINE
          ILINE = 1
        END IF
   10   CONTINUE

      DO 20 J = IROW,NROWS
        IF (LEN.EQ.20) WRITE (LINE(ILINE:ILINE+19),
   +      '(1P,2E10.2)') FACT(APOS)
        IF (LEN.EQ.1) THEN
         IF (FACT(APOS).EQ.CZERO) THEN
            WRITE (LINE(ILINE:ILINE),'(A)') '.'
         ELSE
            WRITE (LINE(ILINE:ILINE),'(A)') '*'
         END IF
        END IF
        APOS = APOS + 1
        IF (J.EQ.IROW+1) THEN
          IF (JPIV.EQ.1) THEN
          IF (LEN.EQ.20) WRITE (LINE(ILINE:ILINE+19),
   +          '(1P,2E10.2)') FACT(APOS2)
          IF (LEN.EQ.1) THEN
            IF (FACT(APOS2).EQ.CZERO) THEN
               WRITE (LINE(ILINE:ILINE),'(A)') '.'
            ELSE
               WRITE (LINE(ILINE:ILINE),'(A)') '*'
            END IF
          END IF
          APOS2 = APOS2 + 1
          END IF
        END IF
        ILINE = ILINE + LEN
        IF (ILINE.GT.80) THEN
          WRITE (MP,'(A)') LINE
          ILINE = 1
        END IF
   20   CONTINUE

      IF (ILINE.GT.1) THEN
        LINE(ILINE:) = ' '
        WRITE (MP,'(A)') LINE
      END IF

   30   CONTINUE

      IWPOS = IWPOS + NROWS
      IF (LEN.EQ.20) WRITE (MP,'(5I16)') (IFACT(K),K=IWPOS,
   +    IWPOS+NCOLS-NROWS-1)
      IF (LEN.EQ.1) WRITE (MP,'(80A1)') (PM(SIGN(1,IFACT(K))),
   +    K=IWPOS,IWPOS+NCOLS-NROWS-1)

      IWPOS = IWPOS + NCOLS - NROWS
      DO 280 IROW = 1,NROWS
      J1 = NROWS
      J2 = NCOLS
      ILINE = 1
      DO 110 J = J1 + 1,J2
        IF (LEN.EQ.20) WRITE (LINE(ILINE:ILINE+19),
   +      '(1P,2E10.2)') FACT(APOS)
        IF (LEN.EQ.1) THEN
         IF (FACT(APOS).EQ.CZERO) THEN
            WRITE (LINE(ILINE:ILINE),'(A)') '.'
         ELSE
            WRITE (LINE(ILINE:ILINE),'(A)') '*'
         END IF
        END IF
        APOS = APOS + 1
        ILINE = ILINE + LEN
        IF (ILINE.GT.80) THEN
          WRITE (MP,'(A)') LINE
          ILINE = 1
        END IF
  110   CONTINUE

      IF (ILINE.GT.1) THEN
        LINE(ILINE:) = ' '
        WRITE (MP,'(A)') LINE
      END IF

  280   CONTINUE
  300 CONTINUE
    END

    SUBROUTINE ME57S(KIND,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
   *            W,LW,ICNTL)
C This subroutine divides a vector by the block diagonal matrix of
C   the matrix factors using factor entries stored in FACT/IFACT
C    by ME57B.
    INTEGER KIND,LFACT
    COMPLEX FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
    COMPLEX W(LW,NRHS),RHS(LRHS,NRHS)
    INTEGER ICNTL(20)
C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C FACT    must be set to hold the real values corresponding to the
C    factors. This must be unchanged since the
C    preceding call to ME57B. It is not altered.
C LFACT   length of array FACT. It is not altered.
C IFACT   holds the integer indexing information for the matrix factors
C    in FACT. This must be unchanged since the preceding call to
C    ME57B. It is not altered.
C LIFACT  length of array IFACT. It is not altered.
C NRHS must be set to number of right-hand sides.
C RHS  on entry, must be set to hold the right hand side modified by
C    the forward substitution operations. On exit, holds the
C    solution vector.
C LHS must be set to the leading dimension of array RHS.
C W    used as workspace to hold the components of the right hand
C    sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C   larger than INFO(21) as returned from ME57B.
C ICNTL Not referenced except:
C   ICNTL(13) Threshold on number of columns in a block for direct
C       addressing using Level 2 and Level 3 BLAS.

C Procedures
    INTRINSIC ABS
    EXTERNAL CGEMM,CTPSV

C
C Local variables.
    INTEGER APOS,APOS2,I,IBLK,II,IPIV,IRHS,IRHS1,
   +      IRHS2,IWPOS,J,JPIV,NCOLS,NROWS
    COMPLEX W1
C APOS  Current position in array FACT.
C APOS2 Current position in array FACT for off-diagonal entry of 2x2
C     pivot.
C I   Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IRHS1 RHS index.
C IRHS2 RHS index.
C IWPOS Position in IFACT of start of current index list.
C J   Temporary DO index
C JPIV  Has the value 1 for the first row of a 2 by 2 pivot and -1 for
C     the second.
C K   Temporary pointer to position in real array.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.
C
    APOS = 1
    APOS2 = IFACT(1)
    IWPOS = 4
C Run through block pivot rows in the reverse order.
    DO 380 IBLK = 1,IFACT(3)

C Find the number of rows and columns in the block.
      NCOLS = IFACT(IWPOS)
      NROWS = IFACT(IWPOS+1)
      IWPOS = IWPOS + 2


      IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN

C Perform operations using direct addressing.

C Multiply by the diagonal matrix (direct addressing)
      DO 10 IPIV = 1,NROWS
        IRHS = ABS(IFACT(IWPOS+IPIV-1))
        DO 9 J = 1,NRHS
          W(IPIV,J) = RHS(IRHS,J)*FACT(APOS)
    9     CONTINUE
        APOS = APOS + (NROWS+1-IPIV)
   10   CONTINUE
      JPIV = 1
      DO 20 IPIV = 1,NROWS
        IRHS = IFACT(IWPOS+IPIV-1)
        IF (IRHS.LT.0) THEN
          IRHS1 = -IFACT(IWPOS+IPIV-1+JPIV)
          IF (KIND.EQ.1 .OR. JPIV.EQ.1) THEN
          DO 18 J = 1,NRHS
            W(IPIV,J) = RHS(IRHS1,J)*FACT(APOS2) + W(IPIV,J)
   18       CONTINUE
          ELSE
          DO 19 J = 1,NRHS
            W(IPIV,J) =
   *          RHS(IRHS1,J)*CONJG(FACT(APOS2)) + W(IPIV,J)
   19       CONTINUE
          ENDIF
          IF (JPIV.EQ.-1) APOS2 = APOS2 + 1
          JPIV = -JPIV
        END IF

   20   CONTINUE

C Reload W back into RHS.
      DO 60 I = 1,NROWS
        II = ABS(IFACT(IWPOS+I-1))
        DO 59 J = 1,NRHS
          RHS(II,J) = W(I,J)
   59     CONTINUE
   60   CONTINUE

      ELSE
C
C Perform operations using indirect addressing.

C Multiply by the diagonal matrix (indirect addressing)
      JPIV = 1
      IF (KIND.EQ.1) THEN
        DO 210 IPIV = 1,NROWS
          IRHS = IFACT(IWPOS+IPIV-1)

          IF (IRHS.GT.0) THEN
C 1 by 1 pivot.
          DO 65 J = 1,NRHS
            RHS(IRHS,J) = RHS(IRHS,J)*FACT(APOS)
   65       CONTINUE
          APOS = APOS + NROWS - IPIV + 1
          ELSE
C 2 by 2 pivot
          IF (JPIV.EQ.1) THEN
            IRHS1 = -IRHS
            IRHS2 = -IFACT(IWPOS+IPIV)
            DO 68 J = 1,NRHS
            W1 = RHS(IRHS1,J)*FACT(APOS) +
   +             RHS(IRHS2,J)*FACT(APOS2)
            RHS(IRHS2,J) = RHS(IRHS1,J)*FACT(APOS2) +
   +                   RHS(IRHS2,J)*FACT(APOS+NROWS-IPIV+1)
            RHS(IRHS1,J) = W1
   68         CONTINUE
            APOS2 = APOS2 + 1
          END IF
          JPIV = -JPIV
          APOS = APOS + NROWS - IPIV + 1
          END IF
  210     CONTINUE
      ELSE
        DO 211 IPIV = 1,NROWS
          IRHS = IFACT(IWPOS+IPIV-1)

          IF (IRHS.GT.0) THEN
C 1 by 1 pivot.
          DO 66 J = 1,NRHS
            RHS(IRHS,J) = RHS(IRHS,J)*FACT(APOS)
   66       CONTINUE
          APOS = APOS + NROWS - IPIV + 1
          ELSE
C 2 by 2 pivot
          IF (JPIV.EQ.1) THEN
            IRHS1 = -IRHS
            IRHS2 = -IFACT(IWPOS+IPIV)
            DO 69 J = 1,NRHS
            W1 = RHS(IRHS1,J)*FACT(APOS) +
   +             RHS(IRHS2,J)*FACT(APOS2)
            RHS(IRHS2,J) = RHS(IRHS1,J)*CONJG(FACT(APOS2)) +
   +                   RHS(IRHS2,J)*FACT(APOS+NROWS-IPIV+1)
            RHS(IRHS1,J) = W1
   69         CONTINUE
            APOS2 = APOS2 + 1
          END IF
          JPIV = -JPIV
          APOS = APOS + NROWS - IPIV + 1
          END IF
  211     CONTINUE
      ENDIF

      END IF

      IWPOS = IWPOS + NCOLS
      APOS = APOS + NROWS*(NCOLS-NROWS)

  380 CONTINUE

    END

    SUBROUTINE ME57T(KIND,N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
   *            W,LW,IW1,ICNTL)
C This subroutine performs backward elimination operations
C   using the factors stored in FACT/IFACT by ME57B.
    INTEGER KIND,N,LFACT
    COMPLEX FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
    COMPLEX W(LW,NRHS),RHS(LRHS,NRHS)
    INTEGER IW1(N),ICNTL(20)
C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C N    must be set to the order of the matrix. It is not altered.
C FACT    must be set to hold the real values corresponding to the
C    factors. This must be unchanged since the
C    preceding call to ME57B. It is not altered.
C LFACT   length of array FACT. It is not altered.
C IFACT   holds the integer indexing information for the matrix factors
C    in FACT. This must be unchanged since the preceding call to
C    ME57B. It is not altered.
C LIFACT  length of array IFACT. It is not altered.
C NRHS must be set to number of right-hand sides.
C RHS  on entry, must be set to hold the right hand side modified by
C    the forward substitution operations. On exit, holds the
C    solution vector.
C LHS must be set to the leading dimension of array RHS.
C W    used as workspace to hold the components of the right hand
C    sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C   larger than INFO(21) as returned from ME57B.
C IW1  on entry IW1(I) (I = 1,NBLK), where  NBLK = IFACT(3) is the
C    number of block pivots, must hold pointers to the beginning of
C    each block pivot in array IFACT, as set by ME57Q\QD.
C    It is not altered.
C ICNTL Not referenced except:
C   ICNTL(13) Threshold on number of columns in a block for direct
C       addressing using Level 2 and Level 3 BLAS.

C Procedures
    INTRINSIC ABS
    EXTERNAL CGEMM,CTPSV

C Constants
    REAL ONE
    PARAMETER (ONE=1.0E0)
    COMPLEX CONE
    PARAMETER (CONE=(1.0E0,0.0E0))
C
C Local variables.
    INTEGER APOS,I,IBLK,II,IPIV,IRHS,
   +      IWPOS,J,J1,J2,K,KK,NCOLS,NROWS
    COMPLEX W1
C APOS  Current position in array FACT.
C I   Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IRHS1 RHS index.
C IRHS2 RHS index.
C IWPOS Position in IFACT of start of current index list.
C J   Temporary DO index
C JPIV  Has the value -1 for the first row of a 2 by 2 pivot and 1 for
C     the second.
C K   Temporary pointer to position in real array.
C J1    Position in IFACT of index of leading entry of row.
C J2    Position in IFACT of index of trailing entry of row.
C K   Temporary variable.
C LROW  Length of current row.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.
C
    APOS = IFACT(1)


C Set IW1 array
    IWPOS = 4
    DO 10 I = 1,IFACT(3)-1
      IW1(I) = IWPOS
      IWPOS = IWPOS + ABS(IFACT(IWPOS))+2
   10 CONTINUE
    IW1(IFACT(3)) = IWPOS

C Run through block pivot rows in the reverse order.
    DO 380 IBLK = IFACT(3),1,-1

C Find the number of rows and columns in the block.
      IWPOS = IW1(IBLK)
      NCOLS = ABS(IFACT(IWPOS))
      NROWS = ABS(IFACT(IWPOS+1))
      APOS = APOS - NROWS* (NCOLS-NROWS)
      IWPOS = IWPOS + 2

      IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN


C Perform operations using direct addressing.

C Load right-hand side into W.
      DO 5 I = 1,NCOLS
        II = ABS(IFACT(IWPOS+I-1))
        DO 3 J = 1,NRHS
          W(I,J) = RHS(II,J)
    3     CONTINUE
    5   CONTINUE


C Treat off-diagonal block (direct addressing)
      K = NCOLS - NROWS
      IF (K.GT.0)  THEN
        IF (KIND.EQ.1) THEN
          CALL CGEMM('T','N',NROWS,NRHS,K,CONE,
   +             FACT(APOS),K,
   +             W(NROWS+1,1),LW,CONE,W,LW)
        ELSE
          CALL CGEMM('C','N',NROWS,NRHS,K,CONE,
   +             FACT(APOS),K,
   +             W(NROWS+1,1),LW,CONE,W,LW)
        ENDIF
      ENDIF

      APOS = APOS-(NROWS*(NROWS+1))/2

C Treat diagonal block (direct addressing)
      IF (KIND.EQ.1) THEN
        DO 22 J = 1,NRHS
          CALL CTPSV('L','T','U',NROWS,FACT(APOS),W(1,J),1)
   22     CONTINUE
      ELSE
        DO 23 J = 1,NRHS
          CALL CTPSV('L','C','U',NROWS,FACT(APOS),W(1,J),1)
   23     CONTINUE
      ENDIF

C Reload W back into RHS.
      DO 60 I = 1,NROWS
        II = ABS(IFACT(IWPOS+I-1))
        DO 59 J = 1,NRHS
          RHS(II,J) = W(I,J)
   59     CONTINUE
   60   CONTINUE

      ELSE
C
C Perform operations using indirect addressing.
      J1 = IWPOS
      J2 = IWPOS + NCOLS - 1

C Treat off-diagonal block (indirect addressing)
      KK = APOS
      J1 = IWPOS + NROWS
      IF (KIND.EQ.1) THEN
        DO 220 IPIV = 1,NROWS
          IRHS = ABS(IFACT(IWPOS+IPIV-1))
          DO 218 II = 1,NRHS
          W1 = RHS(IRHS,II)
          K = KK
          DO 215 J = J1,J2
            W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)),II)
            K = K + 1
  215       CONTINUE
          RHS(IRHS,II) = W1
  218       CONTINUE
          KK = K
  220     CONTINUE
      ELSE
        DO 221 IPIV = 1,NROWS
          IRHS = ABS(IFACT(IWPOS+IPIV-1))
          DO 219 II = 1,NRHS
          W1 = RHS(IRHS,II)
          K = KK
          DO 216 J = J1,J2
            W1 = W1 + CONJG(FACT(K))*RHS(ABS(IFACT(J)),II)
            K = K + 1
  216       CONTINUE
          RHS(IRHS,II) = W1
  219       CONTINUE
          KK = K
  221     CONTINUE
      ENDIF

C Treat diagonal block (indirect addressing)
      J2 = IWPOS + NROWS - 1
      DO 260 IPIV = 1,NROWS
        IRHS = ABS(IFACT(J1-1))
        APOS = APOS - IPIV
        DO 240 II = 1,NRHS
          W1 = RHS(IRHS,II)
          K = APOS + 1
          DO 230 J = J1,J2
          IF (KIND.EQ.1) THEN
            W1 = W1 - FACT(K)*RHS(ABS(IFACT(J)),II)
          ELSE
            W1 = W1 - CONJG(FACT(K))*RHS(ABS(IFACT(J)),II)
          ENDIF
          K = K + 1
  230       CONTINUE
          RHS(IRHS,II) = W1
  240     CONTINUE
        J1 = J1 - 1
  260   CONTINUE

      END IF

  380 CONTINUE


    END

    SUBROUTINE ME57D(KIND,JOB,N,NE,A,IRN,JCN,FACT,LFACT,IFACT,LIFACT,
   *            RHS,X,RESID,WORKC,WORKR,IW,ICNTL,CNTL,INFO,
   *            RINFO)

C This subroutine solves a single system using one or more steps of
C   iterative refinement.
C If ICNTL(9) = 10 (the default), this subroutine performs iterative
C   refinement using the strategy of Arioli, Demmel, and Duff.
C   IF (ICNTL(9) = 1, then one step of iterative refinement is
C   performed.

    INTEGER KIND,JOB,N,NE
    COMPLEX A(NE)
    INTEGER IRN(NE),JCN(NE),LFACT
    COMPLEX FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT)
    COMPLEX RHS(N),X(N),RESID(N),WORKC(N,*)
    REAL WORKR(N,*)
    INTEGER IW(N),ICNTL(20)
    REAL CNTL(5)
    INTEGER INFO(40)
    REAL RINFO(20)

C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C JOB must be set by the user to determine what action is desired by
C   the user.
C   Values of JOB and their effect are:
C IF ICNTL(9)>1, JOB=0 if no estimate of solution in X; JOB=2 if
C      estimate of solution in X.
C IF ICNTL(9)=1, then:
C   0: Solve Ax=b, calculate residual r=b-Ax and exit.
C      (Note that ME57C should be used if solution without residual
C       is required)
C   1: Solve Ax=b, calculate residual r=b-Ax, solve A(dx)=r,
C      update solution and exit.
C If JOB > 1, an estimate of the solution must be input in X.
C   2: Calculate residual r=b-Ax, solve A(dx)=r,
C      update solution and exit.
C If JOB > 2, the residual for this estimate must also be input.
C   3: Solve A(dx)=r, update solution and exit.
C   4: Solve A(dx)=r, update solution, calculate residual for new
C      solution and exit.
C N must be set to the order of the matrix and must be unchanged since
C   the call to ME57B. It is not altered by ME57D.
C NE must be set to the number of entries in the matrix and must be
C   unchanged since the call to ME57A. It is not altered by ME57D.
C A must be set by the user to the values of the matrix as input to
C   ME57A. It is not altered by ME57D.
C IRN,JCN must be set to the row and column indices of the entries in A
C   and must be unchanged since the call to ME57A. They are
C   not altered by ME57D.
C FACT holds information on the factors and must be unchanged since
C   the call to ME57B. It is not altered by ME57D.
C LFACT must be set to the length of FACT. It is not altered by ME57D.
C IFACT holds information on the factors and must be unchanged since
C   the call to ME57B. It is not altered by ME57D.
C LIFACT must be set to the length of IFACT. Not altered by ME57D.
C RHS is a real array of length N that must be set to the right-hand
C   side for the equation being solved. It is not altered by ME57D.
C X is a real array of length N. IF JOB >=2, it must be set on entry to
C   an estimated solution. Otherwise, it need not be set by the user.
C   On exit, the improved solution vector is returned in X.
C RESID is a real array of length N. If JOB > 2, it must be set on entry
C   to the value of the residual for the current solution estimate
C   held in X.  Otherwise, it need not be set by the user on entry.
C   If JOB=0 or 4 or if ICNTL(9)>1, on exit it will hold the residual
C   vector for the equations being solved. If 1<= JOB <= 3, then
C   RESID will hold on exit the last correction vector added to the
C   solution X.
C WORKC is used as a work array.  If ICNTL(9) = 1, it must be of length
C   at least N.  If ICNTL(9)>1 and ICNTL(10)>0, then WORKC must be of
C   length at least 2*N.
C WORKR is used as a work array. If ICNTL(9) = 1, it is not accessed by
C   the subroutine  If ICNTL(9) > 1, it must be of length at least
C   2*N.
C IW is an integer array of length N that is used as a work array if
C   ICNTL(1) > 9.  It is not accessed if ICNTL(9) = 1.
C ICNTL must be set by the user as follows and is not altered.
C   ICNTL(1)  must be set to the stream number for error messages.
C     A value less than 0 suppresses output.
C   ICNTL(2) must be set to the stream number for warning output.
C     A value less than 0 suppresses output.
C   ICNTL(3) must be set to the stream number for monitor output.
C     A value less than 0 suppresses output.
C   ICNTL(4) must be set to the stream number for statistics output.
C     A value less than 0 suppresses output.
C   ICNTL(5) must be set to control the amount of output:
C     0 None.
C     1 Error messages only.
C     2 Error and warning messages.
C     3 As 2, plus scalar parameters and a few entries of array
C       parameters on entry and exit.
C   > 3  As 2, plus all parameters on entry and exit.
C   ICNTL(6:12) Not referenced.
C   ICNTL(9)  Maximum permitted number of steps of iterative
C         refinement.
C   ICNTL(10) Flag to request calculation of error estimate and
C       condition numbers.
C   ICNTL(13) Threshold on number of columns in a block for direct
C       addressing using Level 2 and Level 3 BLAS.
C CNTL must be set by the user as follows and is not altered.
C   CNTL(3) is the required decrease in the scaled residuals required
C       by the Arioli, Demmel, and Duff iteration.
C INFO is an integer array that need not be set by the user.  On exit,
C   a value of INFO(1) equal to zero indicates success. A failure is
C   indicated by a negative value for INFO.
C RINFO is a real array that need not be set by the user. On exit,
C   If ICNTL(9)>1, RINFO is set to information on the matrix and
C   solution including backward errors.

C   .. Local constants ..
    REAL ZERO,ONE
    PARAMETER (ZERO=0.D0,ONE=1.0E0)

C   .. Local variables ..
    REAL COND(2),CTAU,DXMAX,ERROR,OLDOMG(2),OLDOM2,
   *           OMEGA(2),OM2,TAU
    INTEGER I,ICNTLC(20),ITER,J,K,KASE,KK,LDIAG,LP,MP,KEEP71(5)
    LOGICAL LCOND(2)

C
C COND is condition number of system.
C CTAU is set to 1000*machine precision.
C DXMAX used to calculate max norm of solution.
C ERROR used to accumulate error.
C OLDOMG holds value of previous backward errors.
C OLDOM2 holds previous sum of OMEGAs.
C OMEGA used to accumulate backward error.
C OM2 holds sum of OMEGAs.
C TAU is threshold for denominator in scaled residual calculation.
C I  Temporary variable.
C ICNTLC is control array for ME57C/CD.
C ITER is maximum permitted number of steps of iterative refinement.
C J  Temporary variable.
C K  Temporary variable.
C KASE used when calling MF71A.
C KK Temporary variable.
C LDIAG Control for amount of information output.
C LCOND used as switch when calculating condition number.
C LP Stream number for error printing.
C MP Stream number for monitor printing.
C KEEP71 Work array required by MF71.
C

C Procedures
    INTRINSIC MIN
    EXTERNAL ME57C,ME57U,FD15A,MF71A
C EPS is the largest real such that 1+EPS is equal to 1.
    REAL EPS,FD15A

    INFO(1) = 0

C Set local print control variables
    LP = ICNTL(1)
    MP = ICNTL(3)
    LDIAG = ICNTL(5)

C Check input data
    IF (N.LE.0) THEN
      INFO(1) = -1
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I12)')
   +    '**** Error return from ME57D ****  INFO(1) =',INFO(1),
   +    'N has value',N
      INFO(2) = N
      GOTO 500
    ENDIF

    IF (NE.LT.0) THEN
      INFO(1) = -2
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I12)')
   +    '**** Error return from ME57D ****  INFO(1) =',INFO(1),
   +    'NE has value',NE
      INFO(2) = NE
      GOTO 500
    ENDIF

    IF (ICNTL(9).LT.1) THEN
      INFO(1) = -13
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I12)')
   +    '**** Error return from ME57D ****  INFO(1) =',INFO(1),
   +    'ICNTL(9) has value',ICNTL(9)
      INFO(2) = ICNTL(9)
      GOTO 500
    ENDIF

    IF (JOB.LT.0 .OR. JOB.GT.4 .OR. (ICNTL(9).GT.1 .AND.
   *    (JOB.NE.0 .AND. JOB.NE.2)))  THEN
      INFO(1) = -12
      INFO(2) = JOB
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I12)')
   +    '**** Error return from ME57D ****  INFO(1) =',INFO(1),
   +    'JOB has value',JOB
      IF (ICNTL(9).GT.1 .AND. LDIAG.GT.0 .AND. LP.GE.0)
   +    WRITE (LP,'(A,I3)') 'and ICNTL(9) =',ICNTL(9)
      GOTO 500
    ENDIF

C If NE = 0, set variables and return
    IF (NE.EQ.0) THEN
      IF (JOB.NE.3) THEN
      DO 8 I = 1,N
        RESID(I) = ZERO
  8     CONTINUE
      ENDIF
      DO 9 I = 1,N
      X(I) = ZERO
  9   CONTINUE
      INFO(30)=0
      DO 10 I = 6,13
      RINFO(I) = ZERO
 10   CONTINUE
      GO TO 500
    ENDIF

C If requested, print input parameters
    IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
      WRITE (MP,'(/A,A/)')
   +      'Entering iterative refinement solution phase ',
   +      '(ME57D) with ...'
      IF (KIND.EQ.1) THEN
      WRITE(MP,'(A)') 'Matrix is complex symmetric'
      ELSE
      WRITE(MP,'(A)') 'Matrix is Hermitian'
      ENDIF
      WRITE (MP,99999) JOB,N,NE,(ICNTL(I),I=1,5),LFACT,LIFACT,
   +   ICNTL(9),ICNTL(10),ICNTL(13),CNTL(3)
99999 FORMAT(
   +  'JOB     Control for coefficient matrix    =',I12/
   +  'N       Order of matrix             =',I12/
   +  'NE      Number of entries in matrix       =',I12/
   6  'ICNTL(1)  Stream for errors             =',I12/
   7  ' --- (2)  Stream for warnings           =',I12/
   8  ' --- (3)  Stream for monitoring         =',I12/
   9  ' --- (4)  Stream for statistics         =',I12/
   1  ' --- (5)  Level of diagnostic printing      =',I12/
   +  'LFACT   Length of array FACT          =',I12/
   +  'LIFACT    Length of array IFACT         =',I12/
   +  'ICNTL(9)  Number steps iterative refinement   =',I12/
   +  'ICNTL(10) Control for error analysis      =',I12/
   +  'ICNTL(13) Threshold for Level 2 and 3 BLAS    =',I12/
   +  'CNTL(3)   Convergence test for IR         =',1P,E12.4)

C Print out matrix factors.
      CALL ME57U(FACT,LFACT,IFACT,LIFACT,ICNTL)
      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      WRITE(MP,'(/A)') 'Right-hand side'
      WRITE (MP,'((4X, 1P,6E13.3))') (RHS(I),I=1,K)
      IF (K.LT.N) WRITE (MP,'(A)') '   . . .'
    END IF


    DO 15 I=1,5
      ICNTLC(I) = ICNTL(I)
   15 CONTINUE
    ICNTLC(13) = ICNTL(13)
    ICNTLC(15) = ICNTL(15)
C Switch off monitor printing in ME57C/CD
    ICNTLC(3) = -1

    IF (JOB.LE.2) THEN
      IF (JOB .LE. 1) THEN
C No estimate given in X.
      DO 14 I = 1,N
        X(I) = RHS(I)
        RESID(I) = RHS(I)
   14   CONTINUE
C Solve system Ax=b using ME57C/CD
      CALL ME57C(KIND,1,N,FACT,LFACT,IFACT,LIFACT,1,X,N,WORKC,N,
   +          IW,ICNTLC,INFO)
      ELSE
C Estimate of solution was input in X.
      DO 13 I = 1,N
        RESID(I) = RHS(I)
   13   CONTINUE
      ENDIF

      IF (ICNTL(9).EQ.1) THEN
C Compute residual
      DO 16 KK = 1,NE
        I = IRN(KK)
        J = JCN(KK)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 16
        RESID(I) = RESID(I) - A(KK)*X(J)
C Matrix is symmetric
        IF (I.NE.J) THEN
          IF (KIND.EQ.1) THEN
          RESID(J) = RESID(J) - A(KK)*X(I)
          ELSE
          RESID(J) = RESID(J) - CONJG(A(KK))*X(I)
          ENDIF
        ENDIF
   16   CONTINUE
      IF (JOB.EQ.0) GO TO 340
      ELSE
C Calculate values for iterative refinement strategy of Arioli,
C     Demmel and Duff.
      DO 18 I = 1,N
        WORKR(I,1) = ZERO
C Calculate infinity norms of rows of A in WORKR(I,1)
C Sum |a  |, j=1,N (= ||A  ||      )
C     ij         i.  infinity
        WORKR(I,2) = ZERO
   18   CONTINUE
      DO 17 KK = 1,NE
        I = IRN(KK)
        J = JCN(KK)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 17
        RESID(I) = RESID(I) - A(KK)*X(J)
        WORKR(J,1) = WORKR(J,1) + ABS(A(KK)*X(I))
        WORKR(J,2) = WORKR(J,2) + ABS(A(KK))
C Matrix is symmetric
        IF (I.NE.J) THEN
          IF (KIND.EQ.1) THEN
          RESID(J) = RESID(J) - A(KK)*X(I)
          ELSE
          RESID(J) = RESID(J) - CONJG(A(KK))*X(I)
          ENDIF
          WORKR(I,1) = WORKR(I,1) + ABS(A(KK)*X(J))
          WORKR(I,2) = WORKR(I,2) + ABS(A(KK))
        ENDIF
   17   CONTINUE

C Calculate max-norm of solution
      DXMAX = ZERO
      DO 221 I = 1,N
      DXMAX = MAX(DXMAX,ABS(X(I)))
  221   CONTINUE

C Initialize EPS
    EPS = FD15A('E')
C Calculate backward errors OMEGA(1) and OMEGA(2)
C CTAU ... 1000 eps (approx)
      CTAU = 1000.*EPS
C tau is (||A  ||       ||x||   + |b| )*n*1000*epsilon
C        i.  infinity    max   i
      OMEGA(1) = ZERO
      OMEGA(2) = ZERO
      DO 231 I = 1,N
        TAU = (WORKR(I,2)*DXMAX+ABS(RHS(I)))*N*CTAU
        IF ((WORKR(I,1)+ABS(RHS(I))).GT.TAU) THEN
C |Ax-b| /(|A||x| + |b|)
C     i         i
          OMEGA(1) = MAX(OMEGA(1),ABS(RESID(I))/
   +             (WORKR(I,1)+ABS(RHS(I))))
          IW(I) = 1
        ELSE
C TAU will be zero if all zero row in A, for example
          IF (TAU.GT.ZERO) THEN
C |Ax-b| /(|A||x| + ||A  ||      ||x||   )
C     i      i   i.  infinity   max
          OMEGA(2) = MAX(OMEGA(2),ABS(RESID(I))/
   +             (WORKR(I,1)+WORKR(I,2)*DXMAX))
          END IF
          IW(I) = 2
        END IF
  231   CONTINUE
C
C  Stop the calculations if the backward error is small
C
      OM2 = OMEGA(1) + OMEGA(2)
      ITER = 0
C Statement changed because IBM SP held quantities in registers
C       IF ((OM2+ONE).LE.ONE) THEN
      IF (OM2.LE.EPS) THEN
C       write(7,'(A)') 'Go to 270 from after 231'
        GO TO 270
      ENDIF

C Hold current estimate in case needed later
      DO 251 I = 1,N
        WORKC(I,2) = X(I)
  251   CONTINUE
      OLDOMG(1) = OMEGA(1)
      OLDOMG(2) = OMEGA(2)
      OLDOM2 = OM2

      ENDIF

    ENDIF

C At this point JOB >= 1 or ICNTL(9) > 1
C
C Iterative refinement loop
    DO 260 ITER = 1,ICNTL(9)

C     write(7,'(A,I3/A/(5E16.8))') 'Iteration',ITER,
C    +  'Residual',(RESID(I),I=1,N)
C Solve system A(dx) = r using ME57C/CD
      CALL ME57C(KIND,1,N,FACT,LFACT,IFACT,LIFACT,1,RESID,N,WORKC,N,
   +          IW,ICNTLC,INFO)
C     write(7,'(A,I3/A/(5E16.8))') 'Iteration',ITER,
C    +  'Delta x',(RESID(I),I=1,N)

C Update solution
      DO 141 I = 1,N
      X(I) = X(I) + RESID(I)
  141   CONTINUE

C Exit without computing residual
      IF (JOB.LT.4 .AND. ICNTL(9).EQ.1) GO TO 340

C
C Calculate residual using information in A,IRN,JCN
C If ICNTL(9).GT.1 also calculate |A| |x|
C
      IF (ICNTL(9).EQ.1) THEN
C Now JOB = 4
      DO 151 I = 1,N
        RESID(I) = RHS(I)
  151   CONTINUE
      DO 181 KK = 1,NE
        I = IRN(KK)
        J = JCN(KK)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 181
        RESID(I) = RESID(I) - A(KK)*X(J)
C Matrix is symmetric
        IF (I.NE.J) THEN
          IF (KIND.EQ.1) THEN
          RESID(J) = RESID(J) - A(KK)*X(I)
          ELSE
          RESID(J) = RESID(J) - CONJG(A(KK))*X(I)
          ENDIF
        ENDIF
  181   CONTINUE
C Jump because only one step of iterative refinement requested.
      GO TO 340
      ELSE
      DO 153 I = 1,N
C b - Ax
        RESID(I) = RHS(I)
C |A||x|
        WORKR(I,1) = ZERO
  153   CONTINUE
      DO 183 KK = 1,NE
        I = IRN(KK)
        J = JCN(KK)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 183
        RESID(I) = RESID(I) - A(KK)*X(J)
        WORKR(J,1) = WORKR(J,1) + ABS(A(KK)*X(I))
C Matrix is symmetric
        IF (I.NE.J) THEN
          IF (KIND.EQ.1) THEN
          RESID(J) = RESID(J) - A(KK)*X(I)
          ELSE
          RESID(J) = RESID(J) - CONJG(A(KK))*X(I)
          ENDIF
          WORKR(I,1) = WORKR(I,1) + ABS(A(KK)*X(J))
        ENDIF
  183   CONTINUE
      ENDIF

C Calculate max-norm of solution
      DXMAX = ZERO
      DO 220 I = 1,N
      DXMAX = MAX(DXMAX,ABS(X(I)))
  220   CONTINUE

C Calculate omega(1) and omega(2)
C tau is (||A  ||       ||x||   + |b| )*n*1000*epsilon
C        i.  infinity    max   i
      OMEGA(1) = ZERO
      OMEGA(2) = ZERO
      DO 230 I = 1,N
      TAU = (WORKR(I,2)*DXMAX+ABS(RHS(I)))*N*CTAU
      IF ((WORKR(I,1)+ABS(RHS(I))).GT.TAU) THEN
C |Ax-b| /(|A||x| + |b|)
C     i         i
        OMEGA(1) = MAX(OMEGA(1),ABS(RESID(I))/
   +           (WORKR(I,1)+ABS(RHS(I))))
        IW(I) = 1
      ELSE
C TAU will be zero if all zero row in A, for example
        IF (TAU.GT.ZERO) THEN
C |Ax-b| /(|A||x| + ||A  ||      ||x||   )
C     i      i   i.  infinity   max
          OMEGA(2) = MAX(OMEGA(2),ABS(RESID(I))/
   +             (WORKR(I,1)+WORKR(I,2)*DXMAX))
        END IF
        IW(I) = 2
      END IF
  230   CONTINUE
C
C  Stop the calculations if the backward error is small
      OM2 = OMEGA(1) + OMEGA(2)
      IF ((OM2+ONE).LE.ONE) THEN
C       write(7,'(A)') 'Go to 270 from after 230'
      GO TO 270
      ENDIF
C
C  Check the convergence.
C
C     write(7,'(A,2E16.8)') 'OM2,OLDOM2',OM2,OLDOM2
      IF (OM2.GT.OLDOM2*CNTL(3)) THEN
C  Stop if insufficient decrease in omega.
      IF (OM2.GT.OLDOM2) THEN
C Previous estimate was better ... reinstate it.
        OMEGA(1) = OLDOMG(1)
        OMEGA(2) = OLDOMG(2)
        DO 240 I = 1,N
          X(I) = WORKC(I,2)
  240     CONTINUE
      END IF
C       write(7,'(A)') 'Go to 270 from after 240'
      GO TO 270
      ELSE
C Hold current estimate in case needed later
      DO 250 I = 1,N
        WORKC(I,2) = X(I)
  250   CONTINUE
      OLDOMG(1) = OMEGA(1)
      OLDOMG(2) = OMEGA(2)
      OLDOM2 = OM2
      END IF

  260 CONTINUE
C End of iterative refinement loop.

    INFO(1) = -8
    IF (LP.GE.0 .AND. LDIAG.GE.1) WRITE (LP,9170) INFO(1),ICNTL(9)
 9170 FORMAT ('Error return from ME57D/DD because of ','nonconvergence',
   +     ' of iterative refinement'/'Error INFO(1) = ',I2,'  with',
   +     ' ICNTL','(9) = ',I10)

C Set the RINFO parameters
  270 RINFO(6)  = OMEGA(1)
    RINFO(7)  = OMEGA(2)
    RINFO(8) = ZERO
    DO 271 I=1,N
      RINFO(8) = MAX(RINFO(8),WORKR(I,2))
  271 CONTINUE
    RINFO(9) = DXMAX
    RINFO(10) = ZERO
    DO 272 I=1,N
      RINFO(10) = MAX(RINFO(10),ABS(RESID(I)))
  272 CONTINUE
    IF (RINFO(8)*RINFO(9).NE.ZERO)
   *RINFO(10) = RINFO(10)/(RINFO(8)*RINFO(9))
    INFO(30) = ITER

    IF (INFO(1).LT.0) GO TO 340

C Jump if estimate of error not requested.
    IF (ICNTL(10).LE.0) GO TO 340

C
C Calculate condition numbers and estimate of the error.
C
C Condition numbers obtained through use of norm estimation
C   routine MF71A/AD.
C
C Initializations
C
    LCOND(1) = .FALSE.
    LCOND(2) = .FALSE.
    ERROR    = ZERO
    DO 280 I = 1,N
      IF (IW(I).EQ.1) THEN
      WORKR(I,1) = WORKR(I,1) + ABS(RHS(I))
C |A||x| + |b|
      WORKR(I,2) = ZERO
      LCOND(1) = .TRUE.
      ELSE
C |A||x| + ||A  ||      ||x||
C         i.  infinity   max

      WORKR(I,2) = WORKR(I,1) + WORKR(I,2)*DXMAX
      WORKR(I,1) = ZERO
      LCOND(2) = .TRUE.
      END IF
  280 CONTINUE
C
C  Compute the estimate of COND
C
    DO 330 K = 1,2
      IF (LCOND(K)) THEN
C MF71A/AD has its own built in limit of 5 to the number of iterations
C    allowed. It is this limit that will be used to terminate the
C    following loop.
      KASE = 0
      DO 310 KK = 1,40

C MF71A/AD calculates norm of matrix
C We are calculating the infinity norm of INV(A).W

C Initialize WORKC(1,1).  WORKC(1,2) is a work array.
        CALL MF71A(N,KASE,WORKC,COND(K),WORKC(1,2),KEEP71)
C
C  KASE = 0........ Computation completed
C  KASE = 1........ W * INV(TRANSPOSE(A)) * Y
C  KASE = 2........ INV(A) * W * Y
C             W is W/W(*,2) .. Y is WORKC(*,1)
C
        IF (KASE.EQ.0) GO TO 320

        IF (KASE.EQ.1) THEN
C Solve system using ME57C/CD.
C WORKC(1,2) is used as workspace
          CALL ME57C(KIND,1,N,FACT,LFACT,IFACT,LIFACT,1,WORKC,
   *            N,WORKC(1,2),N,IW,ICNTLC,INFO)
          DO 290 I = 1,N
          WORKC(I,1) = WORKR(I,K)*WORKC(I,1)
  290       CONTINUE
        END IF

        IF (KASE.EQ.2) THEN
          DO 300 I = 1,N
          WORKC(I,1) = WORKR(I,K)*WORKC(I,1)
  300       CONTINUE
C Solve system using ME57C/CD.
          CALL ME57C(KIND,1,N,FACT,LFACT,IFACT,LIFACT,1,WORKC,N,
   *            WORKC(1,2),N,IW,ICNTLC,INFO)
        END IF

  310   CONTINUE

C Error return if MF71A does not converge
      INFO(1) = -14
      IF (LP.GE.0 .AND. LDIAG.GE.1) WRITE (LP,9160)
 9160 FORMAT ('Error return from ME57D/DD because of ','error in MF71',
   +     'A/AD'/'Error not calculated')

  320   IF (DXMAX.GT.ZERO) COND(K) = COND(K)/DXMAX
      ERROR = ERROR + OMEGA(K)*COND(K)
      ELSE
      COND(K) = ZERO
      ENDIF

  330 CONTINUE

    RINFO(11)  = COND(1)
    RINFO(12)  = COND(2)
    RINFO(13)  = ERROR

C
C If requested, print output parameters.
 340  IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
      WRITE(MP,99980) INFO(1)
99980 FORMAT (/'Leaving iterative refinement solution phase ',
   +  '(ME57D) with ...'/
   1    'INFO (1)                        =',I12/)
      IF (INFO(1).LT.0) GO TO 500
      IF (ICNTL(9).GT.1) THEN
      WRITE(MP,99981) INFO(30),(RINFO(I),I=6,10)
99981   FORMAT(
   1   'INFO(30)  Number steps iterative ref   =',I10/
   1   'RINFO(6)  Backward errors  (OMEGA(1))  =',1PE10.3/
   2   '-----(7)  Backward errors  (OMEGA(2))  =',1PE10.3/
   3   '-----(8)  Infinity norm of matrix    =',1PE10.3/
   4   '-----(9)  Infinity norm of solution    =',1PE10.3/
   5   '-----(10) Norm of scaled residuals   =',1PE10.3)
      IF (ICNTL(10).GT.0) WRITE(MP,99979) (RINFO(I),I=11,13)
99979     FORMAT (
   1     'RINFO(11) Condition number (COND(1))   =',1PE10.3/
   1     'RINFO(12) Condition number (COND(2))   =',1PE10.3/
   1     'RINFO(13) Error in solution        =',1PE10.3)
      WRITE(MP,'(/A,I10)') 'Residual'
      K=MIN(N,10)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(1P,6E13.3)') (RESID(I),I=1,K)
      IF (K.LT.N) WRITE (MP,'(A)') '   . . .'
      ELSE
C ICNTL(9) = 1
      IF (JOB.GE.1 .AND. JOB.LE.3) THEN
        WRITE(MP,'(/A,I10)') 'Correction to solution'
      ELSE
        WRITE(MP,'(/A,I10)') 'Residual'
      ENDIF
      K=MIN(N,10)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(1P,6E13.3)') (RESID(I),I=1,K)
      IF (K.LT.N) WRITE (MP,'(A)') '   . . .'
      END IF

C Print solution
      K=MIN(N,10)
      IF (LDIAG.GE.4) K = N
      WRITE(MP,'(/A,I10)') 'Solution'
      WRITE (MP,'(1P,6E13.3)') (X(I),I=1,K)
      IF (K.LT.N) WRITE (MP,'(A)') '   . . .'

    END IF

 500  RETURN

    END
    SUBROUTINE ME57E(N,IC,KEEP,FACT,LFACT,NEWFAC,LNEW,
   *            IFACT,LIFACT,NEWIFC,LINEW,INFO)
    INTEGER N,IC,KEEP(*),LFACT,LNEW,LIFACT,LINEW,INFO(40)
    COMPLEX FACT(LFACT),NEWFAC(LNEW)
    INTEGER IFACT(LIFACT),NEWIFC(LINEW)
C Local variables
    INTEGER APOSBB,ASTK,HOLD,I,ISTK,IWPOS,MOVE,NFRONT

C HOLD determines part of keep holding saved variables from ME57O.
    HOLD = N + 3

C Initialize INFO(1) and INFO(2)
    INFO(1) = 0
    INFO(2) = 0

C Test to see whether to map real or integer space
    IF (IC.GE.1) THEN
C Remap of integer space
      IF (LINEW.LE.LIFACT) THEN
      INFO(1) = -7
      INFO(2) = LINEW
      RETURN
      ENDIF
      IWPOS = KEEP(HOLD+7)
      ISTK  = KEEP(HOLD+14)
      NFRONT = KEEP(HOLD+23)
      DO 10 I = 1,IWPOS+NFRONT-1
      NEWIFC(I) = IFACT(I)
   10   CONTINUE
C Move distance
      MOVE = LINEW - LIFACT
      DO 20 I = ISTK+1,LIFACT
      NEWIFC(I+MOVE) = IFACT(I)
   20   CONTINUE
C Reset INPUT, ISTK, PTIRN
      KEEP(HOLD+13) = KEEP(HOLD+13) + MOVE
      KEEP(HOLD+14) = ISTK + MOVE
      KEEP(HOLD+18) = KEEP(HOLD+18) + MOVE
    ENDIF
    IF (IC.NE.1) THEN
C Remap of real space
      IF (LNEW.LE.LFACT) THEN
      INFO(1) = -7
      INFO(2) = LNEW
      RETURN
      ENDIF
C Was .. APOS = KEEP(HOLD+8)
      APOSBB = KEEP(HOLD+9)
      ASTK   = KEEP(HOLD+15)
      DO 60 I = 1, APOSBB-1
      NEWFAC(I) = FACT(I)
   60   CONTINUE
C Move distance
      MOVE = LNEW - LFACT
      DO 70 I = ASTK+1,LFACT
      NEWFAC(I+MOVE) = FACT(I)
   70   CONTINUE
C Reset AINPUT, ASTK, PTRA
      KEEP(HOLD+12) = KEEP(HOLD+12) + MOVE
      KEEP(HOLD+15) = ASTK + MOVE
      KEEP(HOLD+19) = KEEP(HOLD+19) + MOVE
    ENDIF
    RETURN
    END
    SUBROUTINE ME57G(KIND,N,NE,IRN,JCN,IW,IPE,COUNT,FLAG,IWFR,
   *            ICNTL,INFO)

C Given the positions of the entries of a symmetric matrix, construct
C   the sparsity pattern of the whole matrix. Either one of a pair
C   (I,J),(J,I) may be used to represent the pair. Duplicates are
C   ignored.

    INTEGER KIND,N,NE,IRN(NE),JCN(NE),IW(NE*2+N),IPE(N),COUNT(N),
   +      FLAG(N),IWFR,ICNTL(20),INFO(40)
C N must be set to the matrix order. It is not altered.
C NE must be set to the number of entries input. It is not altered.
C IRN(K),K=1,2,...,NE must be set to the row indices of the entries on
C   input.  IRN is not changed.
C JCN(K),K=1,2,...,NE must be set to the column indices of the entries
C   on input.  JCN is not changed.
C IW need not be set on input. On output it contains lists of
C   column indices.
C IPE need not be set on input. On output IPE(I) points to the start of
C   the entry in IW for row I, I=1,2,...,N.
C COUNT need not be set. On output COUNT(I), I=1,2,..,N, contains the
C   number of off-diagonal entries in row I excluding duplicates.
C FLAG is used for workspace to hold flags to permit duplicate entries
C   to be identified quickly.
C IWFR need not be set on input. On output it points to the first
C   unused location in IW.
C ICNTL Warning messages are printed on stream number ICNTL(2) if
C    ICNTL(2).GT.0 and ICNTL(5).GT.1.
C INFO need not be set on input. On output,
C  INFO(1) has one of the values:
C   0 No out-of-range index or duplicate entry found.
C   1 Out-of-range index found.
C   2 Duplicate entry found.
C   3 Out-of-range index found and duplicate entry found.
C  INFO(3) is set to the number of out-of-range entries.
C  INFO(4) is set to the number of off-diagonal duplicate entries.
C
    INTRINSIC MAX,MIN
C
C Local variables
    INTEGER I,J,K,L,LDIAG,LP,WP
C I Row index
C J Column index
C K Position in IRN, JCN, or IW.
C L Position in IW.
C LDIAG Level of diagnostic printing.
C WP Stream for printing warning messages.

C Set LDIAG and WP
    LP = ICNTL(1)
    WP = ICNTL(2)
    LDIAG = ICNTL(5)
    IF (WP.LT.0) LDIAG = 0

C Initialize INFO(1)
    INFO(1) = 0

C Count in INFO(3) the number of out-of-range entries, initialize FLAG,
C    and count in COUNT the numbers of entries in the rows including
C    duplicates.
    INFO(3) = 0
    DO 10 I = 1,N
      FLAG(I) = 0
      COUNT(I) = 0
   10 CONTINUE
    DO 20 K = 1,NE
      I = IRN(K)
      J = JCN(K)
      IF (KIND.NE.1 .AND. I.GT.J) THEN
      INFO(1) = -19
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I8,A)')
   +    '**** Error return from ME57A ****  INFO(1) =',INFO(1),
   +    'Hermitian matrix but entry ',K,' in lower triangle'
      INFO(2) = K
      RETURN
      ENDIF
      IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) THEN
      INFO(3) = INFO(3) + 1
      INFO(1) = 1
      IF (INFO(3).EQ.1 .AND. LDIAG.GT.1) WRITE (WP,'(2A,I2)')
   +      '*** Warning message from subroutine ME57A ***',
   +      ' INFO(1) =',INFO(1)
      IF (INFO(3).LE.10 .AND. LDIAG.GT.1) WRITE (WP,'(3(I10,A))')
   +       K,'th entry (in row',I,' and column',J,') ignored'

      ELSE IF (I.NE.J) THEN
      COUNT(I) = COUNT(I) + 1
      COUNT(J) = COUNT(J) + 1
      END IF

   20 CONTINUE
C
C Accumulate row counts in IPE which is set so that IPE(I) points to
C   position after end of row I.
    IPE(1) = COUNT(1)+1
    DO 30 I = 2,N
      IPE(I) = IPE(I-1) + COUNT(I)
   30 CONTINUE
C
C Run through putting the matrix entries in the right places. IPE is
C   used for holding running pointers and is left holding pointers to
C   row starts.
    DO 40 K = 1,NE
      I = IRN(K)
      J = JCN(K)
      IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N .OR. I.EQ.J) GO TO 40
      IPE(I) = IPE(I) - 1
      IW(IPE(I)) = J
      IPE(J) = IPE(J) - 1
      IW(IPE(J)) = I
   40 CONTINUE
C
C Remove duplicates.
    INFO(4) = 0
C IWFR points to the current position in the compressed set of rows.
    IWFR = 1
C At the start of cycle I of this loop FLAG(J).LT.I for J=1,2,...,N.
C   During the loop FLAG(J) is set to I if an entry in column J is
C   found. This permits duplicates to be recognized quickly.
    DO 60 I = 1,N
      L = IPE(I)
      IPE(I) = IWFR
      DO 50 K = L,L+COUNT(I)-1
      J = IW(K)

      IF (FLAG(J).NE.I) THEN
        FLAG(J) = I
        IW(IWFR) = J
        IWFR = IWFR + 1
      ELSE
C Count duplicates only once.
        IF (I.LT.J) INFO(4) = INFO(4) + 1
      END IF

   50   CONTINUE
C Set COUNT to number without duplicates.
      COUNT(I) = IWFR - IPE(I)
   60 CONTINUE
C
C Test whether duplicates found
    IF (INFO(4).GT.0) THEN
      INFO(1) = INFO(1) + 2
      IF (LDIAG.GT.1 .AND. WP.GE.0) WRITE (WP,'(A/I10,A)')
   +    '*** Warning message from subroutine ME57A ***',INFO(4),
   +    ' off-diagonal duplicate entries found'
    END IF

    END
    SUBROUTINE ME57J(KIND,N,NE,IRN,JCN,PERM,IW,IPE,COUNT,FLAG,IWFR,
   +            ICNTL,INFO)

C Given the positions of the entries of a symmetric matrix and a
C   permutation, construct the sparsity pattern of the upper
C   triangular part of the matrix. Either one of a pair
C   (I,J),(J,I) may be used to represent the pair. Duplicates are
C   ignored.

    INTEGER KIND,N,NE,IRN(NE),JCN(NE),IW(NE+N),IPE(N),COUNT(N),
   +      PERM(N),FLAG(N),IWFR,ICNTL(20),INFO(40)
C N must be set to the matrix order. It is not altered.
C NE must be set to the number of entries input. It is not altered.
C IRN(K),K=1,2,...,NE must be set to the row indices of the entries on
C   input. If IRN(K) or JCN(K) is out of range, IRN(K) is replaced by
C   zero. Otherwise, IRN(K) is not changed.
C JCN(K),K=1,2,...,NE must be set to the column indices of the entries
C   on input. If IRN(K) or JCN(K) is out of range, JCN(K) is replaced
C   by zero. Otherwise, JCN(K) is not changed.
C PERM must be set so that PERM(I) holds the position of variable I
C   in the permuted order.  Its validity as a permutation will have
C   been checked in ME57A/AD.
C IW need not be set on input. On output it contains lists of
C   column indices, each list being headed by its length.
C IPE need not be set on input. On output IPE(I) points to the start of
C   the entry in IW for row I, I=1,2,...,N.
C COUNT need not be set. On output COUNT(I), I=1,2,..,N, contains the
C   number of off-diagonal entries in row I including duplicates.
C   COUNT(0) contains the number of entries with one or both indices
C   out of range.
C FLAG is used for workspace to hold flags to permit duplicate entries
C   to be identified quickly.
C IWFR need not be set on input. On output it points to the first
C   unused location in IW.
C ICNTL Warning messages are printed on stream number ICNTL(2) if
C    ICNTL(2).GT.0 and ICNTL(3).GT.1.
C INFO need not be set on input. On output, INFO(1) has one of the
C   values:
C   0 No out-of-range index or duplicate entry found.
C   1 Out-of-range index found.
C   2 Duplicate entry found.
C   3 Out-of-range index found and duplicate entry found.
C  INFO(3) is set to the number of faulty entries.
C  INFO(4) is set to the number of off-diagonal duplicate entries.

    INTRINSIC MAX,MIN
C
C Local variables
    INTEGER I,J,K,L,LDIAG,LP,WP
C I Row index
C J Column index
C K Position in IRN, JCN, or IW.
C L Position in IW.
C LDIAG Level of monitor printing.
C WP Stream for printing.

C Set LDIAG and WP
    WP = ICNTL(2)
    LP = ICNTL(1)
    LDIAG = ICNTL(5)
    IF (WP.LT.0) LDIAG = 0

C Initialize INFO(1), FLAG, and COUNT.
    INFO(1) = 0
    DO 10 I = 1,N
      FLAG(I) = 0
      COUNT(I) = 0
   10 CONTINUE

C Count in INFO(3) the number of out-of-range entries, initialize FLAG,
C    and count in COUNT the numbers of entries in the rows.
    INFO(3) = 0
    DO 30 K = 1,NE
      I = IRN(K)
      J = JCN(K)
      IF (KIND.NE.1 .AND. I.GT.J) THEN
      INFO(1) = -19
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I8,A)')
   +    '**** Error return from ME57A ****  INFO(1) =',INFO(1),
   *    'Hermitian matrix but entry ',K,' in lower triangle'
      INFO(2) = K
      RETURN
      ENDIF
      IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) THEN
      IRN(K) = 0
      JCN(K) = 0
      INFO(3) = INFO(3) + 1
      INFO(1) = 1
      IF (INFO(3).EQ.1 .AND. LDIAG.GT.1) WRITE (WP,'(2A,I2)')
   +      '*** Warning message from subroutine ME57A ***',
   +      ' INFO(1) =',INFO(1)
      IF (INFO(3).LE.10 .AND. LDIAG.GT.1) WRITE (WP,'(3(I10,A))')
   +      K,'th entry (in row',I,' and column',J,') ignored'

      ELSE IF (PERM(I).LE.PERM(J)) THEN
      COUNT(I) = COUNT(I) + 1
      ELSE
      COUNT(J) = COUNT(J) + 1
      END IF

   30 CONTINUE
C
C Accumulate row counts in IPE ... one added for row length location.
    IPE(1) = COUNT(1) + 1
    DO 40 I = 2,N
      IPE(I) = IPE(I-1) + COUNT(I) + 1
   40 CONTINUE

C Run through putting the matrix entries in the right places. IPE is
C   used for holding running pointers and is left holding pointers to
C   row starts.
    DO 50 K = 1,NE
      I = IRN(K)
      J = JCN(K)
      IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 50
      IF (PERM(I).LE.PERM(J)) THEN
      IW(IPE(I)) = J
      IPE(I) = IPE(I) - 1
      ELSE
      IW(IPE(J)) = I
      IPE(J) = IPE(J) - 1
      END IF
   50 CONTINUE

C Remove duplicates
C IWFR points to the current position in the compressed set of rows.
    IWFR = 1
    INFO(4) = 0

C At the start of cycle I of this loop FLAG(J).LT.I for J=1,2,...,N.
C   During the loop FLAG(J) is set to I if an entry in column J is
C   found. This permits duplicates to be recognized quickly.
    DO 70 I = 1,N
      L = IPE(I)
      IPE(I) = IWFR

      DO 60 K = L + 1,L + COUNT(I)
      J = IW(K)
      IF (FLAG(J).NE.I) THEN
        FLAG(J) = I
        IWFR = IWFR + 1
        IW(IWFR) = J
      ELSE
C Count duplicates only once.
        IF (I.LT.J) INFO(4) = INFO(4) + 1
      END IF
   60   CONTINUE

      IF (IWFR.GT.IPE(I)) THEN
      IW(IPE(I)) = IWFR - IPE(I)
      IWFR = IWFR + 1
      ELSE
      IPE(I) = 0
      END IF

   70 CONTINUE

C Test whether duplicates found
    IF (INFO(4).GT.0) THEN
      INFO(1) = INFO(1) + 2
      IF (LDIAG.GT.1 .AND. WP.GE.0) WRITE (WP,'(A/I10,A)')
   +    '*** Warning message from subroutine ME57A ***',
   +    INFO(4),' off-diagonal duplicate entries found'
    END IF

    END

    SUBROUTINE ME57K(N, IPE, IW, LW, IWFR, PERM, IPS, NV, FLAG,
   *            NCMPA)
    INTEGER N,LW,IWFR,NCMPA
    INTEGER IPE(N)
    INTEGER IW(LW), PERM(N), IPS(N), NV(N), FLAG(N)
C
C Using a given pivotal sequence and a representation of the matrix that
C   includes only non-zeros of the strictly upper-triangular part
C   of the permuted matrix, construct tree pointers.
C
C N must be set to the matrix order. It is not altered.
C IPE(I) must be set to point to the position in IW of the
C   start of row I or have the value zero if row I has no off-
C   diagonal non-zeros. during execution it is used as follows.
C   If variable I is eliminated then IPE(I) points to the list
C   of variables for created element I. If element I is
C   absorbed into newly created element J then IPE(I)=-J.
C IW must be set on entry to hold lists of variables by
C   rows, each list being headed by its length. when a variable
C   is eliminated its list is replaced by a list of variables
C   in the new element.
C LW must be set to the length of IW. It is not altered.
C IWFR must be set to the position in IW of the first free variable.
C   It is revised during execution, continuing to have this meaning.
C PERM(K) must be set to hold the position of variable K in the
C   pivot order. It is not altered.
C IPS(I) need not be set by the user and will be used to hold the
C   inverse permutation to PERM.
C NV need not be set. If variable J has not been eliminated then
C   the last element whose leading variable (variable earliest
C   in the pivot sequence) is J is element NV(J). If element J
C   exists then the last element having the same leading
C   variable is NV(J). In both cases NV(J)=0 if there is no such
C   element. If element J has been merged into a later element
C   then NV(J) is the degree at the time of elimination.
C FLAG is used as workspace for variable flags.
C   FLAG(JS)=ME if JS has been included in the list for ME.
C
    INTEGER I,J,ML,MS,ME,IP,MINJS,IE,KDUMMY,JP
    INTEGER LN,JP1,JS,LWFR,JP2,JE

    EXTERNAL ME57F

C
C Initializations
    DO 10 I=1,N
      FLAG(I) = 0
      NV(I) = 0
      J = PERM(I)
      IPS(J) = I
   10 CONTINUE
    NCMPA = 0
C
C Start of main loop
C
    DO 100 ML=1,N
C ME=MS is the name of the variable eliminated and
C   of the element created in the main loop.
      MS = IPS(ML)
      ME = MS
      FLAG(MS) = ME
C
C Merge row MS with all the elements having MS as leading variable.
C IP points to the start of the new list.
      IP = IWFR
C MINJS is set to the position in the order of the leading variable
C   in the new list.
      MINJS = N
      IE = ME
      DO 70 KDUMMY=1,N
C Search variable list of element IE.
C JP points to the current position in the list being searched.
      JP = IPE(IE)
C LN is the length of the list being searched.
      LN = 0
      IF (JP.LE.0) GO TO 60
      LN = IW(JP)
C
C Search for different variables and add them to list,
C   compressing when necessary
      DO 50 JP1=1,LN
        JP = JP + 1
C Place next variable in JS.
        JS = IW(JP)
C Jump if variable has already been included.
        IF (FLAG(JS).EQ.ME) GO TO 50
        FLAG(JS) = ME
        IF (IWFR.LT.LW) GO TO 40
C Prepare for compressing IW by adjusting pointer to and length of
C   the list for IE to refer to the remaining entries.
        IPE(IE) = JP
        IW(JP) = LN - JP1
C Compress IW.
        CALL ME57F(N, IPE, IW, IP-1, LWFR, NCMPA)
C Copy new list forward
        JP2 = IWFR - 1
        IWFR = LWFR
        IF (IP.GT.JP2) GO TO 30
        DO 20 JP=IP,JP2
          IW(IWFR) = IW(JP)
          IWFR = IWFR + 1
   20     CONTINUE
   30     IP = LWFR
        JP = IPE(IE)
C Add variable JS to new list.
   40     IW(IWFR) = JS
        MINJS = MIN0(MINJS,PERM(JS)+0)
        IWFR = IWFR + 1
   50   CONTINUE
C Record absorption of element IE into new element.
   60   IPE(IE) = -ME
C Pick up next element with leading variable MS.
      JE = NV(IE)
C Store degree of IE.
      NV(IE) = LN + 1
      IE = JE
C Leave loop if there are no more elements.
      IF (IE.EQ.0) GO TO 80
   70   CONTINUE
   80   IF (IWFR.GT.IP) GO TO 90
C Deal with null new element.
      IPE(ME) = 0
      NV(ME) = 1
      GO TO 100
C Link new element with others having same leading variable.
   90   MINJS = IPS(MINJS)
      NV(ME) = NV(MINJS)
      NV(MINJS) = ME
C Move first entry in new list to end to allow room for length at
C   front. Set pointer to front.
      IW(IWFR) = IW(IP)
      IW(IP) = IWFR - IP
      IPE(ME) = IP
      IWFR = IWFR + 1
  100 CONTINUE

    RETURN
    END
C** end of ME57K**

    SUBROUTINE ME57F(N, IPE, IW, LW, IWFR, NCMPA)
C Compress lists held by ME57K in IW and adjust pointers
C   in IPE to correspond.
    INTEGER N,LW,IWFR,NCMPA
    INTEGER IPE(N)
C
    INTEGER   IW(LW)
C N is the matrix order. It is not altered.
C IPE(I) points to the position in IW of the start of list I or is
C   zero if there is no list I. On exit it points to the new position.
C IW holds the lists, each headed by its length. On output the same
C   lists are held, but they are now compressed together.
C LW holds the length of IW. It is not altered.
C IWFR need not be set on entry. On exit it points to the first free
C   location in IW.
C
    INTEGER I,K1,LWFR,IR,K,K2
    NCMPA = NCMPA + 1
C Prepare for compressing by storing the lengths of the
C   lists in IPE and setting the first entry of each list to
C   -(list number).
    DO 10 I=1,N
      K1 = IPE(I)
      IF (K1.LE.0) GO TO 10
      IPE(I) = IW(K1)
      IW(K1) = -I
   10 CONTINUE
C
C Compress
C IWFR points just beyond the end of the compressed file.
C LWFR points just beyond the end of the uncompressed file.
    IWFR = 1
    LWFR = IWFR
    DO 60 IR=1,N
      IF (LWFR.GT.LW) GO TO 70
C Search for the next negative entry.
      DO 20 K=LWFR,LW
      IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
      GO TO 70
C Pick up entry number, store length in new position, set new pointer
C   and prepare to copy list.
   30   I = -IW(K)
      IW(IWFR) = IPE(I)
      IPE(I) = IWFR
      K1 = K + 1
      K2 = K + IW(IWFR)
      IWFR = IWFR + 1
      IF (K1.GT.K2) GO TO 50
C Copy list to new position.
      DO 40 K=K1,K2
      IW(IWFR) = IW(K)
      IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN
    END
C--------------------------------------------------------------------
C        HSL 2000
C      --
C-         Copyright CCLRC Rutherford Appleton Laboratory
C      --
C--------------------------------------------------------------------
    SUBROUTINE ME57L(N, IPE, NV, IPS, NE, NA, NODE, PERM, NSTEPS,
   *            FILS, FRERE, ND, NEMIN, SUBORD)
C
C Tree search
C
C Given son to father tree pointers, reorder so that eldest son has
C   smallest degree and perform depth-first
C   search to find pivot order and number of eliminations
C   and assemblies at each stage.
    INTEGER N, NSTEPS
    INTEGER ND(N)
    INTEGER IPE(N), FILS(N), FRERE(N), SUBORD(N)
C
    INTEGER NV(N), IPS(N), NE(N), NA(N), NODE(N), PERM(N)
    INTEGER NEMIN
C N must be set to the matrix order. It is not altered.
C IPE(I) must be set equal to -(father of node I) or zero if
C    node is a root, if NV(I) > 0. If NV(I) = 0, then I is
C    subordinate variable of a supervariable and -IPE(I) points to
C    principal variable.  It is altered to point to its next
C    younger brother if it has one, but otherwise is not changed.
C NV(I) must be set to zero if variable is a subordinate variable
C    of a supervariable and to the degree otherwise.
C    NV is not altered.
C IPS(I) need not be set. It is used temporarily to hold
C    -(eldest son of node I) if it has one and 0 otherwise. It is
C    finally set to hold the position of node I in the order.
C NE(IS) need not be set. It is set to the number of variables
C    eliminated at stage IS of the elimination.
C NA(IS) need not be set. It is set to the number of elements
C    assembled at stage IS of the elimination.
C NODE (I) need not be set before entry. It is used during the code
C    to hold the number of subordinate variables for variable I and
C    on output it holds
C    the node (in dfs ordering) at which variable I is eliminated.
C    It is also defined for subordinate variables.
C PERM is set to the new permutation after dfs of tree.  PERM(I) is
C    the position of variable I in the pivot order.
C ND(IS) need not be set. It is set to the degree at stage IS of
C   the elimination.
C NSTEPS need not be set. It is set to the number of elimination steps.
C NEMIN is used to control the amalgamation process between
C     a son and its father (if the number of fully summed
C     variables of both nodes is smaller than NEMIN).
C SUBORD(I) need not be set. It holds the first subordinate variable
C     for variable I if I
C     is a principal variable and holds the next subordinate variable
C     if otherwise.  It is zero at the end of the chain.
C
    INTEGER I,IF,IS,NR,NR1,INS,INL,INB,INF,INFS,INSW
    INTEGER K,L,ISON,IN,IFSON,INO
    INTEGER INOS,IB,IL,INT
    INTEGER IPERM

C
C Initialize IPS and NE.
    DO 10 I=1,N
      IPS(I) = 0
      NE(I) = 0
      NODE(I) = 0
      SUBORD(I) = 0
   10 CONTINUE
C
C Set IPS(I) to -(eldest son of node I) and IPE(I) to next younger
C   brother of node I if it has one.
    NR = N + 1
    DO 50 I=1,N
      IF = -IPE(I)
      IF (NV(I).EQ.0) THEN
C I is a subordinate node, principal variable is IF
      IF (SUBORD(IF).NE.0) SUBORD(I) = SUBORD(IF)
      SUBORD(IF) = I
      NODE(IF) = NODE(IF)+1
      ELSE
C Node IF is the father of node I.
      IF (IF.NE.0) THEN
C IS is younger brother of node I.
C IPS(IF) will eventually point to - eldest son of IF.
        IS = -IPS(IF)
        IF (IS.GT.0) IPE(I) = IS
        IPS(IF) = -I
      ELSE
C I is a root node
        NR = NR - 1
        NE(NR) = I
      ENDIF
      ENDIF
   50 CONTINUE
C
C We reorganize the tree so that the eldest son has maximum number of
C variables.  We combine nodes when the number of variables in a son
C is greater than or equal to the number of variables in the father.
C If the eldest son has the maximum number of variables,
C and if a combination is possible, it has to be possible with
C the eldest son.
C
C FILS is just used as workspace during this reorganization and is reset
C afterwards.

    DO 999 I=1,N
     FILS(I) = IPS(I)
 999  CONTINUE

    NR1 = NR
    INS = 0
C Jump if all roots processed.
 1000 IF (NR1.GT.N) GO TO 1151
C Get next root
    INS = NE(NR1)
    NR1 = NR1 + 1
C Depth first search through eldest sons.
 1070 INL = FILS(INS)
    IF (INL.LT.0) THEN
     INS = -INL
     GO TO 1070
    ENDIF
C INS is leaf node.

 1080 IF (IPE(INS).LT.0) THEN
C INS is youngest son otherwise IPE value would be positive.
     INS     = -IPE(INS)
C INS is now the father of the reorganized son so we can
C   clear the pointer to the sons.
     FILS(INS) = 0
C Continue backtracking until we encounter node with younger brother.
     GO TO 1080
    ENDIF

    IF (IPE(INS).EQ.0) THEN
C INS is a root, check for next one.
     INS = 0
     GO TO 1000
    ENDIF
C INB is younger brother of INS.
    INB = IPE(INS)

C?? I think this test is the wrong way round
    IF (NV(INB).GE.NV(INS)) THEN
C?? So reversed
C   IF (NV(INS).GE.NV(INB)) THEN
     INS = INB
C Do depth first search from younger brother
     GO TO 1070
    ENDIF
C
C Exchange INB and INS
C Find previous brother of INS (could be the father)
C then we do depth first search with INS = INB
C
    INF = INB
 1090 INF = IPE(INF)
    IF (INF.GT.0) GO TO 1090
C -INF IS THE FATHER
    INF  = -INF
    INFS = -FILS(INF)
C INFS is eldest son of INF
    IF (INFS.EQ.INS) THEN
C  INS is eldest brother .. a role which INB now assumes
      FILS(INF) = -INB
      IPS(INF)  = -INB
      IPE(INS)  = IPE(INB)
      IPE(INB)  = INS
    ELSE
      INSW = INFS
 1100   INFS = IPE(INSW)
      IF (INFS.NE.INS) THEN
      INSW = INFS
      GO TO 1100
      ENDIF
      IPE(INS) = IPE(INB)
      IPE(INB) = INS
      IPE(INSW)= INB
    ENDIF
      INS    = INB
C Depth first search from moved younger brother
      GO TO 1070
C Set FRERE and FILS
 1151 DO 51 I=1,N
     FRERE(I) = IPE(I)
     FILS(I) = IPS(I)
 51   CONTINUE
C
C Depth-first search.
C IL holds the current tree level. Roots are at level N, their sons
C   are at level N-1, etc.
C IS holds the current elimination stage. We accumulate the number
C   of eliminations at stage is directly in NE(IS). The number of
C   assemblies is accumulated temporarily in NA(IL), for tree
C   level IL, and is transferred to NA(IS) when we reach the
C   appropriate stage IS.
    IS = 1
C I is the current node.
    I = 0
C IPERM is used as pointer to setting permutation vector
    IPERM = 1
    DO 160 K=1,N
      IF (I.GT.0) GO TO 60
C Pick up next root.
C Stop if all roots used (needed because of subordinate variables)
      IF (NR.GT.N) GO TO 161
      I = NE(NR)
      NE(NR) = 0
      NR = NR + 1
      IL = N
      NA(N) = 0
C Go to son for as long as possible, clearing father-son pointers
C   in IPS as each is used and setting NA(IL)=0 for all levels
C   reached.
   60   CONTINUE
      DO 70 L=1,N
      IF (IPS(I).GE.0) GO TO 80
      ISON = -IPS(I)
      IPS(I) = 0
      I = ISON
      IL = IL - 1
      NA(IL) = 0
   70   CONTINUE
   80   CONTINUE
C?? Do we want to expand for subordinate variables
C Record position of node I in the order.
      IPS(I) = K
C Add number of subordinate variables to variable I
      NE(IS) = NE(IS) + NODE(I) + 1
      IF (IL.LT.N) NA(IL+1) = NA(IL+1) + 1
      NA(IS) = NA(IL)
      ND(IS) = NV(I)
      NODE(I) = IS
      PERM(I) = IPERM
      IPERM = IPERM + 1
C Order subordinate variables to node I
      IN = I
  777   IF (SUBORD(IN).EQ.0) GO TO 778
      IN = SUBORD(IN)
      NODE(IN) = IS
      PERM(IN) = IPERM
      IPERM = IPERM + 1
      GO TO 777
C Check for static condensation
  778   IF (NA(IS).NE.1) GO TO 90
      IF (ND(IS-1)-NE(IS-1).EQ.ND(IS)) GO TO 100
C Check for small numbers of eliminations in both last two steps.
   90   IF (NE(IS).GE.NEMIN) GO TO 110
      IF (NA(IS).EQ.0) GO TO 110
      IF (NE(IS-1).GE.NEMIN) GO TO 110

C Combine the last two steps
  100   NA(IS-1) = NA(IS-1) + NA(IS) - 1
      ND(IS-1) = ND(IS) + NE(IS-1)
      NE(IS-1) = NE(IS) + NE(IS-1)
      NE(IS) = 0
      NODE(I) = IS-1
C Find eldest son (IFSON) of node I (IS)
C Note that node I must have a son (node IS-1 is youngest)
      IFSON = -FILS(I)
C Now find youngest son INO (he is node IS-1)
      IN = IFSON
 102    INO = IN
      IN =  FRERE(IN)
      IF (IN.GT.0) GO TO 102
C Cannot be root node .. so points to father
C Merge node IS-1 (INO) into node IS (I)
      NV(INO) = 0
C IPE already set .. was father pointer now principal variable pointer
C Now make subsidiary nodes of INO into subsidiary nodes of I.
C Subordinate nodes of INO become subordinate nodes of I
      IN = I
  888   IF (SUBORD(IN).EQ.0) GO TO 889
      IN = SUBORD(IN)
      NODE(IN) = IS-1
      GO TO 888
  889   SUBORD(IN) = INO
      IN = INO
      IF (SUBORD(IN).EQ.0) GO TO 887
      IN = SUBORD(IN)
      IPE(IN) = -I
  887   CONTINUE

C INOS is eldest son of INO
    INOS = -FILS(INO)

C Find elder brother of node INO
C First check to see if he is only son
    IF (IFSON.EQ.INO) GO TO 107
    IN = IFSON
 105  INS = IN
    IN =  FRERE(IN)
    IF (IN.NE.INO) GO TO 105
C INS is older brother .. make him brother of first son of INO (ie INOS)
C and  make INOS point to I now as father.
C Jump if there is no son of INO
      IF (INOS.EQ.0) THEN
C Elder brother of INO just points to (new) father.
      FRERE(INS) = -I
      GO TO 120
      ELSE
      FRERE(INS) =  INOS
      ENDIF
C INT is youngest brother of INOS.  Make him point to (new) father.
 107    IN = INOS
      IF (IN.EQ.0) GO TO 120
 108    INT = IN
      IN =  FRERE(IN)
      IF (IN.GT.0) GO TO 108
      FRERE(INT) = -I
      GO TO 120
  110   IS = IS + 1
  120   IB = IPE(I)
      IF (IB.GE.0) THEN
C Node I has a younger brother or is a root
      IF (IB.GT.0) NA(IL) = 0
      I = IB
      GO TO 160
      ELSE
C I has no brothers. Go to father of node I
      I = -IB
      IL = IL + 1
      ENDIF
  160 CONTINUE
  161 NSTEPS = IS - 1
    RETURN
    END
    SUBROUTINE ME57M(KIND,N,NE,IRN,JCN,MAP,IRNPRM,
   +            LROW,PERM,COUNT,IDIAG)
C
C This subroutine is called by ME57A/AD and generates the map that
C   reorders the user's input so that the upper triangle of the
C   permuted matrix is held by rows. No check is made for duplicates.
C
    INTEGER KIND,N,NE
C IRNPRM(N+NE) has this dimension to include possibility of no
C    diagonals in input.
    INTEGER IRN(NE),JCN(NE),MAP(NE),IRNPRM(N+NE),LROW(N),PERM(N),
   +      COUNT(N),
   +      IDIAG(N)
C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C N must be set to the matrix order. It is not altered.
C NE must be set to the number of entries input. It is not altered.
C IRN(K) and JCN(K), K=1,2,...,NE must be set to the row and column
C   numbers of the entries. If entry (IRN(K),JCN(K)) lies in the
C   lower triangular part of the permuted matrix, the values of
C   IRN(K) and JCN(K) are interchanged. Otherwise, these arrays are
C   not changed.
C MAP need not be set on input and on return holds the positions of
C   entries when the permuted upper triangle is ordered by rows.
C LROW need not be set. On return, LROW(I),I=1,N holds the number of
C   entries in row I of the permuted matrix.
C PERM(I) must be set to the position of variable I in the
C   pivot order, I=1,2,...,N.
C COUNT is used for workspace. It is set to row counts and then
C   accumulated row counts.
C IDIAG is used for workspace. It is used as pointer to diagonal entry
C   that is first in the row.
C
C Local variables
    INTEGER EXPNE,I,J,K
C I Row index
C J Column index
C K Position in IRN or JCN.

C Accumulate row counts in COUNT, and interchange row and column
C   numbers where necessary.
    DO 10 I = 1,N
C Set to 1 since diagonal will always be present as first entry.
      COUNT(I) = 1
C IDIAG used first as flag to identify duplicate diagonals
      IDIAG(I) = 0
   10 CONTINUE

C EXPNE counts number of entries plus missing diagonals
    EXPNE = NE + N
    DO 20 K = 1,NE

      I = IRN(K)
      J = JCN(K)

      IF (MAX(I,J).GT.N .OR. MIN(I,J).LT.1) THEN
      EXPNE = EXPNE - 1
      GO TO 20
      ENDIF

C Check for duplicate diagonals
      IF (I.EQ.J) THEN
      I = PERM(I)
      IF (IDIAG(I).GE.1) THEN
        COUNT(I) = COUNT(I) + 1
        IDIAG(I) = IDIAG(I) + 1
      ELSE
        IDIAG(I) = 1
        EXPNE = EXPNE - 1
      ENDIF
      GO TO 20
      ENDIF

      IF (PERM(I).LT.PERM(J)) THEN
      I = PERM(I)
      COUNT(I) = COUNT(I) + 1
      ELSE
      J = PERM(J)
      COUNT(J) = COUNT(J) + 1
      END IF

   20 CONTINUE
C
C Store row counts in LROW and accumulate row counts in COUNT.
    LROW(1) = COUNT(1)
    IDIAG(1) = MAX(IDIAG(1),1)
    DO 30 I = 2,N
      LROW(I) = COUNT(I)
C COUNT(I) set to point to position of last entry in row I of permuted
C     upper triangle.
      COUNT(I) = COUNT(I-1) + LROW(I)
      IDIAG(I) = COUNT(I-1) + MAX(IDIAG(I),1)
   30 CONTINUE

C Set diagonal entries in IRNPRM.  This is done separately because some
C   diagonals may not be present in the users input.
    DO 35 I = 1,N
      K = PERM(I)
      IRNPRM(IDIAG(K)) = I
   35 CONTINUE
C
C Run through putting the entries in the right places. COUNT is used for
C   holding running pointers and is left holding pointers to row
C   starts.
C Count number of entries in expanded matrix (allowing for non-input
C   diagonals)
    DO 40 K = 1,NE
      I = IRN(K)
      J = JCN(K)
      IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) THEN
      MAP(K) = 0
      GO TO 40
      ENDIF
      I = PERM(IRN(K))
      J = PERM(JCN(K))
      IF (I.EQ.J) THEN
      MAP(K) = IDIAG(I)
      IRNPRM(IDIAG(I)) = IRN(K)
      IDIAG(I) = IDIAG(I) - 1
      ELSE
      IF (I.GT.J) THEN
C Entry has been mapped to lower triangle
        IF (KIND.EQ.1) THEN
          MAP(K) = COUNT(J)
        ELSE IF (KIND.NE.1) THEN
          MAP(K) = - COUNT(J)
        ENDIF
        IRNPRM(COUNT(J)) = IRN(K)
        COUNT(J) = COUNT(J) - 1
      ELSE
        MAP(K) = COUNT(I)
        IRNPRM(COUNT(I)) = JCN(K)
        COUNT(I) = COUNT(I) - 1
      ENDIF
C       II = MIN(I,J)
C       MAP(K) = COUNT(II)
C       COUNT(I) = COUNT(II) - 1
      ENDIF
   40 CONTINUE

C Set number of entries in expanded matrix
    IDIAG(1) = EXPNE
    RETURN

    END
    SUBROUTINE ME57N(N,LENR,NA,NE,ND,NSTEPS,LSTKI,LSTKR,
   *            INFO,RINFO)
C
C Storage and operation count evaluation.
C
C Evaluate number of operations and space required by factorization
C   using ME57B.  The values given are exact only if no numerical
C   pivoting is performed.
C
C N must be set to the matrix order. It is not altered.
C LENR is number of entries in upper triangle of each row of permuted
C   matrix. It includes diagonal (no duplicates) and duplicates
C   off-diagonal but excludes any out-of-range entries.
C NA,NE,ND must be set to hold, for each tree node, the number of stack
C   elements assembled, the number of eliminations and the size of
C   the assembled front matrix respectively.  They are not altered.
C NSTEPS must be set to hold the number of tree nodes. It is not
C   altered.
C LSTKI is used as a work array by ME57N.
C LSTKR is used as a work array by ME57N.
C
C Counts for operations and storage are accumulated in variables
C   OPS,OPSASS,NRLTOT,NIRTOT,NRLNEC,NIRNEC,NRLADU,NIRADU.
C OPS number of multiplications and additions during factorization.
C OPSASS number of multiplications and additions during assembly.
C NRLADU,NIRADU real and integer storage respectively for the
C   matrix factors.
C NRLTOT,NIRTOT real and integer storage respectively required
C   for the factorization if no compresses are allowed.
C NRLNEC,NIRNEC real and integer storage respectively required for
C   the factorization if compresses are allowed.
C MAXFRT is maximum front size

C   .. Scalar Arguments ..
    INTEGER N,NSTEPS
C   ..
C   .. Array Arguments ..
    INTEGER LENR(N),LSTKI(N),LSTKR(N),NA(NSTEPS),
   +      ND(NSTEPS),NE(NSTEPS),INFO(40)
    REAL RINFO(20)
C   ..
C   .. Local Scalars ..
    INTEGER I,IORG,ISTKI,ISTKR,ITOP,ITREE,JORG,K,
   +      LSTK,NASSR,NELIM,NFR,NSTK,NTOTPV,NZ1,NZ2
    REAL DELIM
C   ..
C   .. Intrinsic Functions ..
    INTRINSIC MAX
C   ..
    REAL OPS,OPSASS
    INTEGER NIRADU,NIRNEC,NIRTOT,NRLADU,NRLNEC,NRLTOT,MAXFRT

C   ..
C   .. Executable Statements ..
C
C Accumulate number of nonzeros with indices in range in NZ1.
C   Duplicates on the diagonal are ignored but NZ1 includes any
C   diagonals not present on input and duplicates off the diagonal.
    NZ1 = 0
    DO 40 I = 1,N
      NZ1 = NZ1 + LENR(I)
   40 CONTINUE
    NZ2 = NZ1
C ISTKR,ISTKI Current number of stack entries in
C   real and integer storage respectively.
C OPS,OPSASS,NRLADU,NIRADU,NIRTOT,NRLTOT,NIRNEC,NRLNEC,NZ2 are defined
C   above.
C NZ2 Current number of original matrix entries not yet processed.
C NTOTPV Current total number of rows eliminated.
C ITOP Current number of elements on the stack.
    ISTKI = 0
    ISTKR = 0
    OPS = 0.0E0
    OPSASS = 0.0E0
    NRLADU = 0
C One location is needed to record the number of blocks actually used.
    NIRADU = 3
C Initialize to what is required in ME57B (as opposed to ME57O).
    NIRTOT = NZ1+N+5
    NRLTOT = NZ1
    NIRNEC = NZ2+N+5
    NRLNEC = NZ2
    NTOTPV = 0
    ITOP = 0
    MAXFRT = 0
C
C Each pass through this loop processes a node of the tree.
    DO 100 ITREE = 1,NSTEPS
      NELIM = NE(ITREE)
      DELIM = NELIM
      NFR = ND(ITREE)
      MAXFRT = MAX(MAXFRT,NFR)
      NSTK = NA(ITREE)
C Adjust storage counts on assembly of current frontal matrix.
      NASSR = NELIM*(NELIM+1)/2 + NFR*NFR
C Data for no compresses so use original number, NZ1.
      NRLTOT = MAX(NRLTOT,NRLADU+NASSR+ISTKR+NZ1)
C Data for compresses so use current number, NZ2.
      NRLNEC = MAX(NRLNEC,NRLADU+NASSR+ISTKR+NZ2)

C Decrease NZ2 by the number of entries in rows being eliminated at
C   this stage.
      DO 70 IORG = 1,NELIM
      JORG = NTOTPV + IORG
      OPSASS = OPSASS + LENR(JORG)
      NZ2 = NZ2 - LENR(JORG)
   70   CONTINUE

      NTOTPV = NTOTPV + NELIM

C Remove elements from the stack.  There are ITOP elements on the
C   stack with the appropriate entries in LSTKR and LSTKI giving
C   the real and integer storage respectively for each stack
C   element.
      DO 80 K = 1,NSTK
      LSTK = LSTKR(ITOP)
      ISTKR = ISTKR - LSTK
      OPSASS = OPSASS + LSTK
      LSTK = LSTKI(ITOP)
      ISTKI = ISTKI - LSTK
      ITOP = ITOP - 1
   80   CONTINUE

C Accumulate nonzeros in factors and number of operations.
      NRLADU = NRLADU + (NELIM*(NELIM+1))/2 + (NFR-NELIM)*NELIM
      NIRADU = NIRADU + 2 + NFR
      OPS = OPS + (DELIM* (12*NFR+6*NFR*NFR - (DELIM+1)*
   +      (6*NFR+6-(2*DELIM+1))))/6 + DELIM

      IF (NFR.GT.NELIM) THEN
C Stack remainder of element.
      ITOP = ITOP + 1
      LSTKR(ITOP) = ((NFR-NELIM)*(NFR-NELIM+1))/2
      LSTKI(ITOP) = NFR - NELIM + 1
      ISTKI = ISTKI + LSTKI(ITOP)
      ISTKR = ISTKR + LSTKR(ITOP)
      ENDIF

C Adjust integer counts to stack elements and allow for next front.
      IF (ITREE.EQ.NSTEPS) THEN
      NIRTOT = MAX(NIRTOT,NIRADU+ISTKI+NZ1)
      NIRNEC = MAX(NIRNEC,NIRADU+ISTKI+NZ2)
      ELSE
      NIRTOT = MAX(NIRTOT,NIRADU+(N-NTOTPV+2)+ISTKI+NZ1)
      NIRNEC = MAX(NIRNEC,NIRADU+(N-NTOTPV+2)+ISTKI+NZ2)
      ENDIF

  100 CONTINUE
C

C Set INFO and RINFO
    INFO(5)   = NRLADU
    INFO(6)   = NIRADU
    INFO(7)   = MAXFRT
    INFO(8)   = NSTEPS
    INFO(9)   = NRLTOT
    INFO(10)  = NIRTOT
    INFO(11)  = NRLNEC
    INFO(12)  = NIRNEC
    RINFO(1)  = OPSASS
    RINFO(2)  = OPS

    RETURN
    END
    SUBROUTINE ME57O(KIND,N,NE,A,LA,IW,LIW,LROW,PERM,NSTEPS,NSTK,
   +            NODE,DIAG,SCHNAB,PPOS,CNTL,ICNTL,INFO,RINFO,
   +            HOLD,BIGA)
C
C Factorization subroutine
C
C This subroutine operates on the input matrix ordered into a tentative
C   pivot order by ME57B and produces the matrices U and inv(D)
C   of the factorization A = (U trans) D U, where D is a block
C   diagonal matrix with blocks of order 1 and 2. Gaussian elimination
C   is used with pivots of order 1 and 2, chosen according to the
C   tentative pivot order unless stability considerations
C   require otherwise.
    INTEGER KIND,N,NE,LA
C SCHNAB is in fact of dimension 5, but must set to * for Fujitsu (and
C   maybe other) compilers on very small matrices.
    COMPLEX A(LA),DIAG(N),SCHNAB(*)
    REAL CNTL(5),RINFO(20),BIGA
    INTEGER LIW,IW(LIW),LROW(N),PERM(N),NSTEPS,NSTK(NSTEPS),
   +      NODE(N),PPOS(N),ICNTL(20),INFO(40),HOLD(40)
C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C N   must be set to the order of the matrix. It is not altered.
C NE  must be set to the number of entries in the upper triangle of
C   the permuted matrix. It is not altered.
C A   must be set so that the upper triangle of the permuted matrix
C   is held by rows in positions LA-NE+1 to LA. Explicit zero
C   diagonals are stored. Duplicate entries are permitted
C   and are summed. During the computation, active frontal matrices
C   are held in A by rows. The working front is held in full form.
C   Stacked elements are held by rows in packed form.
C   On exit, entries 1 to INFO(10) of A hold real information on the
C   factors and should be passed unchanged to ME57C. For each block,
C   the factorized block pivot precedes the rows of the out-of-pivot
C   part held as a rectangular matrix by rows.
C      The factorized pivot has the form:
C          -1  T
C        L  D   L
C
C   where L is unit lower triangular, D is block diagonal with blocks
C   of size 1 or 2. L and the diagonal part of D is held packed by
C   columns and the off-diagonal entries of the 2*2 blocks of D are
C   held from position IW(1).
C LA  length of array A. A value for LA sufficient for the
C   tentative pivot sequence will have been provided by ME57A
C   in variable INFO(11). LA is not altered.
C IW  must be set on input so that IW(I+LIW-NE) holds the column
C   index of the entry in A(I+LA-NE) for I = 1,..., NE.
C   On exit, entries 1 to INFO(16) hold integer information on the
C   factors and should be passed unchanged to ME57C.
C   IW(1) will be set to one greater than the number of entries in the
C   factors.
C   IW(2) points to the end of the factorization.
C   IW(3) will be set to the number of block pivots actually used;
C   this may be different from NSTEPS since numerical considerations
C   may prevent us choosing pivots at some stages.
C   Integer information on each block pivot row follows. For each
C   block pivot row, we have:
C     * no. of columns,
C     * no. of rows,
C     * list of column indices. The column indices for a
C       2x2 pivot are flagged negative.
C   During the computation, the array is used to hold indexing
C   information on stacked elements.  IW stores the number of
C   variables then a list of variables in the element.
C LIW must be set to the length of array IW. A value sufficient for the
C   tentative pivot sequence will have been provided by ME57A in
C   variable INFO(12). LIW is not altered.
C LROW  must be set so that LROW(I) holds the number of entries in row
C   I (in permuted order) of the incoming matrix. LROW is not altered.
C PERM  must be set so that PERM(I) holds the variable that is Ith in
C   the tentative pivot order generated by ME57A. PERM is not
C   altered.
C NSTEPS must be set to the number of nodes in the tree from the
C   analysis. It is the length of array NSTK. Its
C   value will never exceed N. It is not altered.
C NSTK must be set so that NSTK(I) holds the number of
C   stacked elements to be assembled at tree node I.
C NODE must be unchanged since return from ME57A. NODE(I) gives
C   the tree node at which variable I was eliminated in the
C   analysis phase. It is of length N and is not altered.
C DIAG is only accessed if ICNTL(7) is equal to 4.
C   In that case it must be set to the values of the diagonals of
C   matrix.
C SCHNAB is only accessed if ICNTL(7) is equal to 4. It is used to hold
C   parameters for the Schnabel-Eskow modification and max/min values
C   of current diagonal.  Specifically (using notation of S-E):
C   SCHNAB(1) == GAMMA
C   SCHNAB(2) == TAUBAR (in fact root of TAUBAR)
C   SCHNAB(3) == MU
C   SCHNAB(4) == Max entry on diag
C   SCHNAB(5) == Min entry on diag
C This option is only invoked if the matrix is Hermitian
C PPOS  is integer work array of dimension N. If I is a variable in
C   the current front, PPOS(I) is used to indicate its position in the
C   front. For any other uneliminated variable, PPOS(I) is set to N+1.
C CNTL must be set (perhaps by ME57I) so that CNTL(1) holds the
C   pivot threshold and CNTL(2) holds the pivot tolerance.
C ICNTL must be set (perhaps by ME57I).  Entries of ICNTL accessed
C   by ME57O are:
C ICNTL(2) is output unit for warning messages.
C ICNTL(7) is used to control pivoting.  With the default value of 1,
C   1 x 1 and 2 x 2 pivots are used subject to passing a threshold
C   tolerance.  If ICNTL(7) is greater than 1 only 1 x 1 pivots will
C   be used. If ICNTL(7) equal to 2,
C   the subroutine will exit immediately a sign change or zero pivot
C   is detected.  If ICNTL(7) is equal to 3, the subroutine will
C   continue the factorization unless a zero
C   pivot is detected.  If ICNTL(7) is equal to 4, the diagonal of
C   the matrix will be modified so that all pivots are of the same.
C ICNTL(8) is used to control whether, on running out of space, the
C   subroutine exits with an error return (ICNTL(8) = 0), or
C   whether it saves some internal variables so that
C   larger arrays can be allocated and the computation restarted
C   from the point at which it failed.
C ICNTL(11) is the block size used by the Level 3 BLAS (default 32).
C RINFO(3) will be set to the number of floating-point operations
C   required for the assembly.
C RINFO(4) will be set to the number of floating-point operations
C   required for the factorization.  RINFO(5) will be set to the
C   number of extra flops needed for the use of GEMM.
C INFO(1)  holds a diagnostic flag. It need not be set on entry. A zero
C   value on exit indicates success. Possible nonzero values are
C      -3  insufficient storage for A.
C      -4  insufficient storage for IW.
C      -5  zero pivot found when ICNTL(7) = 2 or 3.
C      -6  change in sign of pivots when ICNTL(7) = 2.
C        KIND out of range
C      +4  matrix is singular
C     +10  factorizations pauses because insufficient real space
C     +11  factorizations pauses because insufficient integer space
C INFO(40) is used to accumulate number of reals discarded if
C   elimination continues without restart when space is exhausted.
C

C INFO(32) is set to number of zeros in the triangle of the factors
C INFO(33) is set to number of zeros in the rectangle of the factors
C INFO(34) is set to number of zero columns in rectangle of the factors
C Needed to compute these
    INTEGER ZCOL,RPOS

C
C Constants
    REAL ZERO,HALF,ONE
    PARAMETER (ZERO=0.0E0,HALF=0.5E0,ONE=1.0E0)
    COMPLEX CONE
    PARAMETER (CONE=(1.0E0,0.0E0))
C
C Local variables
    INTEGER AINPUT
    COMPLEX AMULT1,AMULT2
    REAL AMAX
    INTEGER APOS,APOSA,APOSB,APOSBB,APOSBK,APOSC,APOSI,APOSJ,APOSM,
   +      APOS1,APOS2,APOS3,APOS4,ASTK,ATRASH,BLK
    REAL DELTA
    COMPLEX DETPIV
    INTEGER ELT
    REAL FLOPSA,FLOPSB,FLOPSX
    INTEGER I,I1,IASS,IBEG,IELL,IEND,IEXCH,IINPUT,INTSPA,
   +      IORG,IPIV,IPOS,IROW,ISNPIV,ISTK,ISWOP,IWNFS,IWPOS,
   +      J,JAY,JA1,JCOL,JJ,JJJ,JMAX,J1,J2,K,
   +      KB,KBLK,KCT,KR,KROW,K1,K2,L,LASPIV,LDIAG,LIELL,
   +      LP,LPIV, NBSTATIC
    LOGICAL LASTBK,LTWO
    INTEGER MAXFRT
    REAL MAXPIV
    INTEGER NASS,NBLK,NBLOC,NCMPBI,NCMPBR,NEIG,NELL,NFRONT,NIRBDU
    REAL NORMJ
    INTEGER NTWO
C LSTAT is .TRUE. is we are using static pivoting
    LOGICAL SCHUR,LSTAT
    INTEGER MPIV,NPIV,NPOTPV,NRLBDU,NSC1,NST,
   +      NSTACK(2),NSTKAC(2),NTOTPV,
   +      NUMORG,OFFDAG,PHASE,PIVBLK
    COMPLEX PIVOT
    INTEGER PIVSIZ,POSELT,POSPV1,POSPV2,PTRA,PTRIRN,RLSPA,
   +      SIZBLK,SIZC,SIZF,TRLSPA,TINSPA,TOTSTA(2),WP,ZCOUNT
    REAL RMAX,TMAX,TOL,UU,ULOC,UTARG,STCTOL
    COMPLEX SWOP,CMAX

C AINPUT is the position of the first entry of the original matrix
C   reals since the last compress.  It is reset to the current
C   part being processed when ME57P is called.
C AMAX is used to record the largest entry in a row.
C CMAX is the value of the complex number corresponding to AMAX
C AMULT1, AMULT2 are used to hold multipliers.
C   Also used as temporary variables.
C APOS is a pointer to the start of the current front in A.
C APOSA holds the index of the current entry of the matrix A used to
C   form the Schur complement as C = A*B
C APOSB holds the index of the current entry of the matrix B used to
C   form the Schur complement as C = A*B
C APOSBB is a pointer to the start of the buffer before current front
C   in A. It is different from APOS because of the way we store the
C   factors and prevents us overwriting when generating them.
C APOSC holds the index of the current entry of the matrix C used to
C   form the Schur complement as C = A*B
C APOSI is set to the beginning of a row in the working front or just
C    ahead of this position.
C APOSJ is set to the beginning of row JMAX in working front.
C APOSM holds the index of the current entry of the multiplier matrix
C   for the matrix B used to form the Schur complement as C = A*B for
C   the part of the pivot rows outside the pivot block.
C APOS1 is a position in the array A.
C APOS2 is a position in the array A.
C APOS3 is the position in A of the first entry in the copy of U.
C APOS4 is the position in A of the first entry in the Schur complement.
C ASTK indicates position immediately before first stack entry.
C ATRASH is used as limit on space in A being set to zero.
C BLK is DO loop variable for blocks in pivot block.
C DELTA is the amount added to diagonal when in Phase 2 of matrix
C   modification.
C DETPIV is the value of the determinant of the 2x2 pivot or
C   candidate pivot.  It is real when matrix is Hermitian (only).
C ELT is a DO loop variable indicating the element being processed.
C FLOPSA  counts floating-point operations for assembly
C FLOPSB  counts floating-point operations for elimination
C FLOPSX  counts extra flops required by use of GEMM
C I is a DO loop variable.
C I1 is used to hold limit of DO loop.
C IASS is the index of the current tree node.
C IBEG is the position of the beginning of a row in A, used when
C   performing elimination operations on the front.
C IELL Current element being assembled starts in position IELL of IW.
C IEND is the position of the end of a row in A, used when
C   performing elimination operations on the front.
C IEXCH is used to hold the contents of an array entry when doing a
C   swop.
C IINPUT is the position of the first entry of the original matrix
C   integers since the last compress.  It is reset to the current
C   part being processed when ME57P is called. with REAL = .FALSE.
C INTSPA is amount of integer workspace needed to progress to current
C   point in elimination.
C IORG is a DO loop variable indicating which current row from the
C   incoming matrix is being processed.
C IPIV is the current relative position of the pivot search.
C IPOS is used to hold (temporarily) a position in IW.
C IROW is a DO loop variable used when scanning a row of A.
C ISNPIV is +1 if first pivot is positive and -1 if it is negative.
C ISTK points to the bottom of the stack in IW (needed by compress).
C    It indicates position immediately before first stack entry.
C ISWOP is used when swopping two integers.
C IWNFS points to the first free location for a variable that is not
C   fully summed.
C IWPOS points to the first free position for factors in IW.
C J is a temporary variable.
C JA1 is a temporary index.
C JCOL is used as an index into array A.
C JJ and JJJ are Do loop indices.
C JMAX is the relative column index in the front of the largest
C   off-diagonal in the fully summed part of the prospective pivot
C   row.
C J1 and J2 are pointers to the beginning and end of a row segment
C   in the array IW.  They are also used as running pointers in IW.
C K is a temporary variable.
C KBLK Blocked GEMM is performed on a block KBLK by KBLK matrix.
C KCT counts the number of unchecked candidates in a pivot sweep.
C KR is pointer to the current row in the assembled block being
C   tested for a potential pivot.
C KROW is a DO loop variable.
C K1 and K2 are used as running indices for array A.
C KB is block row index for blocked GEMM.
C L is a temporary variable.
C LASPIV is set to value of NPIV at end of previous block of pivots.
C LIELL is the order of the size of the reduced matrix from the front.
C   This is the order of the stacked matrix.
C LPIV is the number of pivots selected in a block pivot.
C LASTBK is flag to indicate when we are processing last block of
C     pivot block.
C LTWO  is logical variable used to indicate if current pivot is a 2 x 2
C     pivot.
C MAXFRT is the maximum front size encountered so far.
C MAXPIV is the largest of two diagonal entries of a 2 x 2 pivot.
C NASS holds the number of fully assembled variables in
C   the newly created element.
C NBLK is the number of block pivots used.
C NBLOC Number of rows of the Schur complement calculated by each GEMM
C   reference. Set to ICNTL(11).
C NCMPBR, NCMPBI are the number of compresses on real and integer space
C   respectively.
C NEIG is number of negative eigenvalues detected.
C NELL is used to hold the current number of son elements.
C NFRONT is the number of variables in the front.
C NIRBDU is number of integer entries in factors.
C NORMJ is used in matrix modification for 1-norm of off-diagonals.
C NTWO is the number of two by two full pivots used.
C SCHUR if set to .TRUE. then the Schur complement will be
C    generated using Level 3 BLAS.
C MPIV is the number of pivots so far chosen in the current block.
C NPIV is the number of pivots so far chosen at the current node.
C NPOTPV is the total number of potential pivots so far.  Variables
C   PERM(1), PERM(2), .... PERM(NPOTPV) are fully assembled.
C NRLBDU is number of real entries in factors.
C NST temporary to hold NSC1+1 if NSC1 > 0, 0 otherwise.
C NSTACK(I), I = 1,2 hold the number of active entries on the
C   real/integer stack.
C NSTKAC(I), I =1,2 hold the number of entries on the real/integer
C   stack and original matrix after a compress.
C NTOTPV is the total number of pivots selected. This is used
C   to determine whether the matrix is singular.
C NUMORG is the number of variables in the tentative pivot from the
C   incoming original rows.
C OFFDAG is the position in A of the off-diagonal entry of a 2 x 2
C      pivot.
C PIVBLK Number of rows of each block when using GEMM in pivot block.
C   Set to minimum of NBLOC and NASS.
C PIVOT is a temporary variable used to hold the value of the current
C      pivot.
C PIVSIZ is order of current pivot (has value 1 or 2).
C POSELT is a pointer into the current element being assembled.
C POSPV1 is the position in A of a 1 x 1 pivot or the first diagonal
C   of a 2 x 2 pivot.
C POSPV2 is the position in A of the second diagonal of a 2 x 2 pivot.
C PTRA  points to the next original row in A.
C PTRIRN points to the next original row in IW.
C RLSPA is amount of real workspace needed to progress to current
C   point in elimination.
C SIZBLK Number of rows in current block when using GEMM in pivot block.
C   Set to minimum of PIVBLK and NASS - NPIV
C SIZC is set to number of rows in remainder of pivot block after
C   blocking is done.
C SIZF is set to number of rows in remainder of pivot row (to NFRONT)
C   after blocking is done.
C TOTSTA(I), I =1,2 hold the number of entries on the stack and
C   original matrix.
C TRLSPA is amount of real workspace needed to progress to current
C   point in elimination if no compresses on the data are performed.
C RMAX is used to record the largest entry in a row.
C SWOP is used when swopping two reals.
C TMAX is used to record the largest entry in a row.
C TOL is the tolerance against which singularity is judged.
C   If static pivoting is used, then TOL is the value for this.
C UU is a local variable used to hold threshold parameter.  Its value is
C   between 0 and 0.5.
C ZCOUNT is number of "zero" rows in current front (used when
C   ICNTL(16) is equal to 1).
C
C Procedures
C ME57P compresses arrays.
C ME57W  adjusts signs in factors, moves the off-diagonal entries of
C    full 2x2 pivots and updates counts.
    REAL FD15A

C?? To identify bug
C   LOGICAL LCASE
C   COMMON /CCASE/LCASE

    INTRINSIC MIN,MAX,ABS,REAL,CONJG
    EXTERNAL CGEMM,FD15A,ME57P,ME57W

C
C Initialization.
    NBLOC = ICNTL(11)
    TOL = CNTL(2)
    LP = ICNTL(1)
    WP = ICNTL(2)
    LDIAG = ICNTL(5)
    INFO(40) = 0
C A local variable UU is used for the threshold parameter, so that
C   CNTL(1) will remain unaltered.
    UU = MIN(CNTL(1),HALF)
    UU = MAX(UU,ZERO)

    LSTAT = .FALSE.
C Check if static pivoting option is on
    IF (CNTL(4).GT.ZERO) THEN
C LSTAT is not now set until number of delayed pivots is CNTL(5)*N
      IF (CNTL(5).EQ.ZERO) LSTAT = .TRUE.
      UTARG = SQRT(UU/CNTL(4))*CNTL(4)
      STCTOL = BIGA*CNTL(4)
C     TOL = STCTOL
    ENDIF

C Action if we are returning in the middle of the factorization
    IF (HOLD(1).GT.0) THEN
      INFO(1) = 0
      NBLK = HOLD(2)
      NTWO = HOLD(3)
      INFO(23) = HOLD(4)
      NCMPBR = 0
      NCMPBI = 0
      NEIG   = HOLD(6)
      MAXFRT = HOLD(7)
C Test compiler by commenting this out
      IWPOS  = HOLD(8)
      APOS   = HOLD(9)
      APOSBB = HOLD(10)
      NSTKAC(1) = HOLD(11)
      NSTKAC(2) = HOLD(12)
      AINPUT  = HOLD(13)
      IINPUT  = HOLD(14)
      ISTK    = HOLD(15)
      ASTK    = HOLD(16)
      INTSPA  = HOLD(17)
      RLSPA   = HOLD(18)
      PTRIRN  = HOLD(19)
      PTRA    = HOLD(20)
      NTOTPV  = HOLD(21)
      NPOTPV  = HOLD(22)
      NUMORG  = HOLD(23)
      NFRONT  = HOLD(24)
      NASS    = HOLD(25)
C     NCOL    = HOLD(26)
      IF (HOLD(1).EQ.1) NELL    = HOLD(27)
      IF (HOLD(1).EQ.2) NPIV    = HOLD(27)
      IASS    = HOLD(28)
      TINSPA  = HOLD(29)
      TRLSPA  = HOLD(30)
      TOTSTA(1) = HOLD(31)
      TOTSTA(2) = HOLD(32)
      NSTACK(1) = HOLD(33)
      NSTACK(2) = HOLD(34)
      INFO(32)  = HOLD(37)
      INFO(33)  = HOLD(38)
      INFO(34)  = HOLD(39)
      NBSTATIC  = HOLD(40)
      IF (ICNTL(7).GT.1) ISNPIV = HOLD(35)
      IF (ICNTL(7).EQ.4) PHASE = HOLD(36)
      IF (HOLD(1).EQ.2) NSC1    = NFRONT-NPIV
      FLOPSA = RINFO(3)
      FLOPSB = RINFO(4)
      FLOPSX = RINFO(5)
      IF (HOLD(1).EQ.1) THEN
C Real arrays expanded
      HOLD(1) = 0
      GO TO 333
      ELSE
      IF (HOLD(1).EQ.3) THEN
C We ran out of space when allocating values for zero pivots at the end
C       of the factorization.  This is most likely to happen when we
C       are running with ICNTL(16) equal to 1 (dropping small entries
C       from front).
        HOLD(1) = 0
        GO TO 555
      ELSE
C Integer arrays expanded
        HOLD(1) = 0
        GO TO 444
      ENDIF
      ENDIF
    ENDIF

C NBSTATIC is the number of modified diagonal entries
    NBSTATIC = 0
C NBLK is the number of block pivots used.
    NBLK = 0
C NTWO is the number of 2 x 2 pivots used.
    NTWO = 0
C NCMPBR, NCMPBI are the number of compresses on real and integer space
C   respectively.
    NCMPBR = 0
    NCMPBI = 0
C FLOPSA is the number of floating-point operations for assembly.
    FLOPSA = ZERO
C FLOPSB is the number of floating-point operations for elimination.
    FLOPSB = ZERO
C FLOPSX  counts extra flops required by use of GEMM
    FLOPSX = ZERO
C NEIG is number of negative eigenvalues detected.
    NEIG = 0
C MAXFRT is the maximum front size encountered so far.
    MAXFRT  = 0
C All relevant INFO and RINFO parameters initialized to zero
C   so they have a valid entry on any error return.
    INFO(1) = 0
    INFO(2) = 0
    INFO(14:29) = 0
    INFO(31:35) = 0
    RINFO(3:5) = ZERO
    RINFO(14:15) = ZERO

C Initialization of array indicating positions of variables in front
    DO 10 I = 1,N
      PPOS(I) = N + 1
   10 CONTINUE
C IWPOS is set to position for first index of first block pivot
    IWPOS = 6
C Set first five entries to dummies to avoid unassigned var in ME57E
    IW(1) = 0
    IW(2) = 0
    IW(3) = 0
    IW(4) = 0
    IW(5) = 0
C APOSBB is a pointer to the next position for storing factors in A
    APOSBB = 1
C Initialize NSTKAC and INTSPA and RLSPA
    NSTACK(1) = 0
    NSTACK(2) = 0
    NSTKAC(1) = NE
    NSTKAC(2) = NE
    TOTSTA(1) = NE
    TOTSTA(2) = NE
    INTSPA = NE+5+N
    RLSPA = NE
    TINSPA = NE+5+N
    TRLSPA = NE
C PTRIRN points to the next original row in IW.
    PTRIRN = LIW - NE + 1
C PTRA  points to the next original row in A.
    PTRA = LA - NE + 1
C ISTK points to the position in IW immediately before the stack.
    ISTK = PTRIRN - 1
C ASTK points to the position in A immediately before the stack.
    ASTK = PTRA - 1
C AINPUT is the position of the first entry of the original matrix
C   reals since the last compress.  It is reset to the current
C   part being processed when ME57P is called with REAL .TRUE.
    AINPUT = PTRA
C IINPUT is the position of the first entry of the original matrix
C   integers since the last compress.  It is reset to the current
C   part being processed when ME57P is called. with REAL = .FALSE.
    IINPUT = PTRIRN
C NTOTPV is the total number of pivots selected.
    NTOTPV = 0
C NPOTPV is the total number of potential pivots so far.
    NPOTPV = 0
C In case we run out of space before the first pivot is chosen, we
C   must initialize ISNPIV.
    IF (ICNTL(7).GT.1) ISNPIV = 0


C Must initialize PHASE here because of ICNTL(8) (OLD) return
C when KIND = 1.
    PHASE = 1
C Calculate diagonal of matrix and store in DIAG
    IF (ICNTL(7).EQ.4 .AND. KIND.NE.1) THEN
      DO 19 I = 1,N
      DIAG(I) = ZERO
   19   CONTINUE
      APOS1 = PTRA-1
      J1 = PTRIRN
      DO 20 I = 1,N
      J2 = J1 + LROW(I) - 1
      DO 25 JJ = J1,J2
        J = IW(JJ)
        APOS1 = APOS1 + 1
        IF (J.EQ.PERM(I)) DIAG(J) = DIAG(J) + REAL(A(APOS1))
   25   CONTINUE
      J1 = J2 + 1
   20   CONTINUE
      SCHNAB(1) = ONE
      SCHNAB(5) = ZERO
      DO 21 I = 1,N
      SCHNAB(1) = MAX(REAL(SCHNAB(1)),REAL(DIAG(I)))
      SCHNAB(5) = MIN(REAL(SCHNAB(5)),REAL(DIAG(I)))
   21   CONTINUE
C Set max entry on diag
      SCHNAB(4) = SCHNAB(1)
C S+E has **2/3 .. Nick wants **1/3
      SCHNAB(2) = FD15A('E')**(1.0/3.0)
      SCHNAB(3) = 0.1
C Initialize RINFO(15) to compute smallest pivot in modified matrix
      RINFO(15) = FD15A('H')
      DELTA   = ZERO
    ENDIF

C   *****************************************************************
C   * Each pass through this main loop performs all the operations  *
C   * associated with one node of the assembly tree.          *
C   *****************************************************************

    IASS = 1
C   DO 2160 IASS = 1,NSTEPS
 2160 CONTINUE

C Find the frontal variables, ordered with the fully summed variables
C   of the incoming rows first, the fully summed rows from previous
C   steps next, followed by the rest in any order.

C NUMORG is the number of variables in the tentative pivot from the
C   incoming original rows.
C Calculate NUMORG and put indices of these fully summed rows in IW.
      NUMORG = 0
      DO 30 I = NPOTPV + 1,N
C J is Ith variable in tentative pivotal sequence.
      J = PERM(I)
C Jump if we have finished with variables in current node.
      IF (ABS(NODE(J)).GT.IASS) GO TO 40
      IW(IWPOS+NUMORG) = J
      NUMORG = NUMORG + 1
      PPOS(J) = NUMORG
   30   CONTINUE

C NASS will be set to the total number of fully assembled variables in
C   the newly created element. First set it to NUMORG.
   40   NASS = NUMORG
C Add indices of fully summed variables of stacked sons to IW.
      NELL = NSTK(IASS)
      IELL = ISTK + 1
      DO 70 ELT = 1,NELL
      DO 50 JJ = IELL + 1,IELL + IW(IELL)
        J = IW(JJ)
        IF (NODE(J).GT.IASS) GO TO 50
C Jump if variable already included.
        IF (PPOS(J).LE.N) GO TO 50
        IW(IWPOS+NASS) = J
        NASS = NASS + 1
        PPOS(J) = NASS
   50   CONTINUE
      IELL = IELL + IW(IELL) + 1
   70   CONTINUE
C IWNFS points to the first free location for a variable that is not
C   fully summed.
      IWNFS = IWPOS + NASS

C Incorporate original rows.
C J1 is the position of the start of the first original row associated
C   with this node of the assembly tree.
      J1 = PTRIRN
      DO 90 IORG = 1,NUMORG
      J2 = J1 + LROW(NPOTPV+IORG) - 1
C Run through index list of original row.
      DO 80 JJ = J1,J2
        J = IW(JJ)
C Jump if variable already included.
        IF (PPOS(J).LE.N) GO TO 80
        IW(IWNFS) = J
        IWNFS = IWNFS + 1
        PPOS(J) = IWNFS - IWPOS
   80   CONTINUE
      J1 = J2 + 1
   90   CONTINUE

C Now incorporate stacked elements.
C J1 is set to beginning
C J2 is set to end
      IELL = ISTK + 1
      DO 170 ELT = 1,NELL
      J1 = IELL+1
      J2 = IELL+IW(IELL)
      DO 150 JJ = J1,J2
        J = IW(JJ)
C Jump if already assembled
        IF (PPOS(J).LE.N) GO TO 150
        IW(IWNFS) = J
        IWNFS = IWNFS + 1
        PPOS(J) = IWNFS - IWPOS
  150   CONTINUE
      IELL = J2 + 1
  170   CONTINUE

C NFRONT is the number of variables in the front.
      NFRONT = IWNFS - IWPOS

C MAXFRT is the largest front size so far encountered.
      MAXFRT = MAX(MAXFRT,NFRONT)

C Set APOS to the position of first entry in frontal matrix.
      IF (INFO(1).NE.-3) THEN
C Buffer space allocated so that triangular part of pivot can be stored
C   without danger of overwrite.
      APOS = APOSBB + (NASS*(NASS+1))/2
      ELSE
      APOS = 1
      END IF
C
C Assemble reals into frontal matrix.
C
C Accumulate real space needed
      RLSPA  = MAX(RLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+NSTKAC(1))
      TRLSPA = MAX(TRLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+TOTSTA(1))

C If necessary, compress A.

  333   IF (APOS+NFRONT*NFRONT-1.GT.ASTK) THEN


      CALL ME57P(A,IW,ASTK,AINPUT,PTRA,.TRUE.)

      NCMPBR = NCMPBR + 1
      IF (APOS+NFRONT*NFRONT-1.GT.ASTK) THEN
        IF (ICNTL(8).NE.0) THEN
C Zero part of A to avoid failure in HSL_MA57
          DO 334 I = APOSBB,ASTK
          A(I) = ZERO
  334       CONTINUE
          HOLD(1) = 1
          HOLD(2) = NBLK
          HOLD(3) = NTWO
          HOLD(4) = INFO(23)
          HOLD(5) = NCMPBI
          HOLD(6) = NEIG
          HOLD(7) = MAXFRT
          HOLD(8) = IWPOS
          HOLD(9) = APOS
          HOLD(10) = APOSBB
          HOLD(11) = NSTKAC(1)
          HOLD(12) = NSTKAC(2)
          HOLD(13) = AINPUT
          HOLD(14) = IINPUT
          HOLD(15) = ISTK
          HOLD(16) = ASTK
          HOLD(17) = INTSPA
          HOLD(18) = RLSPA
          HOLD(19) = PTRIRN
          HOLD(20) = PTRA
          HOLD(21) = NTOTPV
          HOLD(22) = NPOTPV
          HOLD(23) = NUMORG
          HOLD(24) = NFRONT
          HOLD(25) = NASS
C         HOLD(26) = NCOL
          HOLD(27) = NELL
          HOLD(28) = IASS
          HOLD(29) = TINSPA
          HOLD(30) = TRLSPA
          HOLD(31) = TOTSTA(1)
          HOLD(32) = TOTSTA(2)
          HOLD(33) = NSTACK(1)
          HOLD(34) = NSTACK(2)
          IF (ICNTL(7).GT.1) HOLD(35) = ISNPIV
          IF (ICNTL(7).EQ.4) HOLD(36) = PHASE
          HOLD(37) = INFO(32)
          HOLD(38) = INFO(33)
          HOLD(39) = INFO(34)
          RINFO(3) =FLOPSA
          RINFO(4) =FLOPSB
          RINFO(5) =FLOPSX
          HOLD(40) = NBSTATIC
          INFO(35) = HOLD(40)
          INFO(1) = 10
          RETURN
        ELSE
C INFO(40) accumulates number of discards from factors.
          INFO(40) = INFO(40) + APOS - 1
          APOS = 1
          APOSBB = 1
          INFO(1) = -3
          IF (NFRONT*NFRONT.GT.ASTK) THEN
          INFO(17) = MAX(INFO(17),RLSPA)
          IF (ICNTL(7).EQ.4) INFO(17) = MAX(INFO(17),RLSPA + N)
          INFO(2) = LA
          RETURN
          ENDIF
        ENDIF
      ENDIF
      END IF

      ATRASH = APOS + NFRONT*NFRONT - 1
C Zero out appropriate part of A for incoming potential pivot rows.
      DO 210 JJ = APOS,ATRASH
      A(JJ) = ZERO
  210   CONTINUE

C Incorporate reals from original rows.
      J1 = PTRIRN
      DO 230 IORG = 1,NUMORG
C APOSI indicates the position in A just before the beginning of the row
C     being assembled.
      J = PERM(NPOTPV+IORG)
      APOSI = APOS + (PPOS(J)-1)*NFRONT - 1
      J2 = J1 + LROW(NPOTPV+IORG) - 1
      FLOPSA = FLOPSA + J2 - J1 + 1
      DO 220 JJ = J1,J2
        JAY = IW(JJ)
CCC
C Entries always in upper triangle because of ordering.
C Pivot permutations can only affect fully summed variables
        APOS2 = APOSI + PPOS(JAY)
        A(APOS2) = A(APOS2) + A(PTRA)
        PTRA = PTRA + 1
  220   CONTINUE
C       AINPUT = AINPUT + J2 - J1 + 1
      NSTKAC(1) = NSTKAC(1) - J2 + J1 - 1
      J1 = J2 + 1
  230   CONTINUE
C     IINPUT = IINPUT + J2 - PTRIRN
      NSTKAC(2) = NSTKAC(2) - J1 + PTRIRN
      PTRIRN = J1
C Update NPOTPV
      NPOTPV = NPOTPV + NUMORG

C Now assemble reals from stacked elements
C POSELT is a running pointer into that element.
      DO 380 ELT = 1,NELL
      POSELT = ASTK + 1
      LIELL = IW(ISTK+1)
      J1 = ISTK + 2
      J2 = ISTK+1 + LIELL
      FLOPSA = FLOPSA + (LIELL*(LIELL+1))/2
      DO 250 JJ = J1,J2
        J = IW(JJ)
        APOS2 = APOS + (PPOS(J)-1)*NFRONT
        APOS1 = POSELT
        DO 240 JJJ=JJ,J2
          JAY = IW(JJJ)
          IF (PPOS(JAY) .GE. PPOS(J)) THEN
C Upper triangle
          APOS3 = APOS2 + PPOS(JAY) - 1
          A(APOS3) = A(APOS3) + A(APOS1)
          ELSE
          APOS3 = APOS+(PPOS(JAY)-1)*NFRONT+PPOS(J)-1
          IF (KIND.EQ.1) THEN
            A(APOS3) = A(APOS3) + A(APOS1)
          ELSE
            A(APOS3) = A(APOS3) + CONJG(A(APOS1))
          ENDIF
          ENDIF
          APOS1 = APOS1 + 1
  240     CONTINUE
        POSELT = POSELT + LIELL - (JJ-J1)
  250   CONTINUE
C ISTK and ASTK updated to point to posn before next element on stack.
      NSTKAC(2) = NSTKAC(2) - (J2-ISTK)
      NSTACK(2) = NSTACK(2) - (J2-ISTK)
      TOTSTA(2) = TOTSTA(2) - (J2-ISTK)
      ISTK = J2
      ASTK = ASTK + (LIELL*(LIELL+1))/2
      NSTKAC(1) = NSTKAC(1) - (LIELL*(LIELL+1))/2
      NSTACK(1) = NSTACK(1) - (LIELL*(LIELL+1))/2
      TOTSTA(1) = TOTSTA(1) - (LIELL*(LIELL+1))/2
  380   CONTINUE

C       IF (LCASE) THEN
C       write(7,'(/A,I8/A,5I8)') '*** Frontal matrix before step',
C    *    IASS,'NFRONT,NASS,NUMORG,APOS,IWPOS',
C    *      NFRONT,NASS,NUMORG,APOS,IWPOS
C       write(7,'(/A/(10I8))') 'IW array',
C    *    (IW(IWPOS+I-1),I=1,NFRONT)
C       write(7,'(/A)') 'Frontal matrix'
C       DO 1122 J = 1, NFRONT
C       write(7,'(/(6E13.5))')
C    *      (A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT)
C1122   CONTINUE
C       ENDIF

C ******************
C Each time round this loop, we sweep through the remainder
C    of the assembled part of the front looking for pivots.
C ******************

C Set PIVBLK
      PIVBLK = MIN(NBLOC,NASS)
C Set pointer to first entry in block
      APOSBK = APOS
C NPIV is the number of pivots so far selected at the current node.
      NPIV = 0
C Set local value for U
      ULOC = UU

C Each pass through loop processes one block.
      DO 918 BLK = 1,NASS
C Set last block flag
      IF (NPIV+PIVBLK .GE. NASS) THEN
      LASTBK = .TRUE.
      SIZBLK = NASS - NPIV
      ELSE
      LASTBK = .FALSE.
      SIZBLK = PIVBLK
      ENDIF
C Record number of rows processed to date
      LASPIV = NPIV
C MPIV is the number of pivots so far selected in the current block.
      MPIV = 0
C KR is relative position of current pivot candidate
CCCCCCCC
      KR = 0
CCC Set to following to force 2 by 2 pivots in Nocedal examples
C     KR = NUMORG
C KCT is set to one more than the number of unsearched pivot candidates
      KCT = SIZBLK + 1

C**************************
C Loop for pivot searches *
C**************************
  920   CONTINUE
C Increment pointer for circular sweep.
      KR = KR + 1
      KCT = KCT - 1
      IF (KCT.EQ.0) GO TO 930
      IF (KR.GT.SIZBLK) KR = MPIV + 1
C Check to see if no pivot was chosen is complete sweep of pivot block.
C We either take the diagonal entry or the 2 by 2 pivot with the
C   largest fully summed off-diagonal at each stage.
C Note that IPIV is the position within the complete current front.
      IPIV = LASPIV + KR
C APOSI is set to the beginning of row IPIV in working front.
        APOSI = APOS + (IPIV-1)*NFRONT
C Set position and value of potential pivot.
        POSPV1 = APOSI + IPIV - 1
        PIVOT = A(POSPV1)

C Schnabel/Eskow logic
C First check that matrix is Hermitian
   29     IF (KIND.NE.1) THEN
         IF (ICNTL(7).EQ.4) THEN
          IF (PHASE.EQ.2) THEN
C In phase 2 of matrix modification (ICNTL(7) = 4)
C Compute quantity to add to diagonal
C Calculate norm of pivot row
          IF (INFO(27).EQ.0) INFO(27) = NTOTPV + 1
          NORMJ = ZERO
          DO 28 I = POSPV1+1,POSPV1+NFRONT-NPIV-1
            NORMJ = NORMJ + ABS(A(I))
   28       CONTINUE
          DELTA = MAX(ZERO,
   *            - REAL(A(POSPV1))
   *            + MAX(NORMJ,REAL(SCHNAB(2))*REAL(SCHNAB(1))))
          A(POSPV1) = A(POSPV1) + DELTA
          IF (A(POSPV1).EQ.ZERO) GO TO 970
          RINFO(15) = MIN(RINFO(15),REAL(A(POSPV1)))
          DIAG(PERM(NTOTPV+1)) = DELTA
          PIVSIZ = 1
          GO TO 811
          ENDIF
         ENDIF
        ENDIF

        IF (ICNTL(7).GT.1) THEN
C Action if no pivoting requested
          IF (ABS(PIVOT).LE.CNTL(2)) THEN
C Jump if we are not using Schnabel-Eskow modification
          IF (ICNTL(7).LT.4 .OR. KIND.EQ.1) GO TO 970
C We are now in phase 2 of matrix modification (ICNTL(7) = 4)
          PHASE = 2
          GO TO 29
          ENDIF
          IF (NTOTPV.EQ.0) THEN
          IF (REAL(PIVOT).GT.ZERO) ISNPIV = 1
          IF (REAL(PIVOT).LT.ZERO) ISNPIV = -1
          ELSE
          IF (ICNTL(7).EQ.2 .AND. ISNPIV*REAL(PIVOT).LT.ZERO)
   *         GO TO 980
          IF (ICNTL(7).EQ.3 .AND. ISNPIV*REAL(PIVOT).LT.ZERO)
   *         THEN
            INFO(26) = INFO(26) + 1
            ISNPIV = -ISNPIV
          ENDIF
          ENDIF
          IF (ICNTL(7).EQ.4 .AND. KIND.NE.1) THEN
          IF (REAL(PIVOT).GE.REAL(SCHNAB(1))*REAL(SCHNAB(2))
   *          .AND. REAL(SCHNAB(5)) .GE.
   *          -REAL(SCHNAB(3))*REAL(SCHNAB(4))) THEN
C Update and check values of future diagonals
            SCHNAB(5) = ZERO
            SCHNAB(4) = ZERO
            DO 22 I = POSPV1+1,POSPV1+NFRONT-NPIV-1
            J = IW(IWPOS+NPIV+I-POSPV1)
            DIAG(J) = DIAG(J) - REAL(A(I)*A(I)/PIVOT)
            SCHNAB(5) = MIN(REAL(DIAG(J)),REAL(SCHNAB(5)))
            SCHNAB(4) = MAX(REAL(DIAG(J)),REAL(SCHNAB(4)))
            IF (REAL(DIAG(J)).LT.
   *            -REAL(SCHNAB(3))*REAL(SCHNAB(1))) THEN
              PHASE = 2
              GO TO 29
            ENDIF
   22         CONTINUE
            DIAG(PERM(NTOTPV+1)) = ZERO
            RINFO(15) = MIN(RINFO(15),REAL(PIVOT))
          ELSE
            PHASE = 2
            GO TO 29
          ENDIF
          ENDIF
          PIVSIZ = 1
          GO TO 811
        ENDIF

C Numerical pivoting
C Find largest off-diagonal entry in the part of the row in which we
C   seek a pivot.
        AMAX = ZERO
        JMAX = 0
C Split loops in two because only upper triangle is held
C Scan lower triangle by scanning up column of upper triangle.
        DO 110 K = 1, IPIV - NPIV - 1
          IF (ABS(A(POSPV1-K*NFRONT)).GT.AMAX) THEN
          AMAX = ABS(A(POSPV1-K*NFRONT))
          CMAX = A(POSPV1-K*NFRONT)
          JMAX = IPIV - K
          ENDIF
  110     CONTINUE
C Scan upper triangle by scanning along row from first off-diagonal
        DO 111 K =  1, MIN(NASS,LASPIV+PIVBLK) - IPIV
          IF (ABS(A(POSPV1+K)).GT.AMAX) THEN
          AMAX = ABS(A(POSPV1+K))
          CMAX = A(POSPV1+K)
          JMAX = IPIV + K
          ENDIF
  111     CONTINUE
C Do same for the other part.
        RMAX = ZERO

C Restrict partial pivoting check to the fully summed block
C        RMAX = STCTOL
C        DO 112 K = MIN(NASS,LASPIV+PIVBLK)-IPIV+1,NASS
C          RMAX = MAX(RMAX,ABS(A(POSPV1+K)))
C  112     CONTINUE

        DO 112 K = MIN(NASS,LASPIV+PIVBLK)-IPIV+1,NFRONT-IPIV
         RMAX = MAX(RMAX,ABS(A(POSPV1+K)))
 112      CONTINUE

C Action taken if matrix is singular.
        IF (MAX(AMAX,RMAX,ABS(PIVOT)).LE.TOL) THEN
C Skip if all of row is zero.
          GO TO 920
        END IF

C Jump if no nonzero entry in row of pivot block
        IF (MAX(AMAX,ABS(PIVOT)).LE.TOL) GO TO 920
        PIVSIZ = 0
        IF (ABS(PIVOT).GT.ULOC*MAX(RMAX,AMAX)) THEN
          PIVSIZ = 1
          A(POSPV1) = PIVOT
C 1 x 1 pivot is chosen
          GO TO 810
        END IF

C If there is only one remaining fully summed row and column exit.
        IF (NPIV+1.EQ.NASS) THEN
          A(POSPV1) = PIVOT
          GO TO 920
        END IF

C Jump if 2 x 2 candidate is diagonal
        IF (AMAX.LE.TOL) GO TO 920

C Check block pivot of order 2 for stability.
C Find largest entry in row IPIV outwith the pivot.
        IF (RMAX.LT.AMAX) THEN
          RMAX = ZERO
C Split loops in two because only upper triangle is held
C Scan lower triangle by scanning up column of upper triangle.
          DO 113 K = 1, IPIV - NPIV - 1
          IF (IPIV-K.EQ.JMAX) GO TO 113
          RMAX=MAX(RMAX,ABS(A(POSPV1-K*NFRONT)))
  113       CONTINUE
C Scan upper triangle by scanning along row from first off-diagonal
          DO 114 K =  1, NFRONT - IPIV
          IF (IPIV+K.EQ.JMAX) GO TO 114
          RMAX = MAX(RMAX,ABS(A(POSPV1+K)))
  114       CONTINUE
        ENDIF

C APOSJ is set to the beginning of row JMAX in working front.
        APOSJ = APOS + (JMAX-1)*NFRONT
C POSPV2 is the position in A of the second diagonal of a 2 x 2 pivot.
C OFFDAG is the position in A of the off-diagonal of a 2 x 2 pivot.
        POSPV2 = APOSJ + JMAX - 1
        IF (IPIV.GT.JMAX) THEN
          OFFDAG = APOSJ + IPIV - 1
        ELSE
          OFFDAG = APOSI + JMAX - 1
        END IF

C     Find largest entry in row JMAX outwith the pivot.
        TMAX = ZERO
C Split loops in two because only upper triangle is held
C Scan lower triangle by scanning up column of upper triangle.
        DO 115 K = 1, JMAX - NPIV - 1
          IF (JMAX-K.EQ.IPIV) GO TO 115
          TMAX=MAX(TMAX,ABS(A(POSPV2-K*NFRONT)))
  115     CONTINUE
C Scan upper triangle by scanning along row from first off-diagonal
        DO 116 K =  1, NFRONT - JMAX
          IF (JMAX+K.EQ.IPIV) GO TO 116
          TMAX = MAX(TMAX,ABS(A(POSPV2+K)))
  116     CONTINUE


C DETPIV is the value of the determinant of the 2x2 pivot.
        IF (KIND.EQ.1) THEN
          DETPIV = A(POSPV1)*A(POSPV2) - CMAX*CMAX
        ELSE
          DETPIV = A(POSPV1)*A(POSPV2) - CMAX*CONJG(CMAX)
C Make DETPIV real in Hermitian case
          DETPIV = REAL(DETPIV)
        ENDIF

        MAXPIV = MAX(ABS(A(POSPV1)),ABS(A(POSPV2)))
        IF (MAXPIV.EQ.ZERO) MAXPIV = ONE
        IF (ABS(DETPIV)/MAXPIV.LE.TOL) GO TO 920
        PIVSIZ = 2
C Check pivot for stability
C Jump if pivot fails test
C This is componentwise test
        IF ((ABS(A(POSPV2))*RMAX+AMAX*TMAX)*ULOC.GT.
   +      ABS(DETPIV)) GO TO 920
        IF ((ABS(A(POSPV1))*TMAX+AMAX*RMAX)*ULOC.GT.
   +      ABS(DETPIV)) GO TO 920

C 2 x 2 pivot is chosen

C
C       Pivot has been chosen. It has order PIVSIZ.
  810     LPIV = IPIV
        IF (PIVSIZ.EQ.2) LPIV = MIN(IPIV,JMAX)
C Change made at Stephane's suggestion
CCC       KR = MAX(KR,NPIV+PIVSIZ)
        KR = MAX(KR,MPIV+PIVSIZ)
        KCT = SIZBLK - MPIV - PIVSIZ + 1

C The following loop moves the pivot block to the top left
C       hand corner of the uneliminated frontal matrix.
        DO 860 KROW = NPIV,NPIV + PIVSIZ - 1
C We jump if swop is not necessary.
          IF (LPIV.EQ.KROW+1) GO TO 850

C Swop first part of rows (going down columns)
C JA1 is used as running index for row LPIV
          JA1 = APOS + (LPIV-1)
C J1 is used as running index for row KROW+1
          J1 = APOS + KROW
          DO 820 JJ = 1,KROW
          SWOP = A(JA1)
          A(JA1) = A(J1)
          A(J1) = SWOP
          JA1 = JA1 + NFRONT
          J1 = J1 + NFRONT
  820       CONTINUE
C Swop off-diagonal (in position KROW+1, IPIV)
C Only needed in Hermitian case
          IF (KIND.NE.1)  A(JA1) = CONJG(A(JA1))
C Swop middle part of rows (KROW+1 by rows, LPIV by columns)
          JA1 = JA1 + NFRONT
          J1 = J1 + 1
          DO 830 JJ = 1,LPIV - KROW - 2
          SWOP = A(JA1)
          IF (KIND.EQ.1) THEN
            A(JA1) = A(J1)
            A(J1) = SWOP
          ELSE
            A(JA1) = CONJG(A(J1))
            A(J1) = CONJG(SWOP)
          ENDIF
          JA1 = JA1 + NFRONT
          J1 = J1 + 1
  830       CONTINUE
C Swop diagonals
          SWOP = A(APOS+KROW* (NFRONT+1))
          A(APOS+KROW* (NFRONT+1)) = A(JA1)
          A(JA1) = SWOP
C Swop last part of rows
          DO 840 JJ = 1,NFRONT - LPIV
          JA1 = JA1 + 1
          J1 = J1 + 1
          SWOP = A(JA1)
          A(JA1) = A(J1)
          A(J1) = SWOP
  840       CONTINUE
C Swop integer indexing information
          IPOS = IWPOS + KROW
          IEXCH = IWPOS + LPIV - 1
          ISWOP = IW(IPOS)
          IW(IPOS) = IW(IEXCH)
          IW(IEXCH) = ISWOP
C Set LPIV for the swop of the second row of block pivot.
  850       LPIV = MAX(IPIV,JMAX)
  860     CONTINUE

C
C Set POSPV1 and POSPV2 to new position of pivots.
  811     POSPV1 = APOS + NPIV* (NFRONT+1)
        POSPV2 = POSPV1 + NFRONT + 1

        IF (PIVSIZ.EQ.1) THEN
C Perform the elimination using entry A(POSPV1) as pivot.
C We store U and D inverse.
C Later we store D inverse U which is passed to the solution entry.
          FLOPSB = FLOPSB + ONE
          IF (KIND.EQ.1) A(POSPV1) = ONE/A(POSPV1)
C Keep pivots real in Hermitian case
          IF (KIND.NE.1) A(POSPV1) = REAL(ONE/A(POSPV1))
          IF (REAL(A(POSPV1)).LT.ZERO) NEIG = NEIG + 1
          J1 = POSPV1 + 1
          J2 = POSPV1 + NASS - (NPIV+1)
          IBEG = POSPV1 + NFRONT + 1
          IEND = APOS + (NPIV+1)*NFRONT + NFRONT - 1
          DO 880 JJ = J1,J2
C AMULT1 is used to hold the multiplier
          IF (KIND.EQ.1) THEN
            AMULT1 = -A(JJ)*A(POSPV1)
          ELSE
            AMULT1 = -CONJG(A(JJ))*A(POSPV1)
          ENDIF
C Hold original entry for GEMM multiply
          IF (.NOT.LASTBK) THEN
            IF (KIND.EQ.1) THEN
            A(POSPV1+(JJ-J1+1)*NFRONT) = A(JJ)
            ELSE
            A(POSPV1+(JJ-J1+1)*NFRONT) = CONJG(A(JJ))
            ENDIF
          ENDIF
          JCOL = JJ
          FLOPSB = FLOPSB + (IEND-IBEG+1)*2 + 1
          IF (MPIV+JJ-J1+2.GT.PIVBLK) GO TO 871
C          The following special comment forces vectorization on
C             Crays.
CDIR$        IVDEP
          DO 870 IROW = IBEG,IEND
            A(IROW) = A(IROW) + AMULT1*A(JCOL)
            JCOL = JCOL + 1
  870       CONTINUE

CCC
C In Hermitian case make diagonal in position IBEG real
          IF (KIND.NE.1) A(IBEG) = REAL(A(IBEG))

  871       IF (KIND.EQ.1) THEN
            A(JJ) = AMULT1
          ELSE
            A(JJ) = CONJG(AMULT1)
          ENDIF
          IBEG = IBEG + NFRONT + 1
          IEND = IEND + NFRONT
  880       CONTINUE

          NPIV = NPIV + 1
          MPIV = MPIV + 1
          NTOTPV = NTOTPV + 1

C     IF (LCASE) THEN
C       write(7,'(/A,I8/A,7I8)') '*** Frontal matrix at step',IASS,
C    *   'NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ',
C    *    NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ
C       write(7,'(/A/(10I8))') 'IW array',
C    *    (IW(IWPOS+I-1),I=1,NFRONT)
C       write(7,'(/A)') 'Frontal matrix'
C       DO 1123 J = 1, NFRONT
C       write(7,'(/(6E13.5))')
C    *    (A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT)
C1123   CONTINUE
C     ENDIF

          IF (MPIV.EQ.SIZBLK) GO TO 930

        ELSE
C Perform elimination using block pivot of order two.
C Replace block pivot by its inverse.
          OFFDAG = POSPV1 + 1
          FLOPSB = FLOPSB + 6.0

          SWOP = A(POSPV2)
          IF (REAL(DETPIV).LT.ZERO) THEN
          NEIG = NEIG + 1
          ELSE
          IF (REAL(SWOP).LT.ZERO) NEIG = NEIG + 2
          END IF

          IF (KIND.EQ.1) THEN
          A(POSPV2) = A(POSPV1)/DETPIV
          A(POSPV1) = SWOP/DETPIV
          A(OFFDAG) = -A(OFFDAG)/DETPIV
          ELSE
          A(POSPV2) = REAL(A(POSPV1)/DETPIV)
          A(POSPV1) = REAL(SWOP/DETPIV)
          A(OFFDAG) = -A(OFFDAG)/DETPIV
          ENDIF

          J1 = POSPV1 + 2
          J2 = POSPV1 + NASS - (NPIV+1)
CCC NB .. not the following
C         J2 = POSPV1 + NFRONT - (NPIV+1)
          IBEG = POSPV2 + NFRONT + 1
          IEND = APOS + (NPIV+2)*NFRONT + NFRONT - 1
          DO 900 JJ = J1,J2
          K1 = JJ
          K2 = JJ + NFRONT
          IF (KIND.EQ.1) THEN
            AMULT1 = - (A(POSPV1)*A(K1)+A(POSPV1+1)*A(K2))
            AMULT2 = - (A(POSPV1+1)*A(K1)+A(POSPV2)*A(K2))
          ELSE
            AMULT1 = - (A(POSPV1)*CONJG(A(K1))
   *                + CONJG(A(POSPV1+1))*CONJG(A(K2)))
            AMULT2 = - (A(POSPV1+1)*CONJG(A(K1))
   *               + A(POSPV2)*CONJG(A(K2)))
          ENDIF
C Save block in lower triangular part of front
          IF (.NOT.LASTBK) THEN
            IF (KIND.EQ.1) THEN
            A(POSPV1 + (JJ-J1+2)*NFRONT) = A(K1)
            A(POSPV1 + (JJ-J1+2)*NFRONT + 1) = A(K2)
            ELSE
            A(POSPV1 + (JJ-J1+2)*NFRONT) = CONJG(A(K1))
            A(POSPV1 + (JJ-J1+2)*NFRONT + 1) = CONJG(A(K2))
            ENDIF
          ENDIF
          FLOPSB = FLOPSB + (IEND-IBEG+1)*4 + 6
          IF (MPIV+JJ-J1+3.GT.PIVBLK) GO TO 891

C          The following special comment forces vectorization on
C            Crays.
CDIR$        IVDEP
          DO 890 IROW = IBEG,IEND
            A(IROW) = A(IROW) + AMULT1*A(K1) + AMULT2*A(K2)
            K1 = K1 + 1
            K2 = K2 + 1
  890       CONTINUE
  891       IF (KIND.EQ.1) THEN
            A(JJ) = AMULT1
            A(JJ+NFRONT) = AMULT2
          ELSE
            A(JJ) = CONJG(AMULT1)
            A(JJ+NFRONT) = CONJG(AMULT2)
          ENDIF
          IBEG = IBEG + NFRONT + 1
          IEND = IEND + NFRONT
  900       CONTINUE
C Flag column indices of 2 x 2 pivot.
          IPOS = IWPOS + NPIV
          IW(IPOS) = -IW(IPOS)
          IW(IPOS+1) = -IW(IPOS+1)
          NPIV = NPIV + 2
          MPIV = MPIV + 2
          NTOTPV = NTOTPV + 2
          NTWO = NTWO + 1

C      IF (LCASE) THEN
C       write(7,'(/A,I8/A/9I8)') '*** Frontal matrix at step',IASS,
C    *   'NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ,MPIV,SIZBLK',
C    *    NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ,MPIV,SIZBLK
C       write(7,'(/A/(10I8))') 'IW array',
C    *    (IW(IWPOS+I-1),I=1,NFRONT)
C       write(7,'(/A)') 'Frontal matrix'
C       DO 1124 J = 1, NFRONT
C       write(7,'(/(6E13.5))')
C    *    (A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT)
C1124   CONTINUE
C      ENDIF

          IF (MPIV.EQ.SIZBLK) GO TO 930
        END IF

      GO TO 920
C 920   CONTINUE


C Finished pivoting on block BLK ... no update rest of pivot block
C     using GEMM
 930    IF (LASTBK) THEN
      IF (NPIV.EQ.NASS) GO TO 935
      IF (.NOT. LSTAT)  GO TO 935
      ULOC = ULOC/10.0E0
      IF (ULOC.LT.UTARG) THEN
        ULOC = ULOC * 10.0E0
        GO TO 9919
      ENDIF
      KCT = SIZBLK + 1 - MPIV
      GO TO 920
      ENDIF
C Check if any pivots chosen from this block. If not, increase PIVBLK
C     and try again.
      IF (MPIV.EQ.0) THEN
      PIVBLK = 2*PIVBLK
      GO TO 918
      ENDIF

      KBLK = (NASS-(LASPIV+PIVBLK))/PIVBLK
CCC   IF (LCASE) KBLK = 0
      L = NASS - (LASPIV+PIVBLK)
      APOS4 = APOS+(LASPIV+PIVBLK)*(NFRONT+1)
      DO 931 KB = 1,KBLK
      FLOPSX = FLOPSX + PIVBLK*(PIVBLK-1)*MPIV
      IF (KIND.EQ.1) THEN
        CALL CGEMM('N','N',L-(KB-1)*PIVBLK,PIVBLK,MPIV,CONE,
   +           A(APOSBK+PIVBLK*KB),NFRONT,
   +           A(APOSBK+PIVBLK*KB*NFRONT),NFRONT,CONE,
   +           A(APOS4+PIVBLK*(KB-1)*(NFRONT+1)),NFRONT)
      ELSE
        CALL CGEMM('N','N',L-(KB-1)*PIVBLK,PIVBLK,MPIV,CONE,
   +           A(APOSBK+PIVBLK*KB),NFRONT,
   +           A(APOSBK+PIVBLK*KB*NFRONT),NFRONT,CONE,
   +           A(APOS4+PIVBLK*(KB-1)*(NFRONT+1)),NFRONT)
      ENDIF
C And now process the part of the pivot row outside the fs block
      IF (NFRONT.GT.NASS) THEN
        IF (KIND.EQ.1)
   +    CALL CGEMM('N','T',NFRONT-NASS,PIVBLK,MPIV,CONE,
   +           A(APOSBK+NASS-LASPIV),NFRONT,
   +           A(APOSBK+PIVBLK*KB),NFRONT,CONE,
   +           A(APOSBK+KB*NFRONT*PIVBLK+NASS-LASPIV),NFRONT)
        IF (KIND.NE.1)
   +    CALL CGEMM('N','C',NFRONT-NASS,PIVBLK,MPIV,CONE,
   +           A(APOSBK+NASS-LASPIV),NFRONT,
   +           A(APOSBK+PIVBLK*KB),NFRONT,CONE,
   +           A(APOSBK+KB*NFRONT*PIVBLK+NASS-LASPIV),NFRONT)
      ENDIF
  931   CONTINUE

     SIZC = NASS - (KBLK+1)*PIVBLK - LASPIV
     SIZF = NFRONT - (KBLK+1)*PIVBLK - LASPIV
     APOSA = APOSBK + (KBLK+1)*PIVBLK
     DO 934 K = 1,MPIV
C APOSB is in lower triangle and so conjugated already
       APOSB = APOSBK + NFRONT*PIVBLK*(KBLK+1) + K - 1
       APOSM = APOSBK + PIVBLK*(KBLK+1) + (K-1)*NFRONT
       APOSC = APOSBK + PIVBLK*(KBLK+1)*(NFRONT+1)
       DO 933 JJ = 1,SIZC
        IF (KIND.EQ.1) THEN
          DO 932 J = JJ,SIZC
          A(APOSC+J-1) = A(APOSC+J-1) + A(APOSA+J-1)*A(APOSB)
  932       CONTINUE
        ELSE
          DO 937 J = JJ,SIZC
          A(APOSC+J-1) = A(APOSC+J-1)
   *               + A(APOSA+J-1)*A(APOSB)
  937       CONTINUE
        ENDIF
C And now process the part of the pivot row outside the fs block
        IF (KIND.EQ.1) THEN
          DO 936 J = SIZC+1,SIZF
          A(APOSC+J-1) = A(APOSC+J-1) + A(APOSA+J-1)*A(APOSM)
  936       CONTINUE
        ELSE
          DO 938 J = SIZC+1,SIZF
          A(APOSC+J-1) = A(APOSC+J-1)
   *               + A(APOSA+J-1)*CONJG(A(APOSM))
  938       CONTINUE
        ENDIF
        APOSC = APOSC + NFRONT
        APOSB = APOSB + NFRONT
        APOSM = APOSM + 1
  933   CONTINUE
      APOSA = APOSA + NFRONT
  934   CONTINUE

      APOSBK = APOSBK + MPIV*(NFRONT+1)
      LASPIV = NPIV

C     IF (LCASE) THEN
C       write(7,'(/A,I8/A,7I8)') '*** Frontal matrix at step',IASS,
C    *   'NFRONT,NASS,APOS,IWPOS,NPIV',
C    *    NFRONT,NASS,APOS,IWPOS,NPIV
C       write(7,'(/A,2I8)') 'After blocking .. APOSBK,LASPIV',
C    *    APOSBK,LASPIV
C       write(7,'(/A/(10I8))') 'IW array',
C    *    (IW(IWPOS+I-1),I=1,NFRONT)
C       write(7,'(/A)') 'Frontal matrix'
C       DO 1125 J = 1, NFRONT
C       write(7,'(/(6E13.5))')
C    *    (A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT)
C1125   CONTINUE
C     ENDIF

  918   CONTINUE
C End of main elimination loop.

C SCHUR if set to .TRUE. then the Schur complement will be generated
C    using level 3 BLAS at this step so an extra copy of pivot row is
C    made.


C*************************
C   Do static pivoting *
C*************************
 9919    IPIV = LASPIV+MPIV
 9920    IPIV = IPIV + 1
CADD Probably not needed .. use only IPIV
C      IF (KR.GT.SIZBLK) KR = MPIV + 1
C   Note that IPIV is the position within the complete current front.
C      IPIV = LASPIV + KR
C   APOSI is set to the beginning of row IPIV in working front.
       APOSI = APOS + (IPIV-1)*NFRONT
C   Set position and value of potential pivot.
       POSPV1 = APOSI + IPIV - 1
       PIVOT = A(POSPV1)
CADD
C Although these are not needed just now they are kept for when
C we use 2 x 2 static pivots.
CCC      PIVSIZ = 1
CCC      LPIV = IPIV

C Logic has changed so no need to calculate AMAX
C This is code from earlier experiments

CCC      AMAX = ZERO
C Split loops in two because only upper triangle is held
C Scan lower triangle by scanning up column of upper triangle.
CCC      DO 9876 K = 1, IPIV - NPIV - 1
CCC      AMAX = MAX(AMAX,ABS(A(POSPV1-K*NFRONT)))
CCC76    CONTINUE
C Scan upper triangle by scanning along row from first off-diagonal
CCC      DO 9878 K =  1, NFRONT - IPIV
CCC      AMAX = MAX(AMAX,ABS(A(POSPV1+K)))
CCC78    CONTINUE

C   Check size of 1 x 1 pivot and adjust if necessary
C      IF (ABS(A(POSPV1)).LT.CNTL(4)) THEN
C          PIVOT = CNTL(4)
C      IF (ABS(A(POSPV1)).LT.MAX(ULOC*AMAX,STCTOL)) THEN
C          PIVOT = MAX(ULOC*AMAX,STCTOL)
       IF (ABS(A(POSPV1)).LT.STCTOL) THEN
         PIVOT = STCTOL
          IF (REAL(A(POSPV1)) .LT. ZERO) THEN
           A(POSPV1) = -PIVOT
           PIVOT   = -PIVOT
          ELSE
           A(POSPV1) = PIVOT
          ENDIF
          NBSTATIC = NBSTATIC + 1
       ENDIF
C Perform the elimination using entry A(POSPV1) as pivot.
C We store U and D inverse.
C Later we store D inverse U which is passed to the solution entry.
       FLOPSB = FLOPSB + ONE
       A(POSPV1) = ONE/A(POSPV1)
       IF (REAL(A(POSPV1)).LT.ZERO) NEIG = NEIG + 1

       J1 = POSPV1 + 1
       J2 = POSPV1 + NASS - (NPIV+1)
       IBEG = POSPV1 + NFRONT + 1
       IEND = APOSI + 2*NFRONT - 1
       DO 9880 JJ = J1,J2
C AMULT1 is used to hold the multiplier
          IF (KIND.EQ.1) THEN
          AMULT1 = -A(JJ)*A(POSPV1)
          ELSE
          AMULT1 = -CONJG(A(JJ))*A(POSPV1)
          ENDIF
C   Hold original entry for GEMM multiply
          JCOL = JJ
          FLOPSB = FLOPSB + (IEND-IBEG+1)*2 + 1
CADD Not necessary because control is on setting of J2
C         IF (MPIV+JJ-J1+2.GT.PIVBLK) GO TO 9871
C   The following special comment forces vectorization on
C   Crays.
CDIR$        IVDEP
          DO 9870 IROW = IBEG,IEND
           A(IROW) = A(IROW) + AMULT1*A(JCOL)
           JCOL = JCOL + 1
 9870       CONTINUE
          IF (KIND.EQ.1) THEN
          A(JJ) = AMULT1
          ELSE
          A(JJ) = CONJG(AMULT1)
          ENDIF
          IBEG = IBEG + NFRONT + 1
          IEND = IEND + NFRONT
 9880    CONTINUE
       NPIV = NPIV + 1
       MPIV = MPIV + 1
       NTOTPV = NTOTPV + 1

C    IF (LCASE) THEN
C   write(7,'(A)') 'Static pivoting step'
C   write(7,'(/A,I8/A/9I8)') '*** Frontal matrix at step',IASS,
C    *   'NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ,SIZBLK,MPIV',
C    *    NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ,SIZBLK,MPIV
C   write(7,'(/A/(10I8))') 'IW array',
C    *    (IW(IWPOS+I-1),I=1,NFRONT)
C       write(7,'(/A)') 'Frontal matrix'
C       DO 1126 J = 1, NFRONT
C       write(7,'(/(6E13.5))')
C    *    (A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT)
C1126   CONTINUE
C     ENDIF

       IF (MPIV.LT.SIZBLK) GO TO 9920
C*******************************
C End of static pivoting loop  *
C*******************************

  935   SCHUR = (NBLOC.LT.(NFRONT-NASS) .AND. NPIV.GE.NBLOC)
      IF (ICNTL(16).EQ.1) THEN
C Remove "zero" rows from fully summed rows within Schur complement
C ZCOUNT is count of "zero" rows
      ZCOUNT = 0
C APOS4 is beginning of block of fully summed uneliminated variables
      APOS4 = APOS + NPIV*NFRONT + NPIV

C Expand block so that lower triangle is included
C APOSB scans lower triangle by rows
C APOSC sweeps upper triangle by columns
      APOSB = APOS4 + NFRONT
      APOSC = APOS4 + 1
      DO 4444 I = 2,NASS-NPIV
        DO 4443 J = 1,I-1
          IF (KIND.EQ.1) THEN
          A(APOSB) = A(APOSC)
          ELSE
          A(APOSB) = CONJG(A(APOSC))
          ENDIF
          APOSB = APOSB + 1
          APOSC = APOSC + NFRONT
 4443     CONTINUE
        APOSB = APOS4 + NFRONT*I
        APOSC = APOS4 + I
 4444   CONTINUE

C Remove any zero rows by swopping with "first" row so that all zero
C    rows will be swept to beginning of block
C Row I is the row currently being checked
      I = NASS - NPIV
C Also exchange integer information accordingly
 4445   CONTINUE
      IF (ZCOUNT.EQ.I) GO TO 4450
C Check row I for zero
C APOSB is beginning of row I
      APOSB = APOS4 + (I-1)*NFRONT
      DO 4446 J = 1,NFRONT-NPIV
        IF (ABS(A(APOSB+J-1)).GT.TOL) GO TO 4449
 4446   CONTINUE
C Row is all zero
      ZCOUNT = ZCOUNT + 1
C Swop row ZCOUNT with row I
      IF (KIND.EQ.1) THEN
        DO 4447 J = 1,NFRONT-NPIV
          A(APOSB+J-1) = A(APOS4+NFRONT*(ZCOUNT-1)+J-1)
 4447     CONTINUE
      ELSE
C Have to split loop into three parts
        DO 4451 J = 1,I
          A(APOSB+J-1) = A(APOS4+NFRONT*(ZCOUNT-1)+J-1)
 4451     CONTINUE
C Deck does not test this
        DO 4452 J = I+1,ZCOUNT
          A(APOSB+J-1) = CONJG(A(APOS4+NFRONT*(ZCOUNT-1)+J-1))
 4452     CONTINUE
        DO 4453 J = ZCOUNT+1,NFRONT-NPIV
          A(APOSB+J-1) = A(APOS4+NFRONT*(ZCOUNT-1)+J-1)
 4453     CONTINUE
      ENDIF
C Zero row ZCOUNT
      DO 4448 J = 1,NFRONT-NPIV
        A(APOS4+NFRONT*(ZCOUNT-1)+J-1) = ZERO
 4448   CONTINUE
C Swop integers
      ISWOP = IW(IWPOS+NPIV+ZCOUNT-1)
      IW(IWPOS+NPIV+ZCOUNT-1) = IW(IWPOS+NPIV+I-1)
      IW(IWPOS+NPIV+I-1) = ISWOP
      GO TO 4445
 4449   I = I - 1
      GO TO 4445
 4450   CONTINUE
      ELSE
      ZCOUNT = 0
      ENDIF
C Set order of Schur complement (including rows of delayed pivots)
C But not including "zero" rows
      NSC1 = NFRONT - NPIV - ZCOUNT

C     IF (LCASE) THEN
C       write(7,'(/A,I8/A,5I8)') '*** Frontal matrix fter step',IASS,
C    *   'NFRONT',
C    *    NFRONT
C       write(7,'(/A/(10I8))') 'IW array',
C    *    (IW(IWPOS+I-1),I=1,NFRONT)
C       write(7,'(/A)') 'Frontal matrix'
C       DO 1127 J = 1, NFRONT
C       write(7,'(/(6E13.5))')
C    *    (A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT)
C1127   CONTINUE
C    ENDIF

C Accumulate total number of delayed pivots in INFO(23)
C Do not accumulate at last step where matrix is singular
      IF (IASS.NE.NSTEPS) INFO(23) = INFO(23) + NASS - NPIV
C SET LSTAT
      IF (CNTL(4).GT.ZERO .AND. INFO(23).GT.CNTL(5)*N) LSTAT = .TRUE.

C Jump if no Schur complement to form ... just store factors.
      IF (NSC1.EQ.0) GO TO 1830

C Save space for factors and Schur complement (if appropriate).

      IF (.NOT.SCHUR) THEN
C We now compute triangular Schur complement not using BLAS.
C This Schur complement is placed directly on the stack.

C Remove these
      RLSPA = MAX(RLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+
   +              NSTKAC(1))
      TRLSPA = MAX(TRLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+
   +              TOTSTA(1))
C ZCOUNT added because we do store this (for now)
      NSTKAC(1) = NSTKAC(1) + ((NSC1+ZCOUNT+1)*(NSC1+ZCOUNT))/2
      NSTACK(1) = NSTACK(1) + ((NSC1+ZCOUNT+1)*(NSC1+ZCOUNT))/2
      TOTSTA(1) = TOTSTA(1) + ((NSC1+ZCOUNT+1)*(NSC1+ZCOUNT))/2

C Initialize Schur complement
C Copying from the back to avoid overwriting important data
      APOSI = APOS + NFRONT*NFRONT - 1
      DO 1370 JJ = 1,NFRONT-NPIV
        J = APOSI
        DO 1360 JJJ = 1,JJ
          A(ASTK) = A(J)
          ASTK = ASTK - 1
          J = J - 1
 1360     CONTINUE
        APOSI = APOSI - NFRONT
 1370   CONTINUE
C APOS4 is the position in A of the first entry in the Schur complement.
      APOS4 = ASTK + 1


C Perform pivoting operations.
C Initialize variables
      J1 = IWPOS
      LTWO = .FALSE.
      POSPV1 = APOS
      DO 1450 I1 = 1,NPIV
        IF (LTWO) GO TO 1440
        APOSI = APOS + (I1-1)*NFRONT + NASS
        J2 = APOS + NFRONT* (I1-1) + NFRONT - 1
CCC What happens here ??
        APOSC = APOS4 +
   *     ((NASS-NPIV-ZCOUNT)*(2*NFRONT-NPIV-ZCOUNT-NASS+1))/2
C Check to see if current pivot is 1 x 1  or 2 x 2.
        IF (IW(J1).GT.0) THEN
          FLOPSB = FLOPSB + (NFRONT-NASS) +
   *                (NFRONT-NASS)* (NFRONT-NASS+1)
          DO 1410 JJ = APOSI,J2
          IF (KIND.EQ.1) THEN
            AMULT1 = -A(JJ)*A(POSPV1)
          ELSE
            AMULT1 = -CONJG(A(JJ))*A(POSPV1)
          ENDIF
          DO 1400 JJJ = JJ,J2
            A(APOSC) = A(APOSC) + AMULT1*A(JJJ)
            APOSC = APOSC + 1
 1400       CONTINUE
          IF (KIND.EQ.1) THEN
            A(JJ) = AMULT1
          ELSE
            A(JJ) = CONJG(AMULT1)
          ENDIF
 1410       CONTINUE
          J1 = J1 + 1
        ELSE
          POSPV2 = POSPV1 + NFRONT + 1
          OFFDAG = POSPV1 + 1
          FLOPSB = FLOPSB + 6* (NFRONT-NASS) +
   +           2* (NFRONT-NASS)* (NFRONT-NASS+1)
          DO 1430 JJ = APOSI,J2
          IF (KIND.EQ.1) THEN
            AMULT1 = - (A(POSPV1)*A(JJ)+A(OFFDAG)*A(JJ+NFRONT))
            AMULT2 = -A(POSPV2)*A(JJ+NFRONT) - A(OFFDAG)*A(JJ)
          ELSE
            AMULT1 = - CONJG(A(POSPV1)*A(JJ)
   *               + A(OFFDAG)*A(JJ+NFRONT))
            AMULT2 = - CONJG(A(POSPV2)*A(JJ+NFRONT)
   *               + CONJG(A(OFFDAG))*A(JJ))
          ENDIF
          DO 1420 JJJ = JJ,J2
            A(APOSC) = A(APOSC) + AMULT1*A(JJJ) +
   +               AMULT2*A(JJJ+NFRONT)
            APOSC = APOSC + 1
 1420       CONTINUE
          IF (KIND.EQ.1) THEN
            A(JJ) = AMULT1
            A(JJ+NFRONT) = AMULT2
          ELSE
            A(JJ) = CONJG(AMULT1)
            A(JJ+NFRONT) = CONJG(AMULT2)
          ENDIF
 1430       CONTINUE
          J1 = J1 + 2
          POSPV1 = POSPV2
          LTWO = .TRUE.
          GO TO 1450
        END IF

 1440     LTWO = .FALSE.
C Move to next pivot position
        POSPV1 = POSPV1 + NFRONT + 1
 1450   CONTINUE

C       IF (LCASE) THEN
C       write(7,'(A,I8)') 'GEMM not used at stage',IASS
C       write(7,'(A/(6E13.5))') 'Stacking Schur',
C    *       (A(I),I=APOS4,APOS4+(NSC1*(NSC1+1))/2-1)
C       ENDIF

      ELSE

C Action if SCHUR is true, We now use GEMM.

C Since SCHUR is true, copy U,
C   divide factors by D, and generate Schur complement using GEMM.
C   Then compress factors (to upper trapezoidal), and stack
C   half of the Schur complement.

C APOS4 is position in A of first entry of Schur complement to be
C   updated using GEMM.
      APOS4 = APOS+NASS*(NFRONT+1)

C APOS3 is the position in A of the first entry in the copy of U.
      APOS3 = APOS+NASS*NFRONT

C Copy U and divide factors by D
C Conjugating the factors would be good for SOLVE !!!
C Initialize variables
      J1 = IWPOS
      LTWO = .FALSE.
      POSPV1 = APOS
      DO 1490 I = 1,NPIV
        IF (LTWO) GO TO 1480
        APOSI = APOS + (I-1)*NFRONT + NASS
        POSELT = APOS3 + I - 1
C Check to see if current pivot is 1 x 1  or 2 x 2.
        IF (IW(J1).GT.0) THEN
          FLOPSB = FLOPSB + (NFRONT-NASS)
          DO 1460 JJ = APOSI,APOS + NFRONT*I - 1
          IF (KIND.EQ.1) THEN
            A(POSELT) = A(JJ)
          ELSE
            A(POSELT) = CONJG(A(JJ))
          ENDIF
          A(JJ) = -A(JJ)*A(POSPV1)
          POSELT = POSELT + NFRONT
 1460       CONTINUE
          J1 = J1 + 1
        ELSE
          POSPV2 = POSPV1 + NFRONT + 1
          OFFDAG = POSPV1 + 1
          FLOPSB = FLOPSB + 6* (NFRONT-NASS)
          DO 1470 JJ = APOSI,APOS + NFRONT*I - 1
          IF (KIND.EQ.1) THEN
            A(POSELT) = A(JJ)
C           A(POSELT+NFRONT-NASS) = A(JJ+NFRONT)
            A(POSELT+1) = A(JJ+NFRONT)
          ELSE
            A(POSELT) = CONJG(A(JJ))
            A(POSELT+1) = CONJG(A(JJ+NFRONT))
          ENDIF
          A(JJ) = - (A(POSPV1)*A(JJ)+A(OFFDAG)*A(JJ+NFRONT))
          IF (KIND.EQ.1) THEN
          A(JJ+NFRONT) = -A(POSPV2)*A(JJ+NFRONT) -
   +               A(OFFDAG)*A(POSELT)
          ELSE
          A(JJ+NFRONT) = -A(POSPV2)*A(JJ+NFRONT) -
   +               CONJG(A(OFFDAG)*A(POSELT))
          ENDIF
          POSELT = POSELT + NFRONT
 1470       CONTINUE
          J1 = J1 + 2
          POSPV1 = POSPV2
          LTWO = .TRUE.
          GO TO 1490
        END IF

 1480     LTWO = .FALSE.
C Move to next pivot position
        POSPV1 = POSPV1 + NFRONT + 1
 1490   CONTINUE

C Now create Schur complement by using Level 3 BLAS GEMM.
C Jump if Schur complement is null.
C Increment FLOPSB
      FLOPSB = FLOPSB + NPIV* (NFRONT-NASS)**2 +
   *              NPIV* (NFRONT-NASS)
C We divide the multiply into blocks to avoid too many extra
C    computations when using GEMM with a symmetric result.
C Block formed by GEMM has NBLOC rows.
      KBLK = ( NFRONT-NASS)/NBLOC
CCC     IF (LCASE) KBLK = 0
      L =  NFRONT - NASS
      DO 1500 KB = 1,KBLK
C Accumulate extra flops caused by using GEMM
        FLOPSX = FLOPSX + NBLOC* (NBLOC-1)* (NPIV)
        IF (KIND.EQ.1) THEN
          CALL CGEMM('N','N',L-(KB-1)*NBLOC,NBLOC,NPIV,CONE,
   +             A(APOS+NASS+NBLOC*(KB-1)),NFRONT,
   +             A(APOS3+NBLOC*(KB-1)*NFRONT),NFRONT,CONE,
   +             A(APOS4+NBLOC*(NFRONT+1)*(KB-1)),NFRONT)
        ELSE
          CALL CGEMM('N','N',L-(KB-1)*NBLOC,NBLOC,NPIV,CONE,
   +             A(APOS+NASS+NBLOC*(KB-1)),NFRONT,
   +             A(APOS3+NBLOC*(KB-1)*NFRONT),NFRONT,CONE,
   +             A(APOS4+NBLOC*(NFRONT+1)*(KB-1)),NFRONT)
        ENDIF
 1500   CONTINUE

C Calculate the block upper triangular part of the Schur complement.
      DO 1550 I = 1 + KBLK*NBLOC,L
C APOSA holds the index of the current entry of the matrix A used to
C   form the Schur complement as C = A*B
C APOSB holds the index of the current entry of the matrix B used to
C   form the Schur complement as C = A*B
C APOSC holds the index of the current entry of the matrix C used to
C   form the Schur complement as C = A*B
        APOSA = APOS + NASS
        APOSB = APOS3 +(I-1)*NFRONT
        APOSC = APOS4 + (I-1)*NFRONT - 1
        DO 1540 K = 1,NPIV
          DO 1530 J = I,L
          A(APOSC+J) = A(APOSC+J) + A(APOSA+J-1)*A(APOSB)
 1530       CONTINUE
          APOSA = APOSA + NFRONT
          APOSB = APOSB + 1
 1540     CONTINUE
 1550   CONTINUE


C Stack half of Schur complement.

C Stack reals
C Stack in reverse order to avoid need for compresses.
      JA1 = APOS+NFRONT*NFRONT-1
      NSTKAC(1) = NSTKAC(1) + ((NSC1+1)* (NSC1))/2
      NSTACK(1) = NSTACK(1) + ((NSC1+1)* (NSC1))/2
      TOTSTA(1) = TOTSTA(1) + ((NSC1+1)* (NSC1))/2
C Stack by rows
      DO 1710 I = NSC1,1,-1
        DO 1700 JJ = JA1,JA1-(NSC1-I),-1
          A(ASTK) = A(JJ)
          ASTK = ASTK - 1
 1700     CONTINUE
        JA1 = JA1 - NFRONT
 1710   CONTINUE

C       IF (LCASE) THEN
C       write(7,'(A,I8)') 'GEMM used at stage',IASS
C       write(7,'(A/(6E13.5))') 'Stacking Schur',
C    *       (A(I),I=ASTK+1,ASTK+(NSC1*(NSC1+1))/2)
C       ENDIF

C END of SCHUR being true action (started after label 1450)
      END IF

C Stack integers
      NSTKAC(2) = NSTKAC(2) + NSC1 + 1
      NSTACK(2) = NSTACK(2) + NSC1 + 1
      TOTSTA(2) = TOTSTA(2) + NSC1 + 1
C Record space needed to this point
 1830   IF (IASS.EQ.NSTEPS) THEN
      INTSPA = MAX(INTSPA,IWPOS+NFRONT-1+NSTKAC(2))
      TINSPA = MAX(TINSPA,IWPOS+NFRONT-1+TOTSTA(2))
      GO TO 2158
      ELSE
      INTSPA = MAX(INTSPA,IWPOS+NFRONT-1+(N-NTOTPV+2)+NSTKAC(2))
      TINSPA = MAX(TINSPA,IWPOS+NFRONT-1+(N-NTOTPV+2)+TOTSTA(2))
      ENDIF

C Check space and compress if necessary
  444   NST = 0
C +1 for length of stacked entry
      IF (NSC1.GT.0) NST = NSC1 + 1
      IF (IWPOS+NFRONT-1+(N-NTOTPV+2)+NST.GT.ISTK) THEN
C Compress integer storage
      CALL ME57P(A,IW,ISTK,IINPUT,PTRIRN,.FALSE.)
      NCMPBI = NCMPBI + 1

      IF (IWPOS+NFRONT-1+(N-NTOTPV+2)+NST.GT.ISTK) THEN
C Still insufficient space after compress
        IF (ICNTL(8).NE.0) THEN
          HOLD(1) = 2
          HOLD(2) = NBLK
          HOLD(3) = NTWO
          HOLD(4) = INFO(23)
          HOLD(5) = NCMPBI
          HOLD(6) = NEIG
          HOLD(7) = MAXFRT
          HOLD(8) = IWPOS
          HOLD(9) = APOS
          HOLD(10) = APOSBB
          HOLD(11) = NSTKAC(1)
          HOLD(12) = NSTKAC(2)
          HOLD(13) = AINPUT
          HOLD(14) = IINPUT
          HOLD(15) = ISTK
          HOLD(16) = ASTK
          HOLD(17) = INTSPA
          HOLD(18) = RLSPA
          HOLD(19) = PTRIRN
          HOLD(20) = PTRA
          HOLD(21) = NTOTPV
          HOLD(22) = NPOTPV
          HOLD(23) = NUMORG
          HOLD(24) = NFRONT
          HOLD(25) = NASS
C         HOLD(26) = NCOL
          HOLD(27) = NPIV
          HOLD(28) = IASS
          HOLD(29) = TINSPA
          HOLD(30) = TRLSPA
          HOLD(31) = TOTSTA(1)
          HOLD(32) = TOTSTA(2)
          HOLD(33) = NSTACK(1)
          HOLD(34) = NSTACK(2)
          IF (ICNTL(7).GT.1) HOLD(35) = ISNPIV
          IF (ICNTL(7).EQ.4) HOLD(36) = PHASE
          HOLD(37) = INFO(32)
          HOLD(38) = INFO(33)
          HOLD(39) = INFO(34)
          NSC1    = NFRONT-NPIV
          RINFO(3) =FLOPSA
          RINFO(4) =FLOPSB
          RINFO(5) =FLOPSX
          INFO(1) = 11
          HOLD(40) = NBSTATIC
          INFO(35) = HOLD(40)
        ELSE
          INFO(1)  = -4
          INFO(2)  = LIW
          INFO(18) = INTSPA
        ENDIF
        RETURN
      END IF
      END IF

      IF (NSC1.GT.0) THEN
      DO 1720 I = 1,NSC1
        IW(ISTK) = IW(IWPOS+NFRONT-I)
        ISTK = ISTK - 1
 1720   CONTINUE
      IW(ISTK) = NSC1
      ISTK = ISTK - 1
      ENDIF

C    IF (LCASE) THEN
C       write(7,'(/A,I8/A,5I8)') '*** Frontal matrix end step',IASS,
C    *   'NFRONT',
C    *    NFRONT
C       write(7,'(/A/(10I8))') 'IW array',
C    *    (IW(IWPOS+I-1),I=1,NFRONT)
C       write(7,'(/A)') 'Frontal matrix'
C       DO 1128 J = 1, NFRONT
C       write(7,'(/(6E13.5))')
C    *    (A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT)
C1128   CONTINUE
C    ENDIF

C Reset PPOS.
      DO 1840 JJ = IWPOS + NPIV,IWPOS + NFRONT - 1
      J = ABS(IW(JJ))
      PPOS(J) = N + 1
 1840   CONTINUE


C********************************
C    STORE FACTORS
C********************************
C Complete the integer information in the factors
 2158   IF (NPIV.EQ.0) GO TO 2159
      NBLK = NBLK + 1

      IW(IWPOS-2) = NFRONT
      IW(IWPOS-1) = NPIV
      IWPOS = IWPOS + NFRONT + 2

C Store information on the reals for the factors.
C We copy from A(JA1) to A(APOS2) ... the use of buffer space from
C   APOSBB to APOS ensures no overwrites.
      IF (INFO(1).EQ.-3) THEN
      INFO(40) = INFO(40) + (NPIV * (2*NFRONT-NPIV+1))/2
      GO TO 2159
      END IF

      APOS2 = APOSBB
C Store reals from full pivot.
      DO 2130 I = 1,NPIV
C JA1 points to the diagonal
      JA1 = APOS + (I-1)* (NFRONT+1)
      DO 2120 J = I,NPIV
C Store conjugate if matrix Hermitian
        IF (KIND.EQ.1) THEN
          A(APOS2) = A(JA1)
        ELSE
          A(APOS2) = CONJG(A(JA1))
        ENDIF
        IF (A(APOS2).EQ.ZERO) INFO(32) = INFO(32) + 1
        APOS2 = APOS2 + 1
        JA1 = JA1 + 1
 2120   CONTINUE
 2130   CONTINUE
      RPOS = APOS2
C Store rectangle
      DO 2150 I = 1,NPIV
      JA1 = APOS + (I-1)*NFRONT + NPIV
      DO 2140 J = 1,NFRONT - NPIV
        IF (KIND.EQ.1) THEN
          A(APOS2) = A(JA1)
        ELSE
          A(APOS2) = CONJG(A(JA1))
        ENDIF
        APOS2 = APOS2 + 1
        JA1 = JA1 + 1
 2140   CONTINUE
 2150   CONTINUE
C Set APOSBB for next block of factors
      APOSBB = APOS2

C Check rectangle for zeros
      DO 2152 J = 1,NFRONT-NPIV
      APOS2 = RPOS+J-1
      ZCOL = 1
      DO 2151 I = 1,NPIV
        IF (A(APOS2).EQ.ZERO) INFO(33) = INFO(33)+1
        IF (A(APOS2).NE.ZERO) ZCOL = 0
        APOS2 = APOS2 + NFRONT - NPIV
 2151   CONTINUE
      IF (ZCOL.EQ.1) INFO(34) = INFO(34)+1
 2152   CONTINUE

 2159   IASS = IASS + 1
    IF (IASS.LE.NSTEPS) THEN
C Initialize to zero to avoid problem when calling ME57E
      IW(IWPOS-2) = 0
      IW(IWPOS-1) = 0
      GO TO 2160
    ENDIF
C2160 CONTINUE
C
C End of loop on tree nodes.
C

    INFO(35) = NBSTATIC
    IF (INFO(1).EQ.-3) THEN
      INFO(2)  = LA
      INFO(17) = MAX(INFO(17),RLSPA)
CCC Should use ALENB in BD for this?
      IF (ICNTL(7).EQ.4) INFO(17) = MAX(INFO(17),RLSPA + N)
      RETURN
    END IF
    GO TO 1000
 970  INFO(1) = -5
    INFO(2) = NTOTPV + 1
    IF (LDIAG.GT.0 .AND. LP.GE.0)
   *    WRITE(LP,99992) INFO(1),PIVOT,CNTL(2),INFO(2),ICNTL(7)
99992 FORMAT (/'*** Error message from routine ME57B **',
   *     '   INFO(1) = ',I3/'Pivot has value ',2E16.8/'when ',
   *     'CNTL(2) has value ',E16.8/
   *     'at stage',I11,2X,'when ICNTL(7) =',I3)
    RETURN
 980  INFO(1) = -6
    INFO(2) = NTOTPV + 1
    IF (LDIAG.GT.0 .AND. LP.GE.0)
   *    WRITE(LP,99993) INFO(1),INFO(2),ICNTL(7)
99993 FORMAT (/'*** Error message from routine ME57B **',
   *     '   INFO(1) = ',I3/'Change in sign of pivot at stage',
   *     I10,2X,'when ICNTL(7) = ',I3)
    RETURN
 1000 NRLBDU = APOSBB - 1
    NIRBDU = IWPOS - 3
    IF (NTOTPV.NE.N) THEN
      INFO(1) = 4
      IF (LDIAG.GT.0 .AND. WP.GE.0)
   *    WRITE(WP,99994) INFO(1),NTOTPV
99994 FORMAT (/'*** Warning message from routine ME57B **',
   *       '   INFO(1) =',I2/5X, 'Matrix is singular, rank =', I5)

    ENDIF
C Recent change was to remove condition that ICNTL(16) was equal to 1
C More recent change to remove deficiency test.  This change means that
C we now test there is sufficient room to move the off-diagonal entries
C of the two by two pivots.
C 555 IF (NTOTPV.NE.N) THEN
C Check space (by this time there is nothing to compress)
  555 NRLBDU = APOSBB - 1
    NIRBDU = IWPOS - 3
      IF (NIRBDU+3*(N-NTOTPV) .GT. LIW
   +    .OR. NRLBDU+(N-NTOTPV)+NTWO .GT. LA) THEN
C I don't think this can happen ... at least I can't make it happen
C It is left in for "safety" :-)
C Still insufficient space after compress
      IF (ICNTL(8).NE.0) THEN
        HOLD(1) = 3
        HOLD(2) = NBLK
        HOLD(3) = NTWO
        HOLD(4) = INFO(23)
        HOLD(5) = NCMPBI
        HOLD(6) = NEIG
        HOLD(7) = MAXFRT
        HOLD(8) = IWPOS
        HOLD(9) = APOS
        HOLD(10) = APOSBB
        HOLD(11) = NSTKAC(1)
        HOLD(12) = NSTKAC(2)
        HOLD(13) = AINPUT
        HOLD(14) = IINPUT
        HOLD(15) = ISTK
        HOLD(16) = ASTK
        HOLD(17) = INTSPA
        HOLD(18) = RLSPA
        HOLD(19) = PTRIRN
        HOLD(20) = PTRA
        HOLD(21) = NTOTPV
        HOLD(22) = NPOTPV
        HOLD(23) = NUMORG
        HOLD(24) = NFRONT
        HOLD(25) = NASS
C         HOLD(26) = NCOL
        HOLD(27) = NPIV
        HOLD(28) = IASS
        HOLD(29) = TINSPA
        HOLD(30) = TRLSPA
        HOLD(31) = TOTSTA(1)
        HOLD(32) = TOTSTA(2)
        HOLD(33) = NSTACK(1)
        HOLD(34) = NSTACK(2)
        IF (ICNTL(7).GT.1) HOLD(35) = ISNPIV
        IF (ICNTL(7).EQ.4) HOLD(36) = PHASE
        HOLD(37) = INFO(32)
        HOLD(38) = INFO(33)
        HOLD(39) = INFO(34)
        NSC1    = NFRONT-NPIV
        RINFO(3) =FLOPSA
        RINFO(4) =FLOPSB
        RINFO(5) =FLOPSX
        IF (NRLBDU+(N-NTOTPV)+NTWO .GT. LA) INFO(1) = 10
        IF (NIRBDU+3*(N-NTOTPV) .GT. LIW)  INFO(1) = 11
        HOLD(40) = NBSTATIC
        INFO(35) = HOLD(40)
      ELSE
        IF (NIRBDU+3*(N-NTOTPV) .GT. LIW) THEN
          INFO(1)  = -4
          INFO(2)  = LIW
          INFO(18) = MAX(INTSPA,NIRBDU+3*(N-NTOTPV))
        ELSE
          INFO(1)  = -3
          INFO(2) = LA
          INFO(17) = MAX(INFO(17),RLSPA,NRLBDU+(N-NTOTPV)+NTWO)
          IF (ICNTL(7).EQ.4) INFO(17) =
   +      MAX(INFO(17),RLSPA + N,NRLBDU+(N-NTOTPV)+NTWO)
        ENDIF
      ENDIF
      RETURN
      ENDIF
C   ENDIF

C Add explicit entries in factors for zero pivots (now set to 1.0)
C Initialize flag array to identify indices of zero pivots
    IF (NTOTPV.NE.N) THEN
      DO 3331 I = 1,N
      PPOS(I) = 0
 3331   CONTINUE
      IWPOS = 4
      DO 3332 I = 1,NBLK
      NFRONT = IW(IWPOS)
      NPIV = IW(IWPOS+1)
      DO 3330 J = IWPOS+2,IWPOS+NPIV+1
        PPOS(ABS(IW(J))) = 1
 3330   CONTINUE
      IWPOS = IWPOS + NFRONT + 2
 3332   CONTINUE
      K= 0
      DO 3333 I=1,N
      IF (PPOS(I).EQ.0) THEN
        K=K+1
        NBLK = NBLK + 1
        NRLBDU = NRLBDU+1
        A(NRLBDU) = ONE
        IW(NIRBDU+1) = 1
        IW(NIRBDU+2) = 1
        IW(NIRBDU+3) = I
        NIRBDU = NIRBDU+3
      ENDIF
 3333   CONTINUE
    ENDIF

C
    INFO(14) = NRLBDU
C   Move the off-diagonal entries of the 2x2 pivots within full blocks
C   to the end of A.
    IW(1) = NRLBDU + 1
    IW(2) = NRLBDU + NTWO
    INFO(15) = IW(2)
    IW(3) = NBLK
    INFO(31) = NBLK
C   Negate the entries of L, move the off-diagonal entries of the 2x2
C   pivots within full blocks to the end of A and update NRLBDU to
C   correspond
    CALL ME57W(KIND,A,LA,IW,LIW,NRLBDU)
    INFO(16) = NIRBDU
    INFO(18) = INTSPA
    INFO(20) = TINSPA
    INFO(17) = RLSPA
    INFO(19) = TRLSPA
    INFO(21) = MAXFRT
    INFO(22) = NTWO
C   INFO(23)  .. computed as sum of NASS-NPIV
    INFO(24) = NEIG
    INFO(25) = NTOTPV
    INFO(28) = NCMPBR
    INFO(29) = NCMPBI
    RINFO(3) = FLOPSA
    RINFO(4) = FLOPSB
    RINFO(5) = FLOPSX
    IF (INFO(27).GT.0) THEN
      RINFO(14) = ZERO
      DO 332 I = 1,N
      RINFO(14) = MAX(RINFO(14),REAL(DIAG(I)))
 332    CONTINUE
    ENDIF

    RETURN

    END


    SUBROUTINE ME57P(A,IW,J1,J2,ITOP,REAL)
C This subroutine performs a very simple compress (block move).
C   Entries J1+1 to J2-1 (incl.) in A or IW as appropriate are moved
C   to occupy the positions immediately prior to position ITOP.
C A/IW hold the array being compressed.
C J1/J2 define the entries being moved.
C ITOP defines the position immediately after the positions to which
C   J1 to J2 are moved.
C REAL must be set by the user to .TRUE. if the move is on array A,
C   a value of .FALSE. will perform the move on A.
C   .. Scalar Arguments ..
    INTEGER ITOP,J1,J2
    LOGICAL REAL
C   ..
C   .. Array Arguments ..
    COMPLEX A(*)
    INTEGER IW(*)
C   ..
C   .. Local Scalars ..
    INTEGER IPOS,JJ
C   ..
C   .. Executable Statements ..
    IF (J2.EQ.ITOP) GO TO 50
    IPOS = ITOP - 1
    IF (REAL) THEN
      DO 10 JJ = J2-1,J1+1,-1
      A(IPOS) = A(JJ)
      IPOS = IPOS - 1
   10   CONTINUE
    ELSE
      DO 20 JJ = J2-1,J1+1,-1
      IW(IPOS) = IW(JJ)
      IPOS = IPOS - 1
   20   CONTINUE
    ENDIF
    J2 = ITOP
    J1 = IPOS
   50 RETURN
    END
    SUBROUTINE ME57W(KIND,A,LA,IW,LIW,NRLBDU)
C   Negate the entries of L, move the off-diagonal entries of the 2x2
C   pivots within full blocks to the end of A and update NRLBDU to
C   correspond.
    INTEGER KIND,LA,LIW
    COMPLEX A(LA)
    INTEGER IW(LIW)

    INTEGER NRLBDU

C Constants
    REAL ZERO
    PARAMETER (ZERO=0.0E0)
    COMPLEX CZERO
    PARAMETER (CZERO=(0.0E0,0.0E0))
C Local variables
    INTEGER APOS,IBLK,IROW,IWPOS,J,JPIV,NCOLS,NROWS
C APOS  Position in A of current diagonal entry.
C IBLK  Current block.
C IROW  Current row.
C IWPOS Current position in IW.
C J   Do loop variable
C JPIV  Used as a flag so that IPIV is incremented correctly after the
C     use of a 2 by 2 pivot.
C NCOLS Number of columns in the block.
C NROWS Number of rows in the block.

    APOS = 1
    IWPOS = 6
    DO 40 IBLK = 1,IW(3)
      NCOLS = IW(IWPOS-2)
      NROWS = IW(IWPOS-1)
      JPIV = 1
      DO 30 IROW = 1,NROWS
      JPIV = JPIV - 1
      IF (JPIV.EQ.1) GO TO 10
      IF (IW(IWPOS+IROW-1).LT.0) THEN
        JPIV = 2
        NRLBDU = NRLBDU + 1
CC Must reconjugate if KIND not equal to 1
        IF (KIND.EQ.1) THEN
          A(NRLBDU) = A(APOS+1)
        ELSE
          A(NRLBDU) = CONJG(A(APOS+1))
        ENDIF
        A(APOS+1) = ZERO
      END IF

   10   DO 20 J = APOS + 1,APOS + NROWS - IROW
        A(J) = -A(J)
   20   CONTINUE
      APOS = APOS + NROWS - IROW + 1
   30   CONTINUE
C Negate entries in rectangular block (was done earlier by MA47OD)
C     DO 35 J = APOS,APOS+NROWS*(NCOLS-NROWS)-1
C       A(J) = -A(J)
C  35   CONTINUE
      APOS = APOS + NROWS* (NCOLS-NROWS)
      IWPOS = IWPOS + NCOLS + 2
   40 CONTINUE
    END

    SUBROUTINE ME57X(N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
   *            W,LW,IW1,ICNTL)
C This subroutine performs forward elimination using the factors
C   stored in FACT/IFACT by ME57B.
C It is designed for efficiency on one right-hand side.
    INTEGER N,LFACT
    COMPLEX FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT),LRHS,LW
    COMPLEX W(LW),RHS(LRHS)
    INTEGER IW1(N),ICNTL(20)
C N   must be set to the order of the matrix. It is not altered.
C FACT   must be set to hold the real values corresponding to the
C   factors. This must be unchanged since the preceding call to
C   ME57B. It is not altered.
C LFACT  length of array FACT. It is not altered.
C IFACT  holds the integer indexing information for the matrix factors
C   in FACT. This must be unchanged since the preceding call to
C   ME57B. It is not altered.
C LIFACT length of array IFACT. It is not altered.
C RHS on input, must be set to hold the right hand side vector.  On
C   return, it will hold the modified vector following forward
C   elimination.
C LHS must be set to the leading dimension of array RHS.
C W   used as workspace to hold the components of the right hand
C   sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C   larger than INFO(21) as returned from ME57B.
C IW1 need not be set on entry. On exit IW1(I) (I = 1,NBLK), where
C   NBLK = IFACT(3) is the number of block pivots, will
C   hold pointers to the beginning of each block pivot in array IFACT.
C ICNTL Not referenced except:
C   ICNTL(13) Threshold on number of columns in a block for using
C       addressing using Level 2 and Level 3 BLAS.
C

C Procedures
    INTRINSIC ABS
    EXTERNAL CGEMV,CTPSV

C Constant
    REAL ONE
    PARAMETER (ONE=1.0E0)
    COMPLEX CONE
    PARAMETER (CONE=(1.0E0,0.0E0))
C
C Local variables
    INTEGER APOS,I,IBLK,II,IPIV,IRHS,IWPOS,J,J1,J2,K,K1,K2,
   +      NCOLS,NROWS
    COMPLEX W1,W2
C
C APOS  Current position in array FACT.
C I   Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IWPOS Position in IFACT of start of current index list.
C J   Temporary DO index
C K   Temporary pointer to position in real array.
C J1    Position in IFACT of index of leading entry of row.
C J2    Position in IFACT of index of trailing entry of row.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.

    APOS = 1
    IWPOS = 4
    DO 270 IBLK = 1,IFACT(3)

C Find the number of rows and columns in the block.
      IW1(IBLK) = IWPOS
      NCOLS = IFACT(IWPOS)
      NROWS = IFACT(IWPOS+1)
      IWPOS = IWPOS + 2

      IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN


C Perform operations using direct addressing.

C Load appropriate components of right-hand sides into W.
      DO 10 I = 1,NCOLS
        II = ABS(IFACT(IWPOS+I-1))
        W(I) = RHS(II)
   10   CONTINUE

C Treat diagonal block (direct addressing)
      CALL CTPSV('L','N','U',NROWS,FACT(APOS),W,1)
      APOS = APOS + (NROWS* (NROWS+1))/2

C Treat off-diagonal block (direct addressing)
C       IF (NCOLS.GT.NROWS) CALL CGEMM('N','N',NCOLS-NROWS,1,NROWS,
C    +                      CONE,FACT(APOS),NCOLS-NROWS,
C    +                      W,LW,CONE,W(NROWS+1),LW)
      IF (NCOLS.GT.NROWS) THEN
        CALL CGEMV('N',NCOLS-NROWS,NROWS,
   +           CONE,FACT(APOS),NCOLS-NROWS,
   +           W,1,CONE,W(NROWS+1),1)
      ENDIF
      APOS = APOS + NROWS* (NCOLS-NROWS)

C Reload W back into RHS.
      DO 35 I = 1,NCOLS
        II = ABS(IFACT(IWPOS+I-1))
        RHS(II) = W(I)
   35   CONTINUE

      ELSE

C Perform operations using indirect addressing.

      J1 = IWPOS
      J2 = IWPOS + NROWS - 1


C Treat diagonal block (indirect addressing)
      DO 130 IPIV = 1,NROWS
      APOS = APOS + 1
      W1 = RHS(ABS(IFACT(J1)))
      K = APOS
      DO 100 J = J1+1,J2
        IRHS = ABS(IFACT(J))
        RHS(IRHS) = RHS(IRHS) - FACT(K)*W1
        K = K + 1
  100   CONTINUE
      APOS = K
      J1 = J1 + 1
  130   CONTINUE

C Treat off-diagonal block (indirect addressing)
C     J2 = IWPOS + NCOLS - 1
C     DO 136 IPIV = 1,NROWS
C       K = APOS
C       W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)))
C       DO 133 J = J1,J2
C       IRHS = ABS(IFACT(J))
C       RHS(IRHS) = RHS(IRHS) + W1*FACT(K)
C       K = K + 1
C 133   CONTINUE
C       APOS = K
C 136   CONTINUE

C Loop unrolling
      J2 = IWPOS + NCOLS - 1
      DO 136 IPIV = 1,NROWS-1,2
      K1 = APOS
      K2 = APOS+NCOLS-NROWS
      W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)))
      W2 = RHS(ABS(IFACT(IWPOS+IPIV)))
      DO 133 J = J1,J2
        IRHS = ABS(IFACT(J))
        RHS(IRHS) = RHS(IRHS) + W1*FACT(K1) + W2*FACT(K2)
        K1 = K1 + 1
        K2 = K2 + 1
  133   CONTINUE
      APOS = K2
  136   CONTINUE

      IF (MOD(NROWS,2).EQ.1) THEN
      K = APOS
      W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)))
      DO 137 J = J1,J2
        IRHS = ABS(IFACT(J))
        RHS(IRHS) = RHS(IRHS) + W1*FACT(K)
        K = K + 1
  137   CONTINUE
      APOS = K
      ENDIF
    END IF

    IWPOS = IWPOS + NCOLS
  270 CONTINUE

    END


    SUBROUTINE ME57Y(KIND,N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
   *            W,LW,IW1,ICNTL)
C This subroutine performs backward elimination operations
C   using the factors stored in FACT/IFACT by ME57B.
C It is designed for efficiency on one right-hand side.
    INTEGER KIND,N,LFACT
    COMPLEX FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT),LRHS,LW
    COMPLEX W(LW),RHS(LRHS)
    INTEGER IW1(N),ICNTL(20)
C N    must be set to the order of the matrix. It is not altered.
C FACT    must be set to hold the real values corresponding to the
C    factors. This must be unchanged since the
C    preceding call to ME57B. It is not altered.
C LFACT   length of array FACT. It is not altered.
C IFACT   holds the integer indexing information for the matrix factors
C    in FACT. This must be unchanged since the preceding call to
C    ME57B. It is not altered.
C LIFACT  length of array IFACT. It is not altered.
C RHS  on entry, must be set to hold the right hand side modified by
C    the forward substitution operations. On exit, holds the
C    solution vector.
C LHS must be set to the leading dimension of array RHS.
C W    used as workspace to hold the components of the right hand
C    sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C   larger than INFO(21) as returned from ME57B.
C IW1  on entry IW1(I) (I = 1,NBLK), where  NBLK = IFACT(3) is the
C   number of block pivots, must hold pointers to the beginning of
C   each block pivot in array IFACT, as set by ME57X/XD. It is not
C   altered.
C ICNTL Not referenced except:
C   ICNTL(13) Threshold on number of columns in a block for using
C       addressing using Level 2 and Level 3 BLAS.

C Procedures
    INTRINSIC ABS
    EXTERNAL CGEMV,CTPSV

C Constants
    REAL ONE
    PARAMETER (ONE=1.0E0)
    COMPLEX CONE
    PARAMETER (CONE=(1.0E0,0.0E0))
C
C Local variables.
    INTEGER APOS,APOS2,I,IBLK,II,IPIV,IRHS,IRHS1,
   +      IRHS2,IWPOS,J,JPIV,J1,J2,K,K2,LROW,NCOLS,NROWS
    COMPLEX W1,W2
C APOS  Current position in array FACT.
C APOS2 Current position in array FACT for off-diagonal entry of 2x2
C     pivot.
C I   Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IRHS1 RHS index.
C IRHS2 RHS index.
C IWPOS Position in IFACT of start of current index list.
C J   Temporary DO index
C JPIV  Has the value -1 for the first row of a 2 by 2 pivot and 1 for
C     the second.
C K   Temporary pointer to position in real array.
C J1    Position in IFACT of index of leading entry of row.
C J2    Position in IFACT of index of trailing entry of row.
C K   Temporary variable.
C LROW  Length of current row.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.
C

    APOS = IFACT(1)
    APOS2 = IFACT(2)
C Run through block pivot rows in the reverse order.
    DO 380 IBLK = IFACT(3),1,-1

C Find the number of rows and columns in the block.
      IWPOS = IW1(IBLK)
      NCOLS = ABS(IFACT(IWPOS))
      NROWS = ABS(IFACT(IWPOS+1))
      APOS = APOS - NROWS* (NCOLS-NROWS)
      IWPOS = IWPOS + 2

      IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN


C Perform operations using direct addressing.

C Load latter part of right-hand side into W.
      DO 5 I = NROWS + 1,NCOLS
        II = ABS(IFACT(IWPOS+I-1))
        W(I) = RHS(II)
    5   CONTINUE


C Multiply by the diagonal matrix (direct addressing)
      DO 10 IPIV = NROWS,1,-1
        IRHS = ABS(IFACT(IWPOS+IPIV-1))
        APOS = APOS - (NROWS+1-IPIV)
        W(IPIV) = RHS(IRHS)*FACT(APOS)
   10   CONTINUE
      JPIV = -1
      DO 20 IPIV = NROWS,1,-1
        IRHS = IFACT(IWPOS+IPIV-1)
        IF (IRHS.LT.0) THEN
          IRHS1 = -IFACT(IWPOS+IPIV-1+JPIV)
          IF (KIND.EQ.1 .OR. JPIV.EQ.1) THEN
          W(IPIV) = RHS(IRHS1)*FACT(APOS2) + W(IPIV)
          ELSE
          W(IPIV) = RHS(IRHS1)*CONJG(FACT(APOS2)) + W(IPIV)
          ENDIF
          IF (JPIV.EQ.1) APOS2 = APOS2 - 1
          JPIV = -JPIV
        END IF

   20   CONTINUE

C Treat off-diagonal block (direct addressing)
      K = NCOLS - NROWS
      IF (K.GT.0) THEN
        IF (KIND.EQ.1) THEN
          CALL CGEMV('T',K,NROWS,CONE,
   +             FACT(APOS+(NROWS*(NROWS+1))/2),K,
   +             W(NROWS+1),1,CONE,W,1)
        ELSE
          CALL CGEMV('C',K,NROWS,CONE,
   +             FACT(APOS+(NROWS*(NROWS+1))/2),K,
   +             W(NROWS+1),1,CONE,W,1)
        ENDIF
      ENDIF
C       IF (K.GT.0) CALL CGEMM('T','N',NROWS,1,K,CONE,
C    +                 FACT(APOS+(NROWS*(NROWS+1))/2),K,
C    +                 W(NROWS+1),LW,CONE,W,LW)

C Treat diagonal block (direct addressing)
      IF (KIND.EQ.1) CALL CTPSV('L','T','U',NROWS,FACT(APOS),W,1)
      IF (KIND.NE.1) CALL CTPSV('L','C','U',NROWS,FACT(APOS),W,1)
C Reload W back into RHS.
      DO 60 I = 1,NROWS
        II = ABS(IFACT(IWPOS+I-1))
        RHS(II) = W(I)
   60   CONTINUE

      ELSE
C
C Perform operations using indirect addressing.
      J1 = IWPOS
      J2 = IWPOS + NCOLS - 1


C Multiply by the diagonal matrix (indirect addressing)
      JPIV = -1
      IF (KIND.EQ.1) THEN
        DO 210 IPIV = NROWS,1,-1
          IRHS = IFACT(IWPOS+IPIV-1)
          LROW = NROWS + 1 - IPIV

          IF (IRHS.GT.0) THEN
C 1 by 1 pivot.
          APOS = APOS - LROW
          RHS(IRHS) = RHS(IRHS)*FACT(APOS)
          ELSE
C 2 by 2 pivot
          IF (JPIV.EQ.-1) THEN
            IRHS1 = -IFACT(IWPOS+IPIV-2)
            IRHS2 = -IRHS
            APOS = APOS - LROW - LROW - 1
            W1 = RHS(IRHS1)*FACT(APOS) +
   +           RHS(IRHS2)*FACT(APOS2)
            RHS(IRHS2) = RHS(IRHS1)*FACT(APOS2) +
   +               RHS(IRHS2)*FACT(APOS+LROW+1)
            RHS(IRHS1) = W1
            APOS2 = APOS2 - 1
          END IF
          JPIV = -JPIV
          END IF

  210     CONTINUE
      ELSE
        DO 211 IPIV = NROWS,1,-1
          IRHS = IFACT(IWPOS+IPIV-1)
          LROW = NROWS + 1 - IPIV

          IF (IRHS.GT.0) THEN
C 1 by 1 pivot.
          APOS = APOS - LROW
          RHS(IRHS) = RHS(IRHS)*FACT(APOS)
          ELSE
C 2 by 2 pivot
          IF (JPIV.EQ.-1) THEN
            IRHS1 = -IFACT(IWPOS+IPIV-2)
            IRHS2 = -IRHS
            APOS = APOS - LROW - LROW - 1
            W1 = RHS(IRHS1)*FACT(APOS) +
   +           RHS(IRHS2)*FACT(APOS2)
            RHS(IRHS2) = RHS(IRHS1)*CONJG(FACT(APOS2)) +
   +               RHS(IRHS2)*FACT(APOS+LROW+1)
            RHS(IRHS1) = W1
            APOS2 = APOS2 - 1
          END IF
          JPIV = -JPIV
          END IF

  211     CONTINUE
      ENDIF
      APOS = APOS + (NROWS* (NROWS+1))/2

C Treat off-diagonal block (indirect addressing)
C       KK = APOS
C       J1 = IWPOS + NROWS
C       DO 220 IPIV = 1,NROWS
C       IRHS = ABS(IFACT(IWPOS+IPIV-1))
C       W1 = RHS(IRHS)
C       K = KK
C       DO 215 J = J1,J2
C         W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)))
C         K = K + 1
C 215     CONTINUE
C       RHS(IRHS) = W1
C       KK = K
C 220   CONTINUE

C Loop unrolling
      K = APOS
      J1 = IWPOS + NROWS
      DO 220 IPIV = 1,NROWS-1,2
        IRHS = ABS(IFACT(IWPOS+IPIV-1))
        W1 = RHS(IRHS)
        IRHS1 = ABS(IFACT(IWPOS+IPIV))
        W2 = RHS(IRHS1)
        K2 = K+(NCOLS-NROWS)
        DO 215 J = J1,J2
          II = ABS(IFACT(J))
          IF (KIND.EQ.1) THEN
          W1 = W1 + FACT(K)*RHS(II)
          W2 = W2 + FACT(K2)*RHS(II)
          ELSE
          W1 = W1 + CONJG(FACT(K))*RHS(II)
          W2 = W2 + CONJG(FACT(K2))*RHS(II)
          ENDIF
          K = K + 1
          K2 = K2 + 1
  215     CONTINUE
        RHS(IRHS) = W1
        RHS(IRHS1) = W2
        K = K2
  220   CONTINUE

      IF (MOD(NROWS,2).EQ.1) THEN
        IRHS = ABS(IFACT(IWPOS+IPIV-1))
        W1 = RHS(IRHS)
        DO 216 J = J1,J2
          IF (KIND.EQ.1) THEN
          W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)))
          ELSE
          W1 = W1 + CONJG(FACT(K))*RHS(ABS(IFACT(J)))
          ENDIF
          K = K + 1
  216     CONTINUE
        RHS(IRHS) = W1
      ENDIF

C Treat diagonal block (indirect addressing)
      J2 = IWPOS + NROWS - 1
      DO 260 IPIV = 1,NROWS
        IRHS = ABS(IFACT(J1-1))
        APOS = APOS - IPIV
        W1 = RHS(IRHS)
        K = APOS + 1
        DO 230 J = J1,J2
          IF (KIND.EQ.1) THEN
          W1 = W1 - FACT(K)*RHS(ABS(IFACT(J)))
          ELSE
          W1 = W1 - CONJG(FACT(K))*RHS(ABS(IFACT(J)))
          ENDIF
          K = K + 1
  230     CONTINUE
        RHS(IRHS) = W1
        J1 = J1 - 1
  260   CONTINUE

      END IF

  380 CONTINUE

    END
    SUBROUTINE ME57V(KIND,N,NZ,IRN,ICN,IW,LW,IPE,IQ,FLAG,IWFR,
   +           ICNTL,INFO)
C Is identical to subroutine MA27GD.  Internal version for ME57.
C Except for addition of parameter KIND and flag if entry in lower
C   triangle
C
C SORT PRIOR TO CALLING ANALYSIS ROUTINE MA27H/HD (internal ME57
C subroutine ME57H/HD).
C
C GIVEN THE POSITIONS OF THE OFF-DIAGONAL NON-ZEROS OF A SYMMETRIC
C   MATRIX, CONSTRUCT THE SPARSITY PATTERN OF THE OFF-DIAGONAL
C   PART OF THE WHOLE MATRIX (UPPER AND LOWER TRIANGULAR PARTS).
C EITHER ONE OF A PAIR (I,J),(J,I) MAY BE USED TO REPRESENT
C   THE PAIR. DIAGONAL ELEMENTS AND DUPLICATES ARE IGNORED.
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT. IT IS NOT
C   ALTERED.
C IRN(I),I=1,2,...,NZ MUST BE SET TO THE ROW NUMBERS OF THE
C   NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
C   TO IW (SEE DESCRIPTION OF IW).
C ICN(I),I=1,2,...,NZ MUST BE SET TO THE COLUMN NUMBERS OF THE
C   NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
C   TO IW (SEE DESCRIPTION OF IW).
C IW NEED NOT BE SET ON INPUT. ON OUTPUT IT CONTAINS LISTS OF
C   COLUMN INDICES, EACH LIST BEING HEADED BY ITS LENGTH.
C   IRN(1) MAY BE EQUIVALENCED TO IW(1) AND ICN(1) MAY BE
C   EQUIVALENCED TO IW(K), WHERE K.GT.NZ.
C LW MUST BE SET TO THE LENGTH OF IW. IT MUST BE AT LEAST 2*NZ+N.
C   IT IS NOT ALTERED.
C IPE NEED NOT BE SET ON INPUT. ON OUTPUT IPE(I) POINTS TO THE START OF
C   THE ENTRY IN IW FOR ROW I, OR IS ZERO IF THERE IS NO ENTRY.
C IQ NEED NOT BE SET.  ON OUTPUT IQ(I),I=1,N CONTAINS THE NUMBER OF
C   OFF-DIAGONAL N0N-ZEROS IN ROW I INCLUDING DUPLICATES.
C FLAG IS USED FOR WORKSPACE TO HOLD FLAGS TO PERMIT DUPLICATE ENTRIES
C   TO BE IDENTIFIED QUICKLY.
C IWFR NEED NOT BE SET ON INPUT. ON OUTPUT IT POINTS TO THE FIRST
C   UNUSED LOCATION IN IW.
C ICNTL is an INTEGER array of assumed size.
C INFO is an INTEGER array of assumed size.
C
C   .. Scalar Arguments ..
    INTEGER KIND
    INTEGER IWFR,LW,N,NZ
C   ..
C   .. Array Arguments ..
    INTEGER FLAG(N),ICN(*),IPE(N),IQ(N),IRN(*),IW(LW)
    INTEGER ICNTL(*),INFO(*)
C   ..
C   .. Local Scalars ..
    INTEGER I,ID,J,JN,K,K1,K2,L,LAST,LDIAG,LP,LR,N1,NDUP
C   ..
C   .. Executable Statements ..
C
C INITIALIZE INFO(2) AND COUNT IN IPE THE
C   NUMBERS OF NON-ZEROS IN THE ROWS AND MOVE ROW AND COLUMN
C   NUMBERS INTO IW.
    INFO(2) = 0
    DO 10 I = 1,N
      IPE(I) = 0
   10 CONTINUE
    LP = ICNTL(1)
    LDIAG = ICNTL(5)
    LR = NZ
    IF (NZ.EQ.0) GO TO 120
    DO 110 K = 1,NZ
      I = IRN(K)
      J = ICN(K)
      IF (KIND.NE.1 .AND. I.GT.J) THEN
      INFO(1) = -19
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I8,A)')
   +    '**** Error return from ME57A ****  INFO(1) =',INFO(1),
   *    'Hermitian matrix but entry ',K,' in lower triangle'
      INFO(2) = K
      RETURN
      ENDIF
      IF (I.LT.J) THEN
      IF (I.GE.1 .AND. J.LE.N) GO TO 90
      ELSE IF (I.GT.J) THEN
      IF (J.GE.1 .AND. I.LE.N) GO TO 90
      ELSE
      IF (I.GE.1 .AND. I.LE.N) GO TO 80
      END IF
      INFO(2) = INFO(2) + 1
      INFO(1) = 1
      IF (INFO(2).LE.1 .AND. ICNTL(2).GT.0) THEN
      WRITE (ICNTL(2),FMT=60) INFO(1)
      END IF

   60   FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE ME57A',
   +      '  *** INFO(1) =',I2)

      IF (INFO(2).LE.10 .AND. ICNTL(2).GT.0) THEN
      WRITE (ICNTL(2),FMT=70) K,I,J
      END IF

   70   FORMAT (I6,'TH NON-ZERO (IN ROW',I6,' AND COLUMN',I6,
   +       ') IGNORED')

   80   I = 0
      J = 0
      GO TO 100

   90   IPE(I) = IPE(I) + 1
      IPE(J) = IPE(J) + 1
  100   IW(K) = J
      LR = LR + 1
      IW(LR) = I
  110 CONTINUE
C
C ACCUMULATE ROW COUNTS TO GET POINTERS TO ROW STARTS IN BOTH IPE AND IQ
C   AND INITIALIZE FLAG
  120 IQ(1) = 1
    N1 = N - 1
    IF (N1.LE.0) GO TO 140
    DO 130 I = 1,N1
      FLAG(I) = 0
      IF (IPE(I).EQ.0) IPE(I) = -1
      IQ(I+1) = IPE(I) + IQ(I) + 1
      IPE(I) = IQ(I)
  130 CONTINUE
  140 LAST = IPE(N) + IQ(N)
    FLAG(N) = 0
    IF (LR.GE.LAST) GO TO 160
    K1 = LR + 1
    DO 150 K = K1,LAST
      IW(K) = 0
  150 CONTINUE
  160 IPE(N) = IQ(N)
    IWFR = LAST + 1
    IF (NZ.EQ.0) GO TO 230
C
C RUN THROUGH PUTTING THE MATRIX ELEMENTS IN THE RIGHT PLACE
C   BUT WITH SIGNS INVERTED. IQ IS USED FOR HOLDING RUNNING POINTERS
C   AND IS LEFT HOLDING POINTERS TO ROW ENDS.
    DO 220 K = 1,NZ
      J = IW(K)
      IF (J.LE.0) GO TO 220
      L = K
      IW(K) = 0
      DO 210 ID = 1,NZ
      IF (L.GT.NZ) GO TO 170
      L = L + NZ
      GO TO 180

  170   L = L - NZ
  180   I = IW(L)
      IW(L) = 0
      IF (I.LT.J) GO TO 190
      L = IQ(J) + 1
      IQ(J) = L
      JN = IW(L)
      IW(L) = -I
      GO TO 200

  190   L = IQ(I) + 1
      IQ(I) = L
      JN = IW(L)
      IW(L) = -J
  200   J = JN
      IF (J.LE.0) GO TO 220
  210   CONTINUE
  220 CONTINUE
C
C RUN THROUGH RESTORING SIGNS, REMOVING DUPLICATES AND SETTING THE
C   MATE OF EACH NON-ZERO.
C NDUP COUNTS THE NUMBER OF DUPLICATE ELEMENTS.
  230 NDUP = 0
    DO 280 I = 1,N
      K1 = IPE(I) + 1
      K2 = IQ(I)
      IF (K1.LE.K2) GO TO 240
C ROW IS EMPTY. SET POINTER TO ZERO.
      IPE(I) = 0
      IQ(I) = 0
      GO TO 280
C ON ENTRY TO THIS LOOP FLAG(J).LT.I FOR J=1,2,...,N. DURING THE LOOP
C   FLAG(J) IS SET TO I IF A NON-ZERO IN COLUMN J IS FOUND. THIS
C   PERMITS DUPLICATES TO BE RECOGNIZED QUICKLY.
  240   DO 260 K = K1,K2
      J = -IW(K)
      IF (J.LE.0) GO TO 270
      L = IQ(J) + 1
      IQ(J) = L
      IW(L) = I
      IW(K) = J
      IF (FLAG(J).NE.I) GO TO 250
      NDUP = NDUP + 1
      IW(L) = 0
      IW(K) = 0
  250   FLAG(J) = I
  260   CONTINUE
  270   IQ(I) = IQ(I) - IPE(I)
      IF (NDUP.EQ.0) IW(K1-1) = IQ(I)
  280 CONTINUE
    IF (NDUP.EQ.0) GO TO 310
C
C COMPRESS IW TO REMOVE DUMMY ENTRIES CAUSED BY DUPLICATES.
    IWFR = 1
    DO 300 I = 1,N
      K1 = IPE(I) + 1
      IF (K1.EQ.1) GO TO 300
      K2 = IQ(I) + IPE(I)
      L = IWFR
      IPE(I) = IWFR
      IWFR = IWFR + 1
      DO 290 K = K1,K2
      IF (IW(K).EQ.0) GO TO 290
      IW(IWFR) = IW(K)
      IWFR = IWFR + 1
  290   CONTINUE
      IW(L) = IWFR - L - 1
  300 CONTINUE
  310 RETURN

    END
    SUBROUTINE ME57H(N,IPE,IW,LW,IWFR,NV,NXT,LST,IPD,FLAG,IOVFLO,
   +           NCMPA,FRATIO)
C Is identical to subroutine MA27HD.  Internal version for ME57.
C
C ANALYSIS SUBROUTINE
C
C GIVEN REPRESENTATION OF THE WHOLE MATRIX (EXCLUDING DIAGONAL)
C   PERFORM MINIMUM DEGREE ORDERING, CONSTRUCTING TREE POINTERS.
C   IT WORKS WITH SUPERVARIABLES WHICH ARE COLLECTIONS OF ONE OR MORE
C   VARIABLES, STARTING WITH SUPERVARIABLE I CONTAINING VARIABLE I FOR
C   I = 1,2,...,N. ALL VARIABLES IN A SUPERVARIABLE ARE ELIMINATED
C   TOGETHER. EACH SUPERVARIABLE HAS AS NUMERICAL NAME THAT OF ONE
C   OF ITS VARIABLES (ITS PRINCIPAL VARIABLE).
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C IPE(I) MUST BE SET TO POINT TO THE POSITION IN IW OF THE
C   START OF ROW I OR HAVE THE VALUE ZERO IF ROW I HAS NO OFF-
C   DIAGONAL NON-ZEROS. DURING EXECUTION IT IS USED AS FOLLOWS. IF
C   SUPERVARIABLE I IS ABSORBED INTO SUPERVARIABLE J THEN IPE(I)=-J.
C   IF SUPERVARIABLE I IS ELIMINATED THEN IPE(I) EITHER POINTS TO THE
C   LIST OF SUPERVARIABLES FOR CREATED ELEMENT I OR IS ZERO IF
C   THE CREATED ELEMENT IS NULL. IF ELEMENT I
C   IS ABSORBED INTO ELEMENT J THEN IPE(I)=-J.
C IW MUST BE SET ON ENTRY TO HOLD LISTS OF VARIABLES BY
C   ROWS, EACH LIST BEING HEADED BY ITS LENGTH.
C   DURING EXECUTION THESE LISTS ARE REVISED AND HOLD
C   LISTS OF ELEMENTS AND SUPERVARIABLES. THE ELEMENTS
C   ALWAYS HEAD THE LISTS. WHEN A SUPERVARIABLE
C   IS ELIMINATED ITS LIST IS REPLACED BY A LIST OF SUPERVARIABLES
C   IN THE NEW ELEMENT.
C LW MUST BE SET TO THE LENGTH OF IW. IT IS NOT ALTERED.
C IWFR MUST BE SET TO THE POSITION IN IW OF THE FIRST FREE VARIABLE.
C   IT IS REVISED DURING EXECUTION AND CONTINUES TO HAVE THIS MEANING.
C NV(JS) NEED NOT BE SET. DURING EXECUTION IT IS ZERO IF
C   JS IS NOT A PRINCIPAL VARIABLE AND IF IT IS IT HOLDS
C   THE NUMBER OF VARIABLES IN SUPERVARIABLE JS. FOR ELIMINATED
C   VARIABLES IT IS SET TO THE DEGREE AT THE TIME OF ELIMINATION.
C NXT(JS) NEED NOT BE SET. DURING EXECUTION IT IS THE NEXT
C   SUPERVARIABLE HAVING THE SAME DEGREE AS JS, OR ZERO
C   IF IT IS LAST IN ITS LIST.
C LST(JS) NEED NOT BE SET. DURING EXECUTION IT IS THE
C   LAST SUPERVARIABLE HAVING THE SAME DEGREE AS JS OR
C   -(ITS DEGREE) IF IT IS FIRST IN ITS LIST.
C IPD(ID) NEED NOT BE SET. DURING EXECUTION IT
C   IS THE FIRST SUPERVARIABLE WITH DEGREE ID OR ZERO
C   IF THERE ARE NONE.
C FLAG IS USED AS WORKSPACE FOR ELEMENT AND SUPERVARIABLE FLAGS.
C   WHILE THE CODE IS FINDING THE DEGREE OF SUPERVARIABLE IS
C   FLAG HAS THE FOLLOWING VALUES.
C   A) FOR THE CURRENT PIVOT/NEW ELEMENT ME
C       FLAG(ME)=-1
C   B) FOR VARIABLES JS
C       FLAG(JS)=-1 IF JS IS NOT A PRINCIPAL VARIABLE
C       FLAG(JS)=0 IF JS IS A SUPERVARIABLE IN THE NEW ELEMENT
C       FLAG(JS)=NFLG IF JS IS A SUPERVARIABLE NOT IN THE NEW
C           ELEMENT THAT HAS BEEN COUNTED IN THE DEGREE
C           CALCULATION
C       FLAG(JS).GT.NFLG IF JS IS A SUPERVARIABLE NOT IN THE NEW
C           ELEMENT THAT HAS NOT BEEN COUNTED IN THE DEGREE
C           CALCULATION
C   C) FOR ELEMENTS IE
C       FLAG(IE)=-1 IF ELEMENT IE HAS BEEN MERGED INTO ANOTHER
C       FLAG(IE)=-NFLG IF ELEMENT IE HAS BEEN USED IN THE DEGREE
C           CALCULATION FOR IS.
C       FLAG(IE).LT.-NFLG IF ELEMENT IE HAS NOT BEEN USED IN THE
C           DEGREE CALCULATION FOR IS
C IOVFLO should be set to a high legitimate integer.  It is used as a
C      flag.
C NCMPA number of compresses.
C FRATIO is set to ICNTL(14)/100 and is the density of rows regarded as
C   dense.
C
C   .. Scalar Arguments ..
    REAL FRATIO
    INTEGER IWFR,LW,N,IOVFLO,NCMPA
C   ..
C   .. Array Arguments ..
    INTEGER FLAG(N),IPD(N),IPE(N),IW(LW),LST(N),NV(N),NXT(N)
C   ..
C   .. Local Scalars ..
    INTEGER I,ID,IDL,IDN,IE,IP,IS,JP,JP1,JP2,JS,K,K1,K2,KE,KP,KP0,KP1,
   +      KP2,KS,L,LEN,LIMIT,LN,LS,LWFR,MD,ME,ML,MS,NEL,NFLG,NP,
   +      NP0,NS,NVPIV,NVROOT,ROOT
C LIMIT  Limit on number of variables for putting node in root.
C NVROOT Number of variables in the root node
C ROOT   Index of the root node (N+1 if none chosen yet).
C   ..
C   .. External Subroutines ..
    EXTERNAL ME57Z
C   ..
C   .. Intrinsic Functions ..
    INTRINSIC ABS,MIN
C   ..
C If a column of the reduced matrix has relative density greater than
C CNTL(2), it is forced into the root. All such columns are taken to
C have sparsity pattern equal to their merged patterns, so the fill
C and operation counts may be overestimated.
C
C IS,JS,KS,LS,MS,NS ARE USED TO REFER TO SUPERVARIABLES.
C IE,JE,KE ARE USED TO REFER TO ELEMENTS.
C IP,JP,KP,K,NP ARE USED TO POINT TO LISTS OF ELEMENTS.
C   OR SUPERVARIABLES.
C ID IS USED FOR THE DEGREE OF A SUPERVARIABLE.
C MD IS USED FOR THE CURRENT MINIMUM DEGREE.
C IDN IS USED FOR THE NO. OF VARIABLES IN A NEWLY CREATED ELEMENT
C NEL IS USED TO HOLD THE NO. OF VARIABLES THAT HAVE BEEN
C   ELIMINATED.
C ME=MS IS THE NAME OF THE SUPERVARIABLE ELIMINATED AND
C   OF THE ELEMENT CREATED IN THE MAIN LOOP.
C NFLG IS USED FOR THE CURRENT FLAG VALUE IN ARRAY FLAG. IT STARTS
C   WITH THE VALUE IOVFLO AND IS REDUCED BY 1 EACH TIME IT IS USED
C   UNTIL IT HAS THE VALUE 2 WHEN IT IS RESET TO THE VALUE IOVFLO.
C
C   .. Executable Statements ..
C INITIALIZATIONS
    DO 10 I = 1,N
      IPD(I) = 0
      NV(I) = 1
      FLAG(I) = IOVFLO
   10 CONTINUE
    MD = 1
    NCMPA = 0
    NFLG = IOVFLO
    NEL = 0
    ROOT = N+1
    NVROOT = 0
C
C LINK TOGETHER VARIABLES HAVING SAME DEGREE
    DO 30 IS = 1,N
      K = IPE(IS)
      IF (K.LE.0) GO TO 20
      ID = IW(K) + 1
      NS = IPD(ID)
      IF (NS.GT.0) LST(NS) = IS
      NXT(IS) = NS
      IPD(ID) = IS
      LST(IS) = -ID
      GO TO 30
C WE HAVE A VARIABLE THAT CAN BE ELIMINATED AT ONCE BECAUSE THERE IS
C   NO OFF-DIAGONAL NON-ZERO IN ITS ROW.
   20   NEL = NEL + 1
      FLAG(IS) = -1
      NXT(IS) = 0
      LST(IS) = 0
   30 CONTINUE
C
C START OF MAIN LOOP
C
    DO 340 ML = 1,N
C LEAVE LOOP IF ALL VARIABLES HAVE BEEN ELIMINATED.
      IF (NEL+NVROOT+1.GE.N) GO TO 350
C
C FIND NEXT SUPERVARIABLE FOR ELIMINATION.
      DO 40 ID = MD,N
      MS = IPD(ID)
      IF (MS.GT.0) GO TO 50
   40   CONTINUE
   50   MD = ID
C NVPIV HOLDS THE NUMBER OF VARIABLES IN THE PIVOT.
      NVPIV = NV(MS)
C
C REMOVE CHOSEN VARIABLE FROM LINKED LIST
      NS = NXT(MS)
      NXT(MS) = 0
      LST(MS) = 0
      IF (NS.GT.0) LST(NS) = -ID
      IPD(ID) = NS
      ME = MS
      NEL = NEL + NVPIV
C IDN HOLDS THE DEGREE OF THE NEW ELEMENT.
      IDN = 0
C
C RUN THROUGH THE LIST OF THE PIVOTAL SUPERVARIABLE, SETTING TREE
C   POINTERS AND CONSTRUCTING NEW LIST OF SUPERVARIABLES.
C KP IS A POINTER TO THE CURRENT POSITION IN THE OLD LIST.
      KP = IPE(ME)
      FLAG(MS) = -1
C IP POINTS TO THE START OF THE NEW LIST.
      IP = IWFR
C LEN HOLDS THE LENGTH OF THE LIST ASSOCIATED WITH THE PIVOT.
      LEN = IW(KP)
      DO 140 KP1 = 1,LEN
      KP = KP + 1
      KE = IW(KP)
C JUMP IF KE IS AN ELEMENT THAT HAS NOT BEEN MERGED INTO ANOTHER.
      IF (FLAG(KE).LE.-2) GO TO 60
C JUMP IF KE IS AN ELEMENT THAT HAS BEEN MERGED INTO ANOTHER OR IS
C   A SUPERVARIABLE THAT HAS BEEN ELIMINATED.
      IF (FLAG(KE).LE.0) THEN
         IF (IPE(KE).NE.-ROOT) GO TO 140
C KE has been merged into the root
         KE = ROOT
         IF (FLAG(KE).LE.0) GO TO 140
      END IF
C WE HAVE A SUPERVARIABLE. PREPARE TO SEARCH REST OF LIST.
      JP = KP - 1
      LN = LEN - KP1 + 1
      IE = MS
      GO TO 70
C SEARCH VARIABLE LIST OF ELEMENT KE, USING JP AS A POINTER TO IT.
   60   IE = KE
      JP = IPE(IE)
      LN = IW(JP)
C
C SEARCH FOR DIFFERENT SUPERVARIABLES AND ADD THEM TO THE NEW LIST,
C   COMPRESSING WHEN NECESSARY. THIS LOOP IS EXECUTED ONCE FOR
C   EACH ELEMENT IN THE LIST AND ONCE FOR ALL THE SUPERVARIABLES
C   IN THE LIST.
   70   DO 130 JP1 = 1,LN
        JP = JP + 1
        IS = IW(JP)
C JUMP IF IS IS NOT A PRINCIPAL VARIABLE OR HAS ALREADY BEEN COUNTED.
        IF (FLAG(IS).LE.0) THEN
         IF (IPE(IS).EQ.-ROOT) THEN
C IS has been merged into the root
            IS = ROOT
            IW(JP) = ROOT
            IF (FLAG(IS).LE.0) GO TO 130
         ELSE
            GO TO 130
         END IF
        END IF
        FLAG(IS) = 0
        IF (IWFR.LT.LW) GO TO 100
C PREPARE FOR COMPRESSING IW BY ADJUSTING POINTERS AND
C   LENGTHS SO THAT THE LISTS BEING SEARCHED IN THE INNER AND OUTER
C   LOOPS CONTAIN ONLY THE REMAINING ENTRIES.
        IPE(MS) = KP
        IW(KP) = LEN - KP1
        IPE(IE) = JP
        IW(JP) = LN - JP1
C COMPRESS IW
        CALL ME57Z(N,IPE,IW,IP-1,LWFR,NCMPA)
C COPY NEW LIST FORWARD
        JP2 = IWFR - 1
        IWFR = LWFR
        IF (IP.GT.JP2) GO TO 90
        DO 80 JP = IP,JP2
          IW(IWFR) = IW(JP)
          IWFR = IWFR + 1
   80     CONTINUE
C ADJUST POINTERS FOR THE NEW LIST AND THE LISTS BEING SEARCHED.
   90     IP = LWFR
        JP = IPE(IE)
        KP = IPE(ME)
C STORE IS IN NEW LIST.
  100     IW(IWFR) = IS
        IDN = IDN + NV(IS)
        IWFR = IWFR + 1
C REMOVE IS FROM DEGREE LINKED LIST
        LS = LST(IS)
        LST(IS) = 0
        NS = NXT(IS)
        NXT(IS) = 0
        IF (NS.GT.0) LST(NS) = LS
        IF (LS.LT.0) THEN
          LS = -LS
          IPD(LS) = NS
        ELSE IF (LS.GT.0) THEN
          NXT(LS) = NS
        END IF
  130   CONTINUE
C JUMP IF WE HAVE JUST BEEN SEARCHING THE VARIABLES AT THE END OF
C   THE LIST OF THE PIVOT.
      IF (IE.EQ.MS) GO TO 150
C SET TREE POINTER AND FLAG TO INDICATE ELEMENT IE IS ABSORBED INTO
C   NEW ELEMENT ME.
      IPE(IE) = -ME
      FLAG(IE) = -1
  140   CONTINUE

C STORE THE DEGREE OF THE PIVOT.
  150   NV(MS) = IDN + NVPIV
C JUMP IF NEW ELEMENT IS NULL.
      IF (IWFR.EQ.IP) GO TO 330
      K1 = IP
      K2 = IWFR - 1
C
C RUN THROUGH NEW LIST OF SUPERVARIABLES REVISING EACH ASSOCIATED LIST,
C   RECALCULATING DEGREES AND REMOVING DUPLICATES.
      LIMIT = NINT(FRATIO*(N-NEL))
      DO 310 K = K1,K2
      IS = IW(K)
      IF (IS.EQ.ROOT) GO TO 310
      IF (NFLG.GT.2) GO TO 170
C RESET FLAG VALUES TO +/-IOVFLO.
      DO 160 I = 1,N
        IF (FLAG(I).GT.0) FLAG(I) = IOVFLO
        IF (FLAG(I).LE.-2) FLAG(I) = -IOVFLO
  160   CONTINUE
      NFLG = IOVFLO
C REDUCE NFLG BY ONE TO CATER FOR THIS SUPERVARIABLE.
  170   NFLG = NFLG - 1
C BEGIN WITH THE DEGREE OF THE NEW ELEMENT. ITS VARIABLES MUST ALWAYS
C   BE COUNTED DURING THE DEGREE CALCULATION AND THEY ARE ALREADY
C   FLAGGED WITH THE VALUE 0.
      ID = IDN
C RUN THROUGH THE LIST ASSOCIATED WITH SUPERVARIABLE IS
      KP1 = IPE(IS) + 1
C NP POINTS TO THE NEXT ENTRY IN THE REVISED LIST.
      NP = KP1
      KP2 = IW(KP1-1) + KP1 - 1
      DO 220 KP = KP1,KP2
        KE = IW(KP)
C TEST WHETHER KE IS AN ELEMENT, A REDUNDANT ENTRY OR A SUPERVARIABLE.
      IF (FLAG(KE).EQ.-1) THEN
         IF (IPE(KE).NE.-ROOT) GO TO 220
C KE has been merged into the root
         KE = ROOT
         IW(KP) = ROOT
         IF (FLAG(KE).EQ.-1) GO TO 220
      END IF
      IF (FLAG(KE).GE.0) GO TO 230
C SEARCH LIST OF ELEMENT KE, REVISING THE DEGREE WHEN NEW VARIABLES
C   FOUND.
        JP1 = IPE(KE) + 1
        JP2 = IW(JP1-1) + JP1 - 1
        IDL = ID
        DO 190 JP = JP1,JP2
          JS = IW(JP)
C JUMP IF JS HAS ALREADY BEEN COUNTED.
          IF (FLAG(JS).LE.NFLG) GO TO 190
          ID = ID + NV(JS)
          FLAG(JS) = NFLG
  190     CONTINUE
C JUMP IF ONE OR MORE NEW SUPERVARIABLES WERE FOUND.
        IF (ID.GT.IDL) GO TO 210
C CHECK WHETHER EVERY VARIABLE OF ELEMENT KE IS IN NEW ELEMENT ME.
        DO 200 JP = JP1,JP2
          JS = IW(JP)
          IF (FLAG(JS).NE.0) GO TO 210
  200     CONTINUE
C SET TREE POINTER AND FLAG TO INDICATE THAT ELEMENT KE IS ABSORBED
C   INTO NEW ELEMENT ME.
        IPE(KE) = -ME
        FLAG(KE) = -1
        GO TO 220
C STORE ELEMENT KE IN THE REVISED LIST FOR SUPERVARIABLE IS AND FLAG IT.
  210     IW(NP) = KE
        FLAG(KE) = -NFLG
        NP = NP + 1
  220   CONTINUE
      NP0 = NP
      GO TO 250
C TREAT THE REST OF THE LIST ASSOCIATED WITH SUPERVARIABLE IS. IT
C   CONSISTS ENTIRELY OF SUPERVARIABLES.
  230   KP0 = KP
      NP0 = NP
      DO 240 KP = KP0,KP2
        KS = IW(KP)
        IF (FLAG(KS).LE.NFLG) THEN
         IF (IPE(KS).EQ.-ROOT) THEN
            KS = ROOT
            IW(KP) = ROOT
            IF (FLAG(KS).LE.NFLG) GO TO 240
         ELSE
            GO TO 240
         END IF
        END IF
C ADD TO DEGREE, FLAG SUPERVARIABLE KS AND ADD IT TO NEW LIST.
        ID = ID + NV(KS)
        FLAG(KS) = NFLG
        IW(NP) = KS
        NP = NP + 1
  240   CONTINUE
C MOVE FIRST SUPERVARIABLE TO END OF LIST, MOVE FIRST ELEMENT TO END
C   OF ELEMENT PART OF LIST AND ADD NEW ELEMENT TO FRONT OF LIST.
  250   IF (ID.GE.LIMIT) GO TO 295
      IW(NP) = IW(NP0)
      IW(NP0) = IW(KP1)
      IW(KP1) = ME
C STORE THE NEW LENGTH OF THE LIST.
      IW(KP1-1) = NP - KP1 + 1
C
C CHECK WHETHER ROW IS IS IDENTICAL TO ANOTHER BY LOOKING IN LINKED
C   LIST OF SUPERVARIABLES WITH DEGREE ID AT THOSE WHOSE LISTS HAVE
C   FIRST ENTRY ME. NOTE THAT THOSE CONTAINING ME COME FIRST SO THE
C   SEARCH CAN BE TERMINATED WHEN A LIST NOT STARTING WITH ME IS
C   FOUND.
      JS = IPD(ID)
      DO 280 L = 1,N
        IF (JS.LE.0) GO TO 300
        KP1 = IPE(JS) + 1
        IF (IW(KP1).NE.ME) GO TO 300
C JS HAS SAME DEGREE AND IS ACTIVE. CHECK IF IDENTICAL TO IS.
        KP2 = KP1 - 1 + IW(KP1-1)
        DO 260 KP = KP1,KP2
          IE = IW(KP)
C JUMP IF IE IS A SUPERVARIABLE OR AN ELEMENT NOT IN THE LIST OF IS.
          IF (ABS(FLAG(IE)+0).GT.NFLG) GO TO 270
  260     CONTINUE
        GO TO 290

  270     JS = NXT(JS)
  280   CONTINUE
C SUPERVARIABLE AMALGAMATION. ROW IS IS IDENTICAL TO ROW JS.
C REGARD ALL VARIABLES IN THE TWO SUPERVARIABLES AS BEING IN IS. SET
C   TREE POINTER, FLAG AND NV ENTRIES.
  290   IPE(JS) = -IS
      NV(IS) = NV(IS) + NV(JS)
      NV(JS) = 0
      FLAG(JS) = -1
C REPLACE JS BY IS IN LINKED LIST.
      NS = NXT(JS)
      LS = LST(JS)
      IF (NS.GT.0) LST(NS) = IS
      IF (LS.GT.0) NXT(LS) = IS
      LST(IS) = LS
      NXT(IS) = NS
      LST(JS) = 0
      NXT(JS) = 0
      IF (IPD(ID).EQ.JS) IPD(ID) = IS
      GO TO 310
C Treat IS as full. Merge it into the root node.
  295   IF (NVROOT.EQ.0) THEN
        ROOT = IS
        IPE(IS) = 0
      ELSE
        IW(K) = ROOT
        IPE(IS) = -ROOT
        NV(ROOT) = NV(ROOT) + NV(IS)
        NV(IS) = 0
        FLAG(IS) = -1
      END IF
      NVROOT = NV(ROOT)
      GO TO 310
C INSERT IS INTO LINKED LIST OF SUPERVARIABLES OF SAME DEGREE.
  300   NS = IPD(ID)
      IF (NS.GT.0) LST(NS) = IS
      NXT(IS) = NS
      IPD(ID) = IS
      LST(IS) = -ID
      MD = MIN(MD,ID)
  310   CONTINUE
C
C RESET FLAGS FOR SUPERVARIABLES IN NEWLY CREATED ELEMENT AND
C   REMOVE THOSE ABSORBED INTO OTHERS.
      DO 320 K = K1,K2
      IS = IW(K)
      IF (NV(IS).EQ.0) GO TO 320
      FLAG(IS) = NFLG
      IW(IP) = IS
      IP = IP + 1
  320   CONTINUE
      IWFR = K1
      FLAG(ME) = -NFLG
C MOVE FIRST ENTRY TO END TO MAKE ROOM FOR LENGTH.
      IW(IP) = IW(K1)
      IW(K1) = IP - K1
C SET POINTER FOR NEW ELEMENT AND RESET IWFR.
      IPE(ME) = K1
      IWFR = IP + 1
      GO TO 335

  330   IPE(ME) = 0
C
  335   CONTINUE
  340 CONTINUE
C

C Absorb any remaining variables into the root
  350 DO 360 IS = 1,N
      IF(NXT(IS).NE.0 .OR. LST(IS).NE.0) THEN
      IF (NVROOT.EQ.0) THEN
        ROOT = IS
        IPE(IS) = 0
      ELSE
        IPE(IS) = -ROOT
      END IF
      NVROOT = NVROOT + NV(IS)
      NV(IS) = 0
       END IF
  360 CONTINUE
C Link any remaining elements to the root
    DO 370 IE = 1,N
      IF (IPE(IE).GT.0) IPE(IE) = -ROOT
  370 CONTINUE
    IF(NVROOT.GT.0)NV(ROOT)=NVROOT
    END
    SUBROUTINE ME57Z(N,IPE,IW,LW,IWFR,NCMPA)
C Is identical to subroutine MA27UD.  Internal version for ME57.
C COMPRESS LISTS HELD BY MA27H/HD (ME57H/HD) IN IW AND ADJUST POINTERS
C   IN IPE TO CORRESPOND.
C N IS THE MATRIX ORDER. IT IS NOT ALTERED.
C IPE(I) POINTS TO THE POSITION IN IW OF THE START OF LIST I OR IS
C   ZERO IF THERE IS NO LIST I. ON EXIT IT POINTS TO THE NEW POSITION.
C IW HOLDS THE LISTS, EACH HEADED BY ITS LENGTH. ON OUTPUT THE SAME
C   LISTS ARE HELD, BUT THEY ARE NOW COMPRESSED TOGETHER.
C LW HOLDS THE LENGTH OF IW. IT IS NOT ALTERED.
C IWFR NEED NOT BE SET ON ENTRY. ON EXIT IT POINTS TO THE FIRST FREE
C   LOCATION IN IW.
C   ON RETURN IT IS SET TO THE FIRST FREE LOCATION IN IW.
C NCMPA is number of compresses.
C
C   .. Scalar Arguments ..
    INTEGER IWFR,LW,N,NCMPA
C   ..
C   .. Array Arguments ..
    INTEGER IPE(N),IW(LW)
C   ..
C   .. Local Scalars ..
    INTEGER I,IR,K,K1,K2,LWFR
C   ..
C   .. Executable Statements ..
    NCMPA = NCMPA + 1
C PREPARE FOR COMPRESSING BY STORING THE LENGTHS OF THE
C   LISTS IN IPE AND SETTING THE FIRST ENTRY OF EACH LIST TO
C   -(LIST NUMBER).
    DO 10 I = 1,N
      K1 = IPE(I)
      IF (K1.LE.0) GO TO 10
      IPE(I) = IW(K1)
      IW(K1) = -I
   10 CONTINUE
C
C COMPRESS
C IWFR POINTS JUST BEYOND THE END OF THE COMPRESSED FILE.
C LWFR POINTS JUST BEYOND THE END OF THE UNCOMPRESSED FILE.
    IWFR = 1
    LWFR = IWFR
    DO 60 IR = 1,N
      IF (LWFR.GT.LW) GO TO 70
C SEARCH FOR THE NEXT NEGATIVE ENTRY.
      DO 20 K = LWFR,LW
      IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
      GO TO 70
C PICK UP ENTRY NUMBER, STORE LENGTH IN NEW POSITION, SET NEW POINTER
C   AND PREPARE TO COPY LIST.
   30   I = -IW(K)
      IW(IWFR) = IPE(I)
      IPE(I) = IWFR
      K1 = K + 1
      K2 = K + IW(IWFR)
      IWFR = IWFR + 1
      IF (K1.GT.K2) GO TO 50
C COPY LIST TO NEW POSITION.
      DO 40 K = K1,K2
      IW(IWFR) = IW(K)
      IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN

    END
C COPYRIGHT (c) 1999 CCLRC Council for the Central Laboratory
*            of the Research Councils
C
C Version 2.3.1
C For full version history see ChangeLog
C

    SUBROUTINE ME57ID(CNTL, ICNTL)
C
C****************************************************************
C
C  Purpose
C  =======
C
C  The entries of the arrays CNTL and ICNTL control the action of
C  ME57. Default values for the entries are set in this routine.
C
    DOUBLE PRECISION    CNTL(5)
    INTEGER         ICNTL(20)
C
C  Parameters
C  ==========
C
C  The entries of the arrays CNTL and ICNTL control the action of
C   ME57. Default values are set by ME57ID.
C
C   CNTL(1) has default value 0.01 and is used for threshold pivoting
C    Values less than 0.0 will be treated as 0.0 and values greater
C    than 0.5 as 0.5. Values near 0.0 may perhaps give faster
C    factorization times and less entries in the factors but may
C    result in a less stable factorization.  This parameter is only
C    accessed if ICNTL(7) is equal to 1.
C
C   CNTL(2) has default value $10^{-20}$. ME57B/BD will treat any
C    pivot whose modulus is less than CNTL(2) as zero.
C    If ICNTL(16) = 1, then blocks of entries less than CNTL(2) can
C    be discarded during the factorization and the corresponding
C    pivots are placed at the end of the ordering. In this case,
C    a normal value for CNTL(2) could be $10{-12}$.
C
C   CNTL(3) has default value 0.5.  It is used by ME57D/DD to
C    monitor the convergence of the iterative refinement.
C    If the norm of the scaled residuals does not decrease by a
C    factor of at least CNTL(3), convergence is deemed to be too
C    slow and ME57D/DD terminates with INFO(1) set to -8.
C
C   CNTL(4) has default value 0.0. It is used by ME57B/BD to control
C    the static pivoting. If CNTL(4) > 0.0, then
C    small pivots may be replaced by entries of value CNTL(4) so that
C    the factorization could be inaccurate.  If CNTL(5) is also
C    greater than zero, then CNTL(4) is treated as zero (that is
C    uneliminated variables are delayed) until CNTL(5)*N fully summed
C    variables have been delayed.  If static pivots are used, it is
C    recommended that iterative refinement is used when computing the
C    solution.
C
C   CNTL(5)  has default value zero. Static pivoting is invoked if
C    CNTL(4) is greater than zero and the accumulated number of
C    delayed pivots exceeds CNTL(5)*N.
C
C  ICNTL(1) has default value 6.
C   It is the output stream for error messages. If it is negative,
C   these messages are suppressed.
C
C  ICNTL(2) has default value 6.
C   It is the output stream for warning messages. If it is negative,
C   these messages are suppressed.
C
C  ICNTL(3) has default value 6. It is the output stream for monitoring
C   printing. If it is negative, these messages are suppressed.
C
C  ICNTL(4) has default value -1.
C   It is the output stream for the printing of statistics.
C   If it is negative, the statistics are not printed.
C
C  ICNTL(5) is used by ME57 to control printing of error,
C   warning, and monitoring messages. It has default value 2.
C   Possible values are:
C
C    <1     No messages output.
C   1     Only error messages printed.
C   2     Errors and warnings printed.
C   3     Errors and warnings and terse monitoring
C         (only first ten entries of arrays printed).
C    >3     Errors and warnings and all information
C         on input and output parameters printed.
C
C  ICNTL(6) has default value 5 and must be set by the user to 1
C   if the pivot order in KEEP is to be used by ME57AD. For any
C   other value of ICNTL(6), a suitable pivot order will be chosen
C   automatically.  The choices available so far are:
C  ICNTL(6) = 0   AMD using MC47 (with dense row detection disabled)
C  ICNTL(6) = 2   AMD using MC47 (this was previously the MC50 option)
C  ICNTL(6) = 3   MA27 minimum degree ordering
C  ICNTL(6) = 4   METIS_NODEND ordering from MeTiS package.
C  ICNTL(6) = 5   Ordering chosen depending on matrix characteristics.
C             At the moment choices are MC47 or METIS.
C             INFO(36) is set to ordering used.
C  ICNTL(6) > 5   At the moment this is treated as 5 (the default).
C
C  ICNTL(7) is used by ME57BD to control numerical pivoting.  It
C   has default value 1.  Values out of range cause an error return
C   with INFO(1) equal to -10.
C   Possible values are:
C
C  1   Numerical pivoting is performed using the threshold value in
C    CNTL(1).
C
C  2   No pivoting will be performed and an error exit will occur
C    immediately a sign change or a zero is detected among the pivots.
C    This is suitable for cases when A is thought to be definite
C    and is likely to decrease the factorization time while still
C    providing a stable decomposition.
C
C  3   No pivoting will be performed and an error exit will occur if a
C    zero pivot is detected. This is likely to decrease the
C    factorization time, but may be unstable if there is a sign
C    change among the pivots.
C
C  4   No pivoting will be performed but the matrix will be altered
C    so that all pivots are of the same sign.
C
C  ICNTL(8) has default value 0. If ME57BD is called with ICNTL(8) NE 0,
C    then the factorization will discard factors and try
C    to continue the factorization to determine the amount of space
C    needed for a successful factorization.  In this case, a
C    factorization will not have been produced.  If the default value
C    of 0 is used and the factorization stops because of lack of
C    space, the user should reallocate the real or integer space for
C    FACT or IFACT, respectively and reset LFACT or LIFACT
C    appropriately, using ME57ED before recalling ME57BD.
C
C  ICNTL(9) has default value 10.  It corresponds to the maximum number
C    of steps of iterative refinement.
C
C  ICNTL(10) has default value 0. A positive value will return the
c    infinity norm of the input matrix, the computed solution, and
C    the scaled residual in RINFO(5) to RINFO(7), respectively,
C    a backward error estimate in RINFO(8) and RINFO(9), and an
C    estimate of the forward error in RINFO(10).  If ICNTL(10) is
C    negative or zero no estimates are returned.
C
C  ICNTL(11) The block size to be used by the Level 3 BLAS.
C
C  ICNTL(12) Two nodes of the assembly tree are merged only if both
C    involve less than ICNTL(12) eliminations.
C
C  ICNTL(13) Threshold on number of rows in a block for using BLAS2 in
C    ME57CD.
C
C  ICNTL(14) Threshold on number of entries in a row for declaring row
C    full. Set as percentage of N.  So 100 means row must be full.
C
C  ICNTL(15) should be set to 1 (the default) if MC64 scaling is
C   requested.
C
C  ICNTL(16) should be set to 1 if "small" entries are to be removed
C   removed from the frontal matrices.  The default is 0.
C
C  ICNTL(17) to ICNTL(20) are set to zero by ME57ID but are not
C   currently used by ME57.
C
C Local variables
    INTEGER I
    DOUBLE PRECISION ZERO
    PARAMETER (ZERO=0.0D0)
C===============================================
C Default values for variables in control arrays
C===============================================
C   Threshold value for pivoting
    CNTL(1)   = 0.01D0
C Test for zero pivot
    CNTL(2)   = 1.0D-20
C Iterative refinement convergence
    CNTL(3)   = 0.5D0
C Static pivoting control
    CNTL(4) = ZERO
C Control to allow some delayed pivots
    CNTL(5) = ZERO
C   Printing control
    ICNTL(1)  = 6
    ICNTL(2)  = 6
    ICNTL(3)  = 6
    ICNTL(4)  = -1
    ICNTL(5)  = 2
C   Provide pivot order (1=NO)
C   Set to make automatic choice between METIS and MC47
    ICNTL(6)  = 5
C   Pivoting control
    ICNTL(7)  = 1
C   Restart facility
    ICNTL(8)  = 0
C   IR steps
    ICNTL(9)  = 10
C   Error estimates
    ICNTL(10) = 0
C   Blocking for Level 3 BLAS
    ICNTL(11) = 16
C   Node amalgamation parameter (NEMIN)
    ICNTL(12) = 16
C   Switch for use of Level 2 BLAS in solve
    ICNTL(13) = 10
C Flag to indicate threshold will be set to N
    ICNTL(14) = 100
C Flag to invoke MC64 scaling (0 off, 1 on)
    ICNTL(15) = 1
C Flag to invoke dropping small entries from front
C Default is not to drop (set to 1 to drop)
    ICNTL(16) = 0

C Set unused parameters
    DO 110 I=17,20
      ICNTL(I) = 0
  110 CONTINUE

    RETURN
    END


    SUBROUTINE ME57AD(KIND,N,NE,IRN,JCN,LKEEP,KEEP,IWORK,
   *            ICNTL,INFO,RINFO)
C This subroutine is a user-callable driver for the analysis phase of
C   ME57. It performs an ordering, a symbolic factorization and
C   computes information for the numerical factorization by ME57B/BD.
    INTEGER KIND,N,NE,IRN(NE),JCN(NE),IWORK(5*N),LKEEP,KEEP(LKEEP),
   *      ICNTL(20),INFO(40)
    DOUBLE PRECISION RINFO(20)
C
C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C
C  N is an INTEGER variable which must be set by the user to the
C    order n of the matrix A.  It is not altered by the subroutine.
C    Restriction: N > 0.
C
C  NE is an INTEGER variable which must be set by the user to the
C    number of entries being input.  It is not altered by the
C    subroutine. Restriction: NE >= 0.
C
C  IRN and JCN are INTEGER  arrays of length NE. The user
C    must set them so that each off-diagonal nonzero $a_{ij}$ is
C    represented by IRN(k)=i and JCN(k)=j or by IRN(k)=j
C    and JCN(k)=i.  Multiple entries are allowed and any with
C    IRN(k) or JCN(k) out of range are ignored. These arrays will
C    be unaltered by the subroutine.
C
C  IWORK is an INTEGER  array of length 5*N. This need not be set
C    by the user and is used as workspace by ME57AD.
C
C  LKEEP is an INTEGER that must be set to length of array KEEP.
C    Restriction: LKEEP >= 5*N+NE+MAX(N,NE)+42
C
C  KEEP is an INTEGER  array of length LKEEP. It need not be set
C    by the user and must be preserved between a call to ME57AD
C    and subsequent calls to ME57BD and ME57CD.
C    If the user wishes to input
C    the pivot sequence, the position of variable {i} in the pivot
C    order should be placed in KEEP(i), i = 1, 2,..., n and ICNTL(6)
C    should be set to 1.  The subroutine may replace the given order
C    by another that gives the same fill-in pattern and
C    virtually identical numerical results.
C
C ICNTL is an INTEGER array of length 10
C    that contains control parameters and must be set by the user.
C    Default values for the components may be set by a call to
C    ME57I/ID. Details of the control parameters are given in ME57I/ID.
C
C INFO is an INTEGER array of length 40 that need not be set by the
C    user.  On return from ME57AD, a value of zero for INFO(1)
C    indicates that the subroutine has performed successfully.

C RINFO is a REAL (DOUBLE_PRECISION in the D version) array of length 20
C    that need not be set by the
C    user.  This array supplies information on the execution of ME57AD.
C
C**** Still to be updated
C    INFO(1):
C     0  Successful entry.
C    -1 N < 1
C    -2 NE < 0.
C    -9 Invalid permutation supplied in KEEP.
C   -15 LKEEP < 5*N+NE+MAX(N,NE)+42
C    +1 One or more indices out of range.
C    +2 One or more duplicate entries found.
C    +3. Combination of warnings +1 and +2.
C    INFO(2):
C      if INFO(1) = -1, the value input for N.
C      if INFO(1) = -2, the value input for NE.
C      if INFO(1) = -9, index at which error first detected.
C      if INFO(1) = -15, the value input for LKEEP.
C    INFO(3) Number of entries with out-of-range indices.
C    INFO(4) Number of off-diagonal duplicate entries.
C    INFO(5) Forecast number of reals to hold the factorization.
C    INFO(6) Forecast number of integers to hold the factorization.
C    INFO(7) Forecast maximum front size.
C    INFO(8) Number of nodes in the assembly tree.
C    INFO(9) Minimum size for LA of ME57BD (without compress).
C    INFO(10) Minimum size for LIFACT of ME57BD (without compress).
C    INFO(11) Minimum size for LA of ME57BD (with compress).
C    INFO(12) Minimum size for LIFACT of ME57BD (with compress).
C    INFO(13) Number of compresses.
C    INFO(14:40) Not used.

C Procedures
    INTRINSIC MIN
    EXTERNAL ME57GD,MC47ID,MC47BD,ME57VD,ME57HD,ME57JD,ME57KD,
   *       ME57LD,ME57MD,ME57ND
C ME57GD Expand representation to whole matrix and sort.
C MC47BD Approximate Minimum Degree avoiding problems with dense rows.
C ME57VD Is same as MA27GD. Sort for using MA27HD/ME57HD.
C ME57HD Is same as MA27HD. Minimum degree ordering from MA27.
C ME57JD Sort to upper triangular form (pivot sequence given).
C ME57KD Construct tree pointers given output from ME57JD.
C ME57LD Depth-first search of tree.
C ME57MD Construct map.
C ME57ND Calculate storage and operation counts.

C Local variables
    INTEGER I,IL,IN,IPE,IRNPRM,COUNT,FILS,FRERE,HOLD,IFCT,INVP,IPS,
   +      IW,IWFR,K,LDIAG,LP,LW,LROW,MAP,EXPNE,
   +      MP,NCMPA,NEMIN,NODE,NST,NSTEPS,NV,PERM,
   +      IW1,IW2,IW3,IW4,IW5,NSTK,ND,NELIM,NZE,ALENB,
   +      J,JJ,J1,J2,SIZE22,OXO
C Local variables for MA27 ordering call
    INTEGER IF27H
C Local variables for MeTiS call
    INTEGER METOPT(8),METFTN,ICNTL6,INF47(10),ICNT47(10)
    DOUBLE PRECISION ZERO,THRESH,AVNUM,MC47FI,RINF47(10)
    PARAMETER (ZERO=0.0D0)

C I      Temporary DO index.
C IPE    Displacement in array KEEP for array IPE of ME57GD/HD/JD/KD.
C COUNT    Displacement in array KEEP for array COUNT of
C      ME57GD/HD/JD/KD/LD/ND.
C IWFR   First unused location in IFCT(1:LW).
C K      Temporary variable.
C LDIAG    Control for amount of information output.
C LP     Stream number for error messages.
C LW     Length of IFACT when four arrays of length N+1 are excluded.
C LROW   Subscript in array KEEP for array LAST of ME57KD/LD
C      and LROW of ME57MD.
C MAP    Subscript in array KEEP for MAP array.
C MP     Stream number for monitoring.
C NODE   Subscript in array KEEP for array FLAG of ME57KD, and
C      NODE of ME57LD.
C NSTEPS   Number of nodes in the assembly tree.
C NV     Displacement in array IFACT for array NV of ME57GD/HD, PERM
C      of ME57J/JD, IPR of ME57K/KD, and NE of ME57L\LD.
C PERM   Subscript in array KEEP for array PERM of ME57LD/ND.
C SP     Stream number for statistics.
C SIZES    subscript in array KEEP.  On exit,
C      KEEP(SIZES) is number of faulty entries.
C      KEEP(SIZES+1) is the number of nodes in the assembly tree.


C Set local print variables
    LP = ICNTL(1)
    MP = ICNTL(3)
    LDIAG = ICNTL(5)

C Initialize information array.
    DO 10 I = 1,40
      INFO(I) = 0
   10 CONTINUE
    DO 11 I = 1,20
      RINFO(I) = ZERO
   11 CONTINUE

C Check N, NE, and LKEEP for obvious errors
    IF (N.LT.1)  GO TO 20
    IF (NE.LT.0) GO TO 30
    IF (LKEEP.LT.5*N+NE+MAX(N,NE)+42) GO TO 40

    IF (ICNTL(6).EQ.1) THEN
C Check permutation array
      DO 12 I = 1,N
      IWORK(I) = 0
   12   CONTINUE
      DO 14 I=1,N
      K = KEEP(I)
      IF (K.LE.0 .OR. K.GT.N) GO TO 80
      IF (IWORK(K).NE.0) GO TO 80
      IWORK(K) = I
   14   CONTINUE
    ENDIF

C If requested, print input variables.
    IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
      WRITE(MP,'(/A)') 'Entering analysis phase (ME57AD) with ...'
      IF (KIND.EQ.1) THEN
      WRITE(MP,'(A)') 'Matrix is complex symmetric'
      ELSE
      WRITE(MP,'(A)') 'Matrix is Hermitian'
      ENDIF
      WRITE(MP,99980) N,NE,(ICNTL(I),I=1,7),ICNTL(12),ICNTL(15)
99980 FORMAT('N       Order of matrix             =',I12/
   2    'NE      Number of entries             =',I12/
   6    'ICNTL(1)  Stream for errors             =',I12/
   7    ' --- (2)  Stream for warnings           =',I12/
   8    ' --- (3)  Stream for monitoring         =',I12/
   9    ' --- (4)  Stream for statistics         =',I12/
   1    ' --- (5)  Level of diagnostic printing      =',I12/
   2    ' --- (6)  Flag for input pivot order      =',I12/
   2    ' --- (7)  Numerical pivoting control (st est) =',I12/
   2    ' --- (12) Node amalgamation parameter       =',I12/
   2    ' --- (15) Scaling control (storage estimate)  =',I12)
      K = MIN(10,NE)
      IF (LDIAG.GE.4) K = NE
C      WRITE (MP,'(/A/(3(I6,A,2I8,A)))') ' Matrix entries:',
C   +      (I,': (',IRN(I),JCN(I),')',I=1,K)
      IF (K.LT.NE) WRITE (MP,'(A)') '   . . .'

      IF (ICNTL(6).EQ.1) THEN
C Print out permutation array.
      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(A,10I6:/(7X,10I6))') ' KEEP =', (KEEP(I),I=1,K)
      IF (K.LT.N) WRITE (MP,'(7X,A)') '   . . .'
      END IF

    END IF

C Partition IWORK
    IW1 = 1
    IW2 = IW1 + N
    IW3 = IW2 + N
    IW4 = IW3 + N
    IW5 = IW4 + N
    FILS  = IW1
    FRERE = IW2
    ND    = IW3
    NELIM = IW4
    NV    = IW5

C Partition KEEP
    PERM = 1
    NSTEPS = PERM + N
    EXPNE  = NSTEPS + 1
    HOLD   = EXPNE + 1
    LROW = HOLD + 40
    NODE = LROW + N
    NSTK = NODE + N
    MAP  = NSTK + N
    IRNPRM = MAP + MAX(N,NE)
    INVP  = NODE
    IW    = NODE
    IPE   = LROW
    IFCT  = MAP
C This is set for MA57VD/MA57HD (MA27 orderings) only to allow more
C   space for expanded factors
    IF27H = NODE
    IPS   = MAP
    COUNT = NSTK

C Set HOLD(1)
    KEEP(HOLD) = 0

C Sort and order ... generate tree

C Set local value for ICNTL(6)
    ICNTL6 = ICNTL(6)
    IF (ICNTL(6).GT.5) ICNTL6 = 5

    IF (ICNTL6.EQ.4 .OR. ICNTL6.EQ.5) THEN
C MeTiS ordering requested.  Use dummy call to see if it has been
C   installed.
C Set flag for Fortran-style numbering of arrays
      METFTN    = 1
C Set default values for parameters.
      METOPT(1) = 0
C Dummy call with 1 by 1 matrix
      KEEP(IPE)   = 1
      KEEP(IPE+1) = 2
      KEEP(IFCT)  = 1
      CALL METIS_NODEND_ME57(1,KEEP(IPE),KEEP(IFCT),METFTN,METOPT,
   *            KEEP(NSTK),KEEP(PERM))
C Flag set if dummy code for METIS_NODEND has been used
      IF (KEEP(PERM).EQ.-1) THEN
      IF (ICNTL6 .EQ. 4) GO TO 90
C Reset ICNTL6 to use MC47
      ICNTL6 = 2
      ENDIF
    ENDIF
    IF (ICNTL6.NE.1) THEN
C Ordering is to be calculated by program

      CALL MC47ID(ICNT47)

      IF (ICNTL6 .NE. 3) THEN
C ELSE clause (ICNTL6.EQ.3) for MA27 ordering
C MC47 (with dense row detection disabled) used if ICNTL6 equal to 0.
C MC47 used if ICNTL6 equal to 2.
C METIS used if ICNTL6 equal to 4.
C Automatic choice of METIS or MC47 if ICNTL6 equal to 5.

C Sort matrix to obtain complete pattern (upper and lower triangle)
C   but omitting diagonals, duplicates, and out-of-range entries.
C   On exit, sorted matrix is given by IPE(pointers), IFCT (indices),
C   and COUNT (lengths).  IWFR is position after last index in IFCT.
      CALL ME57GD(KIND,N,NE,IRN,JCN,KEEP(IFCT),KEEP(IPE),
   +          KEEP(COUNT),KEEP(IW),IWFR,ICNTL,INFO)
      IF (INFO(1).LT.0) RETURN

      IF (ICNTL6.EQ.5) THEN
C Calculate matrix statistics to determine ordering
        IF (ICNTL(7).EQ.2) THEN
C Action if positive definite option has been chosen.
          AVNUM = FLOAT(IWFR+N-1)/FLOAT(N)
          IF (N.GE.50000) THEN
          ICNTL6 = 4
C         IF (AVNUM.LE.6.0) ICNTL6 = 2
          GO TO 97
          ENDIF
          IF (N.LE.30000) THEN
          ICNTL6 = 2
          IF (AVNUM.GT.100.0) ICNTL6 = 4
          GO TO 97
          ENDIF
          IF (N.GT.30000 .AND. N.LT.50000) THEN
          IF (AVNUM.GT.46.0) THEN
            ICNTL6 = 4
          ELSE
            ICNTL6 = 2
          ENDIF
          GO TO 97
          ENDIF
        ELSE
C Matrix has not been declared positive definite.
          AVNUM = FLOAT(IWFR+N-1)/FLOAT(N)
C Set flag for detection of OXO matrix
          OXO = 0
C Check for KKT and OXO.
C Calculate size of possible  trailing block
          J2 = IWFR - 1
          SIZE22 = 0
          DO 100 J = N,1,-1
          J1 = KEEP(IPE+J-1)
C Note that ME57GD does not sort within order in columns
          DO  99 JJ = J1,J2
            IF (KEEP(IFCT+JJ-1).GT.J) GO TO 101
   99       CONTINUE
          SIZE22 = SIZE22 + 1
          J2 = J1-1
  100       CONTINUE
  101       IF (SIZE22 .GT. 0) THEN
C Check to see if there are no entries in (1,1) block.
          DO 98 I = 1,NE
            IF (IRN(I) .LE. N-SIZE22
   *      .AND. JCN(I) .LE. N-SIZE22) THEN
              AVNUM = FLOAT(IWFR+N-SIZE22-1)/FLOAT(N)
              GO TO 96
            ENDIF
   98       CONTINUE
C The (1,1) block is zero.
          OXO = 1
          AVNUM = FLOAT(IWFR-1)/FLOAT(N)
          ENDIF
   96       IF (N .GE. 100000) THEN
          IF (AVNUM.GT.5.42) THEN
            ICNTL6 = 4
          ELSE
            ICNTL6 = 2
          ENDIF
          GO TO 97
          ENDIF
C Logic for OXO matrices
          IF (OXO.EQ.1) THEN
          IF (FLOAT(N-SIZE22)/FLOAT(SIZE22) .GT. 1.8D0) THEN
            ICNTL6 = 2
          ELSE
            ICNTL6 = 4
          ENDIF
          GO TO 97
          ENDIF
C We can try further simple logic here ... then ...
C Call MC47 to test whether fill-in projected is large
          LW = LKEEP-IFCT+1
          CALL MC47BD(N,LW,KEEP(IPE),IWFR,KEEP(COUNT),
   +            KEEP(IFCT),IWORK(NV),
   +            KEEP(INVP),KEEP(PERM),IWORK(IW1),
   +            IWORK(IW2),IWORK(IW3),IWORK(IW4),
   +            ICNT47,INF47,RINF47)
          INFO(13) = INF47(2)
          ICNTL6 = 2
          NEMIN    = ICNTL(12)
    CALL ME57LD(N,KEEP(IPE),IWORK(NV),KEEP(IPS),IWORK(NELIM),
   +        KEEP(NSTK),KEEP(NODE),KEEP(PERM),
   +        KEEP(NSTEPS),IWORK(FILS),IWORK(FRERE),IWORK(ND),
   +        NEMIN,KEEP(IRNPRM))
          NST = KEEP(NSTEPS)
    CALL ME57MD(KIND,N,NE,IRN,JCN,KEEP(MAP),KEEP(IRNPRM),
   +        KEEP(LROW),KEEP(PERM),
   +        IWORK(IW2),IWORK(IW5))
          KEEP(EXPNE) = IWORK(IW5)
    CALL ME57ND(N,KEEP(LROW),KEEP(NSTK),IWORK(NELIM),
   +        IWORK(ND),NST,IWORK(IW1),IWORK(IW2),
   +        INFO,RINFO)
C Check relative fill-in of MC47
          IF (FLOAT(INFO(5))/FLOAT(NE) .LT. 10.0) THEN
C We will run with the MC47 ordering
          GO TO 93
          ELSE
C Save value of relative fill-in for testing against METIS
          MC47FI = FLOAT(INFO(5))/FLOAT(NE)
C Must test METIS ordering now ... ugh
      CALL ME57GD(KIND,N,NE,IRN,JCN,KEEP(IFCT),KEEP(IPE),
   +          KEEP(COUNT),KEEP(IW),IWFR,ICNTL,INFO)
      IF (INFO(1).LT.0) RETURN
          KEEP(IPE+N) = IWFR
          METFTN    = 1
          METOPT(1) = 0
          IF (N.LT.50) GO TO 92
          DO 91 I = 1,N
            IF ((KEEP(IPE+I)-KEEP(IPE+I-1)) .GT. N/10) THEN
            METOPT(1) = 1
            METOPT(2) = 3
            METOPT(3) = 1
            METOPT(4) = 2
            METOPT(5) = 0
            METOPT(6) = 1
            METOPT(7) = 200
            METOPT(8) = 1
            GO TO 92
            ENDIF
   91       CONTINUE
   92   CALL METIS_NODEND_ME57(N,KEEP(IPE),KEEP(IFCT),METFTN,METOPT,
   *              KEEP(NSTK),KEEP(PERM))
      CALL ME57JD(KIND,N,NE,IRN,JCN,KEEP(PERM),KEEP(IFCT),KEEP(IPE),
   +          KEEP(COUNT),IWORK(IW1),IWFR,ICNTL,INFO)
      IF (INFO(1).LT.0) RETURN
          LW = LKEEP - IFCT + 1
      CALL ME57KD(N,KEEP(IPE),KEEP(IFCT),LW,IWFR,KEEP(PERM),
   +          KEEP(INVP),IWORK(NV),IWORK(IW1),NCMPA)
          INFO(13) = NCMPA
          NEMIN = ICNTL(12)
    CALL ME57LD(N,KEEP(IPE),IWORK(NV),KEEP(IPS),IWORK(NELIM),
   +        KEEP(NSTK),KEEP(NODE),KEEP(PERM),
   +        KEEP(NSTEPS),IWORK(FILS),IWORK(FRERE),IWORK(ND),
   +        NEMIN,KEEP(IRNPRM))
          NST = KEEP(NSTEPS)
    CALL ME57MD(KIND,N,NE,IRN,JCN,KEEP(MAP),KEEP(IRNPRM),
   +        KEEP(LROW),KEEP(PERM),
   +        IWORK(IW2),IWORK(IW5))
          KEEP(EXPNE) = IWORK(IW5)
    CALL ME57ND(N,KEEP(LROW),KEEP(NSTK),IWORK(NELIM),
   +        IWORK(ND),NST,IWORK(IW1),IWORK(IW2),
   +        INFO,RINFO)
          IF (FLOAT(INFO(5))/FLOAT(NE).LT.MC47FI) THEN
            ICNTL6 = 4
            GO TO 93
          ELSE
C Double groan  ... we will run with MC47 after all
            ICNTL6=2
C KEEP(IPE) has been corrupted must reset it.
      CALL ME57GD(KIND,N,NE,IRN,JCN,KEEP(IFCT),KEEP(IPE),
   +          KEEP(COUNT),KEEP(IW),IWFR,ICNTL,INFO)
      IF (INFO(1).LT.0) RETURN
            GO TO 97
          ENDIF
C End of METIS check
          ENDIF
C End of indef case calculation
        ENDIF
C End of logic for ICNTL6 = 5
      ENDIF

   97   IF (ICNTL6.EQ.4) THEN
C Set last pointer in IPE
        KEEP(IPE+N) = IWFR
C Use MeTiS ordering
C Set flag for Fortran-style numbering of arrays
        METFTN    = 1
C This would use only defaults
        METOPT(1) = 0
C Set options for METIS, particularly one for dense columns
C First determine if there are any dense columns
        IF (N.LT.50) GO TO 103
        DO 102 I = 1,N
          IF ((KEEP(IPE+I)-KEEP(IPE+I-1)) .GT. N/10) THEN
C The rest are set to default values
          METOPT(1) = 1
          METOPT(2) = 3
          METOPT(3) = 1
          METOPT(4) = 2
          METOPT(5) = 0
          METOPT(6) = 1
          METOPT(7) = 200
          METOPT(8) = 1
          GO TO 103
          ENDIF
  102     CONTINUE
  103   CALL METIS_NODEND_ME57(N,KEEP(IPE),KEEP(IFCT),METFTN,METOPT,
   *              KEEP(NSTK),KEEP(PERM))
        GO TO 111
      ENDIF


C Obtain ordering using approximate minimum degree ordering.
C Input IPE,IWFR,COUNT,IFCT as from ME57GD.
C Output
C   IPE (- father pointer .. if NV(I) > 0, - subordinate variable
C   pointer if NV(I) = 0)
C   NV(I) for subordinate variables of supervariable, otherwise
C   degree when eliminated.
C   IWFR is set to length required by MC47BD if no compresses.
C   COUNT, IFCT undefined on exit.
C   INVP is inverse permutation and PERM is permutation
C Length of LW set to maximum to avoid compresses in MC47B/BD
      LW = LKEEP-IFCT+1
C ICNTL6 =  0.  MC47 uses code for dealing with dense rows disabled
C In HSL 2002 it was only used in the F90 version.
C ICNTL6 =  2.  MC47 implements code for dealing with dense rows
      IF (ICNTL6 .EQ. 0) ICNT47(4) = -1
      CALL MC47BD(N,LW,KEEP(IPE),IWFR,KEEP(COUNT),
   +          KEEP(IFCT),IWORK(NV),
   +          KEEP(INVP),KEEP(PERM),IWORK(IW1),
   +          IWORK(IW2),IWORK(IW3),IWORK(IW4),
   +          ICNT47,INF47,RINF47)
      INFO(13) = INF47(2)

      ELSE
C End of ICNTL6 .NE. 3

C MA27 ordering being used.  Must insert row lengths in KEEP(IFCT)
C Length of LW set to maximum to avoid compresses.
      LW = LKEEP-IF27H+1
      CALL ME57VD(KIND,N,NE,IRN,JCN,KEEP(IF27H),LW,KEEP(IPE),
   *          IWORK(IW1),IWORK(IW2),IWFR,ICNTL,INFO)
      IF (INFO(1).LT.0) RETURN
C Analyse using minimum degree ordering
      THRESH = FLOAT(ICNTL(14))/100.0
      CALL ME57HD(N,KEEP(IPE),KEEP(IF27H),LW,IWFR,IWORK(NV),
   *          IWORK(IW1),IWORK(IW2),IWORK(IW3),IWORK(IW4),
   +          2139062143,INFO(13),THRESH)
C Set IPE correctly
      DO 110 I = 1,N
        IF (IWORK(NV+I-1).NE.0) GO TO 110
        IN = I
  105     IL = IN
        IN = - KEEP(IPE+IL-1)
        IF (IWORK(NV+IN-1).EQ.0) GO TO 105
C Make subordinate node point to principal node
        KEEP(IPE+I-1) = -IN
  110   CONTINUE
      ENDIF

C End of block for generating ordering
    ENDIF

  111 IF (ICNTL6.EQ.1 .OR. ICNTL6.EQ.4) THEN
C If we have generated ordering using MeTiS then we need to feed
C   permutation as if it were coming from the user as we do not
C   have a tight coupling to MeTiS as for other orderings.

C Sort using given order
      CALL ME57JD(KIND,N,NE,IRN,JCN,KEEP(PERM),KEEP(IFCT),
   +       KEEP(IPE),KEEP(COUNT),IWORK(IW1),IWFR,ICNTL,INFO)
      IF (INFO(1).LT.0) RETURN

C Generating tree using given ordering
C Input:  N,KEEP(IPE),KEEP(IFCT),LW,IWFR,KEEP(PERM)
C Output:  KEEP(IPE),IWORK(NV)
CCC For some reason LW was set to 2*NE which failed when NE < N
C     LW = 2*NE
      LW = LKEEP - IFCT + 1
      CALL ME57KD(N,KEEP(IPE),KEEP(IFCT),LW,IWFR,KEEP(PERM),
   +          KEEP(INVP),IWORK(NV),IWORK(IW1),NCMPA)
      INFO(13) = NCMPA

    END IF


C Perform depth-first search of assembly tree
C Set NEMIN
    NEMIN = ICNTL(12)
C Input  IPE,NV,NEMIN
C Output
C   IPE .. father and younger brother pointer
C   NV  .. unchanged
C   NE/NSTK/ND defined for nodes of tree
C   PERM
C   IPS(I) position of node I in order
C   LROW(I) is size of frontal matrix at node I
    CALL ME57LD(N,KEEP(IPE),IWORK(NV),KEEP(IPS),IWORK(NELIM),
   +        KEEP(NSTK),KEEP(NODE),KEEP(PERM),
   +        KEEP(NSTEPS),IWORK(FILS),IWORK(FRERE),IWORK(ND),
   +        NEMIN,KEEP(IRNPRM))
    NST = KEEP(NSTEPS)

C Construct map for storing the permuted upper triangle by rows.
C Input N,NE,IRN,JCN,PERM
C Output MAP,LROW,IRNPRM
    CALL ME57MD(KIND,N,NE,IRN,JCN,KEEP(MAP),KEEP(IRNPRM),
   +        KEEP(LROW),KEEP(PERM),
   +        IWORK(IW2),IWORK(IW5))

C Set number of entries in expanded input matrix
    KEEP(EXPNE) = IWORK(IW5)

C Evaluate storage and operation counts.
C Input  LROW,NSTK,NELIM,ND
C Output LROW,NSTK (unchanged)
    CALL ME57ND(N,KEEP(LROW),KEEP(NSTK),IWORK(NELIM),
   +        IWORK(ND),NST,IWORK(IW1),IWORK(IW2),
   +        INFO,RINFO)

C Set INFO entry to record ordering used
   93 INFO(36) = ICNTL6
C Add for BIGA
    ALENB    = 1
C Add for Schnabel-Eskow
    IF (ICNTL(7).EQ.4) ALENB = ALENB + N + 5
C Add for scaling
    IF (ICNTL(15).EQ.1) ALENB = ALENB + N

C Allow enough to get started
    INFO(9)  = MAX(INFO(9)+ALENB,ALENB+KEEP(EXPNE)+1)
    INFO(11) = MAX(INFO(11)+ALENB,ALENB+KEEP(EXPNE)+1)
C This is N+5 for starting the factorization, N for first row (max).
    INFO(10) = MAX(INFO(10),KEEP(EXPNE)+N+5)
    INFO(12) = MAX(INFO(12),KEEP(EXPNE)+N+5)

C Needed by ME57B/BD
    IF (ICNTL(15).EQ.1) THEN
      INFO(9) = MAX(INFO(9),ALENB+3*KEEP(EXPNE)+3*N)
      INFO(11) = MAX(INFO(11),ALENB+3*KEEP(EXPNE)+3*N)
C Allow space for integers in computing scaling factors
      INFO(10) = MAX(INFO(10),3*KEEP(EXPNE)+5*N+1)
      INFO(12) = MAX(INFO(12),3*KEEP(EXPNE)+5*N+1)
    ENDIF

C If requested, print parameter values on exit.
    IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
      NZE = KEEP(EXPNE)
      WRITE (MP,99999) INFO(1),NZE,
   *            (INFO(I),I=3,13),INFO(36),(RINFO(I),I=1,2)
99999 FORMAT (/'Leaving analysis phase (ME57AD) with ...'/
   1    'INFO(1)  Error indicator              =',I12/
   2    'Number of entries in matrix incl diagonal   =',I12/
   2    'INFO(3)  Number of out-of-range indices     =',I12/
   2    'INFO(4)  Number of off-diagonal duplicates    =',I12/
   2    'INFO(5)  Forecast real storage for factors    =',I12/
   3    '----(6)  Forecast integer storage for factors =',I12/
   3    '----(7)  Forecast maximum front size      =',I12/
   4    '----(8)  Number of nodes in assembly tree   =',I12/
   5    '----(9)  Size of FACT without compress      =',I12/
   6    '----(10) Size of IFACT without compress     =',I12/
   5    '----(11) Size of FACT with compress       =',I12/
   5    '----(12) Size of IFACT with compress      =',I12/
   5    '----(13) Number of compresses           =',I12/
   5    '----(36) Ordering strategy used by code     =',I12/
   9    'RINFO(1) Forecast additions for assembly    =',1P,D12.5/
   9    'RINFO(2) Forecast ops for elimination       =',1P,D12.5)

      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(/A/(5I12))')  'Permutation array:',
   +            (KEEP(I),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
      WRITE (MP,'(/A/(5I12))')
   +      'Number of entries in rows of permuted matrix:',
   +      (KEEP(LROW+I-1),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
      K = MIN(10,NZE)
      IF (LDIAG.GE.4) K = NZE
      WRITE (MP,'(/A/(5I12))')
   *      'Column indices of permuted matrix:',
   *                 (KEEP(IRNPRM+I-1),I=1,K)
      IF (K.LT.NZE) WRITE (MP,'(16X,A)') '   . . .'
      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(/A/(5I12))')
   +    'Tree nodes at which variables eliminated:',
   +    (KEEP(NODE+I-1),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
      K = MIN(10,NE)
      IF (LDIAG.GE.4) K = NE
      WRITE (MP,'(/A/(5I12))') 'Map array:',
   *                   (KEEP(I),I=MAP,MAP+K-1)
      IF (K.LT.NE) WRITE (MP,'(16X,A)') ' . . .'
    END IF

    RETURN

C Error conditions.
   20 INFO(1) = -1
    INFO(2) = N
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A,I10)')
   +    '**** Error return from ME57AD ****  INFO(1) =',INFO(1),
   +    'N has value ',INFO(2)
    RETURN

   30 INFO(1) = -2
    INFO(2) = NE
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A,I10)')
   +    '**** Error return from ME57AD ****  INFO(1) =',INFO(1),
   +    'NE has value',INFO(2)
     RETURN

   40 INFO(1) = -15
    INFO(2) = LKEEP
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A,I10/A,I10)')
   +    '**** Error return from ME57AD ****  INFO(1) =',INFO(1),
   +    'LKEEP has value    ',INFO(2),
   +    'Should be at least ',5*N+NE+MAX(N,NE)+42
     RETURN

   80 INFO(1) = -9
    INFO(2) = I
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A/A,I10,A)')
   +    '**** Error return from ME57AD ****  INFO(1) =',INFO(1),
   +    'Invalid permutation supplied in KEEP',
   +    'Component',INFO(2),' is faulty'
    RETURN

   90 INFO(1) = -18
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A)')
   +    '**** Error return from ME57AD ****  INFO(1) =',INFO(1),
   +    'MeTiS ordering requested but MeTiS not linked'

    END


C--------------------------------------------------------------------
C        HSL 2000 (2000)
C      --
C-         Copyright Rutherford Appleton Laboratory
C      --
C--------------------------------------------------------------------
    SUBROUTINE ME57BD(KIND, N, NE, A, FACT, LFACT, IFACT, LIFACT,
   * LKEEP, KEEP, PPOS, WORK, ICNTL, CNTL, INFO, RINFO)
C
C Purpose
C =======
C
C
C This subroutine computes the factorization of the matrix input in
C   A using information (in KEEP and IFACT) from ME57AD.
C
    INTEGER KIND,N,NE,LFACT,LIFACT,LKEEP
    COMPLEX*16 A(NE),FACT(LFACT)
    DOUBLE PRECISION WORK(2*NE+3*N)
    DOUBLE PRECISION RINFO(20)
C
C   Control parameters: see description in ME57ID
    DOUBLE PRECISION CNTL(5)
    INTEGER ICNTL(20), IFACT(LIFACT)
    INTEGER   INFO(40), KEEP(LKEEP), PPOS(N)
C
C Parameters
C ==========
C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C
C N is an INTEGER variable which must be set by the user to the
C    order n of the matrix A. It must be unchanged since the
C    last call to ME57AD and is not altered by the
C    subroutine.  Restriction: N > 0
C
C NE is an INTEGER variable which must be set by the user to the
C    number of    entries in the matrix A.  It is not altered by
C    the subroutine.  Restriction: NE >= 0.
C
C A is a COMPLEX (COMPLEX*16 in the D version) array of length NE.
C   It is not altered by the subroutine.
C
C FACT is a COMPLEX (COMPLEX*16 in the D version)
C    array of length LFACT. It need not
C    be set by the user. On exit, entries 1 to INFO(15) of FACT hold
C    the real part of the factors and should be passed unchanged to
C    ME57CD.
C
C LFACT is an INTEGER variable that must be set by the user to
C    the size of array FACT.
C    It should be passed unchanged to ME57CD.
C
C IFACT is an INTEGER array of length LIFACT. It need not
C    be set by the user. On exit, entries 1 to INFO(16) of IFACT hold
C    the integer part of the factors and should be passed unchanged to
C    ME57CD.
C
C LIFACT is an INTEGER variable that must be set by the user to
C    the size of array IFACT.
C    It should be passed unchanged to ME57CD.
C
C LKEEP is an INTEGER that must be set to the length of array KEEP.
C    Restriction: LKEEP >= 5*N+NE+MAX(N,NE)+42.
C
C KEEP is an INTEGER  array of length LKEEP which must be
C     passed unchanged since the last call to ME57AD.  It is not
C     altered by ME57BD.
C
C PPOS is an INTEGER array of length N that is used as workspace.
C
C ICNTL is an INTEGER array of length 20
C        that contains control parameters and must be set
C        by the user.
C        Default values for the components may be set by a call
C        to ME57ID. Details
C        of the control parameters are given in ME57ID
C
C CNTL is a REAL (DOUBLE_PRECISION in the D version) array of length 5
C        that contains control parameters and must be set
C        by the user.
C        Default values for the components may be set by a call
C        to ME57ID. Details
C        of the control parameters are given in ME57ID
C
C RINFO is a REAL (DOUBLE_PRECISION in the D version) array of
C       length 20 which need not be set by the user.
C       We describe below the components of this array modified
C       in the subroutine.
C
C    ______(3)  Number of floating point operations involved
C       during the assembly process
C
C    ______(4)  Number of floating point operations involved
C       during the elimination process
C
C    ______(5)  Number of extra floating point operations caused by
C       use of Level 3 BLAS
C
C
C    .. Error Return ..
C       ============
C
C  A successful return from ME57BD
C  is indicated by a value of INFO(1) positive.
C  Negative values of INFO(1) correspond to
C  error message whereas positive values correspond to
C  warning messages. A negative value of INFO(1) is associated with
C  an error message which will be output on unit ICNTL(1).
C
C
C     .. Local variables ..
C      ===============
C
    INTEGER EXPNE,HOLD,I,IRNPRM,K,LDIAG,LLFACT,LP,LROW,MAP,MM1,MM2,MP
    INTEGER J,JJ,KK,ISCALE,NUM,NE64,IDUP,IMAT,IPT,JLOOP,JNEW,NN,ISING
    INTEGER J1,J2
    INTEGER NSTEPS,NODE,NSTK,PERM,INEW,ALENB

    DOUBLE PRECISION ONE,ZERO,RINF,FD15AD,FCT,SMAX,SMIN
    DOUBLE PRECISION BIGA
    PARAMETER (ONE = 1.0D0, ZERO=0.0D0)

C?? To identify bug
C   LOGICAL LCASE
C   COMMON /CCASE/LCASE

    INTRINSIC MIN

C
C EXPNE is number of entries of original matrix plus any missing
C    diagonals.
C HOLD points to position in IFACT for first entry of array that holds
C    values for restart.
C MM1,MM2 are used to define start point of arrays for matrix
C    modification.  Needed if LFACT < N + 5.
C NSTEPS holds the number of nodes in the tree.

C
C External Subroutines
C ====================
C

    EXTERNAL ME57OD,ME57UD,FD15AD,MC34AD,MC64WD

C Set RINF to largest positive real number (infinity)
    RINF = FD15AD('H')
C Set REPS to smallest number st 1.0+REPS > 1.0
C   REPS = FD15AD('E')

C Set INFO for -3 and -4 error return. Ok that is done on every call
    INFO(17) = 0
    INFO(18) = 0

C
C
C Initialisation of printing controls.
C
    LP   = ICNTL(1)
    MP   = ICNTL(3)
    LDIAG  = ICNTL(5)
C
C??
C   Check if analysis has been effectively performed
C
    IF (N.LE.0)  GO TO 25
    IF (NE.LT.0) GO TO 30
    IF (LKEEP.LT.5*N+NE+MAX(N,NE)+42) GO TO 40
    IF (ICNTL(7).LT.1 .OR. ICNTL(7).GT.4) GO TO 35

C   Partition of array KEEP
    NSTEPS = KEEP(N+1)
    EXPNE  = KEEP(N+2)
    PERM = 1
    HOLD = PERM + N + 2
    LROW = HOLD + 40
    NODE = LROW + N
    NSTK = NODE + N
    MAP  = NSTK + N
    IRNPRM = MAP + MAX(NE,N)

    LLFACT = LFACT - 1

    IF (ICNTL(15).EQ.1) THEN
C Matrix is being scaled using MC64SYM
      ISCALE = LLFACT - N + 1
      LLFACT = ISCALE - 1
    ENDIF

    IF (ICNTL(7).EQ.4) THEN
C Schnabel-Eskow modification being used
C Reserve space in FACT for diagonal entries and controls.
      LLFACT = LLFACT - N - 5
C Set MM1 and MM2 to point to first entries in arrays.
      MM1 = LLFACT+6
      MM2 = LLFACT+1
    ELSE
      MM1 = 1
      MM2 = 1
    ENDIF

C One for BIGA
    ALENB = 1
    IF (ICNTL(7).EQ.4)  ALENB = ALENB + N + 5
    IF (ICNTL(15).EQ.1) ALENB = ALENB + N
C +1 because MAP of o-o-r maps to entry 0
    IF (LLFACT.LT.EXPNE+1)   GO TO 85
C The first five positions and room for a whole row are needed
C at the beginning of the ME57O/OD factorization
    IF (LIFACT.LT.EXPNE+N+5)  GO TO 95
C Check that there is enough space for scaling within ME57B/BD.
    IF (ICNTL(15).EQ.1)  THEN
      IF (LFACT .LT. ALENB + 3*EXPNE  + 3*N) GO TO 85
      IF (LIFACT .LT. 3*EXPNE + 5*N + 1) GO TO 95
    ENDIF

C
C PRINTING OF INPUT PARAMETERS
C*****************************
    IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
      WRITE (MP,99999)
99999 FORMAT (/'Entering factorization phase (ME57BD) with ...')
      IF (KEEP(HOLD).GT.0) WRITE (MP,99998)
99998 FORMAT ('Re-entry call after call to ME57ED')
      IF (KIND.EQ.1) THEN
      WRITE(MP,'(A)') 'Matrix is complex symmetric'
      ELSE
      WRITE(MP,'(A)') 'Matrix is Hermitian'
      ENDIF
      WRITE (MP,99997) N,NE,EXPNE,(ICNTL(I),I=1,5),ICNTL(7),ICNTL(8),
   +       ICNTL(11),ICNTL(15),ICNTL(16),LFACT,LIFACT, NSTEPS,
   +       CNTL(1), CNTL(2), CNTL(4), CNTL(5)
99997 FORMAT ('N     Order of input matrix         =',I12/
   2      'NE    Entries in input matrix         =',I12/
   2      '      Entries in input matrix (inc diags) =',I12/
   6      'ICNTL(1)  Stream for errors           =',I12/
   7      ' --- (2)  Stream for warnings         =',I12/
   8      ' --- (3)  Stream for monitoring         =',I12/
   9      ' --- (4)  Stream for statistics         =',I12/
   1      ' --- (5)  Level of diagnostic printing    =',I12/
   1      ' --- (7)  Numerical pivoting control      =',I12/
   1      ' --- (8)  Restart or discard factors      =',I12/
   1      ' --- (11) Block size for Level 3 BLAS     =',I12/
   1      ' --- (15) Scaling control (1 on)        =',I12/
   1      ' --- (16) Dropping control (1 on)       =',I12/
   4      'LFACT   Size of real working space      =',I12/
   5      'LIFACT  Size of integer working space     =',I12/
   7      '      Number nodes in assembly tree     =',I12/
   9      'CNTL(1) Value of threshold parameter      =',D12.5/
   9      'CNTL(2) Threshold for zero pivot        =',D12.5/
   9      'CNTL(4) Control for value of static pivots  =',D12.5/
   9      'CNTL(5) Control for number delayed pivots   =',D12.5)
      K = MIN(10,NE)
      IF (LDIAG.GE.4) K = NE
      IF (NE.GT.0) THEN
C      WRITE (MP,'(/A/(2(I4,A,1P,2D16.8,A)))') 'Matrix entries:',
C   +   (I,': (',A(I),')',I=1,K)
      IF (K.LT.NE) WRITE (MP,'(A)') '   . . .'
      END IF
      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(/A/(5I12))')  'Permutation array:',
   +            (KEEP(I),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
      WRITE (MP,'(/A/(5I12))')
   +      'Number of entries in rows of permuted matrix:',
   +      (KEEP(LROW+I-1),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
      WRITE (MP,'(/A/(5I12))')
   +    'Tree nodes at which variables eliminated:',
   +    (KEEP(NODE+I-1),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
      K = MIN(10,NSTEPS)
      IF (LDIAG.GE.4) K = NSTEPS
      IF (K.GT.0) WRITE (MP,'(/A/(5I12))')
   +   'Number of assemblies at each tree node:',
   +   (KEEP(NSTK+I-1),I=1,K)
      IF (K.LT.NSTEPS) WRITE (MP,'(16X,A)') ' . . .'
      K = MIN(10,NE)
      IF (LDIAG.GE.4) K = NE
      WRITE (MP,'(/A/(5I12))') 'Map array:',
   *                   (KEEP(I),I=MAP,MAP+K-1)
      IF (K.LT.NE) WRITE (MP,'(16X,A)') ' . . .'
      K = MIN(10,EXPNE)
      IF (LDIAG.GE.4) K = EXPNE
      WRITE (MP,'(/A/(5I12))')
   *      'Column indices of permuted matrix:',
   *                   (KEEP(IRNPRM+I-1),I=1,K)
      IF (K.LT.EXPNE) WRITE (MP,'(16X,A)') '   . . .'
    ENDIF

C Jump if it is reentry
    IF (KEEP(HOLD) .GT. 0) THEN
      BIGA = ABS(FACT(LFACT))
      GO TO 22
    ENDIF

C
C***************************************************
C MAP input nonzeros to appropriate position in FACT
C***************************************************
C
C?? For the moment to handle missing diagonals
    DO 19 K = 1,EXPNE
      FACT(LLFACT-EXPNE+K) = ZERO
   19 CONTINUE

    BIGA = ZERO
    IF (KIND.EQ.1) THEN
      DO 20 K = 1,NE
      BIGA = MAX(BIGA,ABS(A(K)))
      FACT(KEEP(MAP+K-1)+LLFACT-EXPNE) = A(K)
   20   CONTINUE
    ENDIF
    IF (KIND.NE.1) THEN
C Matrix is Hermitian
      DO 201 K = 1,NE
      BIGA = MAX(BIGA,ABS(A(K)))
      KK = KEEP(MAP+K-1)
      IF (KK.GT.0) THEN
        FACT(KK+LLFACT-EXPNE) = A(K)
      ELSE
C Entry is switching triangles
        FACT(-KK+LLFACT-EXPNE) = CONJG(A(K))
      ENDIF
  201   CONTINUE
    ENDIF
    RINFO(18) = BIGA
    DO 21 K = 1,EXPNE
      IFACT(LIFACT-EXPNE+K) = KEEP(IRNPRM+K-1)
   21 CONTINUE
C Invert array PERM
    DO 23 I = 1,N
      PPOS(KEEP(PERM+I-1)) = I
   23 CONTINUE

    IF (ICNTL(15).EQ.1) THEN
C Scaling using MC64.  Matrix must be generated in correct format.

    IPT = 1
    IDUP = IPT+N+1
    IMAT = IDUP+N
    ISING = IMAT + MAX(NE,EXPNE)

C Copy matrix, remove duplicates, and initialize IP array.
    DO 4444 I = 1,N
      IFACT(IDUP+I-1) = 0
 4444 CONTINUE
C Must use new coordinates to keep matrix (half) symmetric
    IFACT(IPT) = 1
    KK = 1
    K = 1
    DO 3333 J = 1,N
      DO 2222 JJ = 1,KEEP(LROW+J-1)
      I = KEEP(PERM+IFACT(LIFACT-EXPNE+K)-1)
      IF (IFACT(IDUP+I-1).GE.IFACT(IPT+J-1)) THEN
C Duplicate
        FACT(IFACT(IDUP+I-1)) =
   &      FACT(IFACT(IDUP+I-1)) + FACT(LLFACT-EXPNE+K)
      ELSE
C Remove explicit zeros
        IF (FACT(LLFACT-EXPNE+K).NE.ZERO) THEN
          IFACT(IDUP+I-1) = KK
          FACT(KK) = FACT(LLFACT-EXPNE+K)
          IFACT(IMAT-1+KK) = I
          KK = KK+1
        ENDIF
      ENDIF
      K = K + 1
 2222   CONTINUE
      IFACT(IPT+J) = KK
 3333 CONTINUE

C Check that Hermitian matrix has reals on diagonal
C Note that this check only applies when scaling is switched on
    IF (KIND.NE.1) THEN
      DO 3111 J = 1,N
      J1 = IFACT(IPT+J-1)
      J2 = IFACT(IPT+J)-1
      DO 2111 JJ = J1,J2
        I = IFACT(IMAT-1+JJ)
        IF (I.EQ.J) THEN
          IF (FACT(JJ) .NE. REAL(FACT(JJ))) THEN
          INFO(1) = -20
          DO 2110 K = 1,N
            IF (KEEP(PERM+K-1).EQ.I) THEN
            INFO(2) = K
            GO TO 2109
            ENDIF
 2110       CONTINUE
 2109 IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I8,A)')
   +   '**** Error return from ME57BD ****  INFO(1) =',INFO(1),
   *   'Hermitian matrix but diagonal entry ',INFO(2),' is complex'
          RETURN
          ENDIF
        ENDIF
 2111   CONTINUE
 3111   CONTINUE
    ENDIF

    DO 3334 I=1,KK-1
      WORK(I) = ABS(FACT(I))
 3334 CONTINUE
C Expand matrix
    CALL MC34AD(N,IFACT(IMAT),IFACT(IPT),.TRUE.,WORK,KEEP(PERM))
    NE64 = IFACT(IPT+N)-1
    DO 75 J = 1,N
      FCT = ZERO
      DO 60 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
C       FACT(K) = ABS(FACT(K))
      IF (WORK(K).GT.FCT) FCT = WORK(K)
   60   CONTINUE
      WORK(NE64+2*N+J) = FCT
      IF (FCT.NE.ZERO) THEN
      FCT = LOG(FCT)
      ELSE
C This can happen if only if column is null so matrix is singular.
      FCT = RINF/N
      ENDIF
      DO 70 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
C Note that zeros have been screened out so that FACT(K) always > 0.
      WORK(K) = FCT - LOG(WORK(K))
   70   CONTINUE
   75 CONTINUE
C Scale matrix
C B = DW(3N+1:3N+NE); IW(1:5N) and DW(1:2N) are workspaces
    CALL MC64WD(N,NE64,IFACT(IPT),IFACT(IMAT),WORK,KEEP(PERM),NUM,
   &   IFACT(IDUP),IFACT(IMAT+NE64),IFACT(IMAT+NE64+N),
   &   IFACT(IMAT+NE64+2*N),IFACT(IMAT+NE64+3*N),
   &   WORK(NE64+1),WORK(NE64+N+1))
    IF (NUM.EQ.N) THEN
      DO 80 J = 1,N
      WORK(NE64+N+J) = WORK(NE64+N+J) - LOG(WORK(NE64+2*N+J))
   80   CONTINUE
C Check size of scaling factors
C   FCT = 0.5*LOG(RINF)
C   DO 86 J = 1,N
C     IF (FACT(NE64+J).LT.FCT .AND. FACT(NE64+N+J).LT.FCT) GO TO 86
C     INF64(1) = 2
C  86 CONTINUE

C Scaling is permuted to scaling on original matrix
      DO 5555 I=1,N
        FACT(ISCALE+PPOS(I)-1) =
   &      SQRT(EXP(WORK(NE64+I)+WORK(NE64+N+I)))
 5555   CONTINUE
      ELSE
C Matrix is singular
C Regenerate PERM and set PPOS to indicate nonsingular block
      K = 0
      DO 3501 I = 1,N
      IF (KEEP(PERM+I-1).LT.0) THEN
        PPOS(I) = -PPOS(I)
        IFACT(ISING+I-1) = 0
      ELSE
        K = K + 1
        IFACT(ISING+I-1) = K
      ENDIF
 3501   CONTINUE
      DO 3502 I = 1,N
      KEEP(PERM+ABS(PPOS(I))-1) = I
 3502   CONTINUE
C Copy matrix, remove duplicates, and initialize IP array.
      DO 3503 I = 1,N
      IFACT(IDUP+I-1) = 0
 3503   CONTINUE
C Must use new coordinates to keep matrix (half) symmetric
      IFACT(IPT) = 1
      KK = 1
      K = 1
      JNEW = 0
      NN = N
      DO 3505 J = 1,N
      IF (PPOS(J).LT.0) THEN
        NN = NN - 1
        K = K + KEEP(LROW+J-1)
        GO TO 3505
      ENDIF
      JNEW = JNEW + 1
      DO 3504 JJ = 1,KEEP(LROW+J-1)
        I = KEEP(PERM+IFACT(LIFACT-EXPNE+K)-1)
        IF (PPOS(I).GT.0) THEN
          IF (IFACT(IDUP+I-1).GE.IFACT(IPT+J-1)) THEN
C Duplicate
          FACT(IFACT(IDUP+I-1)) =
   &        FACT(IFACT(IDUP+I-1)) + FACT(LLFACT-EXPNE+K)
          ELSE
C Remove explicit zeros
          IF (FACT(LLFACT-EXPNE+K).NE.ZERO) THEN
            IFACT(IDUP+I-1) = KK
            FACT(KK) = FACT(LLFACT-EXPNE+K)
            IFACT(IMAT-1+KK) = IFACT(ISING+I-1)
            KK = KK+1
          ENDIF
          ENDIF
        ENDIF
        K = K + 1
 3504   CONTINUE
      IFACT(IPT+JNEW) = KK
 3505   CONTINUE
      DO 3500 I=1,KK-1
      WORK(I) = ABS(FACT(I))
 3500   CONTINUE
    NE64 = IFACT(IPT+NN)-1
      CALL MC34AD(NN,IFACT(IMAT),IFACT(IPT),.TRUE.,WORK,KEEP(PERM))
      NE64 = IFACT(IPT+NN)-1
      DO 3508 J = 1,NN
      FCT = ZERO
      DO 3506 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
C       FACT(K) = ABS(FACT(K))
        IF (WORK(K).GT.FCT) FCT = WORK(K)
 3506   CONTINUE
      WORK(NE64+2*N+J) = FCT
      FCT = LOG(FCT)
      DO 3507 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
        WORK(K) = FCT - LOG(WORK(K))
 3507   CONTINUE
 3508   CONTINUE
      CALL MC64WD(NN,NE64,IFACT(IPT),IFACT(IMAT),WORK,KEEP(PERM),NUM,
   &    IFACT(IDUP),IFACT(IMAT+NE64),IFACT(IMAT+NE64+N),
   &    IFACT(IMAT+NE64+2*N),IFACT(IMAT+NE64+3*N),
   &    WORK(NE64+1),WORK(NE64+N+1))
      DO 3509 J = 1,NN
      WORK(NE64+N+J) = WORK(NE64+N+J) - LOG(WORK(NE64+2*N+J))
 3509   CONTINUE
C Check size of scaling factors
C   FCT = 0.5*LOG(RINF)
C   DO 86 J = 1,N
C     IF (FACT(NE64+J).LT.FCT .AND. FACT(NE64+N+J).LT.FCT) GO TO 86
C     INF64(1) = 2
C  86 CONTINUE

C Scaling is permuted to scaling on original matrix for scale factors
C for nonsingular block
      K=0
C Loop is on new indices
      DO 3510 I=1,N
        IF (PPOS(I).LT.0) THEN
          K = K + 1
          FACT(ISCALE-PPOS(I)-1) = ZERO
        ELSE
          FACT(ISCALE+PPOS(I)-1) =
   &      SQRT(EXP(WORK(NE64+I-K)+WORK(NE64+N+I-K)))
        ENDIF
 3510   CONTINUE
C Compute scaling on nonsingular part
C Remember that PPOS maps from new to original but is flag on new
      DO 3516 I = 1,N
        KEEP(PERM+ABS(PPOS(I))-1) = I
 3516   CONTINUE
C Looping on new indices
      K = 1
      DO 3514 JJ = 1,N
        J = PPOS(JJ)
        IF (J.GT.0) THEN
          DO 3511 JLOOP = 1,KEEP(LROW+JJ-1)
          I = IFACT(LIFACT-EXPNE+K)
C I is original index so have to map to new to do PPOS test
          INEW = KEEP(PERM+I-1)
          IF (PPOS(INEW).LT.0)
   &        FACT(ISCALE+I-1) = MAX(ABS(FACT(ISCALE+I-1)),
   &           ABS(ABS(FACT(LLFACT-EXPNE+K))*FACT(ISCALE+J-1)))
          K = K + 1
 3511       CONTINUE
        ELSE
          DO 3512 JLOOP = 1,KEEP(LROW+JJ-1)
          I = IFACT(LIFACT-EXPNE+K)
C I is original index so have to map to new to do PPOS test
          INEW = KEEP(PERM+I-1)
C Shouldn't happen otherwise nonsingular block not maximum
C Sorry can happen but entry is implicit zero on diagonal
C Note that J is negative
          IF (I .NE. -J)  THEN
          FACT(ISCALE-J-1) =
   &          MAX(ABS(FACT(ISCALE-J-1)),
   &          ABS(ABS(FACT(LLFACT-EXPNE+K))*FACT(ISCALE+I-1)))
          ENDIF
          K = K + 1
 3512       CONTINUE
        ENDIF
 3514   CONTINUE
C Set scaling factors for singular block and reset PPOS
      DO 3513 I = 1,N
        INEW = KEEP(PERM+I-1)
        IF (PPOS(INEW) .LT. 0) THEN
          PPOS(INEW) = - PPOS(INEW)
          IF (FACT(ISCALE+I-1) .EQ. ZERO) THEN
          FACT(ISCALE+I-1) = ONE
          ELSE
          FACT(ISCALE+I-1) = ONE/FACT(ISCALE+I-1)
          ENDIF
        ENDIF
 3513   CONTINUE
      ENDIF
C End of logic for singular matrix

C       DO 8888 I = 1, N
C       FACT(ISCALE+I-1) = ONE
C8888   CONTINUE
      SMAX = ABS(FACT(ISCALE))
      SMIN = ABS(FACT(ISCALE))
      DO 5566 I = 1,N
        SMAX = MAX(SMAX,ABS(FACT(ISCALE+I-1)))
        SMIN = MIN(SMIN,ABS(FACT(ISCALE+I-1)))
 5566   CONTINUE
      RINFO(16) = SMIN
      RINFO(17) = SMAX
C Scale matrix
      K = 1
      BIGA = ZERO
      DO 6666 JJ = 1,N
        J = PPOS(JJ)
        DO 7777 JLOOP = 1,KEEP(LROW+JJ-1)
          I = IFACT(LIFACT-EXPNE+K)
          FACT(LLFACT-EXPNE+K) =
   &      FACT(ISCALE+I-1)*FACT(LLFACT-EXPNE+K)*FACT(ISCALE+J-1)
          BIGA = MAX(BIGA, ABS(FACT(LLFACT-EXPNE+K)))
          K = K + 1
 7777     CONTINUE
 6666   CONTINUE
CPRINT
C Only for testing will remove later
C       DO 6661 I = 1,EXPNE
C       IF (ABS(FACT(LLFACT-EXPNE+I)).GT.ONE+100*REPS) THEN
C         write(LP,*) '**** Matrix not scaled properly'
C         write(LP,*) 'Entry ',I,' has value ',FACT(LLFACT-EXPNE+I)
C         GO TO 6663
C       ENDIF
C6661   CONTINUE
C6663   CONTINUE
C End of scaling
    ELSE
      RINFO(16) = ONE
      RINFO(17) = ONE
    ENDIF
C
C**********************************
C Numerical Factorization
C**********************************
C Work arrays FACT(MM1/MM2), KEEP(PERM), IFACT(1)
    FACT(LFACT) = BIGA
   22 CALL ME57OD(KIND, N, EXPNE, FACT, LLFACT, IFACT, LIFACT,
   *        KEEP(LROW), PPOS,
   *        NSTEPS, KEEP(NSTK), KEEP(NODE), FACT(MM1),
   *        FACT(MM2),
   *        KEEP(PERM),
   *        CNTL, ICNTL,
   *        INFO, RINFO, KEEP(HOLD), BIGA)
    IF (INFO(1).EQ.10 .OR. INFO(1).EQ.11) THEN
      IF (LDIAG.GT.2 .AND. MP.GE.0)  THEN
      IF (INFO(1).EQ.10) WRITE (MP,99982) INFO(1)
99982 FORMAT (/'Leaving factorization phase (ME57BD) with ...'/
   1  'Factorization suspended because of lack of real space'/
   1  'INFO (1) = ',I3)
      IF (INFO(1).EQ.11) WRITE (MP,99983) INFO(1)
99983 FORMAT (/'Leaving factorization phase (ME57BD) with ...'/
   1  'Factorization suspended because of lack of integer space'/
   1  'INFO (1) = ',I3)
      ENDIF
      RETURN
    ENDIF
C Regenerate array PERM
    DO 24 I = 1,N
      KEEP(PERM+PPOS(I)-1) = I
   24 CONTINUE
C Compute INFO(17-20)
      INFO(17) = ALENB + INFO(17)
      INFO(19) = ALENB + INFO(19)
C Allow space for scaling in ME57B/BD
    IF (ICNTL(15).EQ.1) THEN
      INFO(17) = MAX(INFO(17),ALENB + 3*EXPNE+3*N)
      INFO(19) = MAX(INFO(19),ALENB + 3*EXPNE+3*N)
      INFO(18) = MAX(INFO(18),3*EXPNE+5*N+1)
      INFO(20) = MAX(INFO(20),3*EXPNE+5*N+1)
    ENDIF
    IF (INFO(1).EQ.-3) GO TO 85
    IF (INFO(1).EQ.-4) GO TO 95
    IF (INFO(1).LT.0) RETURN
    GO TO 100
C************************
C **** Error returns ****
C************************
   25 INFO(1) = -1
    INFO(2) =  N
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
   +    '**** Error return from ME57BD ****  INFO(1) =',INFO(1),
   +    'N has value ',INFO(2)
    RETURN
   30 INFO(1) = -2
    INFO(2) = NE
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
   +    '**** Error return from ME57BD ****  INFO(1) =',INFO(1),
   +    'NE has value',INFO(2)
    RETURN
   40 INFO(1) = -15
    INFO(2) = LKEEP
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10/A,I10)')
   +    '**** Error return from ME57BD ****  INFO(1) =',INFO(1),
   +    'LKEEP has value    ',INFO(2),
   +    'Should be at least ',5*N+NE+MAX(N,NE)+42
    RETURN
   35 INFO(1) = -10
    INFO(2) = ICNTL(7)
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
   +    '**** Error return from ME57BD ****  INFO(1) =',INFO(1),
   +    'ICNTL(7) has value',ICNTL(7)
    RETURN

   85 INFO(1) = -3
    INFO(2) = LFACT
    INFO(17) = MAX(INFO(17), ALENB + EXPNE + 1)
    IF (ICNTL(15).EQ.1)
   *    INFO(17) = MAX(INFO(17), ALENB + 3*EXPNE + 3*N)
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
   +    '**** Error return from ME57BD ****  INFO(1) =',INFO(1),
   +    'Insufficient real space in FACT, LFACT = ',INFO(2)
    RETURN

   95 INFO(1) = -4
    INFO(2) = LIFACT
    INFO(18) = MAX(INFO(18), EXPNE+N+5)
    IF (ICNTL(15).EQ.1)
   *    INFO(18) = MAX(INFO(18), 3*EXPNE + 5*N + 1)
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
   +    '**** Error return from ME57BD ****  INFO(1) =',INFO(1),
   +    'Insufficient integer space in IFACT, LIFACT = ',INFO(2)
    RETURN

C****************
C Printing section
C****************
 100  IF (LDIAG.LE.2 .OR. MP.LT.0) RETURN
    WRITE (MP,99980) INFO(1), INFO(2),
   *    (INFO(I),I=14,25),INFO(28),INFO(29)
    WRITE (MP,99984) (INFO(I),I=31,35),RINFO(3), RINFO(4),
   *           RINFO(5), RINFO(18)
99980 FORMAT (/'Leaving factorization phase (ME57BD) with ...'/
   1  'INFO (1)                        =',I12/
   2  ' --- (2)                        =',I12/
   3  ' --- (14) Number of entries in factors      =',I12/
   4  ' --- (15) Real storage for factors        =',I12/
   5  ' --- (16) Integer storage for factors       =',I12/
   6  ' --- (17) Min LFACT with compresses       =',I12/
   7  ' --- (18) Min LIFACT with compresses      =',I12/
   8  ' --- (19) Min LFACT without compresses      =',I12/
   9  ' --- (20) Min LIFACT without compresses     =',I12/
   *  ' --- (21) Order of largest frontal matrix   =',I12/
   1  ' --- (22) Number of 2x2 pivots          =',I12/
   2  ' --- (23) Number of delayed pivots        =',I12/
   3  ' --- (24) Number of negative eigenvalues    =',I12/
   4  ' --- (25) Rank of factorization         =',I12/
   5  ' --- (28) Number compresses on real data    =',I12/
   6  ' --- (29) Number compresses on integer data   =',I12)
    IF (ICNTL(15).EQ.1) WRITE (MP,99985) RINFO(16),RINFO(17)
99985 FORMAT (
   1  'RINFO(16) Minimum value of scaling factor   =  ',1PD10.3/
   2  '-----(17) Maximum value of scaling factor   =  ',1PD10.3)
99984 FORMAT (
   7  ' --- (31) Number of block pivots in factors   =',I12/
   7  ' --- (32) Number of zeros factors triangle    =',I12/
   7  ' --- (33) Number of zeros factors rectangle   =',I12/
   7  ' --- (34) Number of zero cols factors rect    =',I12/
   7  ' --- (35) Number of static pivots         =',I12/
   1  'RINFO(3)  Operations during node assembly   =  ',1PD10.3/
   2  '-----(4)  Operations during node elimination  =  ',1PD10.3/
   3  '-----(5)  Extra operations because of BLAS    =  ',1PD10.3/
   3  '-----(18) Largest modulus of entry in matrix  =  ',1PD10.3)
    IF (INFO(27).GT.0) WRITE (MP,99981) INFO(27),RINFO(14),RINFO(15)
99981 FORMAT (/'Matrix modification performed'/
   1  'INFO (27) Step at which matrix first modified =',I12/
   2  'RINFO(14) Maximum value added to diagonal   =  ',1PD10.3/
   2  'RINFO(15) Smallest pivot in modified matrix   =  ',1PD10.3)
C Print out matrix factors from ME57BD.
    CALL ME57UD(FACT,LFACT,IFACT,LIFACT,ICNTL)
C Print scaling factors
    IF (ICNTL(15).NE.1) RETURN
    K = MIN(10,N)
    IF (LDIAG.GE.4) K = N
    WRITE (MP,'(/A/(5D12.5))')  'Scaling factors:',
   +            (REAL(FACT(ISCALE+I-1)),I=1,K)
    IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'

    END
    SUBROUTINE ME57CD(KIND,JOB,N,FACT,LFACT,IFACT,LIFACT,
   *            NRHS,RHS,LRHS,W,
   *            LW,IW1,ICNTL,INFO)
C This subroutine uses the factorisation of the matrix in FACT,IFACT to
C   solve a system of equations.
    INTEGER KIND,JOB,N,LFACT
    COMPLEX*16 FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
    COMPLEX*16 W(LW),RHS(LRHS,NRHS)
    INTEGER IW1(N),ICNTL(20),INFO(40)
CCC
C   LOGICAL LCASE
C   COMMON /CCASE/LCASE
C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C JOB must be set by the user to determine the coefficient matrix
C   of the equations being solved.  If the factorization is
C                  T  T
C         A =  P  L  D  L  P
C   then coefficient matrix is:
C
C   JOB <= 1   A
C               T
C   JOB  = 2   P  L  P
C               T
C   JOB  = 3   P  D  P
C            T  T
C   JOB >= 4   P  L  P
C
C N must be set to the order of the matrix and must be unchanged since
C   the call to ME57BD. It is not altered.
C FACT holds information on the factors and must be unchanged since
C   the call to ME57BD. It is not altered by ME57CD.
C LFACT must be set to the length of FACT. It is not altered.
C IFACT holds information on the factors and must be unchanged since
C   the call to ME57BD. It is not altered by ME57CD.
C LIFACT must be set to the length of IFACT. It is not altered.
C NRHS is the number of right-hand sides being solved for.
C RHS must be set to the right hand sides for the equations being
C   solved. On exit, this array will hold the solutions.
C LHS must be set to the leading dimension of array RHS.
C W is used as a work array.
C LW  must be set to the length of array W.  It must be at least
C   as large as N*NRHS.  (Actually only INFO(21)*NRHS but no way to
C   check this).
C IW1 is used as a work array.
C ICNTL must be set by the user as follows and is not altered.
C   ICNTL(1)  must be set to the stream number for error messages.
C     A value less than 0 suppresses output.
C   ICNTL(2) must be set to the stream number for warning output.
C     A value less than 0 suppresses output.
C   ICNTL(3) must be set to the stream number for monitor output.
C     A value less than 0 suppresses output.
C   ICNTL(4) must be set to the stream number for statistics output.
C     A value less than 0 suppresses output.
C   ICNTL(5) must be set to control the amount of output:
C     0 None.
C     1 Error messages only.
C     2 Error and warning messages.
C     3 As 2, plus scalar parameters and a few entries of array
C       parameters on entry and exit.
C   > 3  As 2, plus all parameters on entry and exit.
C   ICNTL(6:12) Not referenced.
C   ICNTL(13) Threshold on number of columns in a block for direct
C       addressing using Level 2 and Level 3 BLAS.

C Procedures
    INTRINSIC MIN
    EXTERNAL ME57QD,ME57RD,ME57SD,ME57TD,ME57UD,ME57XD,ME57YD

C
C Local variables
    DOUBLE PRECISION SCALE,ONE
    PARAMETER (ONE = 1.0D0)
    INTEGER I,J,K,LDIAG,LLW,LP,MP,ISCALE
C I  Temporary variable.
C J  Temporary variable.
C K  Temporary variable.
C LDIAG Control for amount of information output.
C LP Stream number for error printing.
C MP Stream number for monitor printing.
C

C Set local print control variables
    LP = ICNTL(1)
    MP = ICNTL(3)
    LDIAG = ICNTL(5)

    INFO(1) = 0

C Check input data
    IF (N.LE.0) THEN
      INFO(1) = -1
      INFO(2) = N
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
   +    '**** Error return from ME57CD ****  INFO(1) =',INFO(1),
   +    'N has value',N
      GOTO 500
    ENDIF

    IF (NRHS.LT.1) THEN
      INFO(1) = -16
      INFO(2) = NRHS
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I10,A)')
   +    '**** Error return from ME57CD ****  INFO(1) =',INFO(1),
   +    'value of NRHS =',NRHS,' is less than 1'
      GOTO 500
    ENDIF

    IF (LRHS.LT.N) THEN
      INFO(1) = -11
      INFO(2) = LRHS
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I10,A,I10)')
   +    '**** Error return from ME57CD ****  INFO(1) =',INFO(1),
   +    'value of LRHS =',LRHS,' is less than N=',N
      GOTO 500
    ENDIF

    IF (LW.LT.N*NRHS) THEN
      INFO(1) = -17
      INFO(2) = N*NRHS
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I10,A,I10)')
   +    '**** Error return from ME57CD ****  INFO(1) =',INFO(1),
   +    'value of LW =',LW,' is less than', N*NRHS
      GOTO 500
    ENDIF

C If requested, print input parameters
    IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
      WRITE (MP,'(/A/)') 'Entering solution phase (ME57CD) with ...'
      IF (KIND.EQ.1) THEN
      WRITE(MP,'(A)') 'Matrix is complex symmetric'
      ELSE
      WRITE(MP,'(A)') 'Matrix is Hermitian'
      ENDIF
      WRITE (MP,99999) JOB,N,(ICNTL(I),I=1,5),LFACT,LIFACT,NRHS,
   +       LRHS,LW,ICNTL(13)
99999 FORMAT(
   +    'JOB     Control on coefficient matrix     =',I12/
   +    'N       Order of matrix             =',I12/
   6    'ICNTL(1)  Stream for errors             =',I12/
   7    ' --- (2)  Stream for warnings           =',I12/
   8    ' --- (3)  Stream for monitoring         =',I12/
   9    ' --- (4)  Stream for statistics         =',I12/
   1    ' --- (5)  Level of diagnostic printing      =',I12/
   +    'LFACT   Length of array FACT          =',I12/
   +    'LIFACT    Length of array IFACT         =',I12/
   +    'NRHS    Number of right-hand sides      =',I12/
   +    'LRHS    Leading dimension of RHS array    =',I12/
   +    'LW      Leading dimension of work array   =',I12/
   +    'ICNTL(13) Threshold for Level 2 and 3 BLAS    =',I12)

C Print out matrix factors.
      CALL ME57UD(FACT,LFACT,IFACT,LIFACT,ICNTL)
C Print scaling factors
      IF (ICNTL(15).EQ.1) THEN
      ISCALE = LFACT-N
      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(/A/(5D12.5))')  'Scaling factors:',
   +                (REAL(FACT(ISCALE+I-1)),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
      ENDIF
      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      DO 10 J = 1,NRHS
      WRITE(MP,'(/A,I10)') 'Right-hand side',J
      WRITE (MP,'((1P,6D13.3))') (RHS(I,J),I=1,K)
      IF (K.LT.N) WRITE (MP,'(A)') '   . . .'
   10   CONTINUE
    END IF

    LLW = LW/NRHS


C Scale right-hand side
    IF (ICNTL(15).EQ.1) THEN
      ISCALE = LFACT-N
C   write(7,*) 'N,NRHS,ISCALE', N,NRHS,ISCALE
      DO 5555 I = 1, N
      SCALE = REAL(FACT(ISCALE+I-1))
C   write(7,*) 'SCALE',SCALE
      IF (JOB.GE.4) SCALE = ONE/SCALE
      DO 4444 J = 1, NRHS
C   write(7,*) 'I,J,RHS', I,J,RHS(I,J)
        RHS(I,J) = SCALE*RHS(I,J)
C   write(7,*) 'I,J,RHS', I,J,RHS(I,J)
 4444   CONTINUE
 5555   CONTINUE
    ENDIF

C Forward substitution
    IF (JOB.LE.2) THEN
      IF (NRHS.EQ.1) THEN
      CALL ME57XD(N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
   *          W,LLW,IW1,ICNTL)
CCC
C     IF (LCASE) THEN
C       DO 201 J = 1,NRHS
C       WRITE(7,'(/A,I10)') 'Solution     ',J
C       WRITE (7,'(1P,6D13.3)') (RHS(I,J),I=1,N)
C 201   CONTINUE
C     ENDIF
      ELSE
      CALL ME57QD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
   *          W,LLW,IW1,ICNTL)
CCC
C     IF (LCASE) THEN
C       DO 202 J = 1,NRHS
C       WRITE(7,'(/A,I10)') 'Solution     ',J
C       WRITE (7,'(1P,6D13.3)') (RHS(I,J),I=1,N)
C 202   CONTINUE
C     ENDIF
      ENDIF
      IF (JOB.EQ.2) GO TO 15
C Back substitution.
      IF (NRHS.EQ.1) THEN
      CALL ME57YD(KIND,N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
   *          W,LLW,IW1,ICNTL)
      ELSE
      CALL ME57RD(KIND,N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
   *          W,LLW,IW1,ICNTL)
      ENDIF
    ENDIF
    IF (JOB.EQ.3)
   *  CALL ME57SD(KIND,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
   *          W,LLW,ICNTL)
    IF (JOB.GE.4)
   *  CALL ME57TD(KIND,N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
   *          W,LLW,IW1,ICNTL)

C Scale solution
   15 IF (ICNTL(15).EQ.1) THEN
      ISCALE = LFACT-N
      DO 6666 I = 1, N
      SCALE = REAL(FACT(ISCALE+I-1))
      IF (JOB.EQ.2) SCALE = ONE/SCALE
      DO 7777 J = 1, NRHS
        RHS(I,J) = SCALE*RHS(I,J)
 7777   CONTINUE
 6666   CONTINUE
    ENDIF

C
C If requested, print output parameters.
    IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
      WRITE (MP,'(//A)')
   *     'Leaving solution phase (ME57CD) with ...'
      DO 20 J = 1,NRHS
      WRITE(MP,'(/A,I10)') 'Solution     ',J
      WRITE (MP,'(1P,6D13.3)') (RHS(I,J),I=1,K)
      IF (K.LT.N) WRITE (MP,'(A)') '   . . .'
   20   CONTINUE
    ENDIF

  500 RETURN

    END


    SUBROUTINE ME57QD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
   *            W,LW,IW1,ICNTL)
C This subroutine performs forward elimination using the factors
C   stored in FACT/IFACT by ME57BD.
C   Because of how the factors are stored this subroutine is the
C   same for complex symmetric and Hermitian matrices.
    INTEGER N,LFACT
    COMPLEX*16 FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
    COMPLEX*16 W(LW,NRHS),RHS(LRHS,NRHS)
    INTEGER IW1(N),ICNTL(20)
C N   must be set to the order of the matrix. It is not altered.
C FACT   must be set to hold the real values corresponding to the
C   factors. This must be unchanged since the preceding call to
C   ME57BD. It is not altered.
C LFACT  length of array FACT. It is not altered.
C IFACT  holds the integer indexing information for the matrix factors
C   in FACT. This must be unchanged since the preceding call to
C   ME57BD. It is not altered.
C LIFACT length of array IFACT. It is not altered.
C NRHS must be set to number of right-hand sides.
C RHS on input, must be set to hold the right hand side vector.  On
C   return, it will hold the modified vector following forward
C   elimination.
C LHS must be set to the leading dimension of array RHS.
C W   used as workspace to hold the components of the right hand
C   sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C   larger than INFO(21) as returned from ME57BD.
C IW1 need not be set on entry. On exit IW1(I) (I = 1,NBLK), where
C   NBLK = IFACT(3) is the number of block pivots, will
C   hold pointers to the beginning of each block pivot in array IFACT.
C ICNTL Not referenced except:
C   ICNTL(13) Threshold on number of columns in a block for using
C       addressing using Level 2 and Level 3 BLAS.
C

C Procedures
    INTRINSIC ABS
    EXTERNAL ZGEMM,ZTPSV

C Constant
    DOUBLE PRECISION ONE
    PARAMETER (ONE=1.0D0)
    COMPLEX*16 CONE
    PARAMETER (CONE=(1.0D0,0.0D0))
C
C Local variables
    INTEGER APOS,I,IBLK,II,IPIV,IRHS,IWPOS,J,J1,J2,K,
   +      NCOLS,NROWS
    COMPLEX*16 W1
C
C APOS  Current position in array FACT.
C I   Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IWPOS Position in IFACT of start of current index list.
C J   Temporary DO index
C K   Temporary pointer to position in real array.
C J1    Position in IFACT of index of leading entry of row.
C J2    Position in IFACT of index of trailing entry of row.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.


    APOS = 1
    IWPOS = 4
    DO 270 IBLK = 1,IFACT(3)

C Find the number of rows and columns in the block.
      IW1(IBLK) = IWPOS
      NCOLS = IFACT(IWPOS)
      NROWS = IFACT(IWPOS+1)
      IWPOS = IWPOS + 2

      IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN


C Perform operations using direct addressing.

C Load appropriate components of right-hand sides into W.
      DO 10 I = 1,NCOLS
        II = ABS(IFACT(IWPOS+I-1))
        DO 11 J = 1,NRHS
          W(I,J) = RHS(II,J)
   11     CONTINUE
   10   CONTINUE


C Treat diagonal block (direct addressing)
      DO 12 J = 1,NRHS
        CALL ZTPSV('L','N','U',NROWS,FACT(APOS),W(1,J),1)
   12   CONTINUE
      APOS = APOS + (NROWS* (NROWS+1))/2

C Treat off-diagonal block (direct addressing)
      IF (NCOLS.GT.NROWS) CALL ZGEMM('N','N',NCOLS-NROWS,NRHS,NROWS,
   +                      CONE,FACT(APOS),NCOLS-NROWS,
   +                      W,LW,CONE,W(NROWS+1,1),LW)
      APOS = APOS + NROWS* (NCOLS-NROWS)

C Reload W back into RHS.
      DO 35 I = 1,NCOLS
        II = ABS(IFACT(IWPOS+I-1))
        DO 36 J = 1,NRHS
          RHS(II,J) = W(I,J)
   36     CONTINUE
   35   CONTINUE

      ELSE

C Perform operations using indirect addressing.

      J1 = IWPOS
      J2 = IWPOS + NROWS - 1


C Treat diagonal block (indirect addressing)
      DO 130 IPIV = 1,NROWS
      APOS = APOS + 1
      DO 101 II = 1,NRHS
        W1 = RHS(ABS(IFACT(J1)),II)
        K = APOS
        DO 100 J = J1+1,J2
          IRHS = ABS(IFACT(J))
          RHS(IRHS,II) = RHS(IRHS,II) - FACT(K)*W1
          K = K + 1
  100     CONTINUE
  101   CONTINUE
      APOS = K
      J1 = J1 + 1
  130   CONTINUE

C Treat off-diagonal block (indirect addressing)
      J2 = IWPOS + NCOLS - 1
      DO 136 IPIV = 1,NROWS
      DO 135 II = 1,NRHS
        K = APOS
        W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)),II)
        DO 133 J = J1,J2
          IRHS = ABS(IFACT(J))
          RHS(IRHS,II) = RHS(IRHS,II) + W1*FACT(K)
          K = K + 1
  133     CONTINUE
  135   CONTINUE
      APOS = K
  136   CONTINUE

    END IF

    IWPOS = IWPOS + NCOLS
  270 CONTINUE

    END


    SUBROUTINE ME57RD(KIND,N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
   *            W,LW,IW1,ICNTL)
C This subroutine performs backward elimination operations
C   using the factors stored in FACT/IFACT by ME57BD.
    INTEGER KIND,N,LFACT
    COMPLEX*16 FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
    COMPLEX*16 W(LW,NRHS),RHS(LRHS,NRHS)
    INTEGER IW1(N),ICNTL(20)
C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C N    must be set to the order of the matrix. It is not altered.
C FACT    must be set to hold the real values corresponding to the
C    factors. This must be unchanged since the
C    preceding call to ME57BD. It is not altered.
C LFACT   length of array FACT. It is not altered.
C IFACT   holds the integer indexing information for the matrix factors
C    in FACT. This must be unchanged since the preceding call to
C    ME57BD. It is not altered.
C LIFACT  length of array IFACT. It is not altered.
C NRHS must be set to number of right-hand sides.
C RHS  on entry, must be set to hold the right hand side modified by
C    the forward substitution operations. On exit, holds the
C    solution vector.
C LHS must be set to the leading dimension of array RHS.
C W    used as workspace to hold the components of the right hand
C    sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C   larger than INFO(21) as returned from ME57BD.
C IW1  on entry IW1(I) (I = 1,NBLK), where  NBLK = IFACT(3) is the
C    number of block pivots, must hold pointers to the beginning of
C    each block pivot in array IFACT, as set by ME57Q/QD.
C    It is not altered.
C ICNTL Not referenced except:
C   ICNTL(13) Threshold on number of columns in a block for using
C       addressing using Level 2 and Level 3 BLAS.

C Procedures
    INTRINSIC ABS
    EXTERNAL ZGEMM,ZTPSV

C Constants
    DOUBLE PRECISION ONE
    PARAMETER (ONE=1.0D0)
    COMPLEX*16 CONE
    PARAMETER (CONE=(1.0D0,0.0D0))
C
C Local variables.
    INTEGER APOS,APOS2,I,IBLK,II,IPIV,IRHS,IRHS1,
   +      IRHS2,IWPOS,J,JPIV,J1,J2,K,KK,LROW,NCOLS,NROWS
    COMPLEX*16 W1
C APOS  Current position in array FACT.
C APOS2 Current position in array FACT for off-diagonal entry of 2x2
C     pivot.
C I   Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IRHS1 RHS index.
C IRHS2 RHS index.
C IWPOS Position in IFACT of start of current index list.
C J   Temporary DO index
C JPIV  Has the value -1 for the first row of a 2 by 2 pivot and 1 for
C     the second.
C K   Temporary pointer to position in real array.
C J1    Position in IFACT of index of leading entry of row.
C J2    Position in IFACT of index of trailing entry of row.
C K   Temporary variable.
C LROW  Length of current row.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.
C

    APOS = IFACT(1)
    APOS2 = IFACT(2)
C Run through block pivot rows in the reverse order.
    DO 380 IBLK = IFACT(3),1,-1

C Find the number of rows and columns in the block.
      IWPOS = IW1(IBLK)
      NCOLS = ABS(IFACT(IWPOS))
      NROWS = ABS(IFACT(IWPOS+1))
      APOS = APOS - NROWS* (NCOLS-NROWS)
      IWPOS = IWPOS + 2

      IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN


C Perform operations using direct addressing.

C Load latter part of right-hand side into W.
      DO 5 I = NROWS + 1,NCOLS
        II = ABS(IFACT(IWPOS+I-1))
        DO 3 J = 1,NRHS
          W(I,J) = RHS(II,J)
    3     CONTINUE
    5   CONTINUE


C Multiply by the diagonal matrix (direct addressing)
      DO 10 IPIV = NROWS,1,-1
        IRHS = ABS(IFACT(IWPOS+IPIV-1))
        APOS = APOS - (NROWS+1-IPIV)
        DO 9 J = 1,NRHS
          W(IPIV,J) = RHS(IRHS,J)*FACT(APOS)
    9     CONTINUE
   10   CONTINUE
      JPIV = -1
      DO 20 IPIV = NROWS,1,-1
        IRHS = IFACT(IWPOS+IPIV-1)
        IF (IRHS.LT.0) THEN
          IRHS1 = -IFACT(IWPOS+IPIV-1+JPIV)
          IF (KIND.EQ.1 .OR. JPIV.EQ.1) THEN
          DO 19 J = 1,NRHS
            W(IPIV,J) = RHS(IRHS1,J)*FACT(APOS2) + W(IPIV,J)
   19       CONTINUE
          ELSE
          DO 21 J = 1,NRHS
            W(IPIV,J) = RHS(IRHS1,J)*CONJG(FACT(APOS2))
   *                + W(IPIV,J)
   21       CONTINUE
          ENDIF
          IF (JPIV.EQ.1) APOS2 = APOS2 - 1
          JPIV = -JPIV
        END IF

   20   CONTINUE

C Treat off-diagonal block (direct addressing)
      K = NCOLS - NROWS
      IF (K.GT.0) THEN
        IF (KIND.EQ.1) CALL ZGEMM('T','N',NROWS,NRHS,K,CONE,
   +                 FACT(APOS+(NROWS*(NROWS+1))/2),K,
   +                 W(NROWS+1,1),LW,CONE,W,LW)
        IF (KIND.NE.1) CALL ZGEMM('C','N',NROWS,NRHS,K,CONE,
   +                 FACT(APOS+(NROWS*(NROWS+1))/2),K,
   +                 W(NROWS+1,1),LW,CONE,W,LW)
      ENDIF

C Treat diagonal block (direct addressing)
      IF (KIND.EQ.1) THEN
        DO 22 J = 1,NRHS
          CALL ZTPSV('L','T','U',NROWS,FACT(APOS),W(1,J),1)
   22     CONTINUE
      ELSE
        DO 23 J = 1,NRHS
          CALL ZTPSV('L','C','U',NROWS,FACT(APOS),W(1,J),1)
   23     CONTINUE
      ENDIF
C Reload W back into RHS.
      DO 60 I = 1,NROWS
        II = ABS(IFACT(IWPOS+I-1))
        DO 59 J = 1,NRHS
          RHS(II,J) = W(I,J)
   59     CONTINUE
   60   CONTINUE

      ELSE
C
C Perform operations using indirect addressing.
      J1 = IWPOS
      J2 = IWPOS + NCOLS - 1


C Multiply by the diagonal matrix (indirect addressing)
      JPIV = -1
      IF (KIND.EQ.1) THEN
        DO 210 IPIV = NROWS,1,-1
          IRHS = IFACT(IWPOS+IPIV-1)
          LROW = NROWS + 1 - IPIV
          IF (IRHS.GT.0) THEN
C 1 by 1 pivot.
          APOS = APOS - LROW
          DO 65 J = 1,NRHS
            RHS(IRHS,J) = RHS(IRHS,J)*FACT(APOS)
   65       CONTINUE
          ELSE
C 2 by 2 pivot
          IF (JPIV.EQ.-1) THEN
            IRHS1 = -IFACT(IWPOS+IPIV-2)
            IRHS2 = -IRHS
            APOS = APOS - LROW - LROW - 1
            DO 68 J = 1,NRHS
            W1 = RHS(IRHS1,J)*FACT(APOS) +
   +             RHS(IRHS2,J)*FACT(APOS2)
            RHS(IRHS2,J) = RHS(IRHS1,J)*FACT(APOS2) +
   +                   RHS(IRHS2,J)*FACT(APOS+LROW+1)
            RHS(IRHS1,J) = W1
   68         CONTINUE
            APOS2 = APOS2 - 1
          END IF
          JPIV = -JPIV
          END IF
  210     CONTINUE
      ELSE
        DO 211 IPIV = NROWS,1,-1
          IRHS = IFACT(IWPOS+IPIV-1)
          LROW = NROWS + 1 - IPIV
          IF (IRHS.GT.0) THEN
C 1 by 1 pivot.
          APOS = APOS - LROW
          DO 66 J = 1,NRHS
            RHS(IRHS,J) = RHS(IRHS,J)*FACT(APOS)
   66       CONTINUE
          ELSE
C 2 by 2 pivot
          IF (JPIV.EQ.-1) THEN
            IRHS1 = -IFACT(IWPOS+IPIV-2)
            IRHS2 = -IRHS
            APOS = APOS - LROW - LROW - 1
            DO 69 J = 1,NRHS
            W1 = RHS(IRHS1,J)*FACT(APOS) +
   +             RHS(IRHS2,J)*FACT(APOS2)
            RHS(IRHS2,J) = RHS(IRHS1,J)*CONJG(FACT(APOS2)) +
   +                   RHS(IRHS2,J)*FACT(APOS+LROW+1)
            RHS(IRHS1,J) = W1
   69         CONTINUE
            APOS2 = APOS2 - 1
          END IF
          JPIV = -JPIV
          END IF
  211     CONTINUE
      ENDIF
      APOS = APOS + (NROWS* (NROWS+1))/2

C Treat off-diagonal block (indirect addressing)
      KK = APOS
      J1 = IWPOS + NROWS
      IF (KIND.EQ.1) THEN
        DO 220 IPIV = 1,NROWS
          IRHS = ABS(IFACT(IWPOS+IPIV-1))
          DO 218 II = 1,NRHS
          W1 = RHS(IRHS,II)
          K = KK
          DO 215 J = J1,J2
            W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)),II)
            K = K + 1
  215       CONTINUE
          RHS(IRHS,II) = W1
  218       CONTINUE
          KK = K
  220     CONTINUE
      ELSE
        DO 221 IPIV = 1,NROWS
          IRHS = ABS(IFACT(IWPOS+IPIV-1))
          DO 219 II = 1,NRHS
          W1 = RHS(IRHS,II)
          K = KK
          DO 216 J = J1,J2
            W1 = W1 + CONJG(FACT(K))*RHS(ABS(IFACT(J)),II)
            K = K + 1
  216       CONTINUE
          RHS(IRHS,II) = W1
  219       CONTINUE
          KK = K
  221     CONTINUE
      ENDIF

C Treat diagonal block (indirect addressing)
      J2 = IWPOS + NROWS - 1
      DO 260 IPIV = 1,NROWS
        IRHS = ABS(IFACT(J1-1))
        APOS = APOS - IPIV
        DO 240 II = 1,NRHS
          W1 = RHS(IRHS,II)
          K = APOS + 1
          DO 230 J = J1,J2
          IF (KIND.EQ.1) THEN
            W1 = W1 - FACT(K)*RHS(ABS(IFACT(J)),II)
          ELSE
            W1 = W1 - CONJG(FACT(K))*RHS(ABS(IFACT(J)),II)
          ENDIF
          K = K + 1
  230       CONTINUE
          RHS(IRHS,II) = W1
  240     CONTINUE
        J1 = J1 - 1
  260   CONTINUE

      END IF

  380 CONTINUE

    END

    SUBROUTINE ME57UD(FACT,LFACT,IFACT,LIFACT,ICNTL)
C   Print out matrix factors from ME57BD or a symbolic representation
C   of them.
    INTEGER LFACT
    COMPLEX*16 FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT),ICNTL(20)
C FACT   array holding the reals of the factorization.
C      It is not altered.
C LFACT  length of array FACT. It is not altered.
C IFACT  array holding the integers of the factorization. It is not
C   altered.
C LIFACT length of array IFACT. It is not altered.
C ICNTL is not referenced except:
C   ICNTL(3) must be set to the stream number for diagnostic output.
C     A value less than 1 suppresses output.
C   ICNTL(5) must be set to control the amount of output:
C    <3 None.
C     3 First block only.
C     4 All blocks.
C     5 All blocks, but each entry represented by a single character:
C        + for a positive integer
C        - for a negative integer
C        * for a nonzero entry
C        . for a zero entry

C Procedures
    INTRINSIC MIN,SIGN

C Local variables
    COMPLEX*16 CZERO
    PARAMETER (CZERO = (0.0D0,0.0D0))
    CHARACTER*80 LINE
    INTEGER APOS,APOS2,IBLK,ILINE,IROW,IWPOS,J,JPIV,J1,J2,K,
   +      LDIAG,LEN,MP,NBLK,NCOLS,NROWS
C APOS Current position in FACT.
C APOS2 Position in FACT of next off-diagonal entry of 2x2 pivot.
C ILINE Current position in the line.
C IBLK  Current block.
C IROW  Current row.
C IWPOS Current position in IFACT.
C JPIV  has value 1 only for the first row of a 2x2 pivot.
C J   Column index.
C K   Temporary pointer to position in real array.
C J1    Position of last zero in leading part of row.
C J2    Position of last nonzero in leading part of row.
C K   Temporary DO index.
C LDIAG Control for diagnostic printing.
C LEN   1 if pattern only to be printed and 12 if values to be printed.
C LINE  Character variable in which an output line is built.
C MP    Stream number for warning or diagnostic messages
C NBLK  Number of blocks to be printed.
C NCOLS Number of columns in the block.
C NROWS Number of rows in the block.

    CHARACTER*1 PM(-2:2)
    DATA PM/'*','-','.','+','.'/
    DOUBLE PRECISION ZERO
    PARAMETER (ZERO=0.0D0)

C Initialize MP and LDIAG
    MP = ICNTL(3)
    LDIAG = ICNTL(5)

    APOS2 = IFACT(1)
    NBLK = IFACT(3)
    IF (LDIAG.EQ.3) NBLK = MIN(1,NBLK)
    LEN = 20
    IF (LDIAG.EQ.5) LEN = 1
    IF (LEN.EQ.20) THEN
      IF (NBLK.EQ.IFACT(3)) THEN
      WRITE (MP,'(/A)')
   +    'For each block, the following information is provided:'

      ELSE
      WRITE (MP,'(/A,A)') 'For the first block only,',
   +    ' the following information is provided:'
      END IF

    END IF

    IF (LEN.EQ.20) WRITE (MP,'(A)')
   +    '   1. Block number, number of rows, number of columns',
   +    '   2. List of indices for the pivot, each negated if part of'
   +    ,'    a 2x2 pivot',
   +    '   3. The factorized block pivot',
   +    '    It has the form',
   +    '        -1  T',
   +    '      L  D   L ',
   +    '               -1    T',
   +    '    and is printed as D and L  packed together.',
   +    '   4. List of indices for the non-pivot columns',
   +    '   5. The non-pivot part as rectangular block by rows'

    IWPOS = 4
    APOS = 1

    DO 300 IBLK = 1,NBLK
      NCOLS = IFACT(IWPOS)
      NROWS = IFACT(IWPOS+1)
      IWPOS = IWPOS + 2

      WRITE (MP,'(/4(A,I6))') 'Block pivot',IBLK,' with',NROWS,
   +      ' rows and', NCOLS,' columns'
C      IF (LEN.EQ.20) WRITE (MP,'(6I12)')
C   +               (IFACT(K),K=IWPOS,IWPOS+NROWS-1)
C      IF (LEN.EQ.1) WRITE (MP,'(80A1)') (PM(SIGN(1,IFACT(K))),
C   +    K=IWPOS,IWPOS+NROWS-1)

      JPIV = 0
      DO 30 IROW = 1,NROWS
      IF (JPIV.EQ.1) THEN
        JPIV = 0
      ELSE
        IF (IFACT(IWPOS+IROW-1).LT.0) JPIV = 1
      END IF

      ILINE = 1
      DO 10 J = 1,IROW - 1
        WRITE (LINE(ILINE:ILINE+LEN-1),'(A)') ' '
        ILINE = ILINE + LEN
        IF (ILINE.GT.80) THEN
          WRITE (MP,'(A)') LINE
          ILINE = 1
        END IF
   10   CONTINUE

      DO 20 J = IROW,NROWS
        IF (LEN.EQ.20) WRITE (LINE(ILINE:ILINE+19),
   +      '(1P,2D10.2)') FACT(APOS)
        IF (LEN.EQ.1) THEN
         IF (FACT(APOS).EQ.CZERO) THEN
            WRITE (LINE(ILINE:ILINE),'(A)') '.'
         ELSE
            WRITE (LINE(ILINE:ILINE),'(A)') '*'
         END IF
        END IF
        APOS = APOS + 1
        IF (J.EQ.IROW+1) THEN
          IF (JPIV.EQ.1) THEN
          IF (LEN.EQ.20) WRITE (LINE(ILINE:ILINE+19),
   +          '(1P,2D10.2)') FACT(APOS2)
          IF (LEN.EQ.1) THEN
            IF (FACT(APOS2).EQ.CZERO) THEN
               WRITE (LINE(ILINE:ILINE),'(A)') '.'
            ELSE
               WRITE (LINE(ILINE:ILINE),'(A)') '*'
            END IF
          END IF
          APOS2 = APOS2 + 1
          END IF
        END IF
        ILINE = ILINE + LEN
        IF (ILINE.GT.80) THEN
          WRITE (MP,'(A)') LINE
          ILINE = 1
        END IF
   20   CONTINUE

      IF (ILINE.GT.1) THEN
        LINE(ILINE:) = ' '
        WRITE (MP,'(A)') LINE
      END IF

   30   CONTINUE

      IWPOS = IWPOS + NROWS
      IF (LEN.EQ.20) WRITE (MP,'(5I16)') (IFACT(K),K=IWPOS,
   +    IWPOS+NCOLS-NROWS-1)
      IF (LEN.EQ.1) WRITE (MP,'(80A1)') (PM(SIGN(1,IFACT(K))),
   +    K=IWPOS,IWPOS+NCOLS-NROWS-1)

      IWPOS = IWPOS + NCOLS - NROWS
      DO 280 IROW = 1,NROWS
      J1 = NROWS
      J2 = NCOLS
      ILINE = 1
      DO 110 J = J1 + 1,J2
        IF (LEN.EQ.20) WRITE (LINE(ILINE:ILINE+19),
   +      '(1P,2D10.2)') FACT(APOS)
        IF (LEN.EQ.1) THEN
         IF (FACT(APOS).EQ.CZERO) THEN
            WRITE (LINE(ILINE:ILINE),'(A)') '.'
         ELSE
            WRITE (LINE(ILINE:ILINE),'(A)') '*'
         END IF
        END IF
        APOS = APOS + 1
        ILINE = ILINE + LEN
        IF (ILINE.GT.80) THEN
          WRITE (MP,'(A)') LINE
          ILINE = 1
        END IF
  110   CONTINUE

      IF (ILINE.GT.1) THEN
        LINE(ILINE:) = ' '
        WRITE (MP,'(A)') LINE
      END IF

  280   CONTINUE
  300 CONTINUE
    END

    SUBROUTINE ME57SD(KIND,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
   *            W,LW,ICNTL)
C This subroutine divides a vector by the block diagonal matrix of
C   the matrix factors using factor entries stored in FACT/IFACT
C    by ME57BD.
    INTEGER KIND,LFACT
    COMPLEX*16 FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
    COMPLEX*16 W(LW,NRHS),RHS(LRHS,NRHS)
    INTEGER ICNTL(20)
C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C FACT    must be set to hold the real values corresponding to the
C    factors. This must be unchanged since the
C    preceding call to ME57BD. It is not altered.
C LFACT   length of array FACT. It is not altered.
C IFACT   holds the integer indexing information for the matrix factors
C    in FACT. This must be unchanged since the preceding call to
C    ME57BD. It is not altered.
C LIFACT  length of array IFACT. It is not altered.
C NRHS must be set to number of right-hand sides.
C RHS  on entry, must be set to hold the right hand side modified by
C    the forward substitution operations. On exit, holds the
C    solution vector.
C LHS must be set to the leading dimension of array RHS.
C W    used as workspace to hold the components of the right hand
C    sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C   larger than INFO(21) as returned from ME57BD.
C ICNTL Not referenced except:
C   ICNTL(13) Threshold on number of columns in a block for direct
C       addressing using Level 2 and Level 3 BLAS.

C Procedures
    INTRINSIC ABS
    EXTERNAL ZGEMM,ZTPSV

C
C Local variables.
    INTEGER APOS,APOS2,I,IBLK,II,IPIV,IRHS,IRHS1,
   +      IRHS2,IWPOS,J,JPIV,NCOLS,NROWS
    COMPLEX*16 W1
C APOS  Current position in array FACT.
C APOS2 Current position in array FACT for off-diagonal entry of 2x2
C     pivot.
C I   Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IRHS1 RHS index.
C IRHS2 RHS index.
C IWPOS Position in IFACT of start of current index list.
C J   Temporary DO index
C JPIV  Has the value 1 for the first row of a 2 by 2 pivot and -1 for
C     the second.
C K   Temporary pointer to position in real array.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.
C
    APOS = 1
    APOS2 = IFACT(1)
    IWPOS = 4
C Run through block pivot rows in the reverse order.
    DO 380 IBLK = 1,IFACT(3)

C Find the number of rows and columns in the block.
      NCOLS = IFACT(IWPOS)
      NROWS = IFACT(IWPOS+1)
      IWPOS = IWPOS + 2


      IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN

C Perform operations using direct addressing.

C Multiply by the diagonal matrix (direct addressing)
      DO 10 IPIV = 1,NROWS
        IRHS = ABS(IFACT(IWPOS+IPIV-1))
        DO 9 J = 1,NRHS
          W(IPIV,J) = RHS(IRHS,J)*FACT(APOS)
    9     CONTINUE
        APOS = APOS + (NROWS+1-IPIV)
   10   CONTINUE
      JPIV = 1
      DO 20 IPIV = 1,NROWS
        IRHS = IFACT(IWPOS+IPIV-1)
        IF (IRHS.LT.0) THEN
          IRHS1 = -IFACT(IWPOS+IPIV-1+JPIV)
          IF (KIND.EQ.1 .OR. JPIV.EQ.1) THEN
          DO 18 J = 1,NRHS
            W(IPIV,J) = RHS(IRHS1,J)*FACT(APOS2) + W(IPIV,J)
   18       CONTINUE
          ELSE
          DO 19 J = 1,NRHS
            W(IPIV,J) =
   *          RHS(IRHS1,J)*CONJG(FACT(APOS2)) + W(IPIV,J)
   19       CONTINUE
          ENDIF
          IF (JPIV.EQ.-1) APOS2 = APOS2 + 1
          JPIV = -JPIV
        END IF

   20   CONTINUE

C Reload W back into RHS.
      DO 60 I = 1,NROWS
        II = ABS(IFACT(IWPOS+I-1))
        DO 59 J = 1,NRHS
          RHS(II,J) = W(I,J)
   59     CONTINUE
   60   CONTINUE

      ELSE
C
C Perform operations using indirect addressing.

C Multiply by the diagonal matrix (indirect addressing)
      JPIV = 1
      IF (KIND.EQ.1) THEN
        DO 210 IPIV = 1,NROWS
          IRHS = IFACT(IWPOS+IPIV-1)

          IF (IRHS.GT.0) THEN
C 1 by 1 pivot.
          DO 65 J = 1,NRHS
            RHS(IRHS,J) = RHS(IRHS,J)*FACT(APOS)
   65       CONTINUE
          APOS = APOS + NROWS - IPIV + 1
          ELSE
C 2 by 2 pivot
          IF (JPIV.EQ.1) THEN
            IRHS1 = -IRHS
            IRHS2 = -IFACT(IWPOS+IPIV)
            DO 68 J = 1,NRHS
            W1 = RHS(IRHS1,J)*FACT(APOS) +
   +             RHS(IRHS2,J)*FACT(APOS2)
            RHS(IRHS2,J) = RHS(IRHS1,J)*FACT(APOS2) +
   +                   RHS(IRHS2,J)*FACT(APOS+NROWS-IPIV+1)
            RHS(IRHS1,J) = W1
   68         CONTINUE
            APOS2 = APOS2 + 1
          END IF
          JPIV = -JPIV
          APOS = APOS + NROWS - IPIV + 1
          END IF
  210     CONTINUE
      ELSE
        DO 211 IPIV = 1,NROWS
          IRHS = IFACT(IWPOS+IPIV-1)

          IF (IRHS.GT.0) THEN
C 1 by 1 pivot.
          DO 66 J = 1,NRHS
            RHS(IRHS,J) = RHS(IRHS,J)*FACT(APOS)
   66       CONTINUE
          APOS = APOS + NROWS - IPIV + 1
          ELSE
C 2 by 2 pivot
          IF (JPIV.EQ.1) THEN
            IRHS1 = -IRHS
            IRHS2 = -IFACT(IWPOS+IPIV)
            DO 69 J = 1,NRHS
            W1 = RHS(IRHS1,J)*FACT(APOS) +
   +             RHS(IRHS2,J)*FACT(APOS2)
            RHS(IRHS2,J) = RHS(IRHS1,J)*CONJG(FACT(APOS2)) +
   +                   RHS(IRHS2,J)*FACT(APOS+NROWS-IPIV+1)
            RHS(IRHS1,J) = W1
   69         CONTINUE
            APOS2 = APOS2 + 1
          END IF
          JPIV = -JPIV
          APOS = APOS + NROWS - IPIV + 1
          END IF
  211     CONTINUE
      ENDIF

      END IF

      IWPOS = IWPOS + NCOLS
      APOS = APOS + NROWS*(NCOLS-NROWS)

  380 CONTINUE

    END

    SUBROUTINE ME57TD(KIND,N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
   *            W,LW,IW1,ICNTL)
C This subroutine performs backward elimination operations
C   using the factors stored in FACT/IFACT by ME57BD.
    INTEGER KIND,N,LFACT
    COMPLEX*16 FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
    COMPLEX*16 W(LW,NRHS),RHS(LRHS,NRHS)
    INTEGER IW1(N),ICNTL(20)
C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C N    must be set to the order of the matrix. It is not altered.
C FACT    must be set to hold the real values corresponding to the
C    factors. This must be unchanged since the
C    preceding call to ME57BD. It is not altered.
C LFACT   length of array FACT. It is not altered.
C IFACT   holds the integer indexing information for the matrix factors
C    in FACT. This must be unchanged since the preceding call to
C    ME57BD. It is not altered.
C LIFACT  length of array IFACT. It is not altered.
C NRHS must be set to number of right-hand sides.
C RHS  on entry, must be set to hold the right hand side modified by
C    the forward substitution operations. On exit, holds the
C    solution vector.
C LHS must be set to the leading dimension of array RHS.
C W    used as workspace to hold the components of the right hand
C    sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C   larger than INFO(21) as returned from ME57BD.
C IW1  on entry IW1(I) (I = 1,NBLK), where  NBLK = IFACT(3) is the
C    number of block pivots, must hold pointers to the beginning of
C    each block pivot in array IFACT, as set by ME57Q\QD.
C    It is not altered.
C ICNTL Not referenced except:
C   ICNTL(13) Threshold on number of columns in a block for direct
C       addressing using Level 2 and Level 3 BLAS.

C Procedures
    INTRINSIC ABS
    EXTERNAL ZGEMM,ZTPSV

C Constants
    DOUBLE PRECISION ONE
    PARAMETER (ONE=1.0D0)
    COMPLEX*16 CONE
    PARAMETER (CONE=(1.0D0,0.0D0))
C
C Local variables.
    INTEGER APOS,I,IBLK,II,IPIV,IRHS,
   +      IWPOS,J,J1,J2,K,KK,NCOLS,NROWS
    COMPLEX*16 W1
C APOS  Current position in array FACT.
C I   Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IRHS1 RHS index.
C IRHS2 RHS index.
C IWPOS Position in IFACT of start of current index list.
C J   Temporary DO index
C JPIV  Has the value -1 for the first row of a 2 by 2 pivot and 1 for
C     the second.
C K   Temporary pointer to position in real array.
C J1    Position in IFACT of index of leading entry of row.
C J2    Position in IFACT of index of trailing entry of row.
C K   Temporary variable.
C LROW  Length of current row.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.
C
    APOS = IFACT(1)


C Set IW1 array
    IWPOS = 4
    DO 10 I = 1,IFACT(3)-1
      IW1(I) = IWPOS
      IWPOS = IWPOS + ABS(IFACT(IWPOS))+2
   10 CONTINUE
    IW1(IFACT(3)) = IWPOS

C Run through block pivot rows in the reverse order.
    DO 380 IBLK = IFACT(3),1,-1

C Find the number of rows and columns in the block.
      IWPOS = IW1(IBLK)
      NCOLS = ABS(IFACT(IWPOS))
      NROWS = ABS(IFACT(IWPOS+1))
      APOS = APOS - NROWS* (NCOLS-NROWS)
      IWPOS = IWPOS + 2

      IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN


C Perform operations using direct addressing.

C Load right-hand side into W.
      DO 5 I = 1,NCOLS
        II = ABS(IFACT(IWPOS+I-1))
        DO 3 J = 1,NRHS
          W(I,J) = RHS(II,J)
    3     CONTINUE
    5   CONTINUE


C Treat off-diagonal block (direct addressing)
      K = NCOLS - NROWS
      IF (K.GT.0)  THEN
        IF (KIND.EQ.1) THEN
          CALL ZGEMM('T','N',NROWS,NRHS,K,CONE,
   +             FACT(APOS),K,
   +             W(NROWS+1,1),LW,CONE,W,LW)
        ELSE
          CALL ZGEMM('C','N',NROWS,NRHS,K,CONE,
   +             FACT(APOS),K,
   +             W(NROWS+1,1),LW,CONE,W,LW)
        ENDIF
      ENDIF

      APOS = APOS-(NROWS*(NROWS+1))/2

C Treat diagonal block (direct addressing)
      IF (KIND.EQ.1) THEN
        DO 22 J = 1,NRHS
          CALL ZTPSV('L','T','U',NROWS,FACT(APOS),W(1,J),1)
   22     CONTINUE
      ELSE
        DO 23 J = 1,NRHS
          CALL ZTPSV('L','C','U',NROWS,FACT(APOS),W(1,J),1)
   23     CONTINUE
      ENDIF

C Reload W back into RHS.
      DO 60 I = 1,NROWS
        II = ABS(IFACT(IWPOS+I-1))
        DO 59 J = 1,NRHS
          RHS(II,J) = W(I,J)
   59     CONTINUE
   60   CONTINUE

      ELSE
C
C Perform operations using indirect addressing.
      J1 = IWPOS
      J2 = IWPOS + NCOLS - 1

C Treat off-diagonal block (indirect addressing)
      KK = APOS
      J1 = IWPOS + NROWS
      IF (KIND.EQ.1) THEN
        DO 220 IPIV = 1,NROWS
          IRHS = ABS(IFACT(IWPOS+IPIV-1))
          DO 218 II = 1,NRHS
          W1 = RHS(IRHS,II)
          K = KK
          DO 215 J = J1,J2
            W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)),II)
            K = K + 1
  215       CONTINUE
          RHS(IRHS,II) = W1
  218       CONTINUE
          KK = K
  220     CONTINUE
      ELSE
        DO 221 IPIV = 1,NROWS
          IRHS = ABS(IFACT(IWPOS+IPIV-1))
          DO 219 II = 1,NRHS
          W1 = RHS(IRHS,II)
          K = KK
          DO 216 J = J1,J2
            W1 = W1 + CONJG(FACT(K))*RHS(ABS(IFACT(J)),II)
            K = K + 1
  216       CONTINUE
          RHS(IRHS,II) = W1
  219       CONTINUE
          KK = K
  221     CONTINUE
      ENDIF

C Treat diagonal block (indirect addressing)
      J2 = IWPOS + NROWS - 1
      DO 260 IPIV = 1,NROWS
        IRHS = ABS(IFACT(J1-1))
        APOS = APOS - IPIV
        DO 240 II = 1,NRHS
          W1 = RHS(IRHS,II)
          K = APOS + 1
          DO 230 J = J1,J2
          IF (KIND.EQ.1) THEN
            W1 = W1 - FACT(K)*RHS(ABS(IFACT(J)),II)
          ELSE
            W1 = W1 - CONJG(FACT(K))*RHS(ABS(IFACT(J)),II)
          ENDIF
          K = K + 1
  230       CONTINUE
          RHS(IRHS,II) = W1
  240     CONTINUE
        J1 = J1 - 1
  260   CONTINUE

      END IF

  380 CONTINUE


    END

    SUBROUTINE ME57DD(KIND,JOB,N,NE,A,IRN,JCN,FACT,LFACT,IFACT,LIFACT,
   *            RHS,X,RESID,WORKC,WORKR,IW,ICNTL,CNTL,INFO,
   *            RINFO)

C This subroutine solves a single system using one or more steps of
C   iterative refinement.
C If ICNTL(9) = 10 (the default), this subroutine performs iterative
C   refinement using the strategy of Arioli, Demmel, and Duff.
C   IF (ICNTL(9) = 1, then one step of iterative refinement is
C   performed.

    INTEGER KIND,JOB,N,NE
    COMPLEX*16 A(NE)
    INTEGER IRN(NE),JCN(NE),LFACT
    COMPLEX*16 FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT)
    COMPLEX*16 RHS(N),X(N),RESID(N),WORKC(N,*)
    DOUBLE PRECISION WORKR(N,*)
    INTEGER IW(N),ICNTL(20)
    DOUBLE PRECISION CNTL(5)
    INTEGER INFO(40)
    DOUBLE PRECISION RINFO(20)

C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C JOB must be set by the user to determine what action is desired by
C   the user.
C   Values of JOB and their effect are:
C IF ICNTL(9)>1, JOB=0 if no estimate of solution in X; JOB=2 if
C      estimate of solution in X.
C IF ICNTL(9)=1, then:
C   0: Solve Ax=b, calculate residual r=b-Ax and exit.
C      (Note that ME57CD should be used if solution without residual
C       is required)
C   1: Solve Ax=b, calculate residual r=b-Ax, solve A(dx)=r,
C      update solution and exit.
C If JOB > 1, an estimate of the solution must be input in X.
C   2: Calculate residual r=b-Ax, solve A(dx)=r,
C      update solution and exit.
C If JOB > 2, the residual for this estimate must also be input.
C   3: Solve A(dx)=r, update solution and exit.
C   4: Solve A(dx)=r, update solution, calculate residual for new
C      solution and exit.
C N must be set to the order of the matrix and must be unchanged since
C   the call to ME57BD. It is not altered by ME57DD.
C NE must be set to the number of entries in the matrix and must be
C   unchanged since the call to ME57AD. It is not altered by ME57DD.
C A must be set by the user to the values of the matrix as input to
C   ME57AD. It is not altered by ME57DD.
C IRN,JCN must be set to the row and column indices of the entries in A
C   and must be unchanged since the call to ME57AD. They are
C   not altered by ME57DD.
C FACT holds information on the factors and must be unchanged since
C   the call to ME57BD. It is not altered by ME57DD.
C LFACT must be set to the length of FACT. It is not altered by ME57DD.
C IFACT holds information on the factors and must be unchanged since
C   the call to ME57BD. It is not altered by ME57DD.
C LIFACT must be set to the length of IFACT. Not altered by ME57DD.
C RHS is a real array of length N that must be set to the right-hand
C   side for the equation being solved. It is not altered by ME57DD.
C X is a real array of length N. IF JOB >=2, it must be set on entry to
C   an estimated solution. Otherwise, it need not be set by the user.
C   On exit, the improved solution vector is returned in X.
C RESID is a real array of length N. If JOB > 2, it must be set on entry
C   to the value of the residual for the current solution estimate
C   held in X.  Otherwise, it need not be set by the user on entry.
C   If JOB=0 or 4 or if ICNTL(9)>1, on exit it will hold the residual
C   vector for the equations being solved. If 1<= JOB <= 3, then
C   RESID will hold on exit the last correction vector added to the
C   solution X.
C WORKC is used as a work array.  If ICNTL(9) = 1, it must be of length
C   at least N.  If ICNTL(9)>1 and ICNTL(10)>0, then WORKC must be of
C   length at least 2*N.
C WORKR is used as a work array. If ICNTL(9) = 1, it is not accessed by
C   the subroutine  If ICNTL(9) > 1, it must be of length at least
C   2*N.
C IW is an integer array of length N that is used as a work array if
C   ICNTL(1) > 9.  It is not accessed if ICNTL(9) = 1.
C ICNTL must be set by the user as follows and is not altered.
C   ICNTL(1)  must be set to the stream number for error messages.
C     A value less than 0 suppresses output.
C   ICNTL(2) must be set to the stream number for warning output.
C     A value less than 0 suppresses output.
C   ICNTL(3) must be set to the stream number for monitor output.
C     A value less than 0 suppresses output.
C   ICNTL(4) must be set to the stream number for statistics output.
C     A value less than 0 suppresses output.
C   ICNTL(5) must be set to control the amount of output:
C     0 None.
C     1 Error messages only.
C     2 Error and warning messages.
C     3 As 2, plus scalar parameters and a few entries of array
C       parameters on entry and exit.
C   > 3  As 2, plus all parameters on entry and exit.
C   ICNTL(6:12) Not referenced.
C   ICNTL(9)  Maximum permitted number of steps of iterative
C         refinement.
C   ICNTL(10) Flag to request calculation of error estimate and
C       condition numbers.
C   ICNTL(13) Threshold on number of columns in a block for direct
C       addressing using Level 2 and Level 3 BLAS.
C CNTL must be set by the user as follows and is not altered.
C   CNTL(3) is the required decrease in the scaled residuals required
C       by the Arioli, Demmel, and Duff iteration.
C INFO is an integer array that need not be set by the user.  On exit,
C   a value of INFO(1) equal to zero indicates success. A failure is
C   indicated by a negative value for INFO.
C RINFO is a real array that need not be set by the user. On exit,
C   If ICNTL(9)>1, RINFO is set to information on the matrix and
C   solution including backward errors.

C   .. Local constants ..
    DOUBLE PRECISION ZERO,ONE
    PARAMETER (ZERO=0.D0,ONE=1.0D0)

C   .. Local variables ..
    DOUBLE PRECISION COND(2),CTAU,DXMAX,ERROR,OLDOMG(2),OLDOM2,
   *           OMEGA(2),OM2,TAU
    INTEGER I,ICNTLC(20),ITER,J,K,KASE,KK,LDIAG,LP,MP,KEEP71(5)
    LOGICAL LCOND(2)

C
C COND is condition number of system.
C CTAU is set to 1000*machine precision.
C DXMAX used to calculate max norm of solution.
C ERROR used to accumulate error.
C OLDOMG holds value of previous backward errors.
C OLDOM2 holds previous sum of OMEGAs.
C OMEGA used to accumulate backward error.
C OM2 holds sum of OMEGAs.
C TAU is threshold for denominator in scaled residual calculation.
C I  Temporary variable.
C ICNTLC is control array for ME57C/CD.
C ITER is maximum permitted number of steps of iterative refinement.
C J  Temporary variable.
C K  Temporary variable.
C KASE used when calling MF71AD.
C KK Temporary variable.
C LDIAG Control for amount of information output.
C LCOND used as switch when calculating condition number.
C LP Stream number for error printing.
C MP Stream number for monitor printing.
C KEEP71 Work array required by MF71.
C

C Procedures
    INTRINSIC MIN
    EXTERNAL ME57CD,ME57UD,FD15AD,MF71AD
C EPS is the largest real such that 1+EPS is equal to 1.
    DOUBLE PRECISION EPS,FD15AD

    INFO(1) = 0

C Set local print control variables
    LP = ICNTL(1)
    MP = ICNTL(3)
    LDIAG = ICNTL(5)

C Check input data
    IF (N.LE.0) THEN
      INFO(1) = -1
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I12)')
   +    '**** Error return from ME57DD ****  INFO(1) =',INFO(1),
   +    'N has value',N
      INFO(2) = N
      GOTO 500
    ENDIF

    IF (NE.LT.0) THEN
      INFO(1) = -2
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I12)')
   +    '**** Error return from ME57DD ****  INFO(1) =',INFO(1),
   +    'NE has value',NE
      INFO(2) = NE
      GOTO 500
    ENDIF

    IF (ICNTL(9).LT.1) THEN
      INFO(1) = -13
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I12)')
   +    '**** Error return from ME57DD ****  INFO(1) =',INFO(1),
   +    'ICNTL(9) has value',ICNTL(9)
      INFO(2) = ICNTL(9)
      GOTO 500
    ENDIF

    IF (JOB.LT.0 .OR. JOB.GT.4 .OR. (ICNTL(9).GT.1 .AND.
   *    (JOB.NE.0 .AND. JOB.NE.2)))  THEN
      INFO(1) = -12
      INFO(2) = JOB
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I12)')
   +    '**** Error return from ME57DD ****  INFO(1) =',INFO(1),
   +    'JOB has value',JOB
      IF (ICNTL(9).GT.1 .AND. LDIAG.GT.0 .AND. LP.GE.0)
   +    WRITE (LP,'(A,I3)') 'and ICNTL(9) =',ICNTL(9)
      GOTO 500
    ENDIF

C If NE = 0, set variables and return
    IF (NE.EQ.0) THEN
      IF (JOB.NE.3) THEN
      DO 8 I = 1,N
        RESID(I) = ZERO
  8     CONTINUE
      ENDIF
      DO 9 I = 1,N
      X(I) = ZERO
  9   CONTINUE
      INFO(30)=0
      DO 10 I = 6,13
      RINFO(I) = ZERO
 10   CONTINUE
      GO TO 500
    ENDIF

C If requested, print input parameters
    IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
      WRITE (MP,'(/A,A/)')
   +      'Entering iterative refinement solution phase ',
   +      '(ME57DD) with ...'
      IF (KIND.EQ.1) THEN
      WRITE(MP,'(A)') 'Matrix is complex symmetric'
      ELSE
      WRITE(MP,'(A)') 'Matrix is Hermitian'
      ENDIF
      WRITE (MP,99999) JOB,N,NE,(ICNTL(I),I=1,5),LFACT,LIFACT,
   +   ICNTL(9),ICNTL(10),ICNTL(13),CNTL(3)
99999 FORMAT(
   +  'JOB     Control for coefficient matrix    =',I12/
   +  'N       Order of matrix             =',I12/
   +  'NE      Number of entries in matrix       =',I12/
   6  'ICNTL(1)  Stream for errors             =',I12/
   7  ' --- (2)  Stream for warnings           =',I12/
   8  ' --- (3)  Stream for monitoring         =',I12/
   9  ' --- (4)  Stream for statistics         =',I12/
   1  ' --- (5)  Level of diagnostic printing      =',I12/
   +  'LFACT   Length of array FACT          =',I12/
   +  'LIFACT    Length of array IFACT         =',I12/
   +  'ICNTL(9)  Number steps iterative refinement   =',I12/
   +  'ICNTL(10) Control for error analysis      =',I12/
   +  'ICNTL(13) Threshold for Level 2 and 3 BLAS    =',I12/
   +  'CNTL(3)   Convergence test for IR         =',1P,D12.4)

C Print out matrix factors.
      CALL ME57UD(FACT,LFACT,IFACT,LIFACT,ICNTL)
      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      WRITE(MP,'(/A)') 'Right-hand side'
      WRITE (MP,'((4X, 1P,6D13.3))') (RHS(I),I=1,K)
      IF (K.LT.N) WRITE (MP,'(A)') '   . . .'
    END IF


    DO 15 I=1,5
      ICNTLC(I) = ICNTL(I)
   15 CONTINUE
    ICNTLC(13) = ICNTL(13)
    ICNTLC(15) = ICNTL(15)
C Switch off monitor printing in ME57C/CD
    ICNTLC(3) = -1

    IF (JOB.LE.2) THEN
      IF (JOB .LE. 1) THEN
C No estimate given in X.
      DO 14 I = 1,N
        X(I) = RHS(I)
        RESID(I) = RHS(I)
   14   CONTINUE
C Solve system Ax=b using ME57C/CD
      CALL ME57CD(KIND,1,N,FACT,LFACT,IFACT,LIFACT,1,X,N,WORKC,N,
   +          IW,ICNTLC,INFO)
      ELSE
C Estimate of solution was input in X.
      DO 13 I = 1,N
        RESID(I) = RHS(I)
   13   CONTINUE
      ENDIF

      IF (ICNTL(9).EQ.1) THEN
C Compute residual
      DO 16 KK = 1,NE
        I = IRN(KK)
        J = JCN(KK)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 16
        RESID(I) = RESID(I) - A(KK)*X(J)
C Matrix is symmetric
        IF (I.NE.J) THEN
          IF (KIND.EQ.1) THEN
          RESID(J) = RESID(J) - A(KK)*X(I)
          ELSE
          RESID(J) = RESID(J) - CONJG(A(KK))*X(I)
          ENDIF
        ENDIF
   16   CONTINUE
      IF (JOB.EQ.0) GO TO 340
      ELSE
C Calculate values for iterative refinement strategy of Arioli,
C     Demmel and Duff.
      DO 18 I = 1,N
        WORKR(I,1) = ZERO
C Calculate infinity norms of rows of A in WORKR(I,1)
C Sum |a  |, j=1,N (= ||A  ||      )
C     ij         i.  infinity
        WORKR(I,2) = ZERO
   18   CONTINUE
      DO 17 KK = 1,NE
        I = IRN(KK)
        J = JCN(KK)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 17
        RESID(I) = RESID(I) - A(KK)*X(J)
        WORKR(J,1) = WORKR(J,1) + ABS(A(KK)*X(I))
        WORKR(J,2) = WORKR(J,2) + ABS(A(KK))
C Matrix is symmetric
        IF (I.NE.J) THEN
          IF (KIND.EQ.1) THEN
          RESID(J) = RESID(J) - A(KK)*X(I)
          ELSE
          RESID(J) = RESID(J) - CONJG(A(KK))*X(I)
          ENDIF
          WORKR(I,1) = WORKR(I,1) + ABS(A(KK)*X(J))
          WORKR(I,2) = WORKR(I,2) + ABS(A(KK))
        ENDIF
   17   CONTINUE

C Calculate max-norm of solution
      DXMAX = ZERO
      DO 221 I = 1,N
      DXMAX = MAX(DXMAX,ABS(X(I)))
  221   CONTINUE

C Initialize EPS
    EPS = FD15AD('E')
C Calculate backward errors OMEGA(1) and OMEGA(2)
C CTAU ... 1000 eps (approx)
      CTAU = 1000.*EPS
C tau is (||A  ||       ||x||   + |b| )*n*1000*epsilon
C        i.  infinity    max   i
      OMEGA(1) = ZERO
      OMEGA(2) = ZERO
      DO 231 I = 1,N
        TAU = (WORKR(I,2)*DXMAX+ABS(RHS(I)))*N*CTAU
        IF ((WORKR(I,1)+ABS(RHS(I))).GT.TAU) THEN
C |Ax-b| /(|A||x| + |b|)
C     i         i
          OMEGA(1) = MAX(OMEGA(1),ABS(RESID(I))/
   +             (WORKR(I,1)+ABS(RHS(I))))
          IW(I) = 1
        ELSE
C TAU will be zero if all zero row in A, for example
          IF (TAU.GT.ZERO) THEN
C |Ax-b| /(|A||x| + ||A  ||      ||x||   )
C     i      i   i.  infinity   max
          OMEGA(2) = MAX(OMEGA(2),ABS(RESID(I))/
   +             (WORKR(I,1)+WORKR(I,2)*DXMAX))
          END IF
          IW(I) = 2
        END IF
  231   CONTINUE
C
C  Stop the calculations if the backward error is small
C
      OM2 = OMEGA(1) + OMEGA(2)
      ITER = 0
C Statement changed because IBM SP held quantities in registers
C       IF ((OM2+ONE).LE.ONE) THEN
      IF (OM2.LE.EPS) THEN
C       write(7,'(A)') 'Go to 270 from after 231'
        GO TO 270
      ENDIF

C Hold current estimate in case needed later
      DO 251 I = 1,N
        WORKC(I,2) = X(I)
  251   CONTINUE
      OLDOMG(1) = OMEGA(1)
      OLDOMG(2) = OMEGA(2)
      OLDOM2 = OM2

      ENDIF

    ENDIF

C At this point JOB >= 1 or ICNTL(9) > 1
C
C Iterative refinement loop
    DO 260 ITER = 1,ICNTL(9)

C     write(7,'(A,I3/A/(5D16.8))') 'Iteration',ITER,
C    +  'Residual',(RESID(I),I=1,N)
C Solve system A(dx) = r using ME57C/CD
      CALL ME57CD(KIND,1,N,FACT,LFACT,IFACT,LIFACT,1,RESID,N,WORKC,N,
   +          IW,ICNTLC,INFO)
C     write(7,'(A,I3/A/(5D16.8))') 'Iteration',ITER,
C    +  'Delta x',(RESID(I),I=1,N)

C Update solution
      DO 141 I = 1,N
      X(I) = X(I) + RESID(I)
  141   CONTINUE

C Exit without computing residual
      IF (JOB.LT.4 .AND. ICNTL(9).EQ.1) GO TO 340

C
C Calculate residual using information in A,IRN,JCN
C If ICNTL(9).GT.1 also calculate |A| |x|
C
      IF (ICNTL(9).EQ.1) THEN
C Now JOB = 4
      DO 151 I = 1,N
        RESID(I) = RHS(I)
  151   CONTINUE
      DO 181 KK = 1,NE
        I = IRN(KK)
        J = JCN(KK)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 181
        RESID(I) = RESID(I) - A(KK)*X(J)
C Matrix is symmetric
        IF (I.NE.J) THEN
          IF (KIND.EQ.1) THEN
          RESID(J) = RESID(J) - A(KK)*X(I)
          ELSE
          RESID(J) = RESID(J) - CONJG(A(KK))*X(I)
          ENDIF
        ENDIF
  181   CONTINUE
C Jump because only one step of iterative refinement requested.
      GO TO 340
      ELSE
      DO 153 I = 1,N
C b - Ax
        RESID(I) = RHS(I)
C |A||x|
        WORKR(I,1) = ZERO
  153   CONTINUE
      DO 183 KK = 1,NE
        I = IRN(KK)
        J = JCN(KK)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 183
        RESID(I) = RESID(I) - A(KK)*X(J)
        WORKR(J,1) = WORKR(J,1) + ABS(A(KK)*X(I))
C Matrix is symmetric
        IF (I.NE.J) THEN
          IF (KIND.EQ.1) THEN
          RESID(J) = RESID(J) - A(KK)*X(I)
          ELSE
          RESID(J) = RESID(J) - CONJG(A(KK))*X(I)
          ENDIF
          WORKR(I,1) = WORKR(I,1) + ABS(A(KK)*X(J))
        ENDIF
  183   CONTINUE
      ENDIF

C Calculate max-norm of solution
      DXMAX = ZERO
      DO 220 I = 1,N
      DXMAX = MAX(DXMAX,ABS(X(I)))
  220   CONTINUE

C Calculate omega(1) and omega(2)
C tau is (||A  ||       ||x||   + |b| )*n*1000*epsilon
C        i.  infinity    max   i
      OMEGA(1) = ZERO
      OMEGA(2) = ZERO
      DO 230 I = 1,N
      TAU = (WORKR(I,2)*DXMAX+ABS(RHS(I)))*N*CTAU
      IF ((WORKR(I,1)+ABS(RHS(I))).GT.TAU) THEN
C |Ax-b| /(|A||x| + |b|)
C     i         i
        OMEGA(1) = MAX(OMEGA(1),ABS(RESID(I))/
   +           (WORKR(I,1)+ABS(RHS(I))))
        IW(I) = 1
      ELSE
C TAU will be zero if all zero row in A, for example
        IF (TAU.GT.ZERO) THEN
C |Ax-b| /(|A||x| + ||A  ||      ||x||   )
C     i      i   i.  infinity   max
          OMEGA(2) = MAX(OMEGA(2),ABS(RESID(I))/
   +             (WORKR(I,1)+WORKR(I,2)*DXMAX))
        END IF
        IW(I) = 2
      END IF
  230   CONTINUE
C
C  Stop the calculations if the backward error is small
      OM2 = OMEGA(1) + OMEGA(2)
      IF ((OM2+ONE).LE.ONE) THEN
C       write(7,'(A)') 'Go to 270 from after 230'
      GO TO 270
      ENDIF
C
C  Check the convergence.
C
C     write(7,'(A,2D16.8)') 'OM2,OLDOM2',OM2,OLDOM2
      IF (OM2.GT.OLDOM2*CNTL(3)) THEN
C  Stop if insufficient decrease in omega.
      IF (OM2.GT.OLDOM2) THEN
C Previous estimate was better ... reinstate it.
        OMEGA(1) = OLDOMG(1)
        OMEGA(2) = OLDOMG(2)
        DO 240 I = 1,N
          X(I) = WORKC(I,2)
  240     CONTINUE
      END IF
C       write(7,'(A)') 'Go to 270 from after 240'
      GO TO 270
      ELSE
C Hold current estimate in case needed later
      DO 250 I = 1,N
        WORKC(I,2) = X(I)
  250   CONTINUE
      OLDOMG(1) = OMEGA(1)
      OLDOMG(2) = OMEGA(2)
      OLDOM2 = OM2
      END IF

  260 CONTINUE
C End of iterative refinement loop.

    INFO(1) = -8
    IF (LP.GE.0 .AND. LDIAG.GE.1) WRITE (LP,9170) INFO(1),ICNTL(9)
 9170 FORMAT ('Error return from ME57D/DD because of ','nonconvergence',
   +     ' of iterative refinement'/'Error INFO(1) = ',I2,'  with',
   +     ' ICNTL','(9) = ',I10)

C Set the RINFO parameters
  270 RINFO(6)  = OMEGA(1)
    RINFO(7)  = OMEGA(2)
    RINFO(8) = ZERO
    DO 271 I=1,N
      RINFO(8) = MAX(RINFO(8),WORKR(I,2))
  271 CONTINUE
    RINFO(9) = DXMAX
    RINFO(10) = ZERO
    DO 272 I=1,N
      RINFO(10) = MAX(RINFO(10),ABS(RESID(I)))
  272 CONTINUE
    IF (RINFO(8)*RINFO(9).NE.ZERO)
   *RINFO(10) = RINFO(10)/(RINFO(8)*RINFO(9))
    INFO(30) = ITER

    IF (INFO(1).LT.0) GO TO 340

C Jump if estimate of error not requested.
    IF (ICNTL(10).LE.0) GO TO 340

C
C Calculate condition numbers and estimate of the error.
C
C Condition numbers obtained through use of norm estimation
C   routine MF71A/AD.
C
C Initializations
C
    LCOND(1) = .FALSE.
    LCOND(2) = .FALSE.
    ERROR    = ZERO
    DO 280 I = 1,N
      IF (IW(I).EQ.1) THEN
      WORKR(I,1) = WORKR(I,1) + ABS(RHS(I))
C |A||x| + |b|
      WORKR(I,2) = ZERO
      LCOND(1) = .TRUE.
      ELSE
C |A||x| + ||A  ||      ||x||
C         i.  infinity   max

      WORKR(I,2) = WORKR(I,1) + WORKR(I,2)*DXMAX
      WORKR(I,1) = ZERO
      LCOND(2) = .TRUE.
      END IF
  280 CONTINUE
C
C  Compute the estimate of COND
C
    DO 330 K = 1,2
      IF (LCOND(K)) THEN
C MF71A/AD has its own built in limit of 5 to the number of iterations
C    allowed. It is this limit that will be used to terminate the
C    following loop.
      KASE = 0
      DO 310 KK = 1,40

C MF71A/AD calculates norm of matrix
C We are calculating the infinity norm of INV(A).W

C Initialize WORKC(1,1).  WORKC(1,2) is a work array.
        CALL MF71AD(N,KASE,WORKC,COND(K),WORKC(1,2),KEEP71)
C
C  KASE = 0........ Computation completed
C  KASE = 1........ W * INV(TRANSPOSE(A)) * Y
C  KASE = 2........ INV(A) * W * Y
C             W is W/W(*,2) .. Y is WORKC(*,1)
C
        IF (KASE.EQ.0) GO TO 320

        IF (KASE.EQ.1) THEN
C Solve system using ME57C/CD.
C WORKC(1,2) is used as workspace
          CALL ME57CD(KIND,1,N,FACT,LFACT,IFACT,LIFACT,1,WORKC,
   *            N,WORKC(1,2),N,IW,ICNTLC,INFO)
          DO 290 I = 1,N
          WORKC(I,1) = WORKR(I,K)*WORKC(I,1)
  290       CONTINUE
        END IF

        IF (KASE.EQ.2) THEN
          DO 300 I = 1,N
          WORKC(I,1) = WORKR(I,K)*WORKC(I,1)
  300       CONTINUE
C Solve system using ME57C/CD.
          CALL ME57CD(KIND,1,N,FACT,LFACT,IFACT,LIFACT,1,WORKC,N,
   *            WORKC(1,2),N,IW,ICNTLC,INFO)
        END IF

  310   CONTINUE

C Error return if MF71AD does not converge
      INFO(1) = -14
      IF (LP.GE.0 .AND. LDIAG.GE.1) WRITE (LP,9160)
 9160 FORMAT ('Error return from ME57D/DD because of ','error in MF71',
   +     'A/AD'/'Error not calculated')

  320   IF (DXMAX.GT.ZERO) COND(K) = COND(K)/DXMAX
      ERROR = ERROR + OMEGA(K)*COND(K)
      ELSE
      COND(K) = ZERO
      ENDIF

  330 CONTINUE

    RINFO(11)  = COND(1)
    RINFO(12)  = COND(2)
    RINFO(13)  = ERROR

C
C If requested, print output parameters.
 340  IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
      WRITE(MP,99980) INFO(1)
99980 FORMAT (/'Leaving iterative refinement solution phase ',
   +  '(ME57DD) with ...'/
   1    'INFO (1)                        =',I12/)
      IF (INFO(1).LT.0) GO TO 500
      IF (ICNTL(9).GT.1) THEN
      WRITE(MP,99981) INFO(30),(RINFO(I),I=6,10)
99981   FORMAT(
   1   'INFO(30)  Number steps iterative ref   =',I10/
   1   'RINFO(6)  Backward errors  (OMEGA(1))  =',1PD10.3/
   2   '-----(7)  Backward errors  (OMEGA(2))  =',1PD10.3/
   3   '-----(8)  Infinity norm of matrix    =',1PD10.3/
   4   '-----(9)  Infinity norm of solution    =',1PD10.3/
   5   '-----(10) Norm of scaled residuals   =',1PD10.3)
      IF (ICNTL(10).GT.0) WRITE(MP,99979) (RINFO(I),I=11,13)
99979     FORMAT (
   1     'RINFO(11) Condition number (COND(1))   =',1PD10.3/
   1     'RINFO(12) Condition number (COND(2))   =',1PD10.3/
   1     'RINFO(13) Error in solution        =',1PD10.3)
      WRITE(MP,'(/A,I10)') 'Residual'
      K=MIN(N,10)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(1P,6D13.3)') (RESID(I),I=1,K)
      IF (K.LT.N) WRITE (MP,'(A)') '   . . .'
      ELSE
C ICNTL(9) = 1
      IF (JOB.GE.1 .AND. JOB.LE.3) THEN
        WRITE(MP,'(/A,I10)') 'Correction to solution'
      ELSE
        WRITE(MP,'(/A,I10)') 'Residual'
      ENDIF
      K=MIN(N,10)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(1P,6D13.3)') (RESID(I),I=1,K)
      IF (K.LT.N) WRITE (MP,'(A)') '   . . .'
      END IF

C Print solution
      K=MIN(N,10)
      IF (LDIAG.GE.4) K = N
      WRITE(MP,'(/A,I10)') 'Solution'
      WRITE (MP,'(1P,6D13.3)') (X(I),I=1,K)
      IF (K.LT.N) WRITE (MP,'(A)') '   . . .'

    END IF

 500  RETURN

    END
    SUBROUTINE ME57ED(N,IC,KEEP,FACT,LFACT,NEWFAC,LNEW,
   *            IFACT,LIFACT,NEWIFC,LINEW,INFO)
    INTEGER N,IC,KEEP(*),LFACT,LNEW,LIFACT,LINEW,INFO(40)
    COMPLEX*16 FACT(LFACT),NEWFAC(LNEW)
    INTEGER IFACT(LIFACT),NEWIFC(LINEW)
C Local variables
    INTEGER APOSBB,ASTK,HOLD,I,ISTK,IWPOS,MOVE,NFRONT

C HOLD determines part of keep holding saved variables from ME57OD.
    HOLD = N + 3

C Initialize INFO(1) and INFO(2)
    INFO(1) = 0
    INFO(2) = 0

C Test to see whether to map real or integer space
    IF (IC.GE.1) THEN
C Remap of integer space
      IF (LINEW.LE.LIFACT) THEN
      INFO(1) = -7
      INFO(2) = LINEW
      RETURN
      ENDIF
      IWPOS = KEEP(HOLD+7)
      ISTK  = KEEP(HOLD+14)
      NFRONT = KEEP(HOLD+23)
      DO 10 I = 1,IWPOS+NFRONT-1
      NEWIFC(I) = IFACT(I)
   10   CONTINUE
C Move distance
      MOVE = LINEW - LIFACT
      DO 20 I = ISTK+1,LIFACT
      NEWIFC(I+MOVE) = IFACT(I)
   20   CONTINUE
C Reset INPUT, ISTK, PTIRN
      KEEP(HOLD+13) = KEEP(HOLD+13) + MOVE
      KEEP(HOLD+14) = ISTK + MOVE
      KEEP(HOLD+18) = KEEP(HOLD+18) + MOVE
    ENDIF
    IF (IC.NE.1) THEN
C Remap of real space
      IF (LNEW.LE.LFACT) THEN
      INFO(1) = -7
      INFO(2) = LNEW
      RETURN
      ENDIF
C Was .. APOS = KEEP(HOLD+8)
      APOSBB = KEEP(HOLD+9)
      ASTK   = KEEP(HOLD+15)
      DO 60 I = 1, APOSBB-1
      NEWFAC(I) = FACT(I)
   60   CONTINUE
C Move distance
      MOVE = LNEW - LFACT
      DO 70 I = ASTK+1,LFACT
      NEWFAC(I+MOVE) = FACT(I)
   70   CONTINUE
C Reset AINPUT, ASTK, PTRA
      KEEP(HOLD+12) = KEEP(HOLD+12) + MOVE
      KEEP(HOLD+15) = ASTK + MOVE
      KEEP(HOLD+19) = KEEP(HOLD+19) + MOVE
    ENDIF
    RETURN
    END
    SUBROUTINE ME57GD(KIND,N,NE,IRN,JCN,IW,IPE,COUNT,FLAG,IWFR,
   *            ICNTL,INFO)

C Given the positions of the entries of a symmetric matrix, construct
C   the sparsity pattern of the whole matrix. Either one of a pair
C   (I,J),(J,I) may be used to represent the pair. Duplicates are
C   ignored.

    INTEGER KIND,N,NE,IRN(NE),JCN(NE),IW(NE*2+N),IPE(N),COUNT(N),
   +      FLAG(N),IWFR,ICNTL(20),INFO(40)
C N must be set to the matrix order. It is not altered.
C NE must be set to the number of entries input. It is not altered.
C IRN(K),K=1,2,...,NE must be set to the row indices of the entries on
C   input.  IRN is not changed.
C JCN(K),K=1,2,...,NE must be set to the column indices of the entries
C   on input.  JCN is not changed.
C IW need not be set on input. On output it contains lists of
C   column indices.
C IPE need not be set on input. On output IPE(I) points to the start of
C   the entry in IW for row I, I=1,2,...,N.
C COUNT need not be set. On output COUNT(I), I=1,2,..,N, contains the
C   number of off-diagonal entries in row I excluding duplicates.
C FLAG is used for workspace to hold flags to permit duplicate entries
C   to be identified quickly.
C IWFR need not be set on input. On output it points to the first
C   unused location in IW.
C ICNTL Warning messages are printed on stream number ICNTL(2) if
C    ICNTL(2).GT.0 and ICNTL(5).GT.1.
C INFO need not be set on input. On output,
C  INFO(1) has one of the values:
C   0 No out-of-range index or duplicate entry found.
C   1 Out-of-range index found.
C   2 Duplicate entry found.
C   3 Out-of-range index found and duplicate entry found.
C  INFO(3) is set to the number of out-of-range entries.
C  INFO(4) is set to the number of off-diagonal duplicate entries.
C
    INTRINSIC MAX,MIN
C
C Local variables
    INTEGER I,J,K,L,LDIAG,LP,WP
C I Row index
C J Column index
C K Position in IRN, JCN, or IW.
C L Position in IW.
C LDIAG Level of diagnostic printing.
C WP Stream for printing warning messages.

C Set LDIAG and WP
    LP = ICNTL(1)
    WP = ICNTL(2)
    LDIAG = ICNTL(5)
    IF (WP.LT.0) LDIAG = 0

C Initialize INFO(1)
    INFO(1) = 0

C Count in INFO(3) the number of out-of-range entries, initialize FLAG,
C    and count in COUNT the numbers of entries in the rows including
C    duplicates.
    INFO(3) = 0
    DO 10 I = 1,N
      FLAG(I) = 0
      COUNT(I) = 0
   10 CONTINUE
    DO 20 K = 1,NE
      I = IRN(K)
      J = JCN(K)
      IF (KIND.NE.1 .AND. I.GT.J) THEN
      INFO(1) = -19
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I8,A)')
   +    '**** Error return from ME57AD ****  INFO(1) =',INFO(1),
   +    'Hermitian matrix but entry ',K,' in lower triangle'
      INFO(2) = K
      RETURN
      ENDIF
      IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) THEN
      INFO(3) = INFO(3) + 1
      INFO(1) = 1
      IF (INFO(3).EQ.1 .AND. LDIAG.GT.1) WRITE (WP,'(2A,I2)')
   +      '*** Warning message from subroutine ME57AD ***',
   +      ' INFO(1) =',INFO(1)
      IF (INFO(3).LE.10 .AND. LDIAG.GT.1) WRITE (WP,'(3(I10,A))')
   +       K,'th entry (in row',I,' and column',J,') ignored'

      ELSE IF (I.NE.J) THEN
      COUNT(I) = COUNT(I) + 1
      COUNT(J) = COUNT(J) + 1
      END IF

   20 CONTINUE
C
C Accumulate row counts in IPE which is set so that IPE(I) points to
C   position after end of row I.
    IPE(1) = COUNT(1)+1
    DO 30 I = 2,N
      IPE(I) = IPE(I-1) + COUNT(I)
   30 CONTINUE
C
C Run through putting the matrix entries in the right places. IPE is
C   used for holding running pointers and is left holding pointers to
C   row starts.
    DO 40 K = 1,NE
      I = IRN(K)
      J = JCN(K)
      IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N .OR. I.EQ.J) GO TO 40
      IPE(I) = IPE(I) - 1
      IW(IPE(I)) = J
      IPE(J) = IPE(J) - 1
      IW(IPE(J)) = I
   40 CONTINUE
C
C Remove duplicates.
    INFO(4) = 0
C IWFR points to the current position in the compressed set of rows.
    IWFR = 1
C At the start of cycle I of this loop FLAG(J).LT.I for J=1,2,...,N.
C   During the loop FLAG(J) is set to I if an entry in column J is
C   found. This permits duplicates to be recognized quickly.
    DO 60 I = 1,N
      L = IPE(I)
      IPE(I) = IWFR
      DO 50 K = L,L+COUNT(I)-1
      J = IW(K)

      IF (FLAG(J).NE.I) THEN
        FLAG(J) = I
        IW(IWFR) = J
        IWFR = IWFR + 1
      ELSE
C Count duplicates only once.
        IF (I.LT.J) INFO(4) = INFO(4) + 1
      END IF

   50   CONTINUE
C Set COUNT to number without duplicates.
      COUNT(I) = IWFR - IPE(I)
   60 CONTINUE
C
C Test whether duplicates found
    IF (INFO(4).GT.0) THEN
      INFO(1) = INFO(1) + 2
      IF (LDIAG.GT.1 .AND. WP.GE.0) WRITE (WP,'(A/I10,A)')
   +    '*** Warning message from subroutine ME57AD ***',INFO(4),
   +    ' off-diagonal duplicate entries found'
    END IF

    END
    SUBROUTINE ME57JD(KIND,N,NE,IRN,JCN,PERM,IW,IPE,COUNT,FLAG,IWFR,
   +            ICNTL,INFO)

C Given the positions of the entries of a symmetric matrix and a
C   permutation, construct the sparsity pattern of the upper
C   triangular part of the matrix. Either one of a pair
C   (I,J),(J,I) may be used to represent the pair. Duplicates are
C   ignored.

    INTEGER KIND,N,NE,IRN(NE),JCN(NE),IW(NE+N),IPE(N),COUNT(N),
   +      PERM(N),FLAG(N),IWFR,ICNTL(20),INFO(40)
C N must be set to the matrix order. It is not altered.
C NE must be set to the number of entries input. It is not altered.
C IRN(K),K=1,2,...,NE must be set to the row indices of the entries on
C   input. If IRN(K) or JCN(K) is out of range, IRN(K) is replaced by
C   zero. Otherwise, IRN(K) is not changed.
C JCN(K),K=1,2,...,NE must be set to the column indices of the entries
C   on input. If IRN(K) or JCN(K) is out of range, JCN(K) is replaced
C   by zero. Otherwise, JCN(K) is not changed.
C PERM must be set so that PERM(I) holds the position of variable I
C   in the permuted order.  Its validity as a permutation will have
C   been checked in ME57A/AD.
C IW need not be set on input. On output it contains lists of
C   column indices, each list being headed by its length.
C IPE need not be set on input. On output IPE(I) points to the start of
C   the entry in IW for row I, I=1,2,...,N.
C COUNT need not be set. On output COUNT(I), I=1,2,..,N, contains the
C   number of off-diagonal entries in row I including duplicates.
C   COUNT(0) contains the number of entries with one or both indices
C   out of range.
C FLAG is used for workspace to hold flags to permit duplicate entries
C   to be identified quickly.
C IWFR need not be set on input. On output it points to the first
C   unused location in IW.
C ICNTL Warning messages are printed on stream number ICNTL(2) if
C    ICNTL(2).GT.0 and ICNTL(3).GT.1.
C INFO need not be set on input. On output, INFO(1) has one of the
C   values:
C   0 No out-of-range index or duplicate entry found.
C   1 Out-of-range index found.
C   2 Duplicate entry found.
C   3 Out-of-range index found and duplicate entry found.
C  INFO(3) is set to the number of faulty entries.
C  INFO(4) is set to the number of off-diagonal duplicate entries.

    INTRINSIC MAX,MIN
C
C Local variables
    INTEGER I,J,K,L,LDIAG,LP,WP
C I Row index
C J Column index
C K Position in IRN, JCN, or IW.
C L Position in IW.
C LDIAG Level of monitor printing.
C WP Stream for printing.

C Set LDIAG and WP
    WP = ICNTL(2)
    LP = ICNTL(1)
    LDIAG = ICNTL(5)
    IF (WP.LT.0) LDIAG = 0

C Initialize INFO(1), FLAG, and COUNT.
    INFO(1) = 0
    DO 10 I = 1,N
      FLAG(I) = 0
      COUNT(I) = 0
   10 CONTINUE

C Count in INFO(3) the number of out-of-range entries, initialize FLAG,
C    and count in COUNT the numbers of entries in the rows.
    INFO(3) = 0
    DO 30 K = 1,NE
      I = IRN(K)
      J = JCN(K)
      IF (KIND.NE.1 .AND. I.GT.J) THEN
      INFO(1) = -19
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I8,A)')
   +    '**** Error return from ME57AD ****  INFO(1) =',INFO(1),
   *    'Hermitian matrix but entry ',K,' in lower triangle'
      INFO(2) = K
      RETURN
      ENDIF
      IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) THEN
      IRN(K) = 0
      JCN(K) = 0
      INFO(3) = INFO(3) + 1
      INFO(1) = 1
      IF (INFO(3).EQ.1 .AND. LDIAG.GT.1) WRITE (WP,'(2A,I2)')
   +      '*** Warning message from subroutine ME57AD ***',
   +      ' INFO(1) =',INFO(1)
      IF (INFO(3).LE.10 .AND. LDIAG.GT.1) WRITE (WP,'(3(I10,A))')
   +      K,'th entry (in row',I,' and column',J,') ignored'

      ELSE IF (PERM(I).LE.PERM(J)) THEN
      COUNT(I) = COUNT(I) + 1
      ELSE
      COUNT(J) = COUNT(J) + 1
      END IF

   30 CONTINUE
C
C Accumulate row counts in IPE ... one added for row length location.
    IPE(1) = COUNT(1) + 1
    DO 40 I = 2,N
      IPE(I) = IPE(I-1) + COUNT(I) + 1
   40 CONTINUE

C Run through putting the matrix entries in the right places. IPE is
C   used for holding running pointers and is left holding pointers to
C   row starts.
    DO 50 K = 1,NE
      I = IRN(K)
      J = JCN(K)
      IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 50
      IF (PERM(I).LE.PERM(J)) THEN
      IW(IPE(I)) = J
      IPE(I) = IPE(I) - 1
      ELSE
      IW(IPE(J)) = I
      IPE(J) = IPE(J) - 1
      END IF
   50 CONTINUE

C Remove duplicates
C IWFR points to the current position in the compressed set of rows.
    IWFR = 1
    INFO(4) = 0

C At the start of cycle I of this loop FLAG(J).LT.I for J=1,2,...,N.
C   During the loop FLAG(J) is set to I if an entry in column J is
C   found. This permits duplicates to be recognized quickly.
    DO 70 I = 1,N
      L = IPE(I)
      IPE(I) = IWFR

      DO 60 K = L + 1,L + COUNT(I)
      J = IW(K)
      IF (FLAG(J).NE.I) THEN
        FLAG(J) = I
        IWFR = IWFR + 1
        IW(IWFR) = J
      ELSE
C Count duplicates only once.
        IF (I.LT.J) INFO(4) = INFO(4) + 1
      END IF
   60   CONTINUE

      IF (IWFR.GT.IPE(I)) THEN
      IW(IPE(I)) = IWFR - IPE(I)
      IWFR = IWFR + 1
      ELSE
      IPE(I) = 0
      END IF

   70 CONTINUE

C Test whether duplicates found
    IF (INFO(4).GT.0) THEN
      INFO(1) = INFO(1) + 2
      IF (LDIAG.GT.1 .AND. WP.GE.0) WRITE (WP,'(A/I10,A)')
   +    '*** Warning message from subroutine ME57AD ***',
   +    INFO(4),' off-diagonal duplicate entries found'
    END IF

    END

    SUBROUTINE ME57KD(N, IPE, IW, LW, IWFR, PERM, IPS, NV, FLAG,
   *            NCMPA)
    INTEGER N,LW,IWFR,NCMPA
    INTEGER IPE(N)
    INTEGER IW(LW), PERM(N), IPS(N), NV(N), FLAG(N)
C
C Using a given pivotal sequence and a representation of the matrix that
C   includes only non-zeros of the strictly upper-triangular part
C   of the permuted matrix, construct tree pointers.
C
C N must be set to the matrix order. It is not altered.
C IPE(I) must be set to point to the position in IW of the
C   start of row I or have the value zero if row I has no off-
C   diagonal non-zeros. during execution it is used as follows.
C   If variable I is eliminated then IPE(I) points to the list
C   of variables for created element I. If element I is
C   absorbed into newly created element J then IPE(I)=-J.
C IW must be set on entry to hold lists of variables by
C   rows, each list being headed by its length. when a variable
C   is eliminated its list is replaced by a list of variables
C   in the new element.
C LW must be set to the length of IW. It is not altered.
C IWFR must be set to the position in IW of the first free variable.
C   It is revised during execution, continuing to have this meaning.
C PERM(K) must be set to hold the position of variable K in the
C   pivot order. It is not altered.
C IPS(I) need not be set by the user and will be used to hold the
C   inverse permutation to PERM.
C NV need not be set. If variable J has not been eliminated then
C   the last element whose leading variable (variable earliest
C   in the pivot sequence) is J is element NV(J). If element J
C   exists then the last element having the same leading
C   variable is NV(J). In both cases NV(J)=0 if there is no such
C   element. If element J has been merged into a later element
C   then NV(J) is the degree at the time of elimination.
C FLAG is used as workspace for variable flags.
C   FLAG(JS)=ME if JS has been included in the list for ME.
C
    INTEGER I,J,ML,MS,ME,IP,MINJS,IE,KDUMMY,JP
    INTEGER LN,JP1,JS,LWFR,JP2,JE

    EXTERNAL ME57FD

C
C Initializations
    DO 10 I=1,N
      FLAG(I) = 0
      NV(I) = 0
      J = PERM(I)
      IPS(J) = I
   10 CONTINUE
    NCMPA = 0
C
C Start of main loop
C
    DO 100 ML=1,N
C ME=MS is the name of the variable eliminated and
C   of the element created in the main loop.
      MS = IPS(ML)
      ME = MS
      FLAG(MS) = ME
C
C Merge row MS with all the elements having MS as leading variable.
C IP points to the start of the new list.
      IP = IWFR
C MINJS is set to the position in the order of the leading variable
C   in the new list.
      MINJS = N
      IE = ME
      DO 70 KDUMMY=1,N
C Search variable list of element IE.
C JP points to the current position in the list being searched.
      JP = IPE(IE)
C LN is the length of the list being searched.
      LN = 0
      IF (JP.LE.0) GO TO 60
      LN = IW(JP)
C
C Search for different variables and add them to list,
C   compressing when necessary
      DO 50 JP1=1,LN
        JP = JP + 1
C Place next variable in JS.
        JS = IW(JP)
C Jump if variable has already been included.
        IF (FLAG(JS).EQ.ME) GO TO 50
        FLAG(JS) = ME
        IF (IWFR.LT.LW) GO TO 40
C Prepare for compressing IW by adjusting pointer to and length of
C   the list for IE to refer to the remaining entries.
        IPE(IE) = JP
        IW(JP) = LN - JP1
C Compress IW.
        CALL ME57FD(N, IPE, IW, IP-1, LWFR, NCMPA)
C Copy new list forward
        JP2 = IWFR - 1
        IWFR = LWFR
        IF (IP.GT.JP2) GO TO 30
        DO 20 JP=IP,JP2
          IW(IWFR) = IW(JP)
          IWFR = IWFR + 1
   20     CONTINUE
   30     IP = LWFR
        JP = IPE(IE)
C Add variable JS to new list.
   40     IW(IWFR) = JS
        MINJS = MIN0(MINJS,PERM(JS)+0)
        IWFR = IWFR + 1
   50   CONTINUE
C Record absorption of element IE into new element.
   60   IPE(IE) = -ME
C Pick up next element with leading variable MS.
      JE = NV(IE)
C Store degree of IE.
      NV(IE) = LN + 1
      IE = JE
C Leave loop if there are no more elements.
      IF (IE.EQ.0) GO TO 80
   70   CONTINUE
   80   IF (IWFR.GT.IP) GO TO 90
C Deal with null new element.
      IPE(ME) = 0
      NV(ME) = 1
      GO TO 100
C Link new element with others having same leading variable.
   90   MINJS = IPS(MINJS)
      NV(ME) = NV(MINJS)
      NV(MINJS) = ME
C Move first entry in new list to end to allow room for length at
C   front. Set pointer to front.
      IW(IWFR) = IW(IP)
      IW(IP) = IWFR - IP
      IPE(ME) = IP
      IWFR = IWFR + 1
  100 CONTINUE

    RETURN
    END
C** end of ME57KD**

    SUBROUTINE ME57FD(N, IPE, IW, LW, IWFR, NCMPA)
C Compress lists held by ME57KD in IW and adjust pointers
C   in IPE to correspond.
    INTEGER N,LW,IWFR,NCMPA
    INTEGER IPE(N)
C
    INTEGER   IW(LW)
C N is the matrix order. It is not altered.
C IPE(I) points to the position in IW of the start of list I or is
C   zero if there is no list I. On exit it points to the new position.
C IW holds the lists, each headed by its length. On output the same
C   lists are held, but they are now compressed together.
C LW holds the length of IW. It is not altered.
C IWFR need not be set on entry. On exit it points to the first free
C   location in IW.
C
    INTEGER I,K1,LWFR,IR,K,K2
    NCMPA = NCMPA + 1
C Prepare for compressing by storing the lengths of the
C   lists in IPE and setting the first entry of each list to
C   -(list number).
    DO 10 I=1,N
      K1 = IPE(I)
      IF (K1.LE.0) GO TO 10
      IPE(I) = IW(K1)
      IW(K1) = -I
   10 CONTINUE
C
C Compress
C IWFR points just beyond the end of the compressed file.
C LWFR points just beyond the end of the uncompressed file.
    IWFR = 1
    LWFR = IWFR
    DO 60 IR=1,N
      IF (LWFR.GT.LW) GO TO 70
C Search for the next negative entry.
      DO 20 K=LWFR,LW
      IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
      GO TO 70
C Pick up entry number, store length in new position, set new pointer
C   and prepare to copy list.
   30   I = -IW(K)
      IW(IWFR) = IPE(I)
      IPE(I) = IWFR
      K1 = K + 1
      K2 = K + IW(IWFR)
      IWFR = IWFR + 1
      IF (K1.GT.K2) GO TO 50
C Copy list to new position.
      DO 40 K=K1,K2
      IW(IWFR) = IW(K)
      IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN
    END
C--------------------------------------------------------------------
C        HSL 2000
C      --
C-         Copyright CCLRC Rutherford Appleton Laboratory
C      --
C--------------------------------------------------------------------
    SUBROUTINE ME57LD(N, IPE, NV, IPS, NE, NA, NODE, PERM, NSTEPS,
   *            FILS, FRERE, ND, NEMIN, SUBORD)
C
C Tree search
C
C Given son to father tree pointers, reorder so that eldest son has
C   smallest degree and perform depth-first
C   search to find pivot order and number of eliminations
C   and assemblies at each stage.
    INTEGER N, NSTEPS
    INTEGER ND(N)
    INTEGER IPE(N), FILS(N), FRERE(N), SUBORD(N)
C
    INTEGER NV(N), IPS(N), NE(N), NA(N), NODE(N), PERM(N)
    INTEGER NEMIN
C N must be set to the matrix order. It is not altered.
C IPE(I) must be set equal to -(father of node I) or zero if
C    node is a root, if NV(I) > 0. If NV(I) = 0, then I is
C    subordinate variable of a supervariable and -IPE(I) points to
C    principal variable.  It is altered to point to its next
C    younger brother if it has one, but otherwise is not changed.
C NV(I) must be set to zero if variable is a subordinate variable
C    of a supervariable and to the degree otherwise.
C    NV is not altered.
C IPS(I) need not be set. It is used temporarily to hold
C    -(eldest son of node I) if it has one and 0 otherwise. It is
C    finally set to hold the position of node I in the order.
C NE(IS) need not be set. It is set to the number of variables
C    eliminated at stage IS of the elimination.
C NA(IS) need not be set. It is set to the number of elements
C    assembled at stage IS of the elimination.
C NODE (I) need not be set before entry. It is used during the code
C    to hold the number of subordinate variables for variable I and
C    on output it holds
C    the node (in dfs ordering) at which variable I is eliminated.
C    It is also defined for subordinate variables.
C PERM is set to the new permutation after dfs of tree.  PERM(I) is
C    the position of variable I in the pivot order.
C ND(IS) need not be set. It is set to the degree at stage IS of
C   the elimination.
C NSTEPS need not be set. It is set to the number of elimination steps.
C NEMIN is used to control the amalgamation process between
C     a son and its father (if the number of fully summed
C     variables of both nodes is smaller than NEMIN).
C SUBORD(I) need not be set. It holds the first subordinate variable
C     for variable I if I
C     is a principal variable and holds the next subordinate variable
C     if otherwise.  It is zero at the end of the chain.
C
    INTEGER I,IF,IS,NR,NR1,INS,INL,INB,INF,INFS,INSW
    INTEGER K,L,ISON,IN,IFSON,INO
    INTEGER INOS,IB,IL,INT
    INTEGER IPERM

C
C Initialize IPS and NE.
    DO 10 I=1,N
      IPS(I) = 0
      NE(I) = 0
      NODE(I) = 0
      SUBORD(I) = 0
   10 CONTINUE
C
C Set IPS(I) to -(eldest son of node I) and IPE(I) to next younger
C   brother of node I if it has one.
    NR = N + 1
    DO 50 I=1,N
      IF = -IPE(I)
      IF (NV(I).EQ.0) THEN
C I is a subordinate node, principal variable is IF
      IF (SUBORD(IF).NE.0) SUBORD(I) = SUBORD(IF)
      SUBORD(IF) = I
      NODE(IF) = NODE(IF)+1
      ELSE
C Node IF is the father of node I.
      IF (IF.NE.0) THEN
C IS is younger brother of node I.
C IPS(IF) will eventually point to - eldest son of IF.
        IS = -IPS(IF)
        IF (IS.GT.0) IPE(I) = IS
        IPS(IF) = -I
      ELSE
C I is a root node
        NR = NR - 1
        NE(NR) = I
      ENDIF
      ENDIF
   50 CONTINUE
C
C We reorganize the tree so that the eldest son has maximum number of
C variables.  We combine nodes when the number of variables in a son
C is greater than or equal to the number of variables in the father.
C If the eldest son has the maximum number of variables,
C and if a combination is possible, it has to be possible with
C the eldest son.
C
C FILS is just used as workspace during this reorganization and is reset
C afterwards.

    DO 999 I=1,N
     FILS(I) = IPS(I)
 999  CONTINUE

    NR1 = NR
    INS = 0
C Jump if all roots processed.
 1000 IF (NR1.GT.N) GO TO 1151
C Get next root
    INS = NE(NR1)
    NR1 = NR1 + 1
C Depth first search through eldest sons.
 1070 INL = FILS(INS)
    IF (INL.LT.0) THEN
     INS = -INL
     GO TO 1070
    ENDIF
C INS is leaf node.

 1080 IF (IPE(INS).LT.0) THEN
C INS is youngest son otherwise IPE value would be positive.
     INS     = -IPE(INS)
C INS is now the father of the reorganized son so we can
C   clear the pointer to the sons.
     FILS(INS) = 0
C Continue backtracking until we encounter node with younger brother.
     GO TO 1080
    ENDIF

    IF (IPE(INS).EQ.0) THEN
C INS is a root, check for next one.
     INS = 0
     GO TO 1000
    ENDIF
C INB is younger brother of INS.
    INB = IPE(INS)

C?? I think this test is the wrong way round
    IF (NV(INB).GE.NV(INS)) THEN
C?? So reversed
C   IF (NV(INS).GE.NV(INB)) THEN
     INS = INB
C Do depth first search from younger brother
     GO TO 1070
    ENDIF
C
C Exchange INB and INS
C Find previous brother of INS (could be the father)
C then we do depth first search with INS = INB
C
    INF = INB
 1090 INF = IPE(INF)
    IF (INF.GT.0) GO TO 1090
C -INF IS THE FATHER
    INF  = -INF
    INFS = -FILS(INF)
C INFS is eldest son of INF
    IF (INFS.EQ.INS) THEN
C  INS is eldest brother .. a role which INB now assumes
      FILS(INF) = -INB
      IPS(INF)  = -INB
      IPE(INS)  = IPE(INB)
      IPE(INB)  = INS
    ELSE
      INSW = INFS
 1100   INFS = IPE(INSW)
      IF (INFS.NE.INS) THEN
      INSW = INFS
      GO TO 1100
      ENDIF
      IPE(INS) = IPE(INB)
      IPE(INB) = INS
      IPE(INSW)= INB
    ENDIF
      INS    = INB
C Depth first search from moved younger brother
      GO TO 1070
C Set FRERE and FILS
 1151 DO 51 I=1,N
     FRERE(I) = IPE(I)
     FILS(I) = IPS(I)
 51   CONTINUE
C
C Depth-first search.
C IL holds the current tree level. Roots are at level N, their sons
C   are at level N-1, etc.
C IS holds the current elimination stage. We accumulate the number
C   of eliminations at stage is directly in NE(IS). The number of
C   assemblies is accumulated temporarily in NA(IL), for tree
C   level IL, and is transferred to NA(IS) when we reach the
C   appropriate stage IS.
    IS = 1
C I is the current node.
    I = 0
C IPERM is used as pointer to setting permutation vector
    IPERM = 1
    DO 160 K=1,N
      IF (I.GT.0) GO TO 60
C Pick up next root.
C Stop if all roots used (needed because of subordinate variables)
      IF (NR.GT.N) GO TO 161
      I = NE(NR)
      NE(NR) = 0
      NR = NR + 1
      IL = N
      NA(N) = 0
C Go to son for as long as possible, clearing father-son pointers
C   in IPS as each is used and setting NA(IL)=0 for all levels
C   reached.
   60   CONTINUE
      DO 70 L=1,N
      IF (IPS(I).GE.0) GO TO 80
      ISON = -IPS(I)
      IPS(I) = 0
      I = ISON
      IL = IL - 1
      NA(IL) = 0
   70   CONTINUE
   80   CONTINUE
C?? Do we want to expand for subordinate variables
C Record position of node I in the order.
      IPS(I) = K
C Add number of subordinate variables to variable I
      NE(IS) = NE(IS) + NODE(I) + 1
      IF (IL.LT.N) NA(IL+1) = NA(IL+1) + 1
      NA(IS) = NA(IL)
      ND(IS) = NV(I)
      NODE(I) = IS
      PERM(I) = IPERM
      IPERM = IPERM + 1
C Order subordinate variables to node I
      IN = I
  777   IF (SUBORD(IN).EQ.0) GO TO 778
      IN = SUBORD(IN)
      NODE(IN) = IS
      PERM(IN) = IPERM
      IPERM = IPERM + 1
      GO TO 777
C Check for static condensation
  778   IF (NA(IS).NE.1) GO TO 90
      IF (ND(IS-1)-NE(IS-1).EQ.ND(IS)) GO TO 100
C Check for small numbers of eliminations in both last two steps.
   90   IF (NE(IS).GE.NEMIN) GO TO 110
      IF (NA(IS).EQ.0) GO TO 110
      IF (NE(IS-1).GE.NEMIN) GO TO 110

C Combine the last two steps
  100   NA(IS-1) = NA(IS-1) + NA(IS) - 1
      ND(IS-1) = ND(IS) + NE(IS-1)
      NE(IS-1) = NE(IS) + NE(IS-1)
      NE(IS) = 0
      NODE(I) = IS-1
C Find eldest son (IFSON) of node I (IS)
C Note that node I must have a son (node IS-1 is youngest)
      IFSON = -FILS(I)
C Now find youngest son INO (he is node IS-1)
      IN = IFSON
 102    INO = IN
      IN =  FRERE(IN)
      IF (IN.GT.0) GO TO 102
C Cannot be root node .. so points to father
C Merge node IS-1 (INO) into node IS (I)
      NV(INO) = 0
C IPE already set .. was father pointer now principal variable pointer
C Now make subsidiary nodes of INO into subsidiary nodes of I.
C Subordinate nodes of INO become subordinate nodes of I
      IN = I
  888   IF (SUBORD(IN).EQ.0) GO TO 889
      IN = SUBORD(IN)
      NODE(IN) = IS-1
      GO TO 888
  889   SUBORD(IN) = INO
      IN = INO
      IF (SUBORD(IN).EQ.0) GO TO 887
      IN = SUBORD(IN)
      IPE(IN) = -I
  887   CONTINUE

C INOS is eldest son of INO
    INOS = -FILS(INO)

C Find elder brother of node INO
C First check to see if he is only son
    IF (IFSON.EQ.INO) GO TO 107
    IN = IFSON
 105  INS = IN
    IN =  FRERE(IN)
    IF (IN.NE.INO) GO TO 105
C INS is older brother .. make him brother of first son of INO (ie INOS)
C and  make INOS point to I now as father.
C Jump if there is no son of INO
      IF (INOS.EQ.0) THEN
C Elder brother of INO just points to (new) father.
      FRERE(INS) = -I
      GO TO 120
      ELSE
      FRERE(INS) =  INOS
      ENDIF
C INT is youngest brother of INOS.  Make him point to (new) father.
 107    IN = INOS
      IF (IN.EQ.0) GO TO 120
 108    INT = IN
      IN =  FRERE(IN)
      IF (IN.GT.0) GO TO 108
      FRERE(INT) = -I
      GO TO 120
  110   IS = IS + 1
  120   IB = IPE(I)
      IF (IB.GE.0) THEN
C Node I has a younger brother or is a root
      IF (IB.GT.0) NA(IL) = 0
      I = IB
      GO TO 160
      ELSE
C I has no brothers. Go to father of node I
      I = -IB
      IL = IL + 1
      ENDIF
  160 CONTINUE
  161 NSTEPS = IS - 1
    RETURN
    END
    SUBROUTINE ME57MD(KIND,N,NE,IRN,JCN,MAP,IRNPRM,
   +            LROW,PERM,COUNT,IDIAG)
C
C This subroutine is called by ME57A/AD and generates the map that
C   reorders the user's input so that the upper triangle of the
C   permuted matrix is held by rows. No check is made for duplicates.
C
    INTEGER KIND,N,NE
C IRNPRM(N+NE) has this dimension to include possibility of no
C    diagonals in input.
    INTEGER IRN(NE),JCN(NE),MAP(NE),IRNPRM(N+NE),LROW(N),PERM(N),
   +      COUNT(N),
   +      IDIAG(N)
C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C N must be set to the matrix order. It is not altered.
C NE must be set to the number of entries input. It is not altered.
C IRN(K) and JCN(K), K=1,2,...,NE must be set to the row and column
C   numbers of the entries. If entry (IRN(K),JCN(K)) lies in the
C   lower triangular part of the permuted matrix, the values of
C   IRN(K) and JCN(K) are interchanged. Otherwise, these arrays are
C   not changed.
C MAP need not be set on input and on return holds the positions of
C   entries when the permuted upper triangle is ordered by rows.
C LROW need not be set. On return, LROW(I),I=1,N holds the number of
C   entries in row I of the permuted matrix.
C PERM(I) must be set to the position of variable I in the
C   pivot order, I=1,2,...,N.
C COUNT is used for workspace. It is set to row counts and then
C   accumulated row counts.
C IDIAG is used for workspace. It is used as pointer to diagonal entry
C   that is first in the row.
C
C Local variables
    INTEGER EXPNE,I,J,K
C I Row index
C J Column index
C K Position in IRN or JCN.

C Accumulate row counts in COUNT, and interchange row and column
C   numbers where necessary.
    DO 10 I = 1,N
C Set to 1 since diagonal will always be present as first entry.
      COUNT(I) = 1
C IDIAG used first as flag to identify duplicate diagonals
      IDIAG(I) = 0
   10 CONTINUE

C EXPNE counts number of entries plus missing diagonals
    EXPNE = NE + N
    DO 20 K = 1,NE

      I = IRN(K)
      J = JCN(K)

      IF (MAX(I,J).GT.N .OR. MIN(I,J).LT.1) THEN
      EXPNE = EXPNE - 1
      GO TO 20
      ENDIF

C Check for duplicate diagonals
      IF (I.EQ.J) THEN
      I = PERM(I)
      IF (IDIAG(I).GE.1) THEN
        COUNT(I) = COUNT(I) + 1
        IDIAG(I) = IDIAG(I) + 1
      ELSE
        IDIAG(I) = 1
        EXPNE = EXPNE - 1
      ENDIF
      GO TO 20
      ENDIF

      IF (PERM(I).LT.PERM(J)) THEN
      I = PERM(I)
      COUNT(I) = COUNT(I) + 1
      ELSE
      J = PERM(J)
      COUNT(J) = COUNT(J) + 1
      END IF

   20 CONTINUE
C
C Store row counts in LROW and accumulate row counts in COUNT.
    LROW(1) = COUNT(1)
    IDIAG(1) = MAX(IDIAG(1),1)
    DO 30 I = 2,N
      LROW(I) = COUNT(I)
C COUNT(I) set to point to position of last entry in row I of permuted
C     upper triangle.
      COUNT(I) = COUNT(I-1) + LROW(I)
      IDIAG(I) = COUNT(I-1) + MAX(IDIAG(I),1)
   30 CONTINUE

C Set diagonal entries in IRNPRM.  This is done separately because some
C   diagonals may not be present in the users input.
    DO 35 I = 1,N
      K = PERM(I)
      IRNPRM(IDIAG(K)) = I
   35 CONTINUE
C
C Run through putting the entries in the right places. COUNT is used for
C   holding running pointers and is left holding pointers to row
C   starts.
C Count number of entries in expanded matrix (allowing for non-input
C   diagonals)
    DO 40 K = 1,NE
      I = IRN(K)
      J = JCN(K)
      IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) THEN
      MAP(K) = 0
      GO TO 40
      ENDIF
      I = PERM(IRN(K))
      J = PERM(JCN(K))
      IF (I.EQ.J) THEN
      MAP(K) = IDIAG(I)
      IRNPRM(IDIAG(I)) = IRN(K)
      IDIAG(I) = IDIAG(I) - 1
      ELSE
      IF (I.GT.J) THEN
C Entry has been mapped to lower triangle
        IF (KIND.EQ.1) THEN
          MAP(K) = COUNT(J)
        ELSE IF (KIND.NE.1) THEN
          MAP(K) = - COUNT(J)
        ENDIF
        IRNPRM(COUNT(J)) = IRN(K)
        COUNT(J) = COUNT(J) - 1
      ELSE
        MAP(K) = COUNT(I)
        IRNPRM(COUNT(I)) = JCN(K)
        COUNT(I) = COUNT(I) - 1
      ENDIF
C       II = MIN(I,J)
C       MAP(K) = COUNT(II)
C       COUNT(I) = COUNT(II) - 1
      ENDIF
   40 CONTINUE

C Set number of entries in expanded matrix
    IDIAG(1) = EXPNE
    RETURN

    END
    SUBROUTINE ME57ND(N,LENR,NA,NE,ND,NSTEPS,LSTKI,LSTKR,
   *            INFO,RINFO)
C
C Storage and operation count evaluation.
C
C Evaluate number of operations and space required by factorization
C   using ME57BD.  The values given are exact only if no numerical
C   pivoting is performed.
C
C N must be set to the matrix order. It is not altered.
C LENR is number of entries in upper triangle of each row of permuted
C   matrix. It includes diagonal (no duplicates) and duplicates
C   off-diagonal but excludes any out-of-range entries.
C NA,NE,ND must be set to hold, for each tree node, the number of stack
C   elements assembled, the number of eliminations and the size of
C   the assembled front matrix respectively.  They are not altered.
C NSTEPS must be set to hold the number of tree nodes. It is not
C   altered.
C LSTKI is used as a work array by ME57ND.
C LSTKR is used as a work array by ME57ND.
C
C Counts for operations and storage are accumulated in variables
C   OPS,OPSASS,NRLTOT,NIRTOT,NRLNEC,NIRNEC,NRLADU,NIRADU.
C OPS number of multiplications and additions during factorization.
C OPSASS number of multiplications and additions during assembly.
C NRLADU,NIRADU real and integer storage respectively for the
C   matrix factors.
C NRLTOT,NIRTOT real and integer storage respectively required
C   for the factorization if no compresses are allowed.
C NRLNEC,NIRNEC real and integer storage respectively required for
C   the factorization if compresses are allowed.
C MAXFRT is maximum front size

C   .. Scalar Arguments ..
    INTEGER N,NSTEPS
C   ..
C   .. Array Arguments ..
    INTEGER LENR(N),LSTKI(N),LSTKR(N),NA(NSTEPS),
   +      ND(NSTEPS),NE(NSTEPS),INFO(40)
    DOUBLE PRECISION RINFO(20)
C   ..
C   .. Local Scalars ..
    INTEGER I,IORG,ISTKI,ISTKR,ITOP,ITREE,JORG,K,
   +      LSTK,NASSR,NELIM,NFR,NSTK,NTOTPV,NZ1,NZ2
    DOUBLE PRECISION DELIM
C   ..
C   .. Intrinsic Functions ..
    INTRINSIC MAX
C   ..
    DOUBLE PRECISION OPS,OPSASS
    INTEGER NIRADU,NIRNEC,NIRTOT,NRLADU,NRLNEC,NRLTOT,MAXFRT

C   ..
C   .. Executable Statements ..
C
C Accumulate number of nonzeros with indices in range in NZ1.
C   Duplicates on the diagonal are ignored but NZ1 includes any
C   diagonals not present on input and duplicates off the diagonal.
    NZ1 = 0
    DO 40 I = 1,N
      NZ1 = NZ1 + LENR(I)
   40 CONTINUE
    NZ2 = NZ1
C ISTKR,ISTKI Current number of stack entries in
C   real and integer storage respectively.
C OPS,OPSASS,NRLADU,NIRADU,NIRTOT,NRLTOT,NIRNEC,NRLNEC,NZ2 are defined
C   above.
C NZ2 Current number of original matrix entries not yet processed.
C NTOTPV Current total number of rows eliminated.
C ITOP Current number of elements on the stack.
    ISTKI = 0
    ISTKR = 0
    OPS = 0.0D0
    OPSASS = 0.0D0
    NRLADU = 0
C One location is needed to record the number of blocks actually used.
    NIRADU = 3
C Initialize to what is required in ME57BD (as opposed to ME57OD).
    NIRTOT = NZ1+N+5
    NRLTOT = NZ1
    NIRNEC = NZ2+N+5
    NRLNEC = NZ2
    NTOTPV = 0
    ITOP = 0
    MAXFRT = 0
C
C Each pass through this loop processes a node of the tree.
    DO 100 ITREE = 1,NSTEPS
      NELIM = NE(ITREE)
      DELIM = NELIM
      NFR = ND(ITREE)
      MAXFRT = MAX(MAXFRT,NFR)
      NSTK = NA(ITREE)
C Adjust storage counts on assembly of current frontal matrix.
      NASSR = NELIM*(NELIM+1)/2 + NFR*NFR
C Data for no compresses so use original number, NZ1.
      NRLTOT = MAX(NRLTOT,NRLADU+NASSR+ISTKR+NZ1)
C Data for compresses so use current number, NZ2.
      NRLNEC = MAX(NRLNEC,NRLADU+NASSR+ISTKR+NZ2)

C Decrease NZ2 by the number of entries in rows being eliminated at
C   this stage.
      DO 70 IORG = 1,NELIM
      JORG = NTOTPV + IORG
      OPSASS = OPSASS + LENR(JORG)
      NZ2 = NZ2 - LENR(JORG)
   70   CONTINUE

      NTOTPV = NTOTPV + NELIM

C Remove elements from the stack.  There are ITOP elements on the
C   stack with the appropriate entries in LSTKR and LSTKI giving
C   the real and integer storage respectively for each stack
C   element.
      DO 80 K = 1,NSTK
      LSTK = LSTKR(ITOP)
      ISTKR = ISTKR - LSTK
      OPSASS = OPSASS + LSTK
      LSTK = LSTKI(ITOP)
      ISTKI = ISTKI - LSTK
      ITOP = ITOP - 1
   80   CONTINUE

C Accumulate nonzeros in factors and number of operations.
      NRLADU = NRLADU + (NELIM*(NELIM+1))/2 + (NFR-NELIM)*NELIM
      NIRADU = NIRADU + 2 + NFR
      OPS = OPS + (DELIM* (12*NFR+6*NFR*NFR - (DELIM+1)*
   +      (6*NFR+6-(2*DELIM+1))))/6 + DELIM

      IF (NFR.GT.NELIM) THEN
C Stack remainder of element.
      ITOP = ITOP + 1
      LSTKR(ITOP) = ((NFR-NELIM)*(NFR-NELIM+1))/2
      LSTKI(ITOP) = NFR - NELIM + 1
      ISTKI = ISTKI + LSTKI(ITOP)
      ISTKR = ISTKR + LSTKR(ITOP)
      ENDIF

C Adjust integer counts to stack elements and allow for next front.
      IF (ITREE.EQ.NSTEPS) THEN
      NIRTOT = MAX(NIRTOT,NIRADU+ISTKI+NZ1)
      NIRNEC = MAX(NIRNEC,NIRADU+ISTKI+NZ2)
      ELSE
      NIRTOT = MAX(NIRTOT,NIRADU+(N-NTOTPV+2)+ISTKI+NZ1)
      NIRNEC = MAX(NIRNEC,NIRADU+(N-NTOTPV+2)+ISTKI+NZ2)
      ENDIF

  100 CONTINUE
C

C Set INFO and RINFO
    INFO(5)   = NRLADU
    INFO(6)   = NIRADU
    INFO(7)   = MAXFRT
    INFO(8)   = NSTEPS
    INFO(9)   = NRLTOT
    INFO(10)  = NIRTOT
    INFO(11)  = NRLNEC
    INFO(12)  = NIRNEC
    RINFO(1)  = OPSASS
    RINFO(2)  = OPS

    RETURN
    END
    SUBROUTINE ME57OD(KIND,N,NE,A,LA,IW,LIW,LROW,PERM,NSTEPS,NSTK,
   +            NODE,DIAG,SCHNAB,PPOS,CNTL,ICNTL,INFO,RINFO,
   +            HOLD,BIGA)
C
C Factorization subroutine
C
C This subroutine operates on the input matrix ordered into a tentative
C   pivot order by ME57BD and produces the matrices U and inv(D)
C   of the factorization A = (U trans) D U, where D is a block
C   diagonal matrix with blocks of order 1 and 2. Gaussian elimination
C   is used with pivots of order 1 and 2, chosen according to the
C   tentative pivot order unless stability considerations
C   require otherwise.
    INTEGER KIND,N,NE,LA
C SCHNAB is in fact of dimension 5, but must set to * for Fujitsu (and
C   maybe other) compilers on very small matrices.
    COMPLEX*16 A(LA),DIAG(N),SCHNAB(*)
    DOUBLE PRECISION CNTL(5),RINFO(20),BIGA
    INTEGER LIW,IW(LIW),LROW(N),PERM(N),NSTEPS,NSTK(NSTEPS),
   +      NODE(N),PPOS(N),ICNTL(20),INFO(40),HOLD(40)
C KIND is an INTEGER variable that must be set by the user to indicate
C    whether the matrix is complex symmetric (KIND=1), or complex
C    Hermitian (KIND not equal to 1).
C N   must be set to the order of the matrix. It is not altered.
C NE  must be set to the number of entries in the upper triangle of
C   the permuted matrix. It is not altered.
C A   must be set so that the upper triangle of the permuted matrix
C   is held by rows in positions LA-NE+1 to LA. Explicit zero
C   diagonals are stored. Duplicate entries are permitted
C   and are summed. During the computation, active frontal matrices
C   are held in A by rows. The working front is held in full form.
C   Stacked elements are held by rows in packed form.
C   On exit, entries 1 to INFO(10) of A hold real information on the
C   factors and should be passed unchanged to ME57CD. For each block,
C   the factorized block pivot precedes the rows of the out-of-pivot
C   part held as a rectangular matrix by rows.
C      The factorized pivot has the form:
C          -1  T
C        L  D   L
C
C   where L is unit lower triangular, D is block diagonal with blocks
C   of size 1 or 2. L and the diagonal part of D is held packed by
C   columns and the off-diagonal entries of the 2*2 blocks of D are
C   held from position IW(1).
C LA  length of array A. A value for LA sufficient for the
C   tentative pivot sequence will have been provided by ME57AD
C   in variable INFO(11). LA is not altered.
C IW  must be set on input so that IW(I+LIW-NE) holds the column
C   index of the entry in A(I+LA-NE) for I = 1,..., NE.
C   On exit, entries 1 to INFO(16) hold integer information on the
C   factors and should be passed unchanged to ME57CD.
C   IW(1) will be set to one greater than the number of entries in the
C   factors.
C   IW(2) points to the end of the factorization.
C   IW(3) will be set to the number of block pivots actually used;
C   this may be different from NSTEPS since numerical considerations
C   may prevent us choosing pivots at some stages.
C   Integer information on each block pivot row follows. For each
C   block pivot row, we have:
C     * no. of columns,
C     * no. of rows,
C     * list of column indices. The column indices for a
C       2x2 pivot are flagged negative.
C   During the computation, the array is used to hold indexing
C   information on stacked elements.  IW stores the number of
C   variables then a list of variables in the element.
C LIW must be set to the length of array IW. A value sufficient for the
C   tentative pivot sequence will have been provided by ME57AD in
C   variable INFO(12). LIW is not altered.
C LROW  must be set so that LROW(I) holds the number of entries in row
C   I (in permuted order) of the incoming matrix. LROW is not altered.
C PERM  must be set so that PERM(I) holds the variable that is Ith in
C   the tentative pivot order generated by ME57AD. PERM is not
C   altered.
C NSTEPS must be set to the number of nodes in the tree from the
C   analysis. It is the length of array NSTK. Its
C   value will never exceed N. It is not altered.
C NSTK must be set so that NSTK(I) holds the number of
C   stacked elements to be assembled at tree node I.
C NODE must be unchanged since return from ME57AD. NODE(I) gives
C   the tree node at which variable I was eliminated in the
C   analysis phase. It is of length N and is not altered.
C DIAG is only accessed if ICNTL(7) is equal to 4.
C   In that case it must be set to the values of the diagonals of
C   matrix.
C SCHNAB is only accessed if ICNTL(7) is equal to 4. It is used to hold
C   parameters for the Schnabel-Eskow modification and max/min values
C   of current diagonal.  Specifically (using notation of S-E):
C   SCHNAB(1) == GAMMA
C   SCHNAB(2) == TAUBAR (in fact root of TAUBAR)
C   SCHNAB(3) == MU
C   SCHNAB(4) == Max entry on diag
C   SCHNAB(5) == Min entry on diag
C This option is only invoked if the matrix is Hermitian
C PPOS  is integer work array of dimension N. If I is a variable in
C   the current front, PPOS(I) is used to indicate its position in the
C   front. For any other uneliminated variable, PPOS(I) is set to N+1.
C CNTL must be set (perhaps by ME57ID) so that CNTL(1) holds the
C   pivot threshold and CNTL(2) holds the pivot tolerance.
C ICNTL must be set (perhaps by ME57ID).  Entries of ICNTL accessed
C   by ME57OD are:
C ICNTL(2) is output unit for warning messages.
C ICNTL(7) is used to control pivoting.  With the default value of 1,
C   1 x 1 and 2 x 2 pivots are used subject to passing a threshold
C   tolerance.  If ICNTL(7) is greater than 1 only 1 x 1 pivots will
C   be used. If ICNTL(7) equal to 2,
C   the subroutine will exit immediately a sign change or zero pivot
C   is detected.  If ICNTL(7) is equal to 3, the subroutine will
C   continue the factorization unless a zero
C   pivot is detected.  If ICNTL(7) is equal to 4, the diagonal of
C   the matrix will be modified so that all pivots are of the same.
C ICNTL(8) is used to control whether, on running out of space, the
C   subroutine exits with an error return (ICNTL(8) = 0), or
C   whether it saves some internal variables so that
C   larger arrays can be allocated and the computation restarted
C   from the point at which it failed.
C ICNTL(11) is the block size used by the Level 3 BLAS (default 32).
C RINFO(3) will be set to the number of floating-point operations
C   required for the assembly.
C RINFO(4) will be set to the number of floating-point operations
C   required for the factorization.  RINFO(5) will be set to the
C   number of extra flops needed for the use of GEMM.
C INFO(1)  holds a diagnostic flag. It need not be set on entry. A zero
C   value on exit indicates success. Possible nonzero values are
C      -3  insufficient storage for A.
C      -4  insufficient storage for IW.
C      -5  zero pivot found when ICNTL(7) = 2 or 3.
C      -6  change in sign of pivots when ICNTL(7) = 2.
C        KIND out of range
C      +4  matrix is singular
C     +10  factorizations pauses because insufficient real space
C     +11  factorizations pauses because insufficient integer space
C INFO(40) is used to accumulate number of reals discarded if
C   elimination continues without restart when space is exhausted.
C

C INFO(32) is set to number of zeros in the triangle of the factors
C INFO(33) is set to number of zeros in the rectangle of the factors
C INFO(34) is set to number of zero columns in rectangle of the factors
C Needed to compute these
    INTEGER ZCOL,RPOS

C
C Constants
    DOUBLE PRECISION ZERO,HALF,ONE
    PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
    COMPLEX*16 CONE
    PARAMETER (CONE=(1.0D0,0.0D0))
C
C Local variables
    INTEGER AINPUT
    COMPLEX*16 AMULT1,AMULT2
    DOUBLE PRECISION AMAX
    INTEGER APOS,APOSA,APOSB,APOSBB,APOSBK,APOSC,APOSI,APOSJ,APOSM,
   +      APOS1,APOS2,APOS3,APOS4,ASTK,ATRASH,BLK
    DOUBLE PRECISION DELTA
    COMPLEX*16 DETPIV
    INTEGER ELT
    DOUBLE PRECISION FLOPSA,FLOPSB,FLOPSX
    INTEGER I,I1,IASS,IBEG,IELL,IEND,IEXCH,IINPUT,INTSPA,
   +      IORG,IPIV,IPOS,IROW,ISNPIV,ISTK,ISWOP,IWNFS,IWPOS,
   +      J,JAY,JA1,JCOL,JJ,JJJ,JMAX,J1,J2,K,
   +      KB,KBLK,KCT,KR,KROW,K1,K2,L,LASPIV,LDIAG,LIELL,
   +      LP,LPIV, NBSTATIC
    LOGICAL LASTBK,LTWO
    INTEGER MAXFRT
    DOUBLE PRECISION MAXPIV
    INTEGER NASS,NBLK,NBLOC,NCMPBI,NCMPBR,NEIG,NELL,NFRONT,NIRBDU
    DOUBLE PRECISION NORMJ
    INTEGER NTWO
C LSTAT is .TRUE. is we are using static pivoting
    LOGICAL SCHUR,LSTAT
    INTEGER MPIV,NPIV,NPOTPV,NRLBDU,NSC1,NST,
   +      NSTACK(2),NSTKAC(2),NTOTPV,
   +      NUMORG,OFFDAG,PHASE,PIVBLK
    COMPLEX*16 PIVOT
    INTEGER PIVSIZ,POSELT,POSPV1,POSPV2,PTRA,PTRIRN,RLSPA,
   +      SIZBLK,SIZC,SIZF,TRLSPA,TINSPA,TOTSTA(2),WP,ZCOUNT
    DOUBLE PRECISION RMAX,TMAX,TOL,UU,ULOC,UTARG,STCTOL
    COMPLEX*16 SWOP,CMAX

C AINPUT is the position of the first entry of the original matrix
C   reals since the last compress.  It is reset to the current
C   part being processed when ME57PD is called.
C AMAX is used to record the largest entry in a row.
C CMAX is the value of the complex number corresponding to AMAX
C AMULT1, AMULT2 are used to hold multipliers.
C   Also used as temporary variables.
C APOS is a pointer to the start of the current front in A.
C APOSA holds the index of the current entry of the matrix A used to
C   form the Schur complement as C = A*B
C APOSB holds the index of the current entry of the matrix B used to
C   form the Schur complement as C = A*B
C APOSBB is a pointer to the start of the buffer before current front
C   in A. It is different from APOS because of the way we store the
C   factors and prevents us overwriting when generating them.
C APOSC holds the index of the current entry of the matrix C used to
C   form the Schur complement as C = A*B
C APOSI is set to the beginning of a row in the working front or just
C    ahead of this position.
C APOSJ is set to the beginning of row JMAX in working front.
C APOSM holds the index of the current entry of the multiplier matrix
C   for the matrix B used to form the Schur complement as C = A*B for
C   the part of the pivot rows outside the pivot block.
C APOS1 is a position in the array A.
C APOS2 is a position in the array A.
C APOS3 is the position in A of the first entry in the copy of U.
C APOS4 is the position in A of the first entry in the Schur complement.
C ASTK indicates position immediately before first stack entry.
C ATRASH is used as limit on space in A being set to zero.
C BLK is DO loop variable for blocks in pivot block.
C DELTA is the amount added to diagonal when in Phase 2 of matrix
C   modification.
C DETPIV is the value of the determinant of the 2x2 pivot or
C   candidate pivot.  It is real when matrix is Hermitian (only).
C ELT is a DO loop variable indicating the element being processed.
C FLOPSA  counts floating-point operations for assembly
C FLOPSB  counts floating-point operations for elimination
C FLOPSX  counts extra flops required by use of GEMM
C I is a DO loop variable.
C I1 is used to hold limit of DO loop.
C IASS is the index of the current tree node.
C IBEG is the position of the beginning of a row in A, used when
C   performing elimination operations on the front.
C IELL Current element being assembled starts in position IELL of IW.
C IEND is the position of the end of a row in A, used when
C   performing elimination operations on the front.
C IEXCH is used to hold the contents of an array entry when doing a
C   swop.
C IINPUT is the position of the first entry of the original matrix
C   integers since the last compress.  It is reset to the current
C   part being processed when ME57PD is called. with REAL = .FALSE.
C INTSPA is amount of integer workspace needed to progress to current
C   point in elimination.
C IORG is a DO loop variable indicating which current row from the
C   incoming matrix is being processed.
C IPIV is the current relative position of the pivot search.
C IPOS is used to hold (temporarily) a position in IW.
C IROW is a DO loop variable used when scanning a row of A.
C ISNPIV is +1 if first pivot is positive and -1 if it is negative.
C ISTK points to the bottom of the stack in IW (needed by compress).
C    It indicates position immediately before first stack entry.
C ISWOP is used when swopping two integers.
C IWNFS points to the first free location for a variable that is not
C   fully summed.
C IWPOS points to the first free position for factors in IW.
C J is a temporary variable.
C JA1 is a temporary index.
C JCOL is used as an index into array A.
C JJ and JJJ are Do loop indices.
C JMAX is the relative column index in the front of the largest
C   off-diagonal in the fully summed part of the prospective pivot
C   row.
C J1 and J2 are pointers to the beginning and end of a row segment
C   in the array IW.  They are also used as running pointers in IW.
C K is a temporary variable.
C KBLK Blocked GEMM is performed on a block KBLK by KBLK matrix.
C KCT counts the number of unchecked candidates in a pivot sweep.
C KR is pointer to the current row in the assembled block being
C   tested for a potential pivot.
C KROW is a DO loop variable.
C K1 and K2 are used as running indices for array A.
C KB is block row index for blocked GEMM.
C L is a temporary variable.
C LASPIV is set to value of NPIV at end of previous block of pivots.
C LIELL is the order of the size of the reduced matrix from the front.
C   This is the order of the stacked matrix.
C LPIV is the number of pivots selected in a block pivot.
C LASTBK is flag to indicate when we are processing last block of
C     pivot block.
C LTWO  is logical variable used to indicate if current pivot is a 2 x 2
C     pivot.
C MAXFRT is the maximum front size encountered so far.
C MAXPIV is the largest of two diagonal entries of a 2 x 2 pivot.
C NASS holds the number of fully assembled variables in
C   the newly created element.
C NBLK is the number of block pivots used.
C NBLOC Number of rows of the Schur complement calculated by each GEMM
C   reference. Set to ICNTL(11).
C NCMPBR, NCMPBI are the number of compresses on real and integer space
C   respectively.
C NEIG is number of negative eigenvalues detected.
C NELL is used to hold the current number of son elements.
C NFRONT is the number of variables in the front.
C NIRBDU is number of integer entries in factors.
C NORMJ is used in matrix modification for 1-norm of off-diagonals.
C NTWO is the number of two by two full pivots used.
C SCHUR if set to .TRUE. then the Schur complement will be
C    generated using Level 3 BLAS.
C MPIV is the number of pivots so far chosen in the current block.
C NPIV is the number of pivots so far chosen at the current node.
C NPOTPV is the total number of potential pivots so far.  Variables
C   PERM(1), PERM(2), .... PERM(NPOTPV) are fully assembled.
C NRLBDU is number of real entries in factors.
C NST temporary to hold NSC1+1 if NSC1 > 0, 0 otherwise.
C NSTACK(I), I = 1,2 hold the number of active entries on the
C   real/integer stack.
C NSTKAC(I), I =1,2 hold the number of entries on the real/integer
C   stack and original matrix after a compress.
C NTOTPV is the total number of pivots selected. This is used
C   to determine whether the matrix is singular.
C NUMORG is the number of variables in the tentative pivot from the
C   incoming original rows.
C OFFDAG is the position in A of the off-diagonal entry of a 2 x 2
C      pivot.
C PIVBLK Number of rows of each block when using GEMM in pivot block.
C   Set to minimum of NBLOC and NASS.
C PIVOT is a temporary variable used to hold the value of the current
C      pivot.
C PIVSIZ is order of current pivot (has value 1 or 2).
C POSELT is a pointer into the current element being assembled.
C POSPV1 is the position in A of a 1 x 1 pivot or the first diagonal
C   of a 2 x 2 pivot.
C POSPV2 is the position in A of the second diagonal of a 2 x 2 pivot.
C PTRA  points to the next original row in A.
C PTRIRN points to the next original row in IW.
C RLSPA is amount of real workspace needed to progress to current
C   point in elimination.
C SIZBLK Number of rows in current block when using GEMM in pivot block.
C   Set to minimum of PIVBLK and NASS - NPIV
C SIZC is set to number of rows in remainder of pivot block after
C   blocking is done.
C SIZF is set to number of rows in remainder of pivot row (to NFRONT)
C   after blocking is done.
C TOTSTA(I), I =1,2 hold the number of entries on the stack and
C   original matrix.
C TRLSPA is amount of real workspace needed to progress to current
C   point in elimination if no compresses on the data are performed.
C RMAX is used to record the largest entry in a row.
C SWOP is used when swopping two reals.
C TMAX is used to record the largest entry in a row.
C TOL is the tolerance against which singularity is judged.
C   If static pivoting is used, then TOL is the value for this.
C UU is a local variable used to hold threshold parameter.  Its value is
C   between 0 and 0.5.
C ZCOUNT is number of "zero" rows in current front (used when
C   ICNTL(16) is equal to 1).
C
C Procedures
C ME57PD compresses arrays.
C ME57WD  adjusts signs in factors, moves the off-diagonal entries of
C    full 2x2 pivots and updates counts.
    DOUBLE PRECISION FD15AD

C?? To identify bug
C   LOGICAL LCASE
C   COMMON /CCASE/LCASE

    INTRINSIC MIN,MAX,ABS,REAL,CONJG
    EXTERNAL ZGEMM,FD15AD,ME57PD,ME57WD

C
C Initialization.
    NBLOC = ICNTL(11)
    TOL = CNTL(2)
    LP = ICNTL(1)
    WP = ICNTL(2)
    LDIAG = ICNTL(5)
    INFO(40) = 0
C A local variable UU is used for the threshold parameter, so that
C   CNTL(1) will remain unaltered.
    UU = MIN(CNTL(1),HALF)
    UU = MAX(UU,ZERO)

    LSTAT = .FALSE.
C Check if static pivoting option is on
    IF (CNTL(4).GT.ZERO) THEN
C LSTAT is not now set until number of delayed pivots is CNTL(5)*N
      IF (CNTL(5).EQ.ZERO) LSTAT = .TRUE.
      UTARG = SQRT(UU/CNTL(4))*CNTL(4)
      STCTOL = BIGA*CNTL(4)
C     TOL = STCTOL
    ENDIF

C Action if we are returning in the middle of the factorization
    IF (HOLD(1).GT.0) THEN
      INFO(1) = 0
      NBLK = HOLD(2)
      NTWO = HOLD(3)
      INFO(23) = HOLD(4)
      NCMPBR = 0
      NCMPBI = 0
      NEIG   = HOLD(6)
      MAXFRT = HOLD(7)
C Test compiler by commenting this out
      IWPOS  = HOLD(8)
      APOS   = HOLD(9)
      APOSBB = HOLD(10)
      NSTKAC(1) = HOLD(11)
      NSTKAC(2) = HOLD(12)
      AINPUT  = HOLD(13)
      IINPUT  = HOLD(14)
      ISTK    = HOLD(15)
      ASTK    = HOLD(16)
      INTSPA  = HOLD(17)
      RLSPA   = HOLD(18)
      PTRIRN  = HOLD(19)
      PTRA    = HOLD(20)
      NTOTPV  = HOLD(21)
      NPOTPV  = HOLD(22)
      NUMORG  = HOLD(23)
      NFRONT  = HOLD(24)
      NASS    = HOLD(25)
C     NCOL    = HOLD(26)
      IF (HOLD(1).EQ.1) NELL    = HOLD(27)
      IF (HOLD(1).EQ.2) NPIV    = HOLD(27)
      IASS    = HOLD(28)
      TINSPA  = HOLD(29)
      TRLSPA  = HOLD(30)
      TOTSTA(1) = HOLD(31)
      TOTSTA(2) = HOLD(32)
      NSTACK(1) = HOLD(33)
      NSTACK(2) = HOLD(34)
      INFO(32)  = HOLD(37)
      INFO(33)  = HOLD(38)
      INFO(34)  = HOLD(39)
      NBSTATIC  = HOLD(40)
      IF (ICNTL(7).GT.1) ISNPIV = HOLD(35)
      IF (ICNTL(7).EQ.4) PHASE = HOLD(36)
      IF (HOLD(1).EQ.2) NSC1    = NFRONT-NPIV
      FLOPSA = RINFO(3)
      FLOPSB = RINFO(4)
      FLOPSX = RINFO(5)
      IF (HOLD(1).EQ.1) THEN
C Real arrays expanded
      HOLD(1) = 0
      GO TO 333
      ELSE
      IF (HOLD(1).EQ.3) THEN
C We ran out of space when allocating values for zero pivots at the end
C       of the factorization.  This is most likely to happen when we
C       are running with ICNTL(16) equal to 1 (dropping small entries
C       from front).
        HOLD(1) = 0
        GO TO 555
      ELSE
C Integer arrays expanded
        HOLD(1) = 0
        GO TO 444
      ENDIF
      ENDIF
    ENDIF

C NBSTATIC is the number of modified diagonal entries
    NBSTATIC = 0
C NBLK is the number of block pivots used.
    NBLK = 0
C NTWO is the number of 2 x 2 pivots used.
    NTWO = 0
C NCMPBR, NCMPBI are the number of compresses on real and integer space
C   respectively.
    NCMPBR = 0
    NCMPBI = 0
C FLOPSA is the number of floating-point operations for assembly.
    FLOPSA = ZERO
C FLOPSB is the number of floating-point operations for elimination.
    FLOPSB = ZERO
C FLOPSX  counts extra flops required by use of GEMM
    FLOPSX = ZERO
C NEIG is number of negative eigenvalues detected.
    NEIG = 0
C MAXFRT is the maximum front size encountered so far.
    MAXFRT  = 0
C All relevant INFO and RINFO parameters initialized to zero
C   so they have a valid entry on any error return.
    INFO(1) = 0
    INFO(2) = 0
    INFO(14:29) = 0
    INFO(31:35) = 0
    RINFO(3:5) = ZERO
    RINFO(14:15) = ZERO

C Initialization of array indicating positions of variables in front
    DO 10 I = 1,N
      PPOS(I) = N + 1
   10 CONTINUE
C IWPOS is set to position for first index of first block pivot
    IWPOS = 6
C Set first five entries to dummies to avoid unassigned var in ME57ED
    IW(1) = 0
    IW(2) = 0
    IW(3) = 0
    IW(4) = 0
    IW(5) = 0
C APOSBB is a pointer to the next position for storing factors in A
    APOSBB = 1
C Initialize NSTKAC and INTSPA and RLSPA
    NSTACK(1) = 0
    NSTACK(2) = 0
    NSTKAC(1) = NE
    NSTKAC(2) = NE
    TOTSTA(1) = NE
    TOTSTA(2) = NE
    INTSPA = NE+5+N
    RLSPA = NE
    TINSPA = NE+5+N
    TRLSPA = NE
C PTRIRN points to the next original row in IW.
    PTRIRN = LIW - NE + 1
C PTRA  points to the next original row in A.
    PTRA = LA - NE + 1
C ISTK points to the position in IW immediately before the stack.
    ISTK = PTRIRN - 1
C ASTK points to the position in A immediately before the stack.
    ASTK = PTRA - 1
C AINPUT is the position of the first entry of the original matrix
C   reals since the last compress.  It is reset to the current
C   part being processed when ME57PD is called with REAL .TRUE.
    AINPUT = PTRA
C IINPUT is the position of the first entry of the original matrix
C   integers since the last compress.  It is reset to the current
C   part being processed when ME57PD is called. with REAL = .FALSE.
    IINPUT = PTRIRN
C NTOTPV is the total number of pivots selected.
    NTOTPV = 0
C NPOTPV is the total number of potential pivots so far.
    NPOTPV = 0
C In case we run out of space before the first pivot is chosen, we
C   must initialize ISNPIV.
    IF (ICNTL(7).GT.1) ISNPIV = 0


C Must initialize PHASE here because of ICNTL(8) (OLD) return
C when KIND = 1.
    PHASE = 1
C Calculate diagonal of matrix and store in DIAG
    IF (ICNTL(7).EQ.4 .AND. KIND.NE.1) THEN
      DO 19 I = 1,N
      DIAG(I) = ZERO
   19   CONTINUE
      APOS1 = PTRA-1
      J1 = PTRIRN
      DO 20 I = 1,N
      J2 = J1 + LROW(I) - 1
      DO 25 JJ = J1,J2
        J = IW(JJ)
        APOS1 = APOS1 + 1
        IF (J.EQ.PERM(I)) DIAG(J) = DIAG(J) + REAL(A(APOS1))
   25   CONTINUE
      J1 = J2 + 1
   20   CONTINUE
      SCHNAB(1) = ONE
      SCHNAB(5) = ZERO
      DO 21 I = 1,N
      SCHNAB(1) = MAX(REAL(SCHNAB(1)),REAL(DIAG(I)))
      SCHNAB(5) = MIN(REAL(SCHNAB(5)),REAL(DIAG(I)))
   21   CONTINUE
C Set max entry on diag
      SCHNAB(4) = SCHNAB(1)
C S+E has **2/3 .. Nick wants **1/3
      SCHNAB(2) = FD15AD('E')**(1.0/3.0)
      SCHNAB(3) = 0.1
C Initialize RINFO(15) to compute smallest pivot in modified matrix
      RINFO(15) = FD15AD('H')
      DELTA   = ZERO
    ENDIF

C   *****************************************************************
C   * Each pass through this main loop performs all the operations  *
C   * associated with one node of the assembly tree.          *
C   *****************************************************************

    IASS = 1
C   DO 2160 IASS = 1,NSTEPS
 2160 CONTINUE

C Find the frontal variables, ordered with the fully summed variables
C   of the incoming rows first, the fully summed rows from previous
C   steps next, followed by the rest in any order.

C NUMORG is the number of variables in the tentative pivot from the
C   incoming original rows.
C Calculate NUMORG and put indices of these fully summed rows in IW.
      NUMORG = 0
      DO 30 I = NPOTPV + 1,N
C J is Ith variable in tentative pivotal sequence.
      J = PERM(I)
C Jump if we have finished with variables in current node.
      IF (ABS(NODE(J)).GT.IASS) GO TO 40
      IW(IWPOS+NUMORG) = J
      NUMORG = NUMORG + 1
      PPOS(J) = NUMORG
   30   CONTINUE

C NASS will be set to the total number of fully assembled variables in
C   the newly created element. First set it to NUMORG.
   40   NASS = NUMORG
C Add indices of fully summed variables of stacked sons to IW.
      NELL = NSTK(IASS)
      IELL = ISTK + 1
      DO 70 ELT = 1,NELL
      DO 50 JJ = IELL + 1,IELL + IW(IELL)
        J = IW(JJ)
        IF (NODE(J).GT.IASS) GO TO 50
C Jump if variable already included.
        IF (PPOS(J).LE.N) GO TO 50
        IW(IWPOS+NASS) = J
        NASS = NASS + 1
        PPOS(J) = NASS
   50   CONTINUE
      IELL = IELL + IW(IELL) + 1
   70   CONTINUE
C IWNFS points to the first free location for a variable that is not
C   fully summed.
      IWNFS = IWPOS + NASS

C Incorporate original rows.
C J1 is the position of the start of the first original row associated
C   with this node of the assembly tree.
      J1 = PTRIRN
      DO 90 IORG = 1,NUMORG
      J2 = J1 + LROW(NPOTPV+IORG) - 1
C Run through index list of original row.
      DO 80 JJ = J1,J2
        J = IW(JJ)
C Jump if variable already included.
        IF (PPOS(J).LE.N) GO TO 80
        IW(IWNFS) = J
        IWNFS = IWNFS + 1
        PPOS(J) = IWNFS - IWPOS
   80   CONTINUE
      J1 = J2 + 1
   90   CONTINUE

C Now incorporate stacked elements.
C J1 is set to beginning
C J2 is set to end
      IELL = ISTK + 1
      DO 170 ELT = 1,NELL
      J1 = IELL+1
      J2 = IELL+IW(IELL)
      DO 150 JJ = J1,J2
        J = IW(JJ)
C Jump if already assembled
        IF (PPOS(J).LE.N) GO TO 150
        IW(IWNFS) = J
        IWNFS = IWNFS + 1
        PPOS(J) = IWNFS - IWPOS
  150   CONTINUE
      IELL = J2 + 1
  170   CONTINUE

C NFRONT is the number of variables in the front.
      NFRONT = IWNFS - IWPOS

C MAXFRT is the largest front size so far encountered.
      MAXFRT = MAX(MAXFRT,NFRONT)

C Set APOS to the position of first entry in frontal matrix.
      IF (INFO(1).NE.-3) THEN
C Buffer space allocated so that triangular part of pivot can be stored
C   without danger of overwrite.
      APOS = APOSBB + (NASS*(NASS+1))/2
      ELSE
      APOS = 1
      END IF
C
C Assemble reals into frontal matrix.
C
C Accumulate real space needed
      RLSPA  = MAX(RLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+NSTKAC(1))
      TRLSPA = MAX(TRLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+TOTSTA(1))

C If necessary, compress A.

  333   IF (APOS+NFRONT*NFRONT-1.GT.ASTK) THEN


      CALL ME57PD(A,IW,ASTK,AINPUT,PTRA,.TRUE.)

      NCMPBR = NCMPBR + 1
      IF (APOS+NFRONT*NFRONT-1.GT.ASTK) THEN
        IF (ICNTL(8).NE.0) THEN
C Zero part of A to avoid failure in HSL_MA57
          DO 334 I = APOSBB,ASTK
          A(I) = ZERO
  334       CONTINUE
          HOLD(1) = 1
          HOLD(2) = NBLK
          HOLD(3) = NTWO
          HOLD(4) = INFO(23)
          HOLD(5) = NCMPBI
          HOLD(6) = NEIG
          HOLD(7) = MAXFRT
          HOLD(8) = IWPOS
          HOLD(9) = APOS
          HOLD(10) = APOSBB
          HOLD(11) = NSTKAC(1)
          HOLD(12) = NSTKAC(2)
          HOLD(13) = AINPUT
          HOLD(14) = IINPUT
          HOLD(15) = ISTK
          HOLD(16) = ASTK
          HOLD(17) = INTSPA
          HOLD(18) = RLSPA
          HOLD(19) = PTRIRN
          HOLD(20) = PTRA
          HOLD(21) = NTOTPV
          HOLD(22) = NPOTPV
          HOLD(23) = NUMORG
          HOLD(24) = NFRONT
          HOLD(25) = NASS
C         HOLD(26) = NCOL
          HOLD(27) = NELL
          HOLD(28) = IASS
          HOLD(29) = TINSPA
          HOLD(30) = TRLSPA
          HOLD(31) = TOTSTA(1)
          HOLD(32) = TOTSTA(2)
          HOLD(33) = NSTACK(1)
          HOLD(34) = NSTACK(2)
          IF (ICNTL(7).GT.1) HOLD(35) = ISNPIV
          IF (ICNTL(7).EQ.4) HOLD(36) = PHASE
          HOLD(37) = INFO(32)
          HOLD(38) = INFO(33)
          HOLD(39) = INFO(34)
          RINFO(3) =FLOPSA
          RINFO(4) =FLOPSB
          RINFO(5) =FLOPSX
          HOLD(40) = NBSTATIC
          INFO(35) = HOLD(40)
          INFO(1) = 10
          RETURN
        ELSE
C INFO(40) accumulates number of discards from factors.
          INFO(40) = INFO(40) + APOS - 1
          APOS = 1
          APOSBB = 1
          INFO(1) = -3
          IF (NFRONT*NFRONT.GT.ASTK) THEN
          INFO(17) = MAX(INFO(17),RLSPA)
          IF (ICNTL(7).EQ.4) INFO(17) = MAX(INFO(17),RLSPA + N)
          INFO(2) = LA
          RETURN
          ENDIF
        ENDIF
      ENDIF
      END IF

      ATRASH = APOS + NFRONT*NFRONT - 1
C Zero out appropriate part of A for incoming potential pivot rows.
      DO 210 JJ = APOS,ATRASH
      A(JJ) = ZERO
  210   CONTINUE

C Incorporate reals from original rows.
      J1 = PTRIRN
      DO 230 IORG = 1,NUMORG
C APOSI indicates the position in A just before the beginning of the row
C     being assembled.
      J = PERM(NPOTPV+IORG)
      APOSI = APOS + (PPOS(J)-1)*NFRONT - 1
      J2 = J1 + LROW(NPOTPV+IORG) - 1
      FLOPSA = FLOPSA + J2 - J1 + 1
      DO 220 JJ = J1,J2
        JAY = IW(JJ)
CCC
C Entries always in upper triangle because of ordering.
C Pivot permutations can only affect fully summed variables
        APOS2 = APOSI + PPOS(JAY)
        A(APOS2) = A(APOS2) + A(PTRA)
        PTRA = PTRA + 1
  220   CONTINUE
C       AINPUT = AINPUT + J2 - J1 + 1
      NSTKAC(1) = NSTKAC(1) - J2 + J1 - 1
      J1 = J2 + 1
  230   CONTINUE
C     IINPUT = IINPUT + J2 - PTRIRN
      NSTKAC(2) = NSTKAC(2) - J1 + PTRIRN
      PTRIRN = J1
C Update NPOTPV
      NPOTPV = NPOTPV + NUMORG

C Now assemble reals from stacked elements
C POSELT is a running pointer into that element.
      DO 380 ELT = 1,NELL
      POSELT = ASTK + 1
      LIELL = IW(ISTK+1)
      J1 = ISTK + 2
      J2 = ISTK+1 + LIELL
      FLOPSA = FLOPSA + (LIELL*(LIELL+1))/2
      DO 250 JJ = J1,J2
        J = IW(JJ)
        APOS2 = APOS + (PPOS(J)-1)*NFRONT
        APOS1 = POSELT
        DO 240 JJJ=JJ,J2
          JAY = IW(JJJ)
          IF (PPOS(JAY) .GE. PPOS(J)) THEN
C Upper triangle
          APOS3 = APOS2 + PPOS(JAY) - 1
          A(APOS3) = A(APOS3) + A(APOS1)
          ELSE
          APOS3 = APOS+(PPOS(JAY)-1)*NFRONT+PPOS(J)-1
          IF (KIND.EQ.1) THEN
            A(APOS3) = A(APOS3) + A(APOS1)
          ELSE
            A(APOS3) = A(APOS3) + CONJG(A(APOS1))
          ENDIF
          ENDIF
          APOS1 = APOS1 + 1
  240     CONTINUE
        POSELT = POSELT + LIELL - (JJ-J1)
  250   CONTINUE
C ISTK and ASTK updated to point to posn before next element on stack.
      NSTKAC(2) = NSTKAC(2) - (J2-ISTK)
      NSTACK(2) = NSTACK(2) - (J2-ISTK)
      TOTSTA(2) = TOTSTA(2) - (J2-ISTK)
      ISTK = J2
      ASTK = ASTK + (LIELL*(LIELL+1))/2
      NSTKAC(1) = NSTKAC(1) - (LIELL*(LIELL+1))/2
      NSTACK(1) = NSTACK(1) - (LIELL*(LIELL+1))/2
      TOTSTA(1) = TOTSTA(1) - (LIELL*(LIELL+1))/2
  380   CONTINUE

C       IF (LCASE) THEN
C       write(7,'(/A,I8/A,5I8)') '*** Frontal matrix before step',
C    *    IASS,'NFRONT,NASS,NUMORG,APOS,IWPOS',
C    *      NFRONT,NASS,NUMORG,APOS,IWPOS
C       write(7,'(/A/(10I8))') 'IW array',
C    *    (IW(IWPOS+I-1),I=1,NFRONT)
C       write(7,'(/A)') 'Frontal matrix'
C       DO 1122 J = 1, NFRONT
C       write(7,'(/(6D13.5))')
C    *      (A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT)
C1122   CONTINUE
C       ENDIF

C ******************
C Each time round this loop, we sweep through the remainder
C    of the assembled part of the front looking for pivots.
C ******************

C Set PIVBLK
      PIVBLK = MIN(NBLOC,NASS)
C Set pointer to first entry in block
      APOSBK = APOS
C NPIV is the number of pivots so far selected at the current node.
      NPIV = 0
C Set local value for U
      ULOC = UU

C Each pass through loop processes one block.
      DO 918 BLK = 1,NASS
C Set last block flag
      IF (NPIV+PIVBLK .GE. NASS) THEN
      LASTBK = .TRUE.
      SIZBLK = NASS - NPIV
      ELSE
      LASTBK = .FALSE.
      SIZBLK = PIVBLK
      ENDIF
C Record number of rows processed to date
      LASPIV = NPIV
C MPIV is the number of pivots so far selected in the current block.
      MPIV = 0
C KR is relative position of current pivot candidate
CCCCCCCC
      KR = 0
CCC Set to following to force 2 by 2 pivots in Nocedal examples
C     KR = NUMORG
C KCT is set to one more than the number of unsearched pivot candidates
      KCT = SIZBLK + 1

C**************************
C Loop for pivot searches *
C**************************
  920   CONTINUE
C Increment pointer for circular sweep.
      KR = KR + 1
      KCT = KCT - 1
      IF (KCT.EQ.0) GO TO 930
      IF (KR.GT.SIZBLK) KR = MPIV + 1
C Check to see if no pivot was chosen is complete sweep of pivot block.
C We either take the diagonal entry or the 2 by 2 pivot with the
C   largest fully summed off-diagonal at each stage.
C Note that IPIV is the position within the complete current front.
      IPIV = LASPIV + KR
C APOSI is set to the beginning of row IPIV in working front.
        APOSI = APOS + (IPIV-1)*NFRONT
C Set position and value of potential pivot.
        POSPV1 = APOSI + IPIV - 1
        PIVOT = A(POSPV1)

C Schnabel/Eskow logic
C First check that matrix is Hermitian
   29     IF (KIND.NE.1) THEN
         IF (ICNTL(7).EQ.4) THEN
          IF (PHASE.EQ.2) THEN
C In phase 2 of matrix modification (ICNTL(7) = 4)
C Compute quantity to add to diagonal
C Calculate norm of pivot row
          IF (INFO(27).EQ.0) INFO(27) = NTOTPV + 1
          NORMJ = ZERO
          DO 28 I = POSPV1+1,POSPV1+NFRONT-NPIV-1
            NORMJ = NORMJ + ABS(A(I))
   28       CONTINUE
          DELTA = MAX(ZERO,
   *            - REAL(A(POSPV1))
   *            + MAX(NORMJ,REAL(SCHNAB(2))*REAL(SCHNAB(1))))
          A(POSPV1) = A(POSPV1) + DELTA
          IF (A(POSPV1).EQ.ZERO) GO TO 970
          RINFO(15) = MIN(RINFO(15),REAL(A(POSPV1)))
          DIAG(PERM(NTOTPV+1)) = DELTA
          PIVSIZ = 1
          GO TO 811
          ENDIF
         ENDIF
        ENDIF

        IF (ICNTL(7).GT.1) THEN
C Action if no pivoting requested
          IF (ABS(PIVOT).LE.CNTL(2)) THEN
C Jump if we are not using Schnabel-Eskow modification
          IF (ICNTL(7).LT.4 .OR. KIND.EQ.1) GO TO 970
C We are now in phase 2 of matrix modification (ICNTL(7) = 4)
          PHASE = 2
          GO TO 29
          ENDIF
          IF (NTOTPV.EQ.0) THEN
          IF (REAL(PIVOT).GT.ZERO) ISNPIV = 1
          IF (REAL(PIVOT).LT.ZERO) ISNPIV = -1
          ELSE
          IF (ICNTL(7).EQ.2 .AND. ISNPIV*REAL(PIVOT).LT.ZERO)
   *         GO TO 980
          IF (ICNTL(7).EQ.3 .AND. ISNPIV*REAL(PIVOT).LT.ZERO)
   *         THEN
            INFO(26) = INFO(26) + 1
            ISNPIV = -ISNPIV
          ENDIF
          ENDIF
          IF (ICNTL(7).EQ.4 .AND. KIND.NE.1) THEN
          IF (REAL(PIVOT).GE.REAL(SCHNAB(1))*REAL(SCHNAB(2))
   *          .AND. REAL(SCHNAB(5)) .GE.
   *          -REAL(SCHNAB(3))*REAL(SCHNAB(4))) THEN
C Update and check values of future diagonals
            SCHNAB(5) = ZERO
            SCHNAB(4) = ZERO
            DO 22 I = POSPV1+1,POSPV1+NFRONT-NPIV-1
            J = IW(IWPOS+NPIV+I-POSPV1)
            DIAG(J) = DIAG(J) - REAL(A(I)*A(I)/PIVOT)
            SCHNAB(5) = MIN(REAL(DIAG(J)),REAL(SCHNAB(5)))
            SCHNAB(4) = MAX(REAL(DIAG(J)),REAL(SCHNAB(4)))
            IF (REAL(DIAG(J)).LT.
   *            -REAL(SCHNAB(3))*REAL(SCHNAB(1))) THEN
              PHASE = 2
              GO TO 29
            ENDIF
   22         CONTINUE
            DIAG(PERM(NTOTPV+1)) = ZERO
            RINFO(15) = MIN(RINFO(15),REAL(PIVOT))
          ELSE
            PHASE = 2
            GO TO 29
          ENDIF
          ENDIF
          PIVSIZ = 1
          GO TO 811
        ENDIF

C Numerical pivoting
C Find largest off-diagonal entry in the part of the row in which we
C   seek a pivot.
        AMAX = ZERO
        JMAX = 0
C Split loops in two because only upper triangle is held
C Scan lower triangle by scanning up column of upper triangle.
        DO 110 K = 1, IPIV - NPIV - 1
          IF (ABS(A(POSPV1-K*NFRONT)).GT.AMAX) THEN
          AMAX = ABS(A(POSPV1-K*NFRONT))
          CMAX = A(POSPV1-K*NFRONT)
          JMAX = IPIV - K
          ENDIF
  110     CONTINUE
C Scan upper triangle by scanning along row from first off-diagonal
        DO 111 K =  1, MIN(NASS,LASPIV+PIVBLK) - IPIV
          IF (ABS(A(POSPV1+K)).GT.AMAX) THEN
          AMAX = ABS(A(POSPV1+K))
          CMAX = A(POSPV1+K)
          JMAX = IPIV + K
          ENDIF
  111     CONTINUE
C Do same for the other part.
        RMAX = ZERO

C Restrict partial pivoting check to the fully summed block
C        RMAX = STCTOL
C        DO 112 K = MIN(NASS,LASPIV+PIVBLK)-IPIV+1,NASS
C          RMAX = MAX(RMAX,ABS(A(POSPV1+K)))
C  112     CONTINUE

        DO 112 K = MIN(NASS,LASPIV+PIVBLK)-IPIV+1,NFRONT-IPIV
         RMAX = MAX(RMAX,ABS(A(POSPV1+K)))
 112      CONTINUE

C Action taken if matrix is singular.
        IF (MAX(AMAX,RMAX,ABS(PIVOT)).LE.TOL) THEN
C Skip if all of row is zero.
          GO TO 920
        END IF

C Jump if no nonzero entry in row of pivot block
        IF (MAX(AMAX,ABS(PIVOT)).LE.TOL) GO TO 920
        PIVSIZ = 0
        IF (ABS(PIVOT).GT.ULOC*MAX(RMAX,AMAX)) THEN
          PIVSIZ = 1
          A(POSPV1) = PIVOT
C 1 x 1 pivot is chosen
          GO TO 810
        END IF

C If there is only one remaining fully summed row and column exit.
        IF (NPIV+1.EQ.NASS) THEN
          A(POSPV1) = PIVOT
          GO TO 920
        END IF

C Jump if 2 x 2 candidate is diagonal
        IF (AMAX.LE.TOL) GO TO 920

C Check block pivot of order 2 for stability.
C Find largest entry in row IPIV outwith the pivot.
        IF (RMAX.LT.AMAX) THEN
          RMAX = ZERO
C Split loops in two because only upper triangle is held
C Scan lower triangle by scanning up column of upper triangle.
          DO 113 K = 1, IPIV - NPIV - 1
          IF (IPIV-K.EQ.JMAX) GO TO 113
          RMAX=MAX(RMAX,ABS(A(POSPV1-K*NFRONT)))
  113       CONTINUE
C Scan upper triangle by scanning along row from first off-diagonal
          DO 114 K =  1, NFRONT - IPIV
          IF (IPIV+K.EQ.JMAX) GO TO 114
          RMAX = MAX(RMAX,ABS(A(POSPV1+K)))
  114       CONTINUE
        ENDIF

C APOSJ is set to the beginning of row JMAX in working front.
        APOSJ = APOS + (JMAX-1)*NFRONT
C POSPV2 is the position in A of the second diagonal of a 2 x 2 pivot.
C OFFDAG is the position in A of the off-diagonal of a 2 x 2 pivot.
        POSPV2 = APOSJ + JMAX - 1
        IF (IPIV.GT.JMAX) THEN
          OFFDAG = APOSJ + IPIV - 1
        ELSE
          OFFDAG = APOSI + JMAX - 1
        END IF

C     Find largest entry in row JMAX outwith the pivot.
        TMAX = ZERO
C Split loops in two because only upper triangle is held
C Scan lower triangle by scanning up column of upper triangle.
        DO 115 K = 1, JMAX - NPIV - 1
          IF (JMAX-K.EQ.IPIV) GO TO 115
          TMAX=MAX(TMAX,ABS(A(POSPV2-K*NFRONT)))
  115     CONTINUE
C Scan upper triangle by scanning along row from first off-diagonal
        DO 116 K =  1, NFRONT - JMAX
          IF (JMAX+K.EQ.IPIV) GO TO 116
          TMAX = MAX(TMAX,ABS(A(POSPV2+K)))
  116     CONTINUE


C DETPIV is the value of the determinant of the 2x2 pivot.
        IF (KIND.EQ.1) THEN
          DETPIV = A(POSPV1)*A(POSPV2) - CMAX*CMAX
        ELSE
          DETPIV = A(POSPV1)*A(POSPV2) - CMAX*CONJG(CMAX)
C Make DETPIV real in Hermitian case
          DETPIV = REAL(DETPIV)
        ENDIF

        MAXPIV = MAX(ABS(A(POSPV1)),ABS(A(POSPV2)))
        IF (MAXPIV.EQ.ZERO) MAXPIV = ONE
        IF (ABS(DETPIV)/MAXPIV.LE.TOL) GO TO 920
        PIVSIZ = 2
C Check pivot for stability
C Jump if pivot fails test
C This is componentwise test
        IF ((ABS(A(POSPV2))*RMAX+AMAX*TMAX)*ULOC.GT.
   +      ABS(DETPIV)) GO TO 920
        IF ((ABS(A(POSPV1))*TMAX+AMAX*RMAX)*ULOC.GT.
   +      ABS(DETPIV)) GO TO 920

C 2 x 2 pivot is chosen

C
C       Pivot has been chosen. It has order PIVSIZ.
  810     LPIV = IPIV
        IF (PIVSIZ.EQ.2) LPIV = MIN(IPIV,JMAX)
C Change made at Stephane's suggestion
CCC       KR = MAX(KR,NPIV+PIVSIZ)
        KR = MAX(KR,MPIV+PIVSIZ)
        KCT = SIZBLK - MPIV - PIVSIZ + 1

C The following loop moves the pivot block to the top left
C       hand corner of the uneliminated frontal matrix.
        DO 860 KROW = NPIV,NPIV + PIVSIZ - 1
C We jump if swop is not necessary.
          IF (LPIV.EQ.KROW+1) GO TO 850

C Swop first part of rows (going down columns)
C JA1 is used as running index for row LPIV
          JA1 = APOS + (LPIV-1)
C J1 is used as running index for row KROW+1
          J1 = APOS + KROW
          DO 820 JJ = 1,KROW
          SWOP = A(JA1)
          A(JA1) = A(J1)
          A(J1) = SWOP
          JA1 = JA1 + NFRONT
          J1 = J1 + NFRONT
  820       CONTINUE
C Swop off-diagonal (in position KROW+1, IPIV)
C Only needed in Hermitian case
          IF (KIND.NE.1)  A(JA1) = CONJG(A(JA1))
C Swop middle part of rows (KROW+1 by rows, LPIV by columns)
          JA1 = JA1 + NFRONT
          J1 = J1 + 1
          DO 830 JJ = 1,LPIV - KROW - 2
          SWOP = A(JA1)
          IF (KIND.EQ.1) THEN
            A(JA1) = A(J1)
            A(J1) = SWOP
          ELSE
            A(JA1) = CONJG(A(J1))
            A(J1) = CONJG(SWOP)
          ENDIF
          JA1 = JA1 + NFRONT
          J1 = J1 + 1
  830       CONTINUE
C Swop diagonals
          SWOP = A(APOS+KROW* (NFRONT+1))
          A(APOS+KROW* (NFRONT+1)) = A(JA1)
          A(JA1) = SWOP
C Swop last part of rows
          DO 840 JJ = 1,NFRONT - LPIV
          JA1 = JA1 + 1
          J1 = J1 + 1
          SWOP = A(JA1)
          A(JA1) = A(J1)
          A(J1) = SWOP
  840       CONTINUE
C Swop integer indexing information
          IPOS = IWPOS + KROW
          IEXCH = IWPOS + LPIV - 1
          ISWOP = IW(IPOS)
          IW(IPOS) = IW(IEXCH)
          IW(IEXCH) = ISWOP
C Set LPIV for the swop of the second row of block pivot.
  850       LPIV = MAX(IPIV,JMAX)
  860     CONTINUE

C
C Set POSPV1 and POSPV2 to new position of pivots.
  811     POSPV1 = APOS + NPIV* (NFRONT+1)
        POSPV2 = POSPV1 + NFRONT + 1

        IF (PIVSIZ.EQ.1) THEN
C Perform the elimination using entry A(POSPV1) as pivot.
C We store U and D inverse.
C Later we store D inverse U which is passed to the solution entry.
          FLOPSB = FLOPSB + ONE
          IF (KIND.EQ.1) A(POSPV1) = ONE/A(POSPV1)
C Keep pivots real in Hermitian case
          IF (KIND.NE.1) A(POSPV1) = REAL(ONE/A(POSPV1))
          IF (REAL(A(POSPV1)).LT.ZERO) NEIG = NEIG + 1
          J1 = POSPV1 + 1
          J2 = POSPV1 + NASS - (NPIV+1)
          IBEG = POSPV1 + NFRONT + 1
          IEND = APOS + (NPIV+1)*NFRONT + NFRONT - 1
          DO 880 JJ = J1,J2
C AMULT1 is used to hold the multiplier
          IF (KIND.EQ.1) THEN
            AMULT1 = -A(JJ)*A(POSPV1)
          ELSE
            AMULT1 = -CONJG(A(JJ))*A(POSPV1)
          ENDIF
C Hold original entry for GEMM multiply
          IF (.NOT.LASTBK) THEN
            IF (KIND.EQ.1) THEN
            A(POSPV1+(JJ-J1+1)*NFRONT) = A(JJ)
            ELSE
            A(POSPV1+(JJ-J1+1)*NFRONT) = CONJG(A(JJ))
            ENDIF
          ENDIF
          JCOL = JJ
          FLOPSB = FLOPSB + (IEND-IBEG+1)*2 + 1
          IF (MPIV+JJ-J1+2.GT.PIVBLK) GO TO 871
C          The following special comment forces vectorization on
C             Crays.
CDIR$        IVDEP
          DO 870 IROW = IBEG,IEND
            A(IROW) = A(IROW) + AMULT1*A(JCOL)
            JCOL = JCOL + 1
  870       CONTINUE

CCC
C In Hermitian case make diagonal in position IBEG real
          IF (KIND.NE.1) A(IBEG) = REAL(A(IBEG))

  871       IF (KIND.EQ.1) THEN
            A(JJ) = AMULT1
          ELSE
            A(JJ) = CONJG(AMULT1)
          ENDIF
          IBEG = IBEG + NFRONT + 1
          IEND = IEND + NFRONT
  880       CONTINUE

          NPIV = NPIV + 1
          MPIV = MPIV + 1
          NTOTPV = NTOTPV + 1

C     IF (LCASE) THEN
C       write(7,'(/A,I8/A,7I8)') '*** Frontal matrix at step',IASS,
C    *   'NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ',
C    *    NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ
C       write(7,'(/A/(10I8))') 'IW array',
C    *    (IW(IWPOS+I-1),I=1,NFRONT)
C       write(7,'(/A)') 'Frontal matrix'
C       DO 1123 J = 1, NFRONT
C       write(7,'(/(6D13.5))')
C    *    (A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT)
C1123   CONTINUE
C     ENDIF

          IF (MPIV.EQ.SIZBLK) GO TO 930

        ELSE
C Perform elimination using block pivot of order two.
C Replace block pivot by its inverse.
          OFFDAG = POSPV1 + 1
          FLOPSB = FLOPSB + 6.0

          SWOP = A(POSPV2)
          IF (REAL(DETPIV).LT.ZERO) THEN
          NEIG = NEIG + 1
          ELSE
          IF (REAL(SWOP).LT.ZERO) NEIG = NEIG + 2
          END IF

          IF (KIND.EQ.1) THEN
          A(POSPV2) = A(POSPV1)/DETPIV
          A(POSPV1) = SWOP/DETPIV
          A(OFFDAG) = -A(OFFDAG)/DETPIV
          ELSE
          A(POSPV2) = REAL(A(POSPV1)/DETPIV)
          A(POSPV1) = REAL(SWOP/DETPIV)
          A(OFFDAG) = -A(OFFDAG)/DETPIV
          ENDIF

          J1 = POSPV1 + 2
          J2 = POSPV1 + NASS - (NPIV+1)
CCC NB .. not the following
C         J2 = POSPV1 + NFRONT - (NPIV+1)
          IBEG = POSPV2 + NFRONT + 1
          IEND = APOS + (NPIV+2)*NFRONT + NFRONT - 1
          DO 900 JJ = J1,J2
          K1 = JJ
          K2 = JJ + NFRONT
          IF (KIND.EQ.1) THEN
            AMULT1 = - (A(POSPV1)*A(K1)+A(POSPV1+1)*A(K2))
            AMULT2 = - (A(POSPV1+1)*A(K1)+A(POSPV2)*A(K2))
          ELSE
            AMULT1 = - (A(POSPV1)*CONJG(A(K1))
   *                + CONJG(A(POSPV1+1))*CONJG(A(K2)))
            AMULT2 = - (A(POSPV1+1)*CONJG(A(K1))
   *               + A(POSPV2)*CONJG(A(K2)))
          ENDIF
C Save block in lower triangular part of front
          IF (.NOT.LASTBK) THEN
            IF (KIND.EQ.1) THEN
            A(POSPV1 + (JJ-J1+2)*NFRONT) = A(K1)
            A(POSPV1 + (JJ-J1+2)*NFRONT + 1) = A(K2)
            ELSE
            A(POSPV1 + (JJ-J1+2)*NFRONT) = CONJG(A(K1))
            A(POSPV1 + (JJ-J1+2)*NFRONT + 1) = CONJG(A(K2))
            ENDIF
          ENDIF
          FLOPSB = FLOPSB + (IEND-IBEG+1)*4 + 6
          IF (MPIV+JJ-J1+3.GT.PIVBLK) GO TO 891

C          The following special comment forces vectorization on
C            Crays.
CDIR$        IVDEP
          DO 890 IROW = IBEG,IEND
            A(IROW) = A(IROW) + AMULT1*A(K1) + AMULT2*A(K2)
            K1 = K1 + 1
            K2 = K2 + 1
  890       CONTINUE
  891       IF (KIND.EQ.1) THEN
            A(JJ) = AMULT1
            A(JJ+NFRONT) = AMULT2
          ELSE
            A(JJ) = CONJG(AMULT1)
            A(JJ+NFRONT) = CONJG(AMULT2)
          ENDIF
          IBEG = IBEG + NFRONT + 1
          IEND = IEND + NFRONT
  900       CONTINUE
C Flag column indices of 2 x 2 pivot.
          IPOS = IWPOS + NPIV
          IW(IPOS) = -IW(IPOS)
          IW(IPOS+1) = -IW(IPOS+1)
          NPIV = NPIV + 2
          MPIV = MPIV + 2
          NTOTPV = NTOTPV + 2
          NTWO = NTWO + 1

C      IF (LCASE) THEN
C       write(7,'(/A,I8/A/9I8)') '*** Frontal matrix at step',IASS,
C    *   'NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ,MPIV,SIZBLK',
C    *    NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ,MPIV,SIZBLK
C       write(7,'(/A/(10I8))') 'IW array',
C    *    (IW(IWPOS+I-1),I=1,NFRONT)
C       write(7,'(/A)') 'Frontal matrix'
C       DO 1124 J = 1, NFRONT
C       write(7,'(/(6D13.5))')
C    *    (A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT)
C1124   CONTINUE
C      ENDIF

          IF (MPIV.EQ.SIZBLK) GO TO 930
        END IF

      GO TO 920
C 920   CONTINUE


C Finished pivoting on block BLK ... no update rest of pivot block
C     using GEMM
 930    IF (LASTBK) THEN
      IF (NPIV.EQ.NASS) GO TO 935
      IF (.NOT. LSTAT)  GO TO 935
      ULOC = ULOC/10.0D0
      IF (ULOC.LT.UTARG) THEN
        ULOC = ULOC * 10.0D0
        GO TO 9919
      ENDIF
      KCT = SIZBLK + 1 - MPIV
      GO TO 920
      ENDIF
C Check if any pivots chosen from this block. If not, increase PIVBLK
C     and try again.
      IF (MPIV.EQ.0) THEN
      PIVBLK = 2*PIVBLK
      GO TO 918
      ENDIF

      KBLK = (NASS-(LASPIV+PIVBLK))/PIVBLK
CCC   IF (LCASE) KBLK = 0
      L = NASS - (LASPIV+PIVBLK)
      APOS4 = APOS+(LASPIV+PIVBLK)*(NFRONT+1)
      DO 931 KB = 1,KBLK
      FLOPSX = FLOPSX + PIVBLK*(PIVBLK-1)*MPIV
      IF (KIND.EQ.1) THEN
        CALL ZGEMM('N','N',L-(KB-1)*PIVBLK,PIVBLK,MPIV,CONE,
   +           A(APOSBK+PIVBLK*KB),NFRONT,
   +           A(APOSBK+PIVBLK*KB*NFRONT),NFRONT,CONE,
   +           A(APOS4+PIVBLK*(KB-1)*(NFRONT+1)),NFRONT)
      ELSE
        CALL ZGEMM('N','N',L-(KB-1)*PIVBLK,PIVBLK,MPIV,CONE,
   +           A(APOSBK+PIVBLK*KB),NFRONT,
   +           A(APOSBK+PIVBLK*KB*NFRONT),NFRONT,CONE,
   +           A(APOS4+PIVBLK*(KB-1)*(NFRONT+1)),NFRONT)
      ENDIF
C And now process the part of the pivot row outside the fs block
      IF (NFRONT.GT.NASS) THEN
        IF (KIND.EQ.1)
   +    CALL ZGEMM('N','T',NFRONT-NASS,PIVBLK,MPIV,CONE,
   +           A(APOSBK+NASS-LASPIV),NFRONT,
   +           A(APOSBK+PIVBLK*KB),NFRONT,CONE,
   +           A(APOSBK+KB*NFRONT*PIVBLK+NASS-LASPIV),NFRONT)
        IF (KIND.NE.1)
   +    CALL ZGEMM('N','C',NFRONT-NASS,PIVBLK,MPIV,CONE,
   +           A(APOSBK+NASS-LASPIV),NFRONT,
   +           A(APOSBK+PIVBLK*KB),NFRONT,CONE,
   +           A(APOSBK+KB*NFRONT*PIVBLK+NASS-LASPIV),NFRONT)
      ENDIF
  931   CONTINUE

     SIZC = NASS - (KBLK+1)*PIVBLK - LASPIV
     SIZF = NFRONT - (KBLK+1)*PIVBLK - LASPIV
     APOSA = APOSBK + (KBLK+1)*PIVBLK
     DO 934 K = 1,MPIV
C APOSB is in lower triangle and so conjugated already
       APOSB = APOSBK + NFRONT*PIVBLK*(KBLK+1) + K - 1
       APOSM = APOSBK + PIVBLK*(KBLK+1) + (K-1)*NFRONT
       APOSC = APOSBK + PIVBLK*(KBLK+1)*(NFRONT+1)
       DO 933 JJ = 1,SIZC
        IF (KIND.EQ.1) THEN
          DO 932 J = JJ,SIZC
          A(APOSC+J-1) = A(APOSC+J-1) + A(APOSA+J-1)*A(APOSB)
  932       CONTINUE
        ELSE
          DO 937 J = JJ,SIZC
          A(APOSC+J-1) = A(APOSC+J-1)
   *               + A(APOSA+J-1)*A(APOSB)
  937       CONTINUE
        ENDIF
C And now process the part of the pivot row outside the fs block
        IF (KIND.EQ.1) THEN
          DO 936 J = SIZC+1,SIZF
          A(APOSC+J-1) = A(APOSC+J-1) + A(APOSA+J-1)*A(APOSM)
  936       CONTINUE
        ELSE
          DO 938 J = SIZC+1,SIZF
          A(APOSC+J-1) = A(APOSC+J-1)
   *               + A(APOSA+J-1)*CONJG(A(APOSM))
  938       CONTINUE
        ENDIF
        APOSC = APOSC + NFRONT
        APOSB = APOSB + NFRONT
        APOSM = APOSM + 1
  933   CONTINUE
      APOSA = APOSA + NFRONT
  934   CONTINUE

      APOSBK = APOSBK + MPIV*(NFRONT+1)
      LASPIV = NPIV

C     IF (LCASE) THEN
C       write(7,'(/A,I8/A,7I8)') '*** Frontal matrix at step',IASS,
C    *   'NFRONT,NASS,APOS,IWPOS,NPIV',
C    *    NFRONT,NASS,APOS,IWPOS,NPIV
C       write(7,'(/A,2I8)') 'After blocking .. APOSBK,LASPIV',
C    *    APOSBK,LASPIV
C       write(7,'(/A/(10I8))') 'IW array',
C    *    (IW(IWPOS+I-1),I=1,NFRONT)
C       write(7,'(/A)') 'Frontal matrix'
C       DO 1125 J = 1, NFRONT
C       write(7,'(/(6D13.5))')
C    *    (A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT)
C1125   CONTINUE
C     ENDIF

  918   CONTINUE
C End of main elimination loop.

C SCHUR if set to .TRUE. then the Schur complement will be generated
C    using level 3 BLAS at this step so an extra copy of pivot row is
C    made.


C*************************
C   Do static pivoting *
C*************************
 9919    IPIV = LASPIV+MPIV
 9920    IPIV = IPIV + 1
CADD Probably not needed .. use only IPIV
C      IF (KR.GT.SIZBLK) KR = MPIV + 1
C   Note that IPIV is the position within the complete current front.
C      IPIV = LASPIV + KR
C   APOSI is set to the beginning of row IPIV in working front.
       APOSI = APOS + (IPIV-1)*NFRONT
C   Set position and value of potential pivot.
       POSPV1 = APOSI + IPIV - 1
       PIVOT = A(POSPV1)
CADD
C Although these are not needed just now they are kept for when
C we use 2 x 2 static pivots.
CCC      PIVSIZ = 1
CCC      LPIV = IPIV

C Logic has changed so no need to calculate AMAX
C This is code from earlier experiments

CCC      AMAX = ZERO
C Split loops in two because only upper triangle is held
C Scan lower triangle by scanning up column of upper triangle.
CCC      DO 9876 K = 1, IPIV - NPIV - 1
CCC      AMAX = MAX(AMAX,ABS(A(POSPV1-K*NFRONT)))
CCC76    CONTINUE
C Scan upper triangle by scanning along row from first off-diagonal
CCC      DO 9878 K =  1, NFRONT - IPIV
CCC      AMAX = MAX(AMAX,ABS(A(POSPV1+K)))
CCC78    CONTINUE

C   Check size of 1 x 1 pivot and adjust if necessary
C      IF (ABS(A(POSPV1)).LT.CNTL(4)) THEN
C          PIVOT = CNTL(4)
C      IF (ABS(A(POSPV1)).LT.MAX(ULOC*AMAX,STCTOL)) THEN
C          PIVOT = MAX(ULOC*AMAX,STCTOL)
       IF (ABS(A(POSPV1)).LT.STCTOL) THEN
         PIVOT = STCTOL
          IF (REAL(A(POSPV1)) .LT. ZERO) THEN
           A(POSPV1) = -PIVOT
           PIVOT   = -PIVOT
          ELSE
           A(POSPV1) = PIVOT
          ENDIF
          NBSTATIC = NBSTATIC + 1
       ENDIF
C Perform the elimination using entry A(POSPV1) as pivot.
C We store U and D inverse.
C Later we store D inverse U which is passed to the solution entry.
       FLOPSB = FLOPSB + ONE
       A(POSPV1) = ONE/A(POSPV1)
       IF (REAL(A(POSPV1)).LT.ZERO) NEIG = NEIG + 1

       J1 = POSPV1 + 1
       J2 = POSPV1 + NASS - (NPIV+1)
       IBEG = POSPV1 + NFRONT + 1
       IEND = APOSI + 2*NFRONT - 1
       DO 9880 JJ = J1,J2
C AMULT1 is used to hold the multiplier
          IF (KIND.EQ.1) THEN
          AMULT1 = -A(JJ)*A(POSPV1)
          ELSE
          AMULT1 = -CONJG(A(JJ))*A(POSPV1)
          ENDIF
C   Hold original entry for GEMM multiply
          JCOL = JJ
          FLOPSB = FLOPSB + (IEND-IBEG+1)*2 + 1
CADD Not necessary because control is on setting of J2
C         IF (MPIV+JJ-J1+2.GT.PIVBLK) GO TO 9871
C   The following special comment forces vectorization on
C   Crays.
CDIR$        IVDEP
          DO 9870 IROW = IBEG,IEND
           A(IROW) = A(IROW) + AMULT1*A(JCOL)
           JCOL = JCOL + 1
 9870       CONTINUE
          IF (KIND.EQ.1) THEN
          A(JJ) = AMULT1
          ELSE
          A(JJ) = CONJG(AMULT1)
          ENDIF
          IBEG = IBEG + NFRONT + 1
          IEND = IEND + NFRONT
 9880    CONTINUE
       NPIV = NPIV + 1
       MPIV = MPIV + 1
       NTOTPV = NTOTPV + 1

C    IF (LCASE) THEN
C   write(7,'(A)') 'Static pivoting step'
C   write(7,'(/A,I8/A/9I8)') '*** Frontal matrix at step',IASS,
C    *   'NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ,SIZBLK,MPIV',
C    *    NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ,SIZBLK,MPIV
C   write(7,'(/A/(10I8))') 'IW array',
C    *    (IW(IWPOS+I-1),I=1,NFRONT)
C       write(7,'(/A)') 'Frontal matrix'
C       DO 1126 J = 1, NFRONT
C       write(7,'(/(6D13.5))')
C    *    (A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT)
C1126   CONTINUE
C     ENDIF

       IF (MPIV.LT.SIZBLK) GO TO 9920
C*******************************
C End of static pivoting loop  *
C*******************************

  935   SCHUR = (NBLOC.LT.(NFRONT-NASS) .AND. NPIV.GE.NBLOC)
      IF (ICNTL(16).EQ.1) THEN
C Remove "zero" rows from fully summed rows within Schur complement
C ZCOUNT is count of "zero" rows
      ZCOUNT = 0
C APOS4 is beginning of block of fully summed uneliminated variables
      APOS4 = APOS + NPIV*NFRONT + NPIV

C Expand block so that lower triangle is included
C APOSB scans lower triangle by rows
C APOSC sweeps upper triangle by columns
      APOSB = APOS4 + NFRONT
      APOSC = APOS4 + 1
      DO 4444 I = 2,NASS-NPIV
        DO 4443 J = 1,I-1
          IF (KIND.EQ.1) THEN
          A(APOSB) = A(APOSC)
          ELSE
          A(APOSB) = CONJG(A(APOSC))
          ENDIF
          APOSB = APOSB + 1
          APOSC = APOSC + NFRONT
 4443     CONTINUE
        APOSB = APOS4 + NFRONT*I
        APOSC = APOS4 + I
 4444   CONTINUE

C Remove any zero rows by swopping with "first" row so that all zero
C    rows will be swept to beginning of block
C Row I is the row currently being checked
      I = NASS - NPIV
C Also exchange integer information accordingly
 4445   CONTINUE
      IF (ZCOUNT.EQ.I) GO TO 4450
C Check row I for zero
C APOSB is beginning of row I
      APOSB = APOS4 + (I-1)*NFRONT
      DO 4446 J = 1,NFRONT-NPIV
        IF (ABS(A(APOSB+J-1)).GT.TOL) GO TO 4449
 4446   CONTINUE
C Row is all zero
      ZCOUNT = ZCOUNT + 1
C Swop row ZCOUNT with row I
      IF (KIND.EQ.1) THEN
        DO 4447 J = 1,NFRONT-NPIV
          A(APOSB+J-1) = A(APOS4+NFRONT*(ZCOUNT-1)+J-1)
 4447     CONTINUE
      ELSE
C Have to split loop into three parts
        DO 4451 J = 1,I
          A(APOSB+J-1) = A(APOS4+NFRONT*(ZCOUNT-1)+J-1)
 4451     CONTINUE
C Deck does not test this
        DO 4452 J = I+1,ZCOUNT
          A(APOSB+J-1) = CONJG(A(APOS4+NFRONT*(ZCOUNT-1)+J-1))
 4452     CONTINUE
        DO 4453 J = ZCOUNT+1,NFRONT-NPIV
          A(APOSB+J-1) = A(APOS4+NFRONT*(ZCOUNT-1)+J-1)
 4453     CONTINUE
      ENDIF
C Zero row ZCOUNT
      DO 4448 J = 1,NFRONT-NPIV
        A(APOS4+NFRONT*(ZCOUNT-1)+J-1) = ZERO
 4448   CONTINUE
C Swop integers
      ISWOP = IW(IWPOS+NPIV+ZCOUNT-1)
      IW(IWPOS+NPIV+ZCOUNT-1) = IW(IWPOS+NPIV+I-1)
      IW(IWPOS+NPIV+I-1) = ISWOP
      GO TO 4445
 4449   I = I - 1
      GO TO 4445
 4450   CONTINUE
      ELSE
      ZCOUNT = 0
      ENDIF
C Set order of Schur complement (including rows of delayed pivots)
C But not including "zero" rows
      NSC1 = NFRONT - NPIV - ZCOUNT

C     IF (LCASE) THEN
C       write(7,'(/A,I8/A,5I8)') '*** Frontal matrix fter step',IASS,
C    *   'NFRONT',
C    *    NFRONT
C       write(7,'(/A/(10I8))') 'IW array',
C    *    (IW(IWPOS+I-1),I=1,NFRONT)
C       write(7,'(/A)') 'Frontal matrix'
C       DO 1127 J = 1, NFRONT
C       write(7,'(/(6D13.5))')
C    *    (A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT)
C1127   CONTINUE
C    ENDIF

C Accumulate total number of delayed pivots in INFO(23)
C Do not accumulate at last step where matrix is singular
      IF (IASS.NE.NSTEPS) INFO(23) = INFO(23) + NASS - NPIV
C SET LSTAT
      IF (CNTL(4).GT.ZERO .AND. INFO(23).GT.CNTL(5)*N) LSTAT = .TRUE.

C Jump if no Schur complement to form ... just store factors.
      IF (NSC1.EQ.0) GO TO 1830

C Save space for factors and Schur complement (if appropriate).

      IF (.NOT.SCHUR) THEN
C We now compute triangular Schur complement not using BLAS.
C This Schur complement is placed directly on the stack.

C Remove these
      RLSPA = MAX(RLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+
   +              NSTKAC(1))
      TRLSPA = MAX(TRLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+
   +              TOTSTA(1))
C ZCOUNT added because we do store this (for now)
      NSTKAC(1) = NSTKAC(1) + ((NSC1+ZCOUNT+1)*(NSC1+ZCOUNT))/2
      NSTACK(1) = NSTACK(1) + ((NSC1+ZCOUNT+1)*(NSC1+ZCOUNT))/2
      TOTSTA(1) = TOTSTA(1) + ((NSC1+ZCOUNT+1)*(NSC1+ZCOUNT))/2

C Initialize Schur complement
C Copying from the back to avoid overwriting important data
      APOSI = APOS + NFRONT*NFRONT - 1
      DO 1370 JJ = 1,NFRONT-NPIV
        J = APOSI
        DO 1360 JJJ = 1,JJ
          A(ASTK) = A(J)
          ASTK = ASTK - 1
          J = J - 1
 1360     CONTINUE
        APOSI = APOSI - NFRONT
 1370   CONTINUE
C APOS4 is the position in A of the first entry in the Schur complement.
      APOS4 = ASTK + 1


C Perform pivoting operations.
C Initialize variables
      J1 = IWPOS
      LTWO = .FALSE.
      POSPV1 = APOS
      DO 1450 I1 = 1,NPIV
        IF (LTWO) GO TO 1440
        APOSI = APOS + (I1-1)*NFRONT + NASS
        J2 = APOS + NFRONT* (I1-1) + NFRONT - 1
CCC What happens here ??
        APOSC = APOS4 +
   *     ((NASS-NPIV-ZCOUNT)*(2*NFRONT-NPIV-ZCOUNT-NASS+1))/2
C Check to see if current pivot is 1 x 1  or 2 x 2.
        IF (IW(J1).GT.0) THEN
          FLOPSB = FLOPSB + (NFRONT-NASS) +
   *                (NFRONT-NASS)* (NFRONT-NASS+1)
          DO 1410 JJ = APOSI,J2
          IF (KIND.EQ.1) THEN
            AMULT1 = -A(JJ)*A(POSPV1)
          ELSE
            AMULT1 = -CONJG(A(JJ))*A(POSPV1)
          ENDIF
          DO 1400 JJJ = JJ,J2
            A(APOSC) = A(APOSC) + AMULT1*A(JJJ)
            APOSC = APOSC + 1
 1400       CONTINUE
          IF (KIND.EQ.1) THEN
            A(JJ) = AMULT1
          ELSE
            A(JJ) = CONJG(AMULT1)
          ENDIF
 1410       CONTINUE
          J1 = J1 + 1
        ELSE
          POSPV2 = POSPV1 + NFRONT + 1
          OFFDAG = POSPV1 + 1
          FLOPSB = FLOPSB + 6* (NFRONT-NASS) +
   +           2* (NFRONT-NASS)* (NFRONT-NASS+1)
          DO 1430 JJ = APOSI,J2
          IF (KIND.EQ.1) THEN
            AMULT1 = - (A(POSPV1)*A(JJ)+A(OFFDAG)*A(JJ+NFRONT))
            AMULT2 = -A(POSPV2)*A(JJ+NFRONT) - A(OFFDAG)*A(JJ)
          ELSE
            AMULT1 = - CONJG(A(POSPV1)*A(JJ)
   *               + A(OFFDAG)*A(JJ+NFRONT))
            AMULT2 = - CONJG(A(POSPV2)*A(JJ+NFRONT)
   *               + CONJG(A(OFFDAG))*A(JJ))
          ENDIF
          DO 1420 JJJ = JJ,J2
            A(APOSC) = A(APOSC) + AMULT1*A(JJJ) +
   +               AMULT2*A(JJJ+NFRONT)
            APOSC = APOSC + 1
 1420       CONTINUE
          IF (KIND.EQ.1) THEN
            A(JJ) = AMULT1
            A(JJ+NFRONT) = AMULT2
          ELSE
            A(JJ) = CONJG(AMULT1)
            A(JJ+NFRONT) = CONJG(AMULT2)
          ENDIF
 1430       CONTINUE
          J1 = J1 + 2
          POSPV1 = POSPV2
          LTWO = .TRUE.
          GO TO 1450
        END IF

 1440     LTWO = .FALSE.
C Move to next pivot position
        POSPV1 = POSPV1 + NFRONT + 1
 1450   CONTINUE

C       IF (LCASE) THEN
C       write(7,'(A,I8)') 'GEMM not used at stage',IASS
C       write(7,'(A/(6D13.5))') 'Stacking Schur',
C    *       (A(I),I=APOS4,APOS4+(NSC1*(NSC1+1))/2-1)
C       ENDIF

      ELSE

C Action if SCHUR is true, We now use GEMM.

C Since SCHUR is true, copy U,
C   divide factors by D, and generate Schur complement using GEMM.
C   Then compress factors (to upper trapezoidal), and stack
C   half of the Schur complement.

C APOS4 is position in A of first entry of Schur complement to be
C   updated using GEMM.
      APOS4 = APOS+NASS*(NFRONT+1)

C APOS3 is the position in A of the first entry in the copy of U.
      APOS3 = APOS+NASS*NFRONT

C Copy U and divide factors by D
C Conjugating the factors would be good for SOLVE !!!
C Initialize variables
      J1 = IWPOS
      LTWO = .FALSE.
      POSPV1 = APOS
      DO 1490 I = 1,NPIV
        IF (LTWO) GO TO 1480
        APOSI = APOS + (I-1)*NFRONT + NASS
        POSELT = APOS3 + I - 1
C Check to see if current pivot is 1 x 1  or 2 x 2.
        IF (IW(J1).GT.0) THEN
          FLOPSB = FLOPSB + (NFRONT-NASS)
          DO 1460 JJ = APOSI,APOS + NFRONT*I - 1
          IF (KIND.EQ.1) THEN
            A(POSELT) = A(JJ)
          ELSE
            A(POSELT) = CONJG(A(JJ))
          ENDIF
          A(JJ) = -A(JJ)*A(POSPV1)
          POSELT = POSELT + NFRONT
 1460       CONTINUE
          J1 = J1 + 1
        ELSE
          POSPV2 = POSPV1 + NFRONT + 1
          OFFDAG = POSPV1 + 1
          FLOPSB = FLOPSB + 6* (NFRONT-NASS)
          DO 1470 JJ = APOSI,APOS + NFRONT*I - 1
          IF (KIND.EQ.1) THEN
            A(POSELT) = A(JJ)
C           A(POSELT+NFRONT-NASS) = A(JJ+NFRONT)
            A(POSELT+1) = A(JJ+NFRONT)
          ELSE
            A(POSELT) = CONJG(A(JJ))
            A(POSELT+1) = CONJG(A(JJ+NFRONT))
          ENDIF
          A(JJ) = - (A(POSPV1)*A(JJ)+A(OFFDAG)*A(JJ+NFRONT))
          IF (KIND.EQ.1) THEN
          A(JJ+NFRONT) = -A(POSPV2)*A(JJ+NFRONT) -
   +               A(OFFDAG)*A(POSELT)
          ELSE
          A(JJ+NFRONT) = -A(POSPV2)*A(JJ+NFRONT) -
   +               CONJG(A(OFFDAG)*A(POSELT))
          ENDIF
          POSELT = POSELT + NFRONT
 1470       CONTINUE
          J1 = J1 + 2
          POSPV1 = POSPV2
          LTWO = .TRUE.
          GO TO 1490
        END IF

 1480     LTWO = .FALSE.
C Move to next pivot position
        POSPV1 = POSPV1 + NFRONT + 1
 1490   CONTINUE

C Now create Schur complement by using Level 3 BLAS GEMM.
C Jump if Schur complement is null.
C Increment FLOPSB
      FLOPSB = FLOPSB + NPIV* (NFRONT-NASS)**2 +
   *              NPIV* (NFRONT-NASS)
C We divide the multiply into blocks to avoid too many extra
C    computations when using GEMM with a symmetric result.
C Block formed by GEMM has NBLOC rows.
      KBLK = ( NFRONT-NASS)/NBLOC
CCC     IF (LCASE) KBLK = 0
      L =  NFRONT - NASS
      DO 1500 KB = 1,KBLK
C Accumulate extra flops caused by using GEMM
        FLOPSX = FLOPSX + NBLOC* (NBLOC-1)* (NPIV)
        IF (KIND.EQ.1) THEN
          CALL ZGEMM('N','N',L-(KB-1)*NBLOC,NBLOC,NPIV,CONE,
   +             A(APOS+NASS+NBLOC*(KB-1)),NFRONT,
   +             A(APOS3+NBLOC*(KB-1)*NFRONT),NFRONT,CONE,
   +             A(APOS4+NBLOC*(NFRONT+1)*(KB-1)),NFRONT)
        ELSE
          CALL ZGEMM('N','N',L-(KB-1)*NBLOC,NBLOC,NPIV,CONE,
   +             A(APOS+NASS+NBLOC*(KB-1)),NFRONT,
   +             A(APOS3+NBLOC*(KB-1)*NFRONT),NFRONT,CONE,
   +             A(APOS4+NBLOC*(NFRONT+1)*(KB-1)),NFRONT)
        ENDIF
 1500   CONTINUE

C Calculate the block upper triangular part of the Schur complement.
      DO 1550 I = 1 + KBLK*NBLOC,L
C APOSA holds the index of the current entry of the matrix A used to
C   form the Schur complement as C = A*B
C APOSB holds the index of the current entry of the matrix B used to
C   form the Schur complement as C = A*B
C APOSC holds the index of the current entry of the matrix C used to
C   form the Schur complement as C = A*B
        APOSA = APOS + NASS
        APOSB = APOS3 +(I-1)*NFRONT
        APOSC = APOS4 + (I-1)*NFRONT - 1
        DO 1540 K = 1,NPIV
          DO 1530 J = I,L
          A(APOSC+J) = A(APOSC+J) + A(APOSA+J-1)*A(APOSB)
 1530       CONTINUE
          APOSA = APOSA + NFRONT
          APOSB = APOSB + 1
 1540     CONTINUE
 1550   CONTINUE


C Stack half of Schur complement.

C Stack reals
C Stack in reverse order to avoid need for compresses.
      JA1 = APOS+NFRONT*NFRONT-1
      NSTKAC(1) = NSTKAC(1) + ((NSC1+1)* (NSC1))/2
      NSTACK(1) = NSTACK(1) + ((NSC1+1)* (NSC1))/2
      TOTSTA(1) = TOTSTA(1) + ((NSC1+1)* (NSC1))/2
C Stack by rows
      DO 1710 I = NSC1,1,-1
        DO 1700 JJ = JA1,JA1-(NSC1-I),-1
          A(ASTK) = A(JJ)
          ASTK = ASTK - 1
 1700     CONTINUE
        JA1 = JA1 - NFRONT
 1710   CONTINUE

C       IF (LCASE) THEN
C       write(7,'(A,I8)') 'GEMM used at stage',IASS
C       write(7,'(A/(6D13.5))') 'Stacking Schur',
C    *       (A(I),I=ASTK+1,ASTK+(NSC1*(NSC1+1))/2)
C       ENDIF

C END of SCHUR being true action (started after label 1450)
      END IF

C Stack integers
      NSTKAC(2) = NSTKAC(2) + NSC1 + 1
      NSTACK(2) = NSTACK(2) + NSC1 + 1
      TOTSTA(2) = TOTSTA(2) + NSC1 + 1
C Record space needed to this point
 1830   IF (IASS.EQ.NSTEPS) THEN
      INTSPA = MAX(INTSPA,IWPOS+NFRONT-1+NSTKAC(2))
      TINSPA = MAX(TINSPA,IWPOS+NFRONT-1+TOTSTA(2))
      GO TO 2158
      ELSE
      INTSPA = MAX(INTSPA,IWPOS+NFRONT-1+(N-NTOTPV+2)+NSTKAC(2))
      TINSPA = MAX(TINSPA,IWPOS+NFRONT-1+(N-NTOTPV+2)+TOTSTA(2))
      ENDIF

C Check space and compress if necessary
  444   NST = 0
C +1 for length of stacked entry
      IF (NSC1.GT.0) NST = NSC1 + 1
      IF (IWPOS+NFRONT-1+(N-NTOTPV+2)+NST.GT.ISTK) THEN
C Compress integer storage
      CALL ME57PD(A,IW,ISTK,IINPUT,PTRIRN,.FALSE.)
      NCMPBI = NCMPBI + 1

      IF (IWPOS+NFRONT-1+(N-NTOTPV+2)+NST.GT.ISTK) THEN
C Still insufficient space after compress
        IF (ICNTL(8).NE.0) THEN
          HOLD(1) = 2
          HOLD(2) = NBLK
          HOLD(3) = NTWO
          HOLD(4) = INFO(23)
          HOLD(5) = NCMPBI
          HOLD(6) = NEIG
          HOLD(7) = MAXFRT
          HOLD(8) = IWPOS
          HOLD(9) = APOS
          HOLD(10) = APOSBB
          HOLD(11) = NSTKAC(1)
          HOLD(12) = NSTKAC(2)
          HOLD(13) = AINPUT
          HOLD(14) = IINPUT
          HOLD(15) = ISTK
          HOLD(16) = ASTK
          HOLD(17) = INTSPA
          HOLD(18) = RLSPA
          HOLD(19) = PTRIRN
          HOLD(20) = PTRA
          HOLD(21) = NTOTPV
          HOLD(22) = NPOTPV
          HOLD(23) = NUMORG
          HOLD(24) = NFRONT
          HOLD(25) = NASS
C         HOLD(26) = NCOL
          HOLD(27) = NPIV
          HOLD(28) = IASS
          HOLD(29) = TINSPA
          HOLD(30) = TRLSPA
          HOLD(31) = TOTSTA(1)
          HOLD(32) = TOTSTA(2)
          HOLD(33) = NSTACK(1)
          HOLD(34) = NSTACK(2)
          IF (ICNTL(7).GT.1) HOLD(35) = ISNPIV
          IF (ICNTL(7).EQ.4) HOLD(36) = PHASE
          HOLD(37) = INFO(32)
          HOLD(38) = INFO(33)
          HOLD(39) = INFO(34)
          NSC1    = NFRONT-NPIV
          RINFO(3) =FLOPSA
          RINFO(4) =FLOPSB
          RINFO(5) =FLOPSX
          INFO(1) = 11
          HOLD(40) = NBSTATIC
          INFO(35) = HOLD(40)
        ELSE
          INFO(1)  = -4
          INFO(2)  = LIW
          INFO(18) = INTSPA
        ENDIF
        RETURN
      END IF
      END IF

      IF (NSC1.GT.0) THEN
      DO 1720 I = 1,NSC1
        IW(ISTK) = IW(IWPOS+NFRONT-I)
        ISTK = ISTK - 1
 1720   CONTINUE
      IW(ISTK) = NSC1
      ISTK = ISTK - 1
      ENDIF

C    IF (LCASE) THEN
C       write(7,'(/A,I8/A,5I8)') '*** Frontal matrix end step',IASS,
C    *   'NFRONT',
C    *    NFRONT
C       write(7,'(/A/(10I8))') 'IW array',
C    *    (IW(IWPOS+I-1),I=1,NFRONT)
C       write(7,'(/A)') 'Frontal matrix'
C       DO 1128 J = 1, NFRONT
C       write(7,'(/(6D13.5))')
C    *    (A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT)
C1128   CONTINUE
C    ENDIF

C Reset PPOS.
      DO 1840 JJ = IWPOS + NPIV,IWPOS + NFRONT - 1
      J = ABS(IW(JJ))
      PPOS(J) = N + 1
 1840   CONTINUE


C********************************
C    STORE FACTORS
C********************************
C Complete the integer information in the factors
 2158   IF (NPIV.EQ.0) GO TO 2159
      NBLK = NBLK + 1

      IW(IWPOS-2) = NFRONT
      IW(IWPOS-1) = NPIV
      IWPOS = IWPOS + NFRONT + 2

C Store information on the reals for the factors.
C We copy from A(JA1) to A(APOS2) ... the use of buffer space from
C   APOSBB to APOS ensures no overwrites.
      IF (INFO(1).EQ.-3) THEN
      INFO(40) = INFO(40) + (NPIV * (2*NFRONT-NPIV+1))/2
      GO TO 2159
      END IF

      APOS2 = APOSBB
C Store reals from full pivot.
      DO 2130 I = 1,NPIV
C JA1 points to the diagonal
      JA1 = APOS + (I-1)* (NFRONT+1)
      DO 2120 J = I,NPIV
C Store conjugate if matrix Hermitian
        IF (KIND.EQ.1) THEN
          A(APOS2) = A(JA1)
        ELSE
          A(APOS2) = CONJG(A(JA1))
        ENDIF
        IF (A(APOS2).EQ.ZERO) INFO(32) = INFO(32) + 1
        APOS2 = APOS2 + 1
        JA1 = JA1 + 1
 2120   CONTINUE
 2130   CONTINUE
      RPOS = APOS2
C Store rectangle
      DO 2150 I = 1,NPIV
      JA1 = APOS + (I-1)*NFRONT + NPIV
      DO 2140 J = 1,NFRONT - NPIV
        IF (KIND.EQ.1) THEN
          A(APOS2) = A(JA1)
        ELSE
          A(APOS2) = CONJG(A(JA1))
        ENDIF
        APOS2 = APOS2 + 1
        JA1 = JA1 + 1
 2140   CONTINUE
 2150   CONTINUE
C Set APOSBB for next block of factors
      APOSBB = APOS2

C Check rectangle for zeros
      DO 2152 J = 1,NFRONT-NPIV
      APOS2 = RPOS+J-1
      ZCOL = 1
      DO 2151 I = 1,NPIV
        IF (A(APOS2).EQ.ZERO) INFO(33) = INFO(33)+1
        IF (A(APOS2).NE.ZERO) ZCOL = 0
        APOS2 = APOS2 + NFRONT - NPIV
 2151   CONTINUE
      IF (ZCOL.EQ.1) INFO(34) = INFO(34)+1
 2152   CONTINUE

 2159   IASS = IASS + 1
    IF (IASS.LE.NSTEPS) THEN
C Initialize to zero to avoid problem when calling ME57ED
      IW(IWPOS-2) = 0
      IW(IWPOS-1) = 0
      GO TO 2160
    ENDIF
C2160 CONTINUE
C
C End of loop on tree nodes.
C

    INFO(35) = NBSTATIC
    IF (INFO(1).EQ.-3) THEN
      INFO(2)  = LA
      INFO(17) = MAX(INFO(17),RLSPA)
CCC Should use ALENB in BD for this?
      IF (ICNTL(7).EQ.4) INFO(17) = MAX(INFO(17),RLSPA + N)
      RETURN
    END IF
    GO TO 1000
 970  INFO(1) = -5
    INFO(2) = NTOTPV + 1
    IF (LDIAG.GT.0 .AND. LP.GE.0)
   *    WRITE(LP,99992) INFO(1),PIVOT,CNTL(2),INFO(2),ICNTL(7)
99992 FORMAT (/'*** Error message from routine ME57BD **',
   *     '   INFO(1) = ',I3/'Pivot has value ',2D16.8/'when ',
   *     'CNTL(2) has value ',D16.8/
   *     'at stage',I11,2X,'when ICNTL(7) =',I3)
    RETURN
 980  INFO(1) = -6
    INFO(2) = NTOTPV + 1
    IF (LDIAG.GT.0 .AND. LP.GE.0)
   *    WRITE(LP,99993) INFO(1),INFO(2),ICNTL(7)
99993 FORMAT (/'*** Error message from routine ME57BD **',
   *     '   INFO(1) = ',I3/'Change in sign of pivot at stage',
   *     I10,2X,'when ICNTL(7) = ',I3)
    RETURN
 1000 NRLBDU = APOSBB - 1
    NIRBDU = IWPOS - 3
    IF (NTOTPV.NE.N) THEN
      INFO(1) = 4
      IF (LDIAG.GT.0 .AND. WP.GE.0)
   *    WRITE(WP,99994) INFO(1),NTOTPV
99994 FORMAT (/'*** Warning message from routine ME57BD **',
   *       '   INFO(1) =',I2/5X, 'Matrix is singular, rank =', I5)

    ENDIF
C Recent change was to remove condition that ICNTL(16) was equal to 1
C More recent change to remove deficiency test.  This change means that
C we now test there is sufficient room to move the off-diagonal entries
C of the two by two pivots.
C 555 IF (NTOTPV.NE.N) THEN
C Check space (by this time there is nothing to compress)
  555 NRLBDU = APOSBB - 1
    NIRBDU = IWPOS - 3
    IF (NIRBDU+3*(N-NTOTPV) .GT. LIW
   +    .OR. NRLBDU+(N-NTOTPV)+NTWO .GT. LA) THEN
C I don't think this can happen ... at least I can't make it happen
C It is left in for "safety" :-)
C Still insufficient space after compress
      IF (ICNTL(8).NE.0) THEN
        HOLD(1) = 3
        HOLD(2) = NBLK
        HOLD(3) = NTWO
        HOLD(4) = INFO(23)
        HOLD(5) = NCMPBI
        HOLD(6) = NEIG
        HOLD(7) = MAXFRT
        HOLD(8) = IWPOS
        HOLD(9) = APOS
        HOLD(10) = APOSBB
        HOLD(11) = NSTKAC(1)
        HOLD(12) = NSTKAC(2)
        HOLD(13) = AINPUT
        HOLD(14) = IINPUT
        HOLD(15) = ISTK
        HOLD(16) = ASTK
        HOLD(17) = INTSPA
        HOLD(18) = RLSPA
        HOLD(19) = PTRIRN
        HOLD(20) = PTRA
        HOLD(21) = NTOTPV
        HOLD(22) = NPOTPV
        HOLD(23) = NUMORG
        HOLD(24) = NFRONT
        HOLD(25) = NASS
C         HOLD(26) = NCOL
        HOLD(27) = NPIV
        HOLD(28) = IASS
        HOLD(29) = TINSPA
        HOLD(30) = TRLSPA
        HOLD(31) = TOTSTA(1)
        HOLD(32) = TOTSTA(2)
        HOLD(33) = NSTACK(1)
        HOLD(34) = NSTACK(2)
        IF (ICNTL(7).GT.1) HOLD(35) = ISNPIV
        IF (ICNTL(7).EQ.4) HOLD(36) = PHASE
        HOLD(37) = INFO(32)
        HOLD(38) = INFO(33)
        HOLD(39) = INFO(34)
        NSC1    = NFRONT-NPIV
        RINFO(3) =FLOPSA
        RINFO(4) =FLOPSB
        RINFO(5) =FLOPSX
        IF (NRLBDU+(N-NTOTPV)+NTWO .GT. LA) INFO(1) = 10
        IF (NIRBDU+3*(N-NTOTPV) .GT. LIW)  INFO(1) = 11
        HOLD(40) = NBSTATIC
        INFO(35) = HOLD(40)
      ELSE
        IF (NIRBDU+3*(N-NTOTPV) .GT. LIW) THEN
          INFO(1)  = -4
          INFO(2)  = LIW
          INFO(18) = MAX(INTSPA,NIRBDU+3*(N-NTOTPV))
        ELSE
          INFO(1)  = -3
          INFO(2) = LA
          INFO(17) = MAX(INFO(17),RLSPA,NRLBDU+(N-NTOTPV)+NTWO)
          IF (ICNTL(7).EQ.4) INFO(17) =
   +      MAX(INFO(17),RLSPA + N,NRLBDU+(N-NTOTPV)+NTWO)
        ENDIF
      ENDIF
      RETURN
      ENDIF
C   ENDIF

C Add explicit entries in factors for zero pivots (now set to 1.0)
C Initialize flag array to identify indices of zero pivots
    IF (NTOTPV.NE.N) THEN
      DO 3331 I = 1,N
      PPOS(I) = 0
 3331   CONTINUE
      IWPOS = 4
      DO 3332 I = 1,NBLK
      NFRONT = IW(IWPOS)
      NPIV = IW(IWPOS+1)
      DO 3330 J = IWPOS+2,IWPOS+NPIV+1
        PPOS(ABS(IW(J))) = 1
 3330   CONTINUE
      IWPOS = IWPOS + NFRONT + 2
 3332   CONTINUE
      K= 0
      DO 3333 I=1,N
      IF (PPOS(I).EQ.0) THEN
        K=K+1
        NBLK = NBLK + 1
        NRLBDU = NRLBDU+1
        A(NRLBDU) = ONE
        IW(NIRBDU+1) = 1
        IW(NIRBDU+2) = 1
        IW(NIRBDU+3) = I
        NIRBDU = NIRBDU+3
      ENDIF
 3333   CONTINUE
    ENDIF

C
    INFO(14) = NRLBDU
C   Move the off-diagonal entries of the 2x2 pivots within full blocks
C   to the end of A.
    IW(1) = NRLBDU + 1
    IW(2) = NRLBDU + NTWO
    INFO(15) = IW(2)
    IW(3) = NBLK
    INFO(31) = NBLK
C   Negate the entries of L, move the off-diagonal entries of the 2x2
C   pivots within full blocks to the end of A and update NRLBDU to
C   correspond
    CALL ME57WD(KIND,A,LA,IW,LIW,NRLBDU)
    INFO(16) = NIRBDU
    INFO(18) = INTSPA
    INFO(20) = TINSPA
    INFO(17) = RLSPA
    INFO(19) = TRLSPA
    INFO(21) = MAXFRT
    INFO(22) = NTWO
C   INFO(23)  .. computed as sum of NASS-NPIV
    INFO(24) = NEIG
    INFO(25) = NTOTPV
    INFO(28) = NCMPBR
    INFO(29) = NCMPBI
    RINFO(3) = FLOPSA
    RINFO(4) = FLOPSB
    RINFO(5) = FLOPSX
    IF (INFO(27).GT.0) THEN
      RINFO(14) = ZERO
      DO 332 I = 1,N
      RINFO(14) = MAX(RINFO(14),REAL(DIAG(I)))
 332    CONTINUE
    ENDIF

    RETURN

    END


    SUBROUTINE ME57PD(A,IW,J1,J2,ITOP,REAL)
C This subroutine performs a very simple compress (block move).
C   Entries J1+1 to J2-1 (incl.) in A or IW as appropriate are moved
C   to occupy the positions immediately prior to position ITOP.
C A/IW hold the array being compressed.
C J1/J2 define the entries being moved.
C ITOP defines the position immediately after the positions to which
C   J1 to J2 are moved.
C REAL must be set by the user to .TRUE. if the move is on array A,
C   a value of .FALSE. will perform the move on A.
C   .. Scalar Arguments ..
    INTEGER ITOP,J1,J2
    LOGICAL REAL
C   ..
C   .. Array Arguments ..
    COMPLEX*16 A(*)
    INTEGER IW(*)
C   ..
C   .. Local Scalars ..
    INTEGER IPOS,JJ
C   ..
C   .. Executable Statements ..
    IF (J2.EQ.ITOP) GO TO 50
    IPOS = ITOP - 1
    IF (REAL) THEN
      DO 10 JJ = J2-1,J1+1,-1
      A(IPOS) = A(JJ)
      IPOS = IPOS - 1
   10   CONTINUE
    ELSE
      DO 20 JJ = J2-1,J1+1,-1
      IW(IPOS) = IW(JJ)
      IPOS = IPOS - 1
   20   CONTINUE
    ENDIF
    J2 = ITOP
    J1 = IPOS
   50 RETURN
    END
    SUBROUTINE ME57WD(KIND,A,LA,IW,LIW,NRLBDU)
C   Negate the entries of L, move the off-diagonal entries of the 2x2
C   pivots within full blocks to the end of A and update NRLBDU to
C   correspond.
    INTEGER KIND,LA,LIW
    COMPLEX*16 A(LA)
    INTEGER IW(LIW)

    INTEGER NRLBDU

C Constants
    DOUBLE PRECISION ZERO
    PARAMETER (ZERO=0.0D0)
    COMPLEX*16 CZERO
    PARAMETER (CZERO=(0.0D0,0.0D0))
C Local variables
    INTEGER APOS,IBLK,IROW,IWPOS,J,JPIV,NCOLS,NROWS
C APOS  Position in A of current diagonal entry.
C IBLK  Current block.
C IROW  Current row.
C IWPOS Current position in IW.
C J   Do loop variable
C JPIV  Used as a flag so that IPIV is incremented correctly after the
C     use of a 2 by 2 pivot.
C NCOLS Number of columns in the block.
C NROWS Number of rows in the block.

    APOS = 1
    IWPOS = 6
    DO 40 IBLK = 1,IW(3)
      NCOLS = IW(IWPOS-2)
      NROWS = IW(IWPOS-1)
      JPIV = 1
      DO 30 IROW = 1,NROWS
      JPIV = JPIV - 1
      IF (JPIV.EQ.1) GO TO 10
      IF (IW(IWPOS+IROW-1).LT.0) THEN
        JPIV = 2
        NRLBDU = NRLBDU + 1
CC Must reconjugate if KIND not equal to 1
        IF (KIND.EQ.1) THEN
          A(NRLBDU) = A(APOS+1)
        ELSE
          A(NRLBDU) = CONJG(A(APOS+1))
        ENDIF
        A(APOS+1) = ZERO
      END IF

   10   DO 20 J = APOS + 1,APOS + NROWS - IROW
        A(J) = -A(J)
   20   CONTINUE
      APOS = APOS + NROWS - IROW + 1
   30   CONTINUE
C Negate entries in rectangular block (was done earlier by MA47OD)
C     DO 35 J = APOS,APOS+NROWS*(NCOLS-NROWS)-1
C       A(J) = -A(J)
C  35   CONTINUE
      APOS = APOS + NROWS* (NCOLS-NROWS)
      IWPOS = IWPOS + NCOLS + 2
   40 CONTINUE
    END

    SUBROUTINE ME57XD(N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
   *            W,LW,IW1,ICNTL)
C This subroutine performs forward elimination using the factors
C   stored in FACT/IFACT by ME57BD.
C It is designed for efficiency on one right-hand side.
    INTEGER N,LFACT
    COMPLEX*16 FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT),LRHS,LW
    COMPLEX*16 W(LW),RHS(LRHS)
    INTEGER IW1(N),ICNTL(20)
C N   must be set to the order of the matrix. It is not altered.
C FACT   must be set to hold the real values corresponding to the
C   factors. This must be unchanged since the preceding call to
C   ME57BD. It is not altered.
C LFACT  length of array FACT. It is not altered.
C IFACT  holds the integer indexing information for the matrix factors
C   in FACT. This must be unchanged since the preceding call to
C   ME57BD. It is not altered.
C LIFACT length of array IFACT. It is not altered.
C RHS on input, must be set to hold the right hand side vector.  On
C   return, it will hold the modified vector following forward
C   elimination.
C LHS must be set to the leading dimension of array RHS.
C W   used as workspace to hold the components of the right hand
C   sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C   larger than INFO(21) as returned from ME57BD.
C IW1 need not be set on entry. On exit IW1(I) (I = 1,NBLK), where
C   NBLK = IFACT(3) is the number of block pivots, will
C   hold pointers to the beginning of each block pivot in array IFACT.
C ICNTL Not referenced except:
C   ICNTL(13) Threshold on number of columns in a block for using
C       addressing using Level 2 and Level 3 BLAS.
C

C Procedures
    INTRINSIC ABS
    EXTERNAL ZGEMV,ZTPSV

C Constant
    DOUBLE PRECISION ONE
    PARAMETER (ONE=1.0D0)
    COMPLEX*16 CONE
    PARAMETER (CONE=(1.0D0,0.0D0))
C
C Local variables
    INTEGER APOS,I,IBLK,II,IPIV,IRHS,IWPOS,J,J1,J2,K,K1,K2,
   +      NCOLS,NROWS
    COMPLEX*16 W1,W2
C
C APOS  Current position in array FACT.
C I   Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IWPOS Position in IFACT of start of current index list.
C J   Temporary DO index
C K   Temporary pointer to position in real array.
C J1    Position in IFACT of index of leading entry of row.
C J2    Position in IFACT of index of trailing entry of row.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.

    APOS = 1
    IWPOS = 4
    DO 270 IBLK = 1,IFACT(3)

C Find the number of rows and columns in the block.
      IW1(IBLK) = IWPOS
      NCOLS = IFACT(IWPOS)
      NROWS = IFACT(IWPOS+1)
      IWPOS = IWPOS + 2

      IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN


C Perform operations using direct addressing.

C Load appropriate components of right-hand sides into W.
      DO 10 I = 1,NCOLS
        II = ABS(IFACT(IWPOS+I-1))
        W(I) = RHS(II)
   10   CONTINUE

C Treat diagonal block (direct addressing)
      CALL ZTPSV('L','N','U',NROWS,FACT(APOS),W,1)
      APOS = APOS + (NROWS* (NROWS+1))/2

C Treat off-diagonal block (direct addressing)
C       IF (NCOLS.GT.NROWS) CALL ZGEMM('N','N',NCOLS-NROWS,1,NROWS,
C    +                      CONE,FACT(APOS),NCOLS-NROWS,
C    +                      W,LW,CONE,W(NROWS+1),LW)
      IF (NCOLS.GT.NROWS) THEN
        CALL ZGEMV('N',NCOLS-NROWS,NROWS,
   +           CONE,FACT(APOS),NCOLS-NROWS,
   +           W,1,CONE,W(NROWS+1),1)
      ENDIF
      APOS = APOS + NROWS* (NCOLS-NROWS)

C Reload W back into RHS.
      DO 35 I = 1,NCOLS
        II = ABS(IFACT(IWPOS+I-1))
        RHS(II) = W(I)
   35   CONTINUE

      ELSE

C Perform operations using indirect addressing.

      J1 = IWPOS
      J2 = IWPOS + NROWS - 1


C Treat diagonal block (indirect addressing)
      DO 130 IPIV = 1,NROWS
      APOS = APOS + 1
      W1 = RHS(ABS(IFACT(J1)))
      K = APOS
      DO 100 J = J1+1,J2
        IRHS = ABS(IFACT(J))
        RHS(IRHS) = RHS(IRHS) - FACT(K)*W1
        K = K + 1
  100   CONTINUE
      APOS = K
      J1 = J1 + 1
  130   CONTINUE

C Treat off-diagonal block (indirect addressing)
C     J2 = IWPOS + NCOLS - 1
C     DO 136 IPIV = 1,NROWS
C       K = APOS
C       W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)))
C       DO 133 J = J1,J2
C       IRHS = ABS(IFACT(J))
C       RHS(IRHS) = RHS(IRHS) + W1*FACT(K)
C       K = K + 1
C 133   CONTINUE
C       APOS = K
C 136   CONTINUE

C Loop unrolling
      J2 = IWPOS + NCOLS - 1
      DO 136 IPIV = 1,NROWS-1,2
      K1 = APOS
      K2 = APOS+NCOLS-NROWS
      W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)))
      W2 = RHS(ABS(IFACT(IWPOS+IPIV)))
      DO 133 J = J1,J2
        IRHS = ABS(IFACT(J))
        RHS(IRHS) = RHS(IRHS) + W1*FACT(K1) + W2*FACT(K2)
        K1 = K1 + 1
        K2 = K2 + 1
  133   CONTINUE
      APOS = K2
  136   CONTINUE

      IF (MOD(NROWS,2).EQ.1) THEN
      K = APOS
      W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)))
      DO 137 J = J1,J2
        IRHS = ABS(IFACT(J))
        RHS(IRHS) = RHS(IRHS) + W1*FACT(K)
        K = K + 1
  137   CONTINUE
      APOS = K
      ENDIF
    END IF

    IWPOS = IWPOS + NCOLS
  270 CONTINUE

    END


    SUBROUTINE ME57YD(KIND,N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
   *            W,LW,IW1,ICNTL)
C This subroutine performs backward elimination operations
C   using the factors stored in FACT/IFACT by ME57BD.
C It is designed for efficiency on one right-hand side.
    INTEGER KIND,N,LFACT
    COMPLEX*16 FACT(LFACT)
    INTEGER LIFACT,IFACT(LIFACT),LRHS,LW
    COMPLEX*16 W(LW),RHS(LRHS)
    INTEGER IW1(N),ICNTL(20)
C N    must be set to the order of the matrix. It is not altered.
C FACT    must be set to hold the real values corresponding to the
C    factors. This must be unchanged since the
C    preceding call to ME57BD. It is not altered.
C LFACT   length of array FACT. It is not altered.
C IFACT   holds the integer indexing information for the matrix factors
C    in FACT. This must be unchanged since the preceding call to
C    ME57BD. It is not altered.
C LIFACT  length of array IFACT. It is not altered.
C RHS  on entry, must be set to hold the right hand side modified by
C    the forward substitution operations. On exit, holds the
C    solution vector.
C LHS must be set to the leading dimension of array RHS.
C W    used as workspace to hold the components of the right hand
C    sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C   larger than INFO(21) as returned from ME57BD.
C IW1  on entry IW1(I) (I = 1,NBLK), where  NBLK = IFACT(3) is the
C   number of block pivots, must hold pointers to the beginning of
C   each block pivot in array IFACT, as set by ME57X/XD. It is not
C   altered.
C ICNTL Not referenced except:
C   ICNTL(13) Threshold on number of columns in a block for using
C       addressing using Level 2 and Level 3 BLAS.

C Procedures
    INTRINSIC ABS
    EXTERNAL ZGEMV,ZTPSV

C Constants
    DOUBLE PRECISION ONE
    PARAMETER (ONE=1.0D0)
    COMPLEX*16 CONE
    PARAMETER (CONE=(1.0D0,0.0D0))
C
C Local variables.
    INTEGER APOS,APOS2,I,IBLK,II,IPIV,IRHS,IRHS1,
   +      IRHS2,IWPOS,J,JPIV,J1,J2,K,K2,LROW,NCOLS,NROWS
    COMPLEX*16 W1,W2
C APOS  Current position in array FACT.
C APOS2 Current position in array FACT for off-diagonal entry of 2x2
C     pivot.
C I   Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IRHS1 RHS index.
C IRHS2 RHS index.
C IWPOS Position in IFACT of start of current index list.
C J   Temporary DO index
C JPIV  Has the value -1 for the first row of a 2 by 2 pivot and 1 for
C     the second.
C K   Temporary pointer to position in real array.
C J1    Position in IFACT of index of leading entry of row.
C J2    Position in IFACT of index of trailing entry of row.
C K   Temporary variable.
C LROW  Length of current row.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.
C

    APOS = IFACT(1)
    APOS2 = IFACT(2)
C Run through block pivot rows in the reverse order.
    DO 380 IBLK = IFACT(3),1,-1

C Find the number of rows and columns in the block.
      IWPOS = IW1(IBLK)
      NCOLS = ABS(IFACT(IWPOS))
      NROWS = ABS(IFACT(IWPOS+1))
      APOS = APOS - NROWS* (NCOLS-NROWS)
      IWPOS = IWPOS + 2

      IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN


C Perform operations using direct addressing.

C Load latter part of right-hand side into W.
      DO 5 I = NROWS + 1,NCOLS
        II = ABS(IFACT(IWPOS+I-1))
        W(I) = RHS(II)
    5   CONTINUE


C Multiply by the diagonal matrix (direct addressing)
      DO 10 IPIV = NROWS,1,-1
        IRHS = ABS(IFACT(IWPOS+IPIV-1))
        APOS = APOS - (NROWS+1-IPIV)
        W(IPIV) = RHS(IRHS)*FACT(APOS)
   10   CONTINUE
      JPIV = -1
      DO 20 IPIV = NROWS,1,-1
        IRHS = IFACT(IWPOS+IPIV-1)
        IF (IRHS.LT.0) THEN
          IRHS1 = -IFACT(IWPOS+IPIV-1+JPIV)
          IF (KIND.EQ.1 .OR. JPIV.EQ.1) THEN
          W(IPIV) = RHS(IRHS1)*FACT(APOS2) + W(IPIV)
          ELSE
          W(IPIV) = RHS(IRHS1)*CONJG(FACT(APOS2)) + W(IPIV)
          ENDIF
          IF (JPIV.EQ.1) APOS2 = APOS2 - 1
          JPIV = -JPIV
        END IF

   20   CONTINUE

C Treat off-diagonal block (direct addressing)
      K = NCOLS - NROWS
      IF (K.GT.0) THEN
        IF (KIND.EQ.1) THEN
          CALL ZGEMV('T',K,NROWS,CONE,
   +             FACT(APOS+(NROWS*(NROWS+1))/2),K,
   +             W(NROWS+1),1,CONE,W,1)
        ELSE
          CALL ZGEMV('C',K,NROWS,CONE,
   +             FACT(APOS+(NROWS*(NROWS+1))/2),K,
   +             W(NROWS+1),1,CONE,W,1)
        ENDIF
      ENDIF
C       IF (K.GT.0) CALL ZGEMM('T','N',NROWS,1,K,CONE,
C    +                 FACT(APOS+(NROWS*(NROWS+1))/2),K,
C    +                 W(NROWS+1),LW,CONE,W,LW)

C Treat diagonal block (direct addressing)
      IF (KIND.EQ.1) CALL ZTPSV('L','T','U',NROWS,FACT(APOS),W,1)
      IF (KIND.NE.1) CALL ZTPSV('L','C','U',NROWS,FACT(APOS),W,1)
C Reload W back into RHS.
      DO 60 I = 1,NROWS
        II = ABS(IFACT(IWPOS+I-1))
        RHS(II) = W(I)
   60   CONTINUE

      ELSE
C
C Perform operations using indirect addressing.
      J1 = IWPOS
      J2 = IWPOS + NCOLS - 1


C Multiply by the diagonal matrix (indirect addressing)
      JPIV = -1
      IF (KIND.EQ.1) THEN
        DO 210 IPIV = NROWS,1,-1
          IRHS = IFACT(IWPOS+IPIV-1)
          LROW = NROWS + 1 - IPIV

          IF (IRHS.GT.0) THEN
C 1 by 1 pivot.
          APOS = APOS - LROW
          RHS(IRHS) = RHS(IRHS)*FACT(APOS)
          ELSE
C 2 by 2 pivot
          IF (JPIV.EQ.-1) THEN
            IRHS1 = -IFACT(IWPOS+IPIV-2)
            IRHS2 = -IRHS
            APOS = APOS - LROW - LROW - 1
            W1 = RHS(IRHS1)*FACT(APOS) +
   +           RHS(IRHS2)*FACT(APOS2)
            RHS(IRHS2) = RHS(IRHS1)*FACT(APOS2) +
   +               RHS(IRHS2)*FACT(APOS+LROW+1)
            RHS(IRHS1) = W1
            APOS2 = APOS2 - 1
          END IF
          JPIV = -JPIV
          END IF

  210     CONTINUE
      ELSE
        DO 211 IPIV = NROWS,1,-1
          IRHS = IFACT(IWPOS+IPIV-1)
          LROW = NROWS + 1 - IPIV

          IF (IRHS.GT.0) THEN
C 1 by 1 pivot.
          APOS = APOS - LROW
          RHS(IRHS) = RHS(IRHS)*FACT(APOS)
          ELSE
C 2 by 2 pivot
          IF (JPIV.EQ.-1) THEN
            IRHS1 = -IFACT(IWPOS+IPIV-2)
            IRHS2 = -IRHS
            APOS = APOS - LROW - LROW - 1
            W1 = RHS(IRHS1)*FACT(APOS) +
   +           RHS(IRHS2)*FACT(APOS2)
            RHS(IRHS2) = RHS(IRHS1)*CONJG(FACT(APOS2)) +
   +               RHS(IRHS2)*FACT(APOS+LROW+1)
            RHS(IRHS1) = W1
            APOS2 = APOS2 - 1
          END IF
          JPIV = -JPIV
          END IF

  211     CONTINUE
      ENDIF
      APOS = APOS + (NROWS* (NROWS+1))/2

C Treat off-diagonal block (indirect addressing)
C       KK = APOS
C       J1 = IWPOS + NROWS
C       DO 220 IPIV = 1,NROWS
C       IRHS = ABS(IFACT(IWPOS+IPIV-1))
C       W1 = RHS(IRHS)
C       K = KK
C       DO 215 J = J1,J2
C         W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)))
C         K = K + 1
C 215     CONTINUE
C       RHS(IRHS) = W1
C       KK = K
C 220   CONTINUE

C Loop unrolling
      K = APOS
      J1 = IWPOS + NROWS
      DO 220 IPIV = 1,NROWS-1,2
        IRHS = ABS(IFACT(IWPOS+IPIV-1))
        W1 = RHS(IRHS)
        IRHS1 = ABS(IFACT(IWPOS+IPIV))
        W2 = RHS(IRHS1)
        K2 = K+(NCOLS-NROWS)
        DO 215 J = J1,J2
          II = ABS(IFACT(J))
          IF (KIND.EQ.1) THEN
          W1 = W1 + FACT(K)*RHS(II)
          W2 = W2 + FACT(K2)*RHS(II)
          ELSE
          W1 = W1 + CONJG(FACT(K))*RHS(II)
          W2 = W2 + CONJG(FACT(K2))*RHS(II)
          ENDIF
          K = K + 1
          K2 = K2 + 1
  215     CONTINUE
        RHS(IRHS) = W1
        RHS(IRHS1) = W2
        K = K2
  220   CONTINUE

      IF (MOD(NROWS,2).EQ.1) THEN
        IRHS = ABS(IFACT(IWPOS+IPIV-1))
        W1 = RHS(IRHS)
        DO 216 J = J1,J2
          IF (KIND.EQ.1) THEN
          W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)))
          ELSE
          W1 = W1 + CONJG(FACT(K))*RHS(ABS(IFACT(J)))
          ENDIF
          K = K + 1
  216     CONTINUE
        RHS(IRHS) = W1
      ENDIF

C Treat diagonal block (indirect addressing)
      J2 = IWPOS + NROWS - 1
      DO 260 IPIV = 1,NROWS
        IRHS = ABS(IFACT(J1-1))
        APOS = APOS - IPIV
        W1 = RHS(IRHS)
        K = APOS + 1
        DO 230 J = J1,J2
          IF (KIND.EQ.1) THEN
          W1 = W1 - FACT(K)*RHS(ABS(IFACT(J)))
          ELSE
          W1 = W1 - CONJG(FACT(K))*RHS(ABS(IFACT(J)))
          ENDIF
          K = K + 1
  230     CONTINUE
        RHS(IRHS) = W1
        J1 = J1 - 1
  260   CONTINUE

      END IF

  380 CONTINUE

    END
    SUBROUTINE ME57VD(KIND,N,NZ,IRN,ICN,IW,LW,IPE,IQ,FLAG,IWFR,
   +           ICNTL,INFO)
C Is identical to subroutine MA27GD.  Internal version for ME57.
C Except for addition of parameter KIND and flag if entry in lower
C   triangle
C
C SORT PRIOR TO CALLING ANALYSIS ROUTINE MA27H/HD (internal ME57
C subroutine ME57H/HD).
C
C GIVEN THE POSITIONS OF THE OFF-DIAGONAL NON-ZEROS OF A SYMMETRIC
C   MATRIX, CONSTRUCT THE SPARSITY PATTERN OF THE OFF-DIAGONAL
C   PART OF THE WHOLE MATRIX (UPPER AND LOWER TRIANGULAR PARTS).
C EITHER ONE OF A PAIR (I,J),(J,I) MAY BE USED TO REPRESENT
C   THE PAIR. DIAGONAL ELEMENTS AND DUPLICATES ARE IGNORED.
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT. IT IS NOT
C   ALTERED.
C IRN(I),I=1,2,...,NZ MUST BE SET TO THE ROW NUMBERS OF THE
C   NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
C   TO IW (SEE DESCRIPTION OF IW).
C ICN(I),I=1,2,...,NZ MUST BE SET TO THE COLUMN NUMBERS OF THE
C   NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
C   TO IW (SEE DESCRIPTION OF IW).
C IW NEED NOT BE SET ON INPUT. ON OUTPUT IT CONTAINS LISTS OF
C   COLUMN INDICES, EACH LIST BEING HEADED BY ITS LENGTH.
C   IRN(1) MAY BE EQUIVALENCED TO IW(1) AND ICN(1) MAY BE
C   EQUIVALENCED TO IW(K), WHERE K.GT.NZ.
C LW MUST BE SET TO THE LENGTH OF IW. IT MUST BE AT LEAST 2*NZ+N.
C   IT IS NOT ALTERED.
C IPE NEED NOT BE SET ON INPUT. ON OUTPUT IPE(I) POINTS TO THE START OF
C   THE ENTRY IN IW FOR ROW I, OR IS ZERO IF THERE IS NO ENTRY.
C IQ NEED NOT BE SET.  ON OUTPUT IQ(I),I=1,N CONTAINS THE NUMBER OF
C   OFF-DIAGONAL N0N-ZEROS IN ROW I INCLUDING DUPLICATES.
C FLAG IS USED FOR WORKSPACE TO HOLD FLAGS TO PERMIT DUPLICATE ENTRIES
C   TO BE IDENTIFIED QUICKLY.
C IWFR NEED NOT BE SET ON INPUT. ON OUTPUT IT POINTS TO THE FIRST
C   UNUSED LOCATION IN IW.
C ICNTL is an INTEGER array of assumed size.
C INFO is an INTEGER array of assumed size.
C
C   .. Scalar Arguments ..
    INTEGER KIND
    INTEGER IWFR,LW,N,NZ
C   ..
C   .. Array Arguments ..
    INTEGER FLAG(N),ICN(*),IPE(N),IQ(N),IRN(*),IW(LW)
    INTEGER ICNTL(*),INFO(*)
C   ..
C   .. Local Scalars ..
    INTEGER I,ID,J,JN,K,K1,K2,L,LAST,LDIAG,LP,LR,N1,NDUP
C   ..
C   .. Executable Statements ..
C
C INITIALIZE INFO(2) AND COUNT IN IPE THE
C   NUMBERS OF NON-ZEROS IN THE ROWS AND MOVE ROW AND COLUMN
C   NUMBERS INTO IW.
    INFO(2) = 0
    DO 10 I = 1,N
      IPE(I) = 0
   10 CONTINUE
    LP = ICNTL(1)
    LDIAG = ICNTL(5)
    LR = NZ
    IF (NZ.EQ.0) GO TO 120
    DO 110 K = 1,NZ
      I = IRN(K)
      J = ICN(K)
      IF (KIND.NE.1 .AND. I.GT.J) THEN
      INFO(1) = -19
    IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I8,A)')
   +    '**** Error return from ME57AD ****  INFO(1) =',INFO(1),
   *    'Hermitian matrix but entry ',K,' in lower triangle'
      INFO(2) = K
      RETURN
      ENDIF
      IF (I.LT.J) THEN
      IF (I.GE.1 .AND. J.LE.N) GO TO 90
      ELSE IF (I.GT.J) THEN
      IF (J.GE.1 .AND. I.LE.N) GO TO 90
      ELSE
      IF (I.GE.1 .AND. I.LE.N) GO TO 80
      END IF
      INFO(2) = INFO(2) + 1
      INFO(1) = 1
      IF (INFO(2).LE.1 .AND. ICNTL(2).GT.0) THEN
      WRITE (ICNTL(2),FMT=60) INFO(1)
      END IF

   60   FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE ME57AD',
   +      '  *** INFO(1) =',I2)

      IF (INFO(2).LE.10 .AND. ICNTL(2).GT.0) THEN
      WRITE (ICNTL(2),FMT=70) K,I,J
      END IF

   70   FORMAT (I6,'TH NON-ZERO (IN ROW',I6,' AND COLUMN',I6,
   +       ') IGNORED')

   80   I = 0
      J = 0
      GO TO 100

   90   IPE(I) = IPE(I) + 1
      IPE(J) = IPE(J) + 1
  100   IW(K) = J
      LR = LR + 1
      IW(LR) = I
  110 CONTINUE
C
C ACCUMULATE ROW COUNTS TO GET POINTERS TO ROW STARTS IN BOTH IPE AND IQ
C   AND INITIALIZE FLAG
  120 IQ(1) = 1
    N1 = N - 1
    IF (N1.LE.0) GO TO 140
    DO 130 I = 1,N1
      FLAG(I) = 0
      IF (IPE(I).EQ.0) IPE(I) = -1
      IQ(I+1) = IPE(I) + IQ(I) + 1
      IPE(I) = IQ(I)
  130 CONTINUE
  140 LAST = IPE(N) + IQ(N)
    FLAG(N) = 0
    IF (LR.GE.LAST) GO TO 160
    K1 = LR + 1
    DO 150 K = K1,LAST
      IW(K) = 0
  150 CONTINUE
  160 IPE(N) = IQ(N)
    IWFR = LAST + 1
    IF (NZ.EQ.0) GO TO 230
C
C RUN THROUGH PUTTING THE MATRIX ELEMENTS IN THE RIGHT PLACE
C   BUT WITH SIGNS INVERTED. IQ IS USED FOR HOLDING RUNNING POINTERS
C   AND IS LEFT HOLDING POINTERS TO ROW ENDS.
    DO 220 K = 1,NZ
      J = IW(K)
      IF (J.LE.0) GO TO 220
      L = K
      IW(K) = 0
      DO 210 ID = 1,NZ
      IF (L.GT.NZ) GO TO 170
      L = L + NZ
      GO TO 180

  170   L = L - NZ
  180   I = IW(L)
      IW(L) = 0
      IF (I.LT.J) GO TO 190
      L = IQ(J) + 1
      IQ(J) = L
      JN = IW(L)
      IW(L) = -I
      GO TO 200

  190   L = IQ(I) + 1
      IQ(I) = L
      JN = IW(L)
      IW(L) = -J
  200   J = JN
      IF (J.LE.0) GO TO 220
  210   CONTINUE
  220 CONTINUE
C
C RUN THROUGH RESTORING SIGNS, REMOVING DUPLICATES AND SETTING THE
C   MATE OF EACH NON-ZERO.
C NDUP COUNTS THE NUMBER OF DUPLICATE ELEMENTS.
  230 NDUP = 0
    DO 280 I = 1,N
      K1 = IPE(I) + 1
      K2 = IQ(I)
      IF (K1.LE.K2) GO TO 240
C ROW IS EMPTY. SET POINTER TO ZERO.
      IPE(I) = 0
      IQ(I) = 0
      GO TO 280
C ON ENTRY TO THIS LOOP FLAG(J).LT.I FOR J=1,2,...,N. DURING THE LOOP
C   FLAG(J) IS SET TO I IF A NON-ZERO IN COLUMN J IS FOUND. THIS
C   PERMITS DUPLICATES TO BE RECOGNIZED QUICKLY.
  240   DO 260 K = K1,K2
      J = -IW(K)
      IF (J.LE.0) GO TO 270
      L = IQ(J) + 1
      IQ(J) = L
      IW(L) = I
      IW(K) = J
      IF (FLAG(J).NE.I) GO TO 250
      NDUP = NDUP + 1
      IW(L) = 0
      IW(K) = 0
  250   FLAG(J) = I
  260   CONTINUE
  270   IQ(I) = IQ(I) - IPE(I)
      IF (NDUP.EQ.0) IW(K1-1) = IQ(I)
  280 CONTINUE
    IF (NDUP.EQ.0) GO TO 310
C
C COMPRESS IW TO REMOVE DUMMY ENTRIES CAUSED BY DUPLICATES.
    IWFR = 1
    DO 300 I = 1,N
      K1 = IPE(I) + 1
      IF (K1.EQ.1) GO TO 300
      K2 = IQ(I) + IPE(I)
      L = IWFR
      IPE(I) = IWFR
      IWFR = IWFR + 1
      DO 290 K = K1,K2
      IF (IW(K).EQ.0) GO TO 290
      IW(IWFR) = IW(K)
      IWFR = IWFR + 1
  290   CONTINUE
      IW(L) = IWFR - L - 1
  300 CONTINUE
  310 RETURN

    END
    SUBROUTINE ME57HD(N,IPE,IW,LW,IWFR,NV,NXT,LST,IPD,FLAG,IOVFLO,
   +           NCMPA,FRATIO)
C Was identical to subroutine ME27HD.  Internal version for ME57.
C   Changes made in September 2009 because of bug in compress control
C   found by Nick.
C
C ANALYSIS SUBROUTINE
C
C GIVEN REPRESENTATION OF THE WHOLE MATRIX (EXCLUDING DIAGONAL)
C   PERFORM MINIMUM DEGREE ORDERING, CONSTRUCTING TREE POINTERS.
C   IT WORKS WITH SUPERVARIABLES WHICH ARE COLLECTIONS OF ONE OR MORE
C   VARIABLES, STARTING WITH SUPERVARIABLE I CONTAINING VARIABLE I FOR
C   I = 1,2,...,N. ALL VARIABLES IN A SUPERVARIABLE ARE ELIMINATED
C   TOGETHER. EACH SUPERVARIABLE HAS AS NUMERICAL NAME THAT OF ONE
C   OF ITS VARIABLES (ITS PRINCIPAL VARIABLE).
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C IPE(I) MUST BE SET TO POINT TO THE POSITION IN IW OF THE
C   START OF ROW I OR HAVE THE VALUE ZERO IF ROW I HAS NO OFF-
C   DIAGONAL NON-ZEROS. DURING EXECUTION IT IS USED AS FOLLOWS. IF
C   SUPERVARIABLE I IS ABSORBED INTO SUPERVARIABLE J THEN IPE(I)=-J.
C   IF SUPERVARIABLE I IS ELIMINATED THEN IPE(I) EITHER POINTS TO THE
C   LIST OF SUPERVARIABLES FOR CREATED ELEMENT I OR IS ZERO IF
C   THE CREATED ELEMENT IS NULL. IF ELEMENT I
C   IS ABSORBED INTO ELEMENT J THEN IPE(I)=-J.
C IW MUST BE SET ON ENTRY TO HOLD LISTS OF VARIABLES BY
C   ROWS, EACH LIST BEING HEADED BY ITS LENGTH.
C   DURING EXECUTION THESE LISTS ARE REVISED AND HOLD
C   LISTS OF ELEMENTS AND SUPERVARIABLES. THE ELEMENTS
C   ALWAYS HEAD THE LISTS. WHEN A SUPERVARIABLE
C   IS ELIMINATED ITS LIST IS REPLACED BY A LIST OF SUPERVARIABLES
C   IN THE NEW ELEMENT.
C LW MUST BE SET TO THE LENGTH OF IW. IT IS NOT ALTERED.
C IWFR MUST BE SET TO THE POSITION IN IW OF THE FIRST FREE VARIABLE.
C   IT IS REVISED DURING EXECUTION AND CONTINUES TO HAVE THIS MEANING.
C NV(JS) NEED NOT BE SET. DURING EXECUTION IT IS ZERO IF
C   JS IS NOT A PRINCIPAL VARIABLE AND IF IT IS IT HOLDS
C   THE NUMBER OF VARIABLES IN SUPERVARIABLE JS. FOR ELIMINATED
C   VARIABLES IT IS SET TO THE DEGREE AT THE TIME OF ELIMINATION.
C NXT(JS) NEED NOT BE SET. DURING EXECUTION IT IS THE NEXT
C   SUPERVARIABLE HAVING THE SAME DEGREE AS JS, OR ZERO
C   IF IT IS LAST IN ITS LIST.
C LST(JS) NEED NOT BE SET. DURING EXECUTION IT IS THE
C   LAST SUPERVARIABLE HAVING THE SAME DEGREE AS JS OR
C   -(ITS DEGREE) IF IT IS FIRST IN ITS LIST.
C IPD(ID) NEED NOT BE SET. DURING EXECUTION IT
C   IS THE FIRST SUPERVARIABLE WITH DEGREE ID OR ZERO
C   IF THERE ARE NONE.
C FLAG IS USED AS WORKSPACE FOR ELEMENT AND SUPERVARIABLE FLAGS.
C   WHILE THE CODE IS FINDING THE DEGREE OF SUPERVARIABLE IS
C   FLAG HAS THE FOLLOWING VALUES.
C   A) FOR THE CURRENT PIVOT/NEW ELEMENT ME
C       FLAG(ME)=-1
C   B) FOR VARIABLES JS
C       FLAG(JS)=-1 IF JS IS NOT A PRINCIPAL VARIABLE
C       FLAG(JS)=0 IF JS IS A SUPERVARIABLE IN THE NEW ELEMENT
C       FLAG(JS)=NFLG IF JS IS A SUPERVARIABLE NOT IN THE NEW
C           ELEMENT THAT HAS BEEN COUNTED IN THE DEGREE
C           CALCULATION
C       FLAG(JS).GT.NFLG IF JS IS A SUPERVARIABLE NOT IN THE NEW
C           ELEMENT THAT HAS NOT BEEN COUNTED IN THE DEGREE
C           CALCULATION
C   C) FOR ELEMENTS IE
C       FLAG(IE)=-1 IF ELEMENT IE HAS BEEN MERGED INTO ANOTHER
C       FLAG(IE)=-NFLG IF ELEMENT IE HAS BEEN USED IN THE DEGREE
C           CALCULATION FOR IS.
C       FLAG(IE).LT.-NFLG IF ELEMENT IE HAS NOT BEEN USED IN THE
C           DEGREE CALCULATION FOR IS
C IOVFLO should be set to a high legitimate integer.  It is used as a
C      flag.
C NCMPA number of compresses.
C FRATIO is set to ICNTL(14)/100 and is the density of rows regarded as
C   dense.
C
C   .. Scalar Arguments ..
    DOUBLE PRECISION FRATIO
    INTEGER IWFR,LW,N,IOVFLO,NCMPA
C   ..
C   .. Array Arguments ..
    INTEGER FLAG(N),IPD(N),IPE(N),IW(LW),LST(N),NV(N),NXT(N)
C   ..
C   .. Local Scalars ..
    INTEGER I,ID,IDL,IDN,IE,IP,IS,JP,JP1,JP2,JS,K,K1,K2,KE,KP,KP0,KP1,
   +      KP2,KS,L,LEN,LIMIT,LN,LS,LWFR,MD,ME,ML,MS,NEL,NFLG,NP,
   +      NP0,NS,NVPIV,NVROOT,ROOT
C LIMIT  Limit on number of variables for putting node in root.
C NVROOT Number of variables in the root node
C ROOT   Index of the root node (N+1 if none chosen yet).
C   ..
C   .. External Subroutines ..
    EXTERNAL ME57ZD
C   ..
C   .. Intrinsic Functions ..
    INTRINSIC ABS,MIN
C   ..
C If a column of the reduced matrix has relative density greater than
C CNTL(2), it is forced into the root. All such columns are taken to
C have sparsity pattern equal to their merged patterns, so the fill
C and operation counts may be overestimated.
C
C IS,JS,KS,LS,MS,NS are used to refer to supervariables.
C IE,JE,KE are used to refer to elements.
C IP,JP,KP,K,NP are used to point to lists of elements
C   or supervariables.
C ID is used for the degree of a supervariable.
C MD is used for the current minimum degree.
C IDN is used for the no. of variables in a newly created element
C NEL is used to hold the no. of variables that have been
C   eliminated.
C ME=MS is the name of the supervariable eliminated and
C   of the element created in the main loop.
C NFLG is used for the current flag value in array FLAG. It starts
C   with the value IOVFLO and is reduced by 1 each time it is used
C   until it has the value 2 when it is reset to the value IOVFLO.
C
C   .. Executable Statements ..
C Initializations
    DO 10 I = 1,N
      IPD(I) = 0
      NV(I) = 1
      FLAG(I) = IOVFLO
   10 CONTINUE
    MD = 1
    NCMPA = 0
    NFLG = IOVFLO
    NEL = 0
    ROOT = N+1
    NVROOT = 0
C
C Link together variables having same degree
    DO 30 IS = 1,N
      K = IPE(IS)
      IF (K.GT.0) THEN
      ID = IW(K) + 1
      NS = IPD(ID)
      IF (NS.GT.0) LST(NS) = IS
      NXT(IS) = NS
      IPD(ID) = IS
      LST(IS) = -ID
      ELSE
C We have a variable that can be eliminated at once because there is
C   no off-diagonal nonzero in its row.
      NEL = NEL + 1
      FLAG(IS) = -1
      NXT(IS) = 0
      LST(IS) = 0
      ENDIF
   30 CONTINUE

C
C Start of main loop
C
    DO 340 ML = 1,N

C Leave loop if all variables have been eliminated.
      IF (NEL+NVROOT+1.GE.N) GO TO 350
C
C Find next supervariable for elimination.
      DO 40 ID = MD,N
      MS = IPD(ID)
      IF (MS.GT.0) GO TO 50
   40   CONTINUE
   50   MD = ID
C Nvpiv holds the number of variables in the pivot.
      NVPIV = NV(MS)
C
C Remove chosen variable from linked list
      NS = NXT(MS)
      NXT(MS) = 0
      LST(MS) = 0
      IF (NS.GT.0) LST(NS) = -ID
      IPD(ID) = NS
      ME = MS
      NEL = NEL + NVPIV
C IDN holds the degree of the new element.
      IDN = 0
C
C Run through the list of the pivotal supervariable, setting tree
C   pointers and constructing new list of supervariables.
C KP is a pointer to the current position in the old list.
      KP = IPE(ME)
      FLAG(MS) = -1
C IP points to the start of the new list.
      IP = IWFR
C LEN holds the length of the list associated with the pivot.
      LEN = IW(KP)
      DO 140 KP1 = 1,LEN
      KP = KP + 1
      KE = IW(KP)
C Jump if KE is an element that has not been merged into another.
      IF (FLAG(KE).LE.-2) GO TO 60
C Jump if KE is an element that has been merged into another or is
C   a supervariable that has been eliminated.
      IF (FLAG(KE).LE.0) THEN
         IF (IPE(KE).NE.-ROOT) GO TO 140
C KE has been merged into the root
         KE = ROOT
         IF (FLAG(KE).LE.0) GO TO 140
      END IF
C We have a supervariable. Prepare to search rest of list.
      JP = KP - 1
      LN = LEN - KP1 + 1
      IE = MS
      GO TO 70
C Search variable list of element KE, using JP as a pointer to it.
   60   IE = KE
      JP = IPE(IE)
      LN = IW(JP)
C
C Search for different supervariables and add them to the new list,
C   compressing when necessary. This loop is executed once for
C   each element in the list and once for all the supervariables
C   in the list.
   70   DO 130 JP1 = 1,LN
        JP = JP + 1
        IS = IW(JP)
C Jump if IS is not a principal variable or has already been counted.
        IF (FLAG(IS).LE.0) THEN
         IF (IPE(IS).EQ.-ROOT) THEN
C IS has been merged into the root
            IS = ROOT
            IW(JP) = ROOT
            IF (FLAG(IS).LE.0) GO TO 130
         ELSE
            GO TO 130
         END IF
        END IF
        FLAG(IS) = 0

C To fix Nick bug need to add one here to store (eventually) length
C   of new row
        IF (IWFR .GE. LW-1) THEN
C Logic was previously as below
CCC       IF (IWFR.LT.LW) GO TO 100
C Prepare for compressing IW by adjusting pointers and
C   lengths so that the lists being searched in the inner and outer
C   loops contain only the remaining entries.
          IPE(MS) = KP
          IW(KP) = LEN - KP1
          IPE(IE) = JP
          IW(JP) = LN - JP1
C Compress IW
          CALL ME57ZD(N,IPE,IW,IP-1,LWFR,NCMPA)
C Copy new list forward
          JP2 = IWFR - 1
          IWFR = LWFR
          IF (IP.GT.JP2) GO TO 90
          DO 80 JP = IP,JP2
          IW(IWFR) = IW(JP)
          IWFR = IWFR + 1
   80       CONTINUE
C Adjust pointers for the new list and the lists being searched.
   90       IP = LWFR
          JP = IPE(IE)
          KP = IPE(ME)
        ENDIF

C Store IS in new list.
        IW(IWFR) = IS
        IDN = IDN + NV(IS)
        IWFR = IWFR + 1
C Remove IS from degree linked list
        LS = LST(IS)
        LST(IS) = 0
        NS = NXT(IS)
        NXT(IS) = 0
        IF (NS.GT.0) LST(NS) = LS
        IF (LS.LT.0) THEN
          LS = -LS
          IPD(LS) = NS
        ELSE IF (LS.GT.0) THEN
          NXT(LS) = NS
        END IF
  130   CONTINUE
C Jump if we have just been searching the variables at the end of
C   the list of the pivot.
      IF (IE.EQ.MS) GO TO 150
C Set tree pointer and flag to indicate element IE is absorbed into
C   new element ME.
      IPE(IE) = -ME
      FLAG(IE) = -1
  140   CONTINUE

C Store the degree of the pivot.
  150   NV(MS) = IDN + NVPIV

C Jump if new element is null.
      IF (IWFR.EQ.IP) THEN
      IPE(ME) = 0
      GO TO 340
      ENDIF

      K1 = IP
      K2 = IWFR - 1
C
C Run through new list of supervariables revising each associated list,
C   recalculating degrees and removing duplicates.
      LIMIT = NINT(FRATIO*(N-NEL))
      DO 310 K = K1,K2
      IS = IW(K)
      IF (IS.EQ.ROOT) GO TO 310
      IF (NFLG.GT.2) GO TO 170
C Reset FLAG values to +/-IOVFLO.
      DO 160 I = 1,N
        IF (FLAG(I).GT.0) FLAG(I) = IOVFLO
        IF (FLAG(I).LE.-2) FLAG(I) = -IOVFLO
  160   CONTINUE
      NFLG = IOVFLO
C Reduce NFLG by one to cater for this supervariable.
  170   NFLG = NFLG - 1
C Begin with the degree of the new element. Its variables must always
C   be counted during the degree calculation and they are already
C   flagged with the value 0.
      ID = IDN
C Run through the list associated with supervariable IS
      KP1 = IPE(IS) + 1
C NP points to the next entry in the revised list.
      NP = KP1
      KP2 = IW(KP1-1) + KP1 - 1
      DO 220 KP = KP1,KP2
        KE = IW(KP)
C Test whether KE is an element, a redundant entry or a supervariable.
        IF (FLAG(KE).EQ.-1) THEN
          IF (IPE(KE).NE.-ROOT) GO TO 220
C KE has been merged into the root
          KE = ROOT
          IW(KP) = ROOT
          IF (FLAG(KE).EQ.-1) GO TO 220
        END IF
        IF (FLAG(KE).GE.0) GO TO 230
C Search list of element KE, revising the degree when new variables
C   found.
        JP1 = IPE(KE) + 1
        JP2 = IW(JP1-1) + JP1 - 1
        IDL = ID
        DO 190 JP = JP1,JP2
          JS = IW(JP)
C Jump if JS has already been counted.
          IF (FLAG(JS).LE.NFLG) GO TO 190
          ID = ID + NV(JS)
          FLAG(JS) = NFLG
  190     CONTINUE
C Jump if one or more new supervariables were found.
        IF (ID.GT.IDL) GO TO 210
C Check whether every variable of element KE is in new element ME.
        DO 200 JP = JP1,JP2
          JS = IW(JP)
          IF (FLAG(JS).NE.0) GO TO 210
  200     CONTINUE
C Set tree pointer and FLAG to indicate that element KE is absorbed
C   into new element ME.
        IPE(KE) = -ME
        FLAG(KE) = -1
        GO TO 220
C Store element KE in the revised list for supervariable IS and flag it.
  210     IW(NP) = KE
        FLAG(KE) = -NFLG
        NP = NP + 1
  220   CONTINUE
      NP0 = NP
      GO TO 250
C Treat the rest of the list associated with supervariable IS. It
C   consists entirely of supervariables.
  230   KP0 = KP
      NP0 = NP
      DO 240 KP = KP0,KP2
        KS = IW(KP)
        IF (FLAG(KS).LE.NFLG) THEN
         IF (IPE(KS).EQ.-ROOT) THEN
            KS = ROOT
            IW(KP) = ROOT
            IF (FLAG(KS).LE.NFLG) GO TO 240
         ELSE
            GO TO 240
         END IF
        END IF
C Add to degree, flag supervariable KS and add it to new list.
        ID = ID + NV(KS)
        FLAG(KS) = NFLG
        IW(NP) = KS
        NP = NP + 1
  240   CONTINUE
C Move first supervariable to end of list, move first element to end
C   of element part of list and add new element to front of list.
  250   IF (ID.GE.LIMIT) GO TO 295
      IW(NP) = IW(NP0)
      IW(NP0) = IW(KP1)
      IW(KP1) = ME
C Store the new length of the list.
      IW(KP1-1) = NP - KP1 + 1
C
C Check whether row is is identical to another by looking in linked
C   list of supervariables with degree ID at those whose lists have
C   first entry ME. Note that those containing ME come first so the
C   search can be terminated when a list not starting with ME is
C   found.
      JS = IPD(ID)
      DO 280 L = 1,N
        IF (JS.LE.0) GO TO 300
        KP1 = IPE(JS) + 1
        IF (IW(KP1).NE.ME) GO TO 300
C JS has same degree and is active. Check if identical to IS.
        KP2 = KP1 - 1 + IW(KP1-1)
        DO 260 KP = KP1,KP2
          IE = IW(KP)
C Jump if IE is a supervariable or an element not in the list of IS.
          IF (ABS(FLAG(IE)+0).GT.NFLG) GO TO 270
  260     CONTINUE
        GO TO 290

  270     JS = NXT(JS)
  280   CONTINUE
C Supervariable amalgamation. Row IS is identical to row JS.
C Regard all variables in the two supervariables as being in IS. Set
C   tree pointer, FLAG and NV entries.
  290   IPE(JS) = -IS
      NV(IS) = NV(IS) + NV(JS)
      NV(JS) = 0
      FLAG(JS) = -1
C Replace JS by IS in linked list.
      NS = NXT(JS)
      LS = LST(JS)
      IF (NS.GT.0) LST(NS) = IS
      IF (LS.GT.0) NXT(LS) = IS
      LST(IS) = LS
      NXT(IS) = NS
      LST(JS) = 0
      NXT(JS) = 0
      IF (IPD(ID).EQ.JS) IPD(ID) = IS
      GO TO 310
C Treat IS as full. Merge it into the root node.
  295   IF (NVROOT.EQ.0) THEN
        ROOT = IS
        IPE(IS) = 0
      ELSE
        IW(K) = ROOT
        IPE(IS) = -ROOT
        NV(ROOT) = NV(ROOT) + NV(IS)
        NV(IS) = 0
        FLAG(IS) = -1
      END IF
      NVROOT = NV(ROOT)
      GO TO 310
C Insert IS into linked list of supervariables of same degree.
  300   NS = IPD(ID)
      IF (NS.GT.0) LST(NS) = IS
      NXT(IS) = NS
      IPD(ID) = IS
      LST(IS) = -ID
      MD = MIN(MD,ID)
  310   CONTINUE

C
C Reset flags for supervariables in newly created element and
C   remove those absorbed into others.
      DO 320 K = K1,K2
      IS = IW(K)
      IF (NV(IS).EQ.0) GO TO 320
      FLAG(IS) = NFLG
      IW(IP) = IS
      IP = IP + 1
  320   CONTINUE

      FLAG(ME) = -NFLG
C Move first entry to end to make room for length.
      IW(IP) = IW(K1)
      IW(K1) = IP - K1
C Set pointer for new element and reset IWFR.
      IPE(ME) = K1
      IWFR = IP + 1

C  End of main loop
  340 CONTINUE
C

C Absorb any remaining variables into the root
  350 DO 360 IS = 1,N
      IF(NXT(IS).NE.0 .OR. LST(IS).NE.0) THEN
      IF (NVROOT.EQ.0) THEN
        ROOT = IS
        IPE(IS) = 0
      ELSE
        IPE(IS) = -ROOT
      END IF
      NVROOT = NVROOT + NV(IS)
      NV(IS) = 0
       END IF
  360 CONTINUE
C Link any remaining elements to the root
    DO 370 IE = 1,N
      IF (IPE(IE).GT.0) IPE(IE) = -ROOT
  370 CONTINUE
    IF(NVROOT.GT.0)NV(ROOT)=NVROOT
    END
    SUBROUTINE ME57ZD(N,IPE,IW,LW,IWFR,NCMPA)
C Is identical to subroutine ME27UD.  Internal version for ME57.
C COMPRESS LISTS HELD BY ME27H/HD (ME57H/HD) IN IW AND ADJUST POINTERS
C   IN IPE TO CORRESPOND.
C N IS THE MATRIX ORDER. IT IS NOT ALTERED.
C IPE(I) POINTS TO THE POSITION IN IW OF THE START OF LIST I OR IS
C   ZERO IF THERE IS NO LIST I. ON EXIT IT POINTS TO THE NEW POSITION.
C IW HOLDS THE LISTS, EACH HEADED BY ITS LENGTH. ON OUTPUT THE SAME
C   LISTS ARE HELD, BUT THEY ARE NOW COMPRESSED TOGETHER.
C LW HOLDS THE LENGTH OF IW. IT IS NOT ALTERED.
C IWFR NEED NOT BE SET ON ENTRY. ON EXIT IT POINTS TO THE FIRST FREE
C   LOCATION IN IW.
C   ON RETURN IT IS SET TO THE FIRST FREE LOCATION IN IW.
C NCMPA is number of compresses.
C
C   .. Scalar Arguments ..
    INTEGER IWFR,LW,N,NCMPA
C   ..
C   .. Array Arguments ..
    INTEGER IPE(N),IW(LW)
C   ..
C   .. Local Scalars ..
    INTEGER I,IR,K,K1,K2,LWFR
C   ..
C   .. Executable Statements ..
    NCMPA = NCMPA + 1
C PREPARE FOR COMPRESSING BY STORING THE LENGTHS OF THE
C   LISTS IN IPE AND SETTING THE FIRST ENTRY OF EACH LIST TO
C   -(LIST NUMBER).
    DO 10 I = 1,N
      K1 = IPE(I)
      IF (K1.LE.0) GO TO 10
      IPE(I) = IW(K1)
      IW(K1) = -I
   10 CONTINUE
C
C COMPRESS
C IWFR POINTS JUST BEYOND THE END OF THE COMPRESSED FILE.
C LWFR POINTS JUST BEYOND THE END OF THE UNCOMPRESSED FILE.
    IWFR = 1
    LWFR = IWFR
    DO 60 IR = 1,N
      IF (LWFR.GT.LW) GO TO 70
C SEARCH FOR THE NEXT NEGATIVE ENTRY.
      DO 20 K = LWFR,LW
      IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
      GO TO 70
C PICK UP ENTRY NUMBER, STORE LENGTH IN NEW POSITION, SET NEW POINTER
C   AND PREPARE TO COPY LIST.
   30   I = -IW(K)
      IW(IWFR) = IPE(I)
      IPE(I) = IWFR
      K1 = K + 1
      K2 = K + IW(IWFR)
      IWFR = IWFR + 1
      IF (K1.GT.K2) GO TO 50
C COPY LIST TO NEW POSITION.
      DO 40 K = K1,K2
      IW(IWFR) = IW(K)
      IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN

    END
* COPYRIGHT (c) 1993 Council for the Central Laboratory
*            of the Research Councils

C Original date 29 Jan 2001
C 29 January 2001. Modified from MC49 to be threadsafe.

C 12th July 2004 Version 1.0.0. Version numbering added.
C 28 February 2008. Version 1.0.1. Comments flowed to column 72.
C 21 September 2009. Version 1.0.2. Minor change to documentation.

    SUBROUTINE MC59A(ICNTL,NC,NR,NE,IRN,LJCN,JCN,LA,A,LIP,IP,
   &            LIW,IW,INFO)
C
C To sort the sparsity pattern of a matrix to an ordering by columns.
C There is an option for ordering the entries within each column by
C increasing row indices and an option for checking the user-supplied
C matrix entries for indices which are out-of-range or duplicated.
C
C ICNTL:  INTEGER array of length 10. Intent(IN). Used to specify
C       control parameters for the subroutine.
C ICNTL(1): indicates whether the user-supplied matrix entries are to
C       be checked for duplicates, and out-of-range indices.
C       Note  simple checks are always performed.
C       ICNTL(1) = 0, data checking performed.
C       Otherwise, no data checking.
C ICNTL(2): indicates the ordering requested.
C       ICNTL(2) = 0, input is by rows and columns in arbitrary
C       order and the output is sorted by columns.
C       ICNTL(2) = 1, the output is also row ordered
C       within each column.
C       ICNTL(2) = 2, the input is already ordered by
C       columns and is to be row ordered within each column.
C       Values outside the range 0 to 2 are flagged as an error.
C ICNTL(3): indicates whether matrix entries are also being ordered.
C       ICNTL(3) = 0, matrix entries are ordered.
C       Otherwise, only the sparsity pattern is ordered
C       and the array A is not accessed by the routine.
C ICNTL(4): the unit number of the device to
C       which error messages are sent. Error messages
C       can be suppressed by setting ICNTL(4) < 0.
C ICNTL(5): the unit number of the device to
C       which warning messages are sent. Warning
C       messages can be suppressed by setting ICNTL(5) < 0.
C ICNTL(6)  indicates whether matrix symmetric. If unsymmetric, ICNTL(6)
C       must be set to 0.
C       If ICNTL(6) = -1 or 1, symmetric and only the lower
C       triangular part of the reordered matrix is returned.
C       If ICNTL(6) = -2 or 2, Hermitian and only the lower
C       triangular part of the reordered matrix is returned.
C       If error checks are performed (ICNTL(1) = 0)
C       and ICNTL(6)> 1 or 2, the values of duplicate
C       entries are added together; if ICNTL(6) < -1 or -2, the
C       value of the first occurrence of the entry is used.
C ICNTL(7) to ICNTL(10) are not currently accessed by the routine.
C
C NC:    INTEGER variable. Intent(IN). Must be set by the user
C      to the number of columns in the matrix.
C NR:    INTEGER variable. Intent(IN). Must be set by the user
C      to the number of rows in the matrix.
C NE:    INTEGER variable. Intent(IN). Must be set by the user
C      to the number of entries in the matrix.
C IRN: INTEGER array of length NE. Intent (INOUT). Must be set by the
C        user to hold the row indices of the entries in the matrix.
C      If ICNTL(2).NE.2, the entries may be in any order.
C      If ICNTL(2).EQ.2, the entries in column J must be in
C        positions IP(J) to IP(J+1)-1 of IRN. On exit, the row
C        indices are reordered so that the entries of a single
C        column are contiguous with column J preceding column J+1, J
C        = 1, 2, ..., NC-1, with no space between columns.
C      If ICNTL(2).EQ.0, the order within each column is arbitrary;
C        if ICNTL(2) = 1 or 2, the order within each column is by
C        increasing row indices.
C LJCN:    INTEGER variable. Intent(IN). Defines length array
C JCN:   INTEGER array of length LJCN. Intent (INOUT).
C      If ICNTL(2) = 0 or 1, JCN(K) must be set by the user
C      to the column index of the entry
C      whose row index is held in IRN(K), K = 1, 2, ..., NE.
C      On exit, the contents of this array  will have been altered.
C      If ICNTL(2) = 2, the array is not accessed.
C LA:    INTEGER variable. Intent(IN). Defines length of array
C      A.
C A:     is a REAL (DOUBLE PRECISION in the D version, INTEGER in
C      the I version, COMPLEX in the C version,
C      or COMPLEX"*"16 in the Z version) array of length LA.
C      Intent(INOUT).
C      If ICNTL(3).EQ.0, A(K) must be set by the user to
C      hold the value of the entry with row index IRN(K),
C      K = 1, 2, ..., NE. On exit, the array will have been
C      permuted in the same way as the array IRN.
C      If ICNTL(3).NE.0, the array is not accessed.
C LIP:   INTEGER variable. Intent(IN). Defines length of array
C      IP.
C IP:    INTEGER array of length LIP. Intent(INOUT). IP
C      need only be set by the user if ICNTL(2) = 2.
C      In this case, IP(J) holds the position in
C      the array IRN of the first entry in column J, J = 1, 2,
C      ..., NC, and IP(NC+1) is one greater than the number of
C      entries in the matrix.
C      In all cases, the array IP will have this meaning on exit
C      from the subroutine and is altered when ICNTL(2) = 2 only
C      when ICNTL(1) =  0 and there are out-of-range
C      indices or duplicates.
C LIW:   INTEGER variable. Intent(IN). Defines length of array
C      IW.
C IW:    INTEGER array of length LIW. Intent(OUT). Used by the
C      routine as workspace.
C INFO:    INTEGER array of length 10.  Intent(OUT). On exit,
C      a negative value of INFO(1) is used to signal a fatal
C      error in the input data, a positive value of INFO(1)
C      indicates that a warning has been issued, and a
C      zero value is used to indicate a successful call.
C      In cases of error, further information is held in INFO(2).
C      For warnings, further information is
C      provided in INFO(3) to INFO(6).  INFO(7) to INFO(10) are not
C      currently used and are set to zero.
C      Possible nonzero values of INFO(1):
C       -1 -  The restriction ICNTL(2) = 0, 1, or 2 violated.
C         Value of ICNTL(2) is given by INFO(2).
C       -2 -  NC.LE.0. Value of NC is given by INFO(2).
C       -3 -  Error in NR. Value of NR is given by INFO(2).
C       -4 -  NE.LE.0. Value of NE is given by INFO(2).
C       -5 -  LJCN too small. Min. value of LJCN is given by INFO(2).
C       -6 -  LA too small. Min. value of LA is given by INFO(2).
C       -7 -  LIW too small. Value of LIW is given by INFO(2).
C       -8 -  LIP too small. Value of LIP is given by INFO(2).
C       -9 -  The entries of IP not monotonic increasing.
C      -10 -  For each I, IRN(I) or JCN(I) out-of-range.
C      -11 -  ICNTL(6) is out of range.
C       +1 -  One or more duplicated entries. One copy of
C         each such entry is kept and, if ICNTL(3) = 0 and
C         ICNTL(6).GE.0, the values of these entries are
C         added together. If  ICNTL(3) = 0 and ICNTL(6).LT.0,
C         the value of the first occurrence of the entry is used.
C         Initially INFO(3) is set to zero. If an entry appears
C         k times, INFO(3) is incremented by k-1 and INFO(6)
C         is set to the revised number of entries in the
C         matrix.
C       +2 - One or more of the entries in IRN out-of-range. These
C         entries are removed by the routine.`INFO(4) is set to
C         the number of entries which were out-of-range and
C         INFO(6) is set to the revised number of entries in the
C         matrix.
C       +4 - One or more of the entries in JCN out-of-range. These
C         entries are removed by the routine. INFO(5) is set to
C         the number of entries which were out-of-range and
C         INFO(6) is set to the revised number of entries in the
C         matrix. Positive values of INFO(1) are summed so that
C         the user can identify all warnings.
C
C   .. Scalar Arguments ..
    INTEGER LA,LIP,LIW,LJCN,NC,NE,NR
C   ..
C   .. Array Arguments ..
    REAL A(LA)
    INTEGER ICNTL(10),IP(LIP),INFO(10),IRN(NE),IW(LIW),JCN(LJCN)
C   ..
C   .. Local Scalars ..
    INTEGER I,ICNTL1,ICNTL2,ICNTL3,ICNTL6,LAA
    INTEGER IDUP,IOUT,IUP,JOUT,LP,MP,KNE,PART
    LOGICAL LCHECK
C   ..
C   .. External Subroutines ..
    EXTERNAL MC59B,MC59C,MC59D,MC59E,MC59F
C   ..
C   .. Intrinsic Functions ..
    INTRINSIC MAX
C   ..
C   .. Executable Statements ..

C Initialise
    DO 10 I = 1,10
       INFO(I) = 0
   10 CONTINUE

    ICNTL1 = ICNTL(1)
    ICNTL2 = ICNTL(2)
    ICNTL3 = ICNTL(3)
    ICNTL6 = ICNTL(6)
    LCHECK = (ICNTL1.EQ.0)
C Streams for errors/warnings
    LP = ICNTL(4)
    MP = ICNTL(5)

C  Check the input data
    IF (ICNTL2.GT.2 .OR. ICNTL2.LT.0) THEN
       INFO(1) = -1
       INFO(2) = ICNTL2
       IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9010) ICNTL2
       END IF
       GO TO 70
    END IF

    IF (ICNTL6.GT.2 .OR. ICNTL6.LT.-2) THEN
       INFO(1) = -11
       INFO(2) = ICNTL6
       IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9150) ICNTL6
       END IF
       GO TO 70
    END IF
C For real matrices, symmetric = Hermitian so only
C have to distinguish between unsymmetric (ICNTL6 = 0) and
C symmetric (ICNTL6.ne.0)

    IF (NC.LT.1) THEN
      INFO(1) = -2
      INFO(2) = NC
      IF (LP.GT.0) THEN
      WRITE (LP,FMT=9000) INFO(1)
      WRITE (LP,FMT=9020) NC
      END IF
      GO TO 70
    END IF

    IF (NR.LT.1) THEN
      INFO(1) = -3
      INFO(2) = NR
      IF (LP.GT.0) THEN
      WRITE (LP,FMT=9000) INFO(1)
      WRITE (LP,FMT=9030) NR
      END IF
      GO TO 70
    END IF

    IF (ICNTL6.NE.0 .AND. NR.NE.NC) THEN
      INFO(1) = -3
      INFO(2) = NR
      IF (LP.GT.0) THEN
      WRITE (LP,FMT=9000) INFO(1)
      WRITE (LP,FMT=9035) NC,NR
      END IF
      GO TO 70
    END IF

    IF (NE.LT.1) THEN
      INFO(1) = -4
      INFO(2) = NE
      IF (LP.GT.0) THEN
      WRITE (LP,FMT=9000) INFO(1)
      WRITE (LP,FMT=9040) NE
      END IF
      GO TO 70
    END IF

    IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.1) THEN
      IF (LJCN.LT.NE) THEN
      INFO(1) = -5
      INFO(2) = NE
      END IF
    ELSE
      IF (LJCN.LT.1) THEN
      INFO(1) = -5
      INFO(2) = 1
      END IF
    END IF
    IF (INFO(1).EQ.-5) THEN
       IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9050) LJCN,INFO(2)
       END IF
       GO TO 70
    END IF

    IF (ICNTL3.EQ.0) THEN
      IF (LA.LT.NE) THEN
      INFO(1) = -6
      INFO(2) = NE
      END IF
    ELSE
      IF (LA.LT.1) THEN
      INFO(1) = -6
      INFO(2) = 1
      END IF
    END IF
    IF (INFO(1).EQ.-6) THEN
       IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9060) LA,INFO(2)
       END IF
       GO TO 70
    END IF

    IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.2) THEN
      IF (LIP.LT.NC+1) THEN
      INFO(1) = -7
      INFO(2) = NC+1
      END IF
    ELSE IF (LIP.LT.MAX(NR,NC)+1) THEN
      INFO(1) = -7
      INFO(2) = MAX(NR,NC)+1
    END IF
    IF (INFO(1).EQ.-7) THEN
      IF (LP.GT.0) THEN
      WRITE (LP,FMT=9000) INFO(1)
      WRITE (LP,FMT=9065) LIP,INFO(2)
      END IF
      GO TO 70
    END IF

C Check workspace sufficient
    IF (LIW.LT.MAX(NR,NC)+1) THEN
      INFO(1) = -8
      INFO(2) = MAX(NR,NC)+1
      IF (LP.GT.0) THEN
      WRITE (LP,FMT=9000) INFO(1)
      WRITE (LP,FMT=9070) LIW,INFO(2)
      END IF
      GO TO 70
    END IF

    LAA = NE
    IF (ICNTL3.NE.0) LAA = 1
C Initialise counts of number of out-of-range entries and duplicates
    IOUT = 0
    JOUT = 0
    IDUP = 0
    IUP = 0

C PART is used by MC59B to indicate if upper or lower or
C all of matrix is required.
C PART =  0 : unsymmetric case, whole matrix wanted
C PART =  1 : symmetric case, lower triangular part of matrix wanted
C PART = -1 : symmetric case, upper triangular part of matrix wanted
    PART = 0
    IF (ICNTL6.NE.0) PART = 1

    IF (ICNTL2.EQ.0) THEN

C Order directly by columns
C On exit from MC59B, KNE holds number of entries in matrix
C after removal of out-of-range entries. If no data checking, KNE = NE.
      CALL MC59B(LCHECK,PART,NC,NR,NE,IRN,JCN,LAA,A,IP,IW,
   +          IOUT,JOUT,KNE)
C Return if ALL entries out-of-range.
      IF (KNE.EQ.0) GO TO 50

C Check for duplicates
      IF (LCHECK) CALL MC59E(NC,NR,NE,IRN,LIP,IP,LAA,A,IW,IDUP,
   &                KNE,ICNTL6)

    ELSE IF (ICNTL2.EQ.1) THEN

C First order by rows.
C Interchanged roles of IRN and JCN, so set PART = -1
C if matrix is symmetric case
      IF (ICNTL6.NE.0) PART = -1
      CALL MC59B(LCHECK,PART,NR,NC,NE,JCN,IRN,LAA,A,IW,IP,
   +          JOUT,IOUT,KNE)
C Return if ALL entries out-of-range.
      IF (KNE.EQ.0) GO TO 50

C At this point, JCN and IW hold column indices and row pointers
C Optionally, check for duplicates.
      IF (LCHECK) CALL MC59E(NR,NC,NE,JCN,NR+1,IW,LAA,A,IP,
   &                IDUP,KNE,ICNTL6)

C Now order by columns and by rows within each column
      CALL MC59C(NC,NR,KNE,IRN,JCN,LAA,A,IP,IW)

    ELSE IF (ICNTL2.EQ.2) THEN
C Input is using IP, IRN.
C Optionally check for duplicates and remove out-of-range entries
      IF (LCHECK) THEN
      CALL MC59F(NC,NR,NE,IRN,NC+1,IP,LAA,A,LIW,IW,IDUP,
   +          IOUT,IUP,KNE,ICNTL6,INFO)
C Return if IP not monotonic.
      IF (INFO(1).EQ.-9) GO TO 40
C Return if ALL entries out-of-range.
      IF (KNE.EQ.0) GO TO 50
      ELSE
       KNE = NE
      END IF

C  Order by rows within each column
      CALL MC59D(NC,KNE,IRN,IP,LAA,A)

    END IF

    INFO(3) = IDUP
    INFO(4) = IOUT
    INFO(5) = JOUT
    INFO(6) = KNE
    INFO(7) = IUP

C Set warning flag if out-of-range /duplicates found
    IF (IDUP.GT.0) INFO(1) = INFO(1) + 1
    IF (IOUT.GT.0) INFO(1) = INFO(1) + 2
    IF (JOUT.GT.0) INFO(1) = INFO(1) + 4
    IF (INFO(1).GT.0 .AND. MP.GT.0) THEN
      WRITE (MP,FMT=9080) INFO(1)
      IF (IOUT.GT.0) WRITE (MP,FMT=9090) IOUT
      IF (JOUT.GT.0) WRITE (MP,FMT=9110) JOUT
      IF (IDUP.GT.0) WRITE (MP,FMT=9100) IDUP
      IF (IUP.GT.0)  WRITE (MP,FMT=9130) IUP
    END IF
    GO TO 70

   40 INFO(3) = IDUP
    INFO(4) = IOUT
    INFO(7) = IUP
    IF (LP.GT.0) THEN
      WRITE (LP,FMT=9000) INFO(1)
      WRITE (LP,FMT=9140)
    END IF
    GO TO 70

   50 INFO(1) = -10
    INFO(4) = IOUT
    INFO(5) = JOUT
    INFO(2) = IOUT + JOUT
    IF (LP.GT.0) THEN
      WRITE (LP,FMT=9000) INFO(1)
      WRITE (LP,FMT=9120)
    END IF
   70 RETURN

 9000 FORMAT (/,' *** Error return from MC59A *** INFO(1) = ',I3)
 9010 FORMAT (1X,'ICNTL(2) = ',I2,' is out of range')
 9020 FORMAT (1X,'NC = ',I6,' is out of range')
 9030 FORMAT (1X,'NR = ',I6,' is out of range')
 9035 FORMAT (1X,'Symmetric case. NC = ',I6,' but NR = ',I6)
 9040 FORMAT (1X,'NE = ',I10,' is out of range')
 9050 FORMAT (1X,'Increase LJCN from ',I10,' to at least ',I10)
 9060 FORMAT (1X,'Increase LA from ',I10,' to at least ',I10)
 9065 FORMAT (1X,'Increase LIP from ',I8,' to at least ',I10)
 9070 FORMAT (1X,'Increase LIW from ',I8,' to at least ',I10)
 9080 FORMAT (/,' *** Warning message from MC59A *** INFO(1) = ',I3)
 9090 FORMAT (1X,I8,' entries in IRN supplied by the user were ',
   +     /,'     out of range and were ignored by the routine')
 9100 FORMAT (1X,I8,' duplicate entries were supplied by the user')
 9110 FORMAT (1X,I8,' entries in JCN supplied by the user were ',
   +     /,'     out of range and were ignored by the routine')
 9120 FORMAT (1X,'All entries out of range')
 9130 FORMAT (1X,I8,' of these entries were in the upper triangular ',
   +     /,'     part of matrix')
 9140 FORMAT (1X,'Entries in IP are not monotonic increasing')
 9150 FORMAT (1X,'ICNTL(6) = ',I2,' is out of range')
    END
C***********************************************************************
    SUBROUTINE MC59B(LCHECK,PART,NC,NR,NE,IRN,JCN,LA,A,IP,IW,IOUT,
   +            JOUT,KNE)
C
C   To sort a sparse matrix from arbitrary order to
C   column order, unordered within each column. Optionally
C   checks for out-of-range entries in IRN,JCN.
C
C LCHECK - logical variable. Intent(IN). If true, check
C      for out-of-range indices.
C PART -   integer variable. Intent(IN)
C PART =  0 : unsymmetric case, whole matrix wanted
C PART =  1 : symmetric case, lower triangular part of matrix wanted
C         (ie IRN(K) .ge. JCN(K) on exit)
C PART = -1 : symmetric case, upper triangular part of matrix wanted
C         (ie IRN(K) .le. JCN(K) on exit)
C   NC - integer variable. Intent(IN)
C    - on entry must be set to the number of columns in the matrix
C   NR - integer variable. Intent(IN)
C    - on entry must be set to the number of rows in the matrix
C   NE - integer variable. Intent(IN)
C    - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(INOUT)
C    - on entry set to contain the row indices of the nonzeros
C      in arbitrary order.
C    - on exit, the entries in IRN are reordered so that the row
C      indices for column 1 precede those for column 2 and so on,
C      but the order within columns is arbitrary.
C  JCN - integer array of length NE. Intent(INOUT)
C    - on entry set to contain the column indices of the nonzeros
C    - JCN(K) must be the column index of
C      the entry in IRN(K)
C    - on exit, JCN(K) is the column index for the entry with
C      row index IRN(K) (K=1,...,NE).
C  LA  - integer variable which defines the length of the array A.
C      Intent(IN)
C   A  - real (double precision/complex/complex*16) array of length LA
C      Intent(INOUT)
C    - if LA > 1, the array must be of length NE, and A(K)
C      must be set to the value of the entry in (IRN(K), JCN(K));
C      on exit A is reordered in the same way as IRN
C    - if LA = 1, the array is not accessed
C  IP  - integer array of length NC+1. Intent(INOUT)
C    - not set on entry
C    - on exit, IP(J) contains the position in IRN (and A) of the
C      first entry in column J (J=1,...,NC)
C    - IP(NC+1) is set to NE+1
C  IW  - integer array of length NC+1.  Intent(INOUT)
C    - the array is used as workspace
C    - on exit IW(I) = IP(I) (so IW(I) points to the beginning
C      of column I).
C IOUT - integer variable. Intent(OUT). On exit, holds number
C      of entries in IRN found to be out-of-range
C JOUT - integer variable. Intent(OUT). On exit, holds number
C      of entries in JCN found to be out-of-range
C  KNE - integer variable. Intent(OUT). On exit, holds number
C      of entries in matrix after removal of out-of-range entries.
C      If no data checking, KNE = NE.

C   .. Scalar Arguments ..
    INTEGER LA,NC,NE,NR,IOUT,JOUT,KNE,PART
    LOGICAL LCHECK
C   ..
C   .. Array Arguments ..
    REAL A(LA)
    INTEGER IP(NC+1),IRN(NE),IW(NC+1),JCN(NE)
C   ..
C   .. Local Scalars ..
    REAL ACE,ACEP
    INTEGER I,ICE,ICEP,J,JCE,JCEP,K,L,LOC
C   ..
C   .. Executable Statements ..

C Initialise IW
    DO 10 J = 1,NC + 1
      IW(J) = 0
   10 CONTINUE

    KNE = 0
    IOUT = 0
    JOUT = 0
C Count the number of entries in each column and store in IW.
C We also allow checks for out-of-range indices
    IF (LCHECK) THEN
C Check data.
C Treat case of pattern only separately.
      IF (LA.GT.1) THEN
      IF (PART.EQ.0) THEN
C Unsymmetric
        DO 20 K = 1,NE
          I = IRN(K)
          J = JCN(K)
          IF (I.GT.NR .OR. I.LT.1) THEN
          IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
          IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
          ELSE IF (J.GT.NC .OR. J.LT.1) THEN
          JOUT = JOUT + 1
          ELSE
          KNE = KNE + 1
          IRN(KNE) = I
          JCN(KNE) = J
          A(KNE) = A(K)
          IW(J) = IW(J) + 1
          END IF
   20     CONTINUE
      ELSE IF (PART.EQ.1) THEN
C Symmetric, lower triangle
        DO 21 K = 1,NE
          I = IRN(K)
          J = JCN(K)
          IF (I.GT.NR .OR. I.LT.1) THEN
          IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
          IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
          ELSE IF (J.GT.NC .OR. J.LT.1) THEN
          JOUT = JOUT + 1
          ELSE
          KNE = KNE + 1
C Lower triangle ... swap if necessary
          IF (I.LT.J) THEN
            IRN(KNE) = J
            JCN(KNE) = I
            IW(I) = IW(I) + 1
          ELSE
            IRN(KNE) = I
            JCN(KNE) = J
            IW(J) = IW(J) + 1
          END IF
          A(KNE) = A(K)
          END IF
   21     CONTINUE
      ELSE IF (PART.EQ.-1) THEN
C Symmetric, upper triangle
        DO 22 K = 1,NE
          I = IRN(K)
          J = JCN(K)
          IF (I.GT.NR .OR. I.LT.1) THEN
          IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
          IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
          ELSE IF (J.GT.NC .OR. J.LT.1) THEN
          JOUT = JOUT + 1
          ELSE
          KNE = KNE + 1
C Upper triangle ... swap if necessary
          IF (I.GT.J) THEN
            IRN(KNE) = J
            JCN(KNE) = I
            IW(I) = IW(I) + 1
          ELSE
            IRN(KNE) = I
            JCN(KNE) = J
            IW(J) = IW(J) + 1
          END IF
          A(KNE) = A(K)
          END IF
   22     CONTINUE
      END IF
      ELSE
C Pattern only
      IF (PART.EQ.0) THEN
        DO 25 K = 1,NE
          I = IRN(K)
          J = JCN(K)
          IF (I.GT.NR .OR. I.LT.1) THEN
          IOUT = IOUT + 1
          IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
          ELSE IF (J.GT.NC .OR. J.LT.1) THEN
          JOUT = JOUT + 1
          ELSE
          KNE = KNE + 1
          IRN(KNE) = I
          JCN(KNE) = J
          IW(J) = IW(J) + 1
          END IF
   25     CONTINUE
      ELSE IF (PART.EQ.1) THEN
        DO 26 K = 1,NE
          I = IRN(K)
          J = JCN(K)
          IF (I.GT.NR .OR. I.LT.1) THEN
          IOUT = IOUT + 1
          IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
          ELSE IF (J.GT.NC .OR. J.LT.1) THEN
          JOUT = JOUT + 1
          ELSE
          KNE = KNE + 1
C Lower triangle ... swap if necessary
          IF (I.LT.J) THEN
            IRN(KNE) = J
            JCN(KNE) = I
            IW(I) = IW(I) + 1
          ELSE
            IRN(KNE) = I
            JCN(KNE) = J
            IW(J) = IW(J) + 1
          END IF
          END IF
   26     CONTINUE
      ELSE IF (PART.EQ.-1) THEN
        DO 27 K = 1,NE
          I = IRN(K)
          J = JCN(K)
          IF (I.GT.NR .OR. I.LT.1) THEN
          IOUT = IOUT + 1
          IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
          ELSE IF (J.GT.NC .OR. J.LT.1) THEN
          JOUT = JOUT + 1
          ELSE
          KNE = KNE + 1
C Upper triangle ... swap if necessary
          IF (I.GT.J) THEN
            IRN(KNE) = J
            JCN(KNE) = I
            IW(I) = IW(I) + 1
          ELSE
            IRN(KNE) = I
            JCN(KNE) = J
            IW(J) = IW(J) + 1
          END IF
          END IF
   27     CONTINUE
      END IF
      END IF
C Return if ALL entries out-of-range.
      IF (KNE.EQ.0) GO TO 130

    ELSE

C No checks
      KNE = NE
      IF (PART.EQ.0) THEN
      DO 30 K = 1,NE
        J = JCN(K)
        IW(J) = IW(J) + 1
   30   CONTINUE
      ELSE IF (PART.EQ.1) THEN
      DO 35 K = 1,NE
        I = IRN(K)
        J = JCN(K)
C Lower triangle ... swap if necessary
        IF (I.LT.J) THEN
         IRN(K) = J
         JCN(K) = I
         IW(I) = IW(I) + 1
        ELSE
          IW(J) = IW(J) + 1
        END IF
   35   CONTINUE
      ELSE IF (PART.EQ.-1) THEN
      DO 36 K = 1,NE
        I = IRN(K)
        J = JCN(K)
C Upper triangle ... swap if necessary
        IF (I.GT.J) THEN
         IRN(K) = J
         JCN(K) = I
         IW(I) = IW(I) + 1
        ELSE
          IW(J) = IW(J) + 1
        END IF
   36   CONTINUE
      END IF
    END IF

C KNE is now the number of nonzero entries in matrix.

C Put into IP and IW the positions where each column
C would begin in a compressed collection with the columns
C in natural order.

    IP(1) = 1
    DO 37 J = 2,NC + 1
      IP(J) = IW(J-1) + IP(J-1)
      IW(J-1) = IP(J-1)
   37 CONTINUE

C Reorder the elements into column order.
C Fill in each column from the front, and as a new entry is placed
C in column K increase the pointer IW(K) by one.

    IF (LA.EQ.1) THEN
C Pattern only
      DO 70 L = 1,NC
      DO 60 K = IW(L),IP(L+1) - 1
        ICE = IRN(K)
        JCE = JCN(K)
        DO 40 J = 1,NE
          IF (JCE.EQ.L) GO TO 50
          LOC = IW(JCE)
          JCEP = JCN(LOC)
          ICEP = IRN(LOC)
          IW(JCE) = LOC + 1
          JCN(LOC) = JCE
          IRN(LOC) = ICE
          JCE = JCEP
          ICE = ICEP
   40     CONTINUE
   50     JCN(K) = JCE
        IRN(K) = ICE
   60   CONTINUE
   70   CONTINUE
    ELSE

      DO 120 L = 1,NC
      DO 110 K = IW(L),IP(L+1) - 1
        ICE = IRN(K)
        JCE = JCN(K)
        ACE = A(K)
        DO 90 J = 1,NE
          IF (JCE.EQ.L) GO TO 100
          LOC = IW(JCE)
          JCEP = JCN(LOC)
          ICEP = IRN(LOC)
          IW(JCE) = LOC + 1
          JCN(LOC) = JCE
          IRN(LOC) = ICE
          JCE = JCEP
          ICE = ICEP
          ACEP = A(LOC)
          A(LOC) = ACE
          ACE = ACEP
   90     CONTINUE
  100     JCN(K) = JCE
        IRN(K) = ICE
        A(K) = ACE
  110   CONTINUE
  120   CONTINUE
    END IF

  130 CONTINUE

    RETURN
    END
C
C**********************************************************
    SUBROUTINE MC59C(NC,NR,NE,IRN,JCN,LA,A,IP,IW)
C
C   To sort a sparse matrix stored by rows,
C   unordered within each row, to ordering by columns, with
C   ordering by rows within each column.
C
C   NC - integer variable. Intent(IN)
C    - on entry must be set to the number of columns in the matrix
C   NR - integer variable. Intent(IN)
C    - on entry must be set to the number of rows in the matrix
C  NE - integer variable. Intent(IN)
C    - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(OUT).
C    - not set on entry.
C    - on exit,  IRN holds row indices with the row
C      indices for column 1 preceding those for column 2 and so on,
C      with ordering by rows within each column.
C  JCN - integer array of length NE. Intent(INOUT)
C    - on entry set to contain the column indices of the nonzeros
C      with indices for column 1 preceding those for column 2
C      and so on, with the order within columns is arbitrary.
C    - on exit, contents destroyed.
C  LA  - integer variable which defines the length of the array A.
C      Intent(IN)
C   A  - real (double precision/complex/complex*16) array of length LA
C      Intent(INOUT)
C    - if LA > 1, the array must be of length NE, and A(K)
C      must be set to the value of the entry in JCN(K);
C      on exit A, A(K) holds the value of the entry in IRN(K).
C    - if LA = 1, the array is not accessed
C  IP  - integer array of length NC+1. Intent(INOUT)
C    - not set on entry
C    - on exit, IP(J) contains the position in IRN (and A) of the
C      first entry in column J (J=1,...,NC)
C    - IP(NC+1) is set to NE+1
C  IW  - integer array of length NR+1.  Intent(IN)
C    - on entry, must be set on entry so that IW(J) points to the
C      position in JCN of the first entry in row J, J=1,...,NR, and
C      IW(NR+1) must be set to NE+1
C
C   .. Scalar Arguments ..
    INTEGER LA,NC,NE,NR
C   ..
C   .. Array Arguments ..
    REAL A(LA)
    INTEGER IP(NC+1),IRN(NE),IW(NR+1),JCN(NE)
C   ..
C   .. Local Scalars ..
    REAL ACE,ACEP
    INTEGER I,ICE,ICEP,J,J1,J2,K,L,LOC,LOCP
C   ..
C   .. Executable Statements ..

C  Count the number of entries in each column

    DO 10 J = 1,NC
      IP(J) = 0
   10 CONTINUE

    IF (LA.GT.1) THEN

      DO 20 K = 1,NE
      I = JCN(K)
      IP(I) = IP(I) + 1
      IRN(K) = JCN(K)
   20   CONTINUE
      IP(NC+1) = NE + 1

C  Set IP so that IP(I) points to the first entry in column I+1

      IP(1) = IP(1) + 1
      DO 30 J = 2,NC
      IP(J) = IP(J) + IP(J-1)
   30   CONTINUE

      DO 50 I = NR,1,-1
      J1 = IW(I)
      J2 = IW(I+1) - 1
      DO 40 J = J1,J2
        K = IRN(J)
        L = IP(K) - 1
        JCN(J) = L
        IRN(J) = I
        IP(K) = L
   40   CONTINUE
   50   CONTINUE
      IP(NC+1) = NE + 1
      DO 70 J = 1,NE
      LOC = JCN(J)
      IF (LOC.EQ.0) GO TO 70
      ICE = IRN(J)
      ACE = A(J)
      JCN(J) = 0
      DO 60 K = 1,NE
        LOCP = JCN(LOC)
        ICEP = IRN(LOC)
        ACEP = A(LOC)
        JCN(LOC) = 0
        IRN(LOC) = ICE
        A(LOC) = ACE
        IF (LOCP.EQ.0) GO TO 70
        ICE = ICEP
        ACE = ACEP
        LOC = LOCP
   60   CONTINUE
   70   CONTINUE
    ELSE

C Pattern only

C  Count the number of entries in each column

      DO 90 K = 1,NE
      I = JCN(K)
      IP(I) = IP(I) + 1
   90   CONTINUE
      IP(NC+1) = NE + 1

C  Set IP so that IP(I) points to the first entry in column I+1

      IP(1) = IP(1) + 1
      DO 100 J = 2,NC
      IP(J) = IP(J) + IP(J-1)
  100   CONTINUE

      DO 120 I = NR,1,-1
      J1 = IW(I)
      J2 = IW(I+1) - 1
      DO 110 J = J1,J2
        K = JCN(J)
        L = IP(K) - 1
        IRN(L) = I
        IP(K) = L
  110   CONTINUE
  120   CONTINUE

    END IF

    RETURN
    END

C**********************************************************

    SUBROUTINE MC59D(NC,NE,IRN,IP,LA,A)
C
C To sort from arbitrary order within each column to order
C by increasing row index. Note: this is taken from MC20B/BD.
C
C   NC - integer variable. Intent(IN)
C    - on entry must be set to the number of columns in the matrix
C   NE - integer variable. Intent(IN)
C    - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(INOUT)
C    - on entry set to contain the row indices of the nonzeros
C      ordered so that the row
C      indices for column 1 precede those for column 2 and so on,
C      but the order within columns is arbitrary.
C      On exit, the order within each column is by increasing
C      row indices.
C   LA - integer variable which defines the length of the array A.
C      Intent(IN)
C    A - real (double precision/complex/complex*16) array of length LA
C      Intent(INOUT)
C    - if LA > 1, the array must be of length NE, and A(K)
C      must be set to the value of the entry in IRN(K);
C      on exit A is reordered in the same way as IRN
C    - if LA = 1, the array is not accessed
C  IP  - integer array of length NC. Intent(IN)
C    - on entry, IP(J) contains the position in IRN (and A) of the
C      first entry in column J (J=1,...,NC)
C   . .
C   .. Scalar Arguments ..
    INTEGER LA,NC,NE
C   ..
C   .. Array Arguments ..
    REAL A(LA)
    INTEGER IRN(NE),IP(NC)
C   ..
C   .. Local Scalars ..
    REAL ACE
    INTEGER ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
C   ..
C   .. Intrinsic Functions ..
    INTRINSIC ABS
C   ..
C   .. Executable Statements ..

C Jump if pattern only.
    IF (LA.GT.1) THEN
      KMAX = NE
      DO 50 JJ = 1,NC
      J = NC + 1 - JJ
      KLO = IP(J) + 1
      IF (KLO.GT.KMAX) GO TO 40
      KOR = KMAX
      DO 30 KDUMMY = KLO,KMAX
C Items KOR, KOR+1, .... ,KMAX are in order
        ACE = A(KOR-1)
        ICE = IRN(KOR-1)
        DO 10 K = KOR,KMAX
          IK = IRN(K)
          IF (ABS(ICE).LE.ABS(IK)) GO TO 20
          IRN(K-1) = IK
          A(K-1) = A(K)
   10     CONTINUE
        K = KMAX + 1
   20     IRN(K-1) = ICE
        A(K-1) = ACE
        KOR = KOR - 1
   30   CONTINUE
C Next column
   40   KMAX = KLO - 2
   50   CONTINUE
    ELSE

C Pattern only.
      KMAX = NE
      DO 150 JJ = 1,NC
      J = NC + 1 - JJ
      KLO = IP(J) + 1
      IF (KLO.GT.KMAX) GO TO 140
      KOR = KMAX
      DO 130 KDUMMY = KLO,KMAX
C Items KOR, KOR+1, .... ,KMAX are in order
        ICE = IRN(KOR-1)
        DO 110 K = KOR,KMAX
          IK = IRN(K)
          IF (ABS(ICE).LE.ABS(IK)) GO TO 120
          IRN(K-1) = IK
  110     CONTINUE
        K = KMAX + 1
  120     IRN(K-1) = ICE
        KOR = KOR - 1
  130   CONTINUE
C Next column
  140   KMAX = KLO - 2
  150   CONTINUE
    END IF
    END
C***********************************************************************

    SUBROUTINE MC59E(NC,NR,NE,IRN,LIP,IP,LA,A,IW,IDUP,KNE,ICNTL6)

C Checks IRN for duplicate entries.
C On exit, IDUP holds number of duplicates found and KNE is number
C of entries in matrix after removal of duplicates
C   . .
C   .. Scalar Arguments ..
    INTEGER ICNTL6,IDUP,KNE,LIP,LA,NC,NR,NE
C   ..
C   .. Array Arguments ..
    REAL A(LA)
    INTEGER IRN(NE),IP(LIP),IW(NR)
C   ..
C   .. Local Scalars ..
    INTEGER I,J,K,KSTART,KSTOP,NZJ

    IDUP = 0
    KNE = 0
C Initialise IW
    DO 10 I = 1,NR
      IW(I) = 0
   10 CONTINUE

    KSTART = IP(1)
    IF (LA.GT.1) THEN
C Matrix entries considered
      NZJ = 0
      DO 30 J = 1,NC
      KSTOP = IP(J+1)
      IP(J+1) = IP(J)
      DO 20 K = KSTART,KSTOP - 1
        I = IRN(K)
        IF (IW(I).LE.NZJ) THEN
          KNE = KNE + 1
          IRN(KNE) = I
          A(KNE) = A(K)
          IP(J+1) = IP(J+1) + 1
          IW(I) = KNE
        ELSE
C We have a duplicate in column J
          IDUP = IDUP + 1
C If requested, sum duplicates
          IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
        END IF
   20   CONTINUE
      KSTART = KSTOP
      NZJ = KNE
   30   CONTINUE

    ELSE

C Pattern only
      DO 50 J = 1,NC
      KSTOP = IP(J+1)
      IP(J+1) = IP(J)
      DO 40 K = KSTART,KSTOP - 1
        I = IRN(K)
        IF (IW(I).LT.J) THEN
          KNE = KNE + 1
          IRN(KNE) = I
          IP(J+1) = IP(J+1) + 1
          IW(I) = J
        ELSE
C  We have a duplicate in column J
          IDUP = IDUP + 1
        END IF
   40   CONTINUE
      KSTART = KSTOP
   50   CONTINUE
    END IF

    RETURN
    END
C***********************************************************************

    SUBROUTINE MC59F(NC,NR,NE,IRN,LIP,IP,LA,A,LIW,IW,IDUP,IOUT,
   +            IUP,KNE,ICNTL6,INFO)

C Checks IRN for duplicate and out-of-range entries.
C For symmetric matrix, also checks NO entries lie in upper triangle.
C Also checks IP is monotonic.
C On exit:
C IDUP holds number of duplicates found
C IOUT holds number of out-of-range entries
C For symmetric matrix, IUP holds number of entries in upper
C triangular part.
C KNE holds number of entries in matrix after removal of
C out-of-range and duplicate entries.
C Note: this is similar to MC59E except it also checks IP is
C monotonic and removes out-of-range entries in IRN.
C   . .
C   .. Scalar Arguments ..
    INTEGER ICNTL6,IDUP,IOUT,IUP,KNE,LA,LIP,LIW,NC,NR,NE
C   ..
C   .. Array Arguments ..
    REAL A(LA)
    INTEGER IRN(NE),IP(LIP),IW(LIW),INFO(2)
C   ..
C   .. Local Scalars ..
    INTEGER I,J,K,KSTART,KSTOP,NZJ,LOWER

    IDUP = 0
    IOUT = 0
    IUP = 0
    KNE = 0
C Initialise IW
    DO 10 I = 1,NR
      IW(I) = 0
   10 CONTINUE

    KSTART = IP(1)
    LOWER = 1
    IF (LA.GT.1) THEN
      NZJ = 0
      DO 30 J = 1,NC
C In symmetric case, entries out-of-range if they lie
C in upper triangular part.
      IF (ICNTL6.NE.0) LOWER = J
      KSTOP = IP(J+1)
      IF (KSTART.GT.KSTOP) THEN
        INFO(1) = -9
        INFO(2) = J
        RETURN
      END IF
      IP(J+1) = IP(J)
      DO 20 K = KSTART,KSTOP - 1
        I = IRN(K)
C Check for out-of-range
        IF (I.GT.NR .OR. I.LT.LOWER) THEN
          IOUT = IOUT + 1
C In symmetric case, check if entry is out-of-range because
C it lies in upper triangular part.
          IF (ICNTL6.NE.0 .AND. I.LT.J) IUP = IUP + 1
        ELSE IF (IW(I).LE.NZJ) THEN
          KNE = KNE + 1
          IRN(KNE) = I
          A(KNE) = A(K)
          IP(J+1) = IP(J+1) + 1
          IW(I) = KNE
        ELSE
C  We have a duplicate in column J
          IDUP = IDUP + 1
C If requested, sum duplicates
          IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
        END IF
   20   CONTINUE
      KSTART = KSTOP
      NZJ = KNE
   30   CONTINUE

    ELSE

C Pattern only
      DO 50 J = 1,NC
C In symmetric case, entries out-of-range if lie
C in upper triangular part.
      IF (ICNTL6.NE.0) LOWER = J
      KSTOP = IP(J+1)
      IF (KSTART.GT.KSTOP) THEN
        INFO(1) = -9
        INFO(2) = J
        RETURN
      END IF
      IP(J+1) = IP(J)
      DO  40 K = KSTART,KSTOP - 1
        I = IRN(K)
C Check for out-of-range
        IF (I.GT.NR .OR. I.LT.LOWER) THEN
          IOUT = IOUT + 1
          IF (ICNTL6.NE.0 .AND. I.GT.1) IUP = IUP + 1
        ELSE IF (IW(I).LT.J) THEN
          KNE = KNE + 1
          IRN(KNE) = I
          IP(J+1) = IP(J+1) + 1
          IW(I) = J
        ELSE
C  We have a duplicate in column J
          IDUP = IDUP + 1
        END IF
   40   CONTINUE
      KSTART = KSTOP
   50   CONTINUE
    END IF

    RETURN
    END

* COPYRIGHT (c) 1993 Council for the Central Laboratory
*            of the Research Councils

C Original date 29 Jan 2001
C 29 January 2001. Modified from MC49 to be threadsafe.

C 12th July 2004 Version 1.0.0. Version numbering added.
C 28 February 2008. Version 1.0.1. Comments flowed to column 72.
C 21 September 2009. Version 1.0.2. Minor change to documentation.

    SUBROUTINE MC59AD(ICNTL,NC,NR,NE,IRN,LJCN,JCN,LA,A,LIP,IP,
   &            LIW,IW,INFO)
C
C To sort the sparsity pattern of a matrix to an ordering by columns.
C There is an option for ordering the entries within each column by
C increasing row indices and an option for checking the user-supplied
C matrix entries for indices which are out-of-range or duplicated.
C
C ICNTL:  INTEGER array of length 10. Intent(IN). Used to specify
C       control parameters for the subroutine.
C ICNTL(1): indicates whether the user-supplied matrix entries are to
C       be checked for duplicates, and out-of-range indices.
C       Note  simple checks are always performed.
C       ICNTL(1) = 0, data checking performed.
C       Otherwise, no data checking.
C ICNTL(2): indicates the ordering requested.
C       ICNTL(2) = 0, input is by rows and columns in arbitrary
C       order and the output is sorted by columns.
C       ICNTL(2) = 1, the output is also row ordered
C       within each column.
C       ICNTL(2) = 2, the input is already ordered by
C       columns and is to be row ordered within each column.
C       Values outside the range 0 to 2 are flagged as an error.
C ICNTL(3): indicates whether matrix entries are also being ordered.
C       ICNTL(3) = 0, matrix entries are ordered.
C       Otherwise, only the sparsity pattern is ordered
C       and the array A is not accessed by the routine.
C ICNTL(4): the unit number of the device to
C       which error messages are sent. Error messages
C       can be suppressed by setting ICNTL(4) < 0.
C ICNTL(5): the unit number of the device to
C       which warning messages are sent. Warning
C       messages can be suppressed by setting ICNTL(5) < 0.
C ICNTL(6)  indicates whether matrix symmetric. If unsymmetric, ICNTL(6)
C       must be set to 0.
C       If ICNTL(6) = -1 or 1, symmetric and only the lower
C       triangular part of the reordered matrix is returned.
C       If ICNTL(6) = -2 or 2, Hermitian and only the lower
C       triangular part of the reordered matrix is returned.
C       If error checks are performed (ICNTL(1) = 0)
C       and ICNTL(6)> 1 or 2, the values of duplicate
C       entries are added together; if ICNTL(6) < -1 or -2, the
C       value of the first occurrence of the entry is used.
C ICNTL(7) to ICNTL(10) are not currently accessed by the routine.
C
C NC:    INTEGER variable. Intent(IN). Must be set by the user
C      to the number of columns in the matrix.
C NR:    INTEGER variable. Intent(IN). Must be set by the user
C      to the number of rows in the matrix.
C NE:    INTEGER variable. Intent(IN). Must be set by the user
C      to the number of entries in the matrix.
C IRN: INTEGER array of length NE. Intent (INOUT). Must be set by the
C        user to hold the row indices of the entries in the matrix.
C      If ICNTL(2).NE.2, the entries may be in any order.
C      If ICNTL(2).EQ.2, the entries in column J must be in
C        positions IP(J) to IP(J+1)-1 of IRN. On exit, the row
C        indices are reordered so that the entries of a single
C        column are contiguous with column J preceding column J+1, J
C        = 1, 2, ..., NC-1, with no space between columns.
C      If ICNTL(2).EQ.0, the order within each column is arbitrary;
C        if ICNTL(2) = 1 or 2, the order within each column is by
C        increasing row indices.
C LJCN:    INTEGER variable. Intent(IN). Defines length array
C JCN:   INTEGER array of length LJCN. Intent (INOUT).
C      If ICNTL(2) = 0 or 1, JCN(K) must be set by the user
C      to the column index of the entry
C      whose row index is held in IRN(K), K = 1, 2, ..., NE.
C      On exit, the contents of this array  will have been altered.
C      If ICNTL(2) = 2, the array is not accessed.
C LA:    INTEGER variable. Intent(IN). Defines length of array
C      A.
C A:     is a REAL (DOUBLE PRECISION in the D version, INTEGER in
C      the I version, COMPLEX in the C version,
C      or COMPLEX"*"16 in the Z version) array of length LA.
C      Intent(INOUT).
C      If ICNTL(3).EQ.0, A(K) must be set by the user to
C      hold the value of the entry with row index IRN(K),
C      K = 1, 2, ..., NE. On exit, the array will have been
C      permuted in the same way as the array IRN.
C      If ICNTL(3).NE.0, the array is not accessed.
C LIP:   INTEGER variable. Intent(IN). Defines length of array
C      IP.
C IP:    INTEGER array of length LIP. Intent(INOUT). IP
C      need only be set by the user if ICNTL(2) = 2.
C      In this case, IP(J) holds the position in
C      the array IRN of the first entry in column J, J = 1, 2,
C      ..., NC, and IP(NC+1) is one greater than the number of
C      entries in the matrix.
C      In all cases, the array IP will have this meaning on exit
C      from the subroutine and is altered when ICNTL(2) = 2 only
C      when ICNTL(1) =  0 and there are out-of-range
C      indices or duplicates.
C LIW:   INTEGER variable. Intent(IN). Defines length of array
C      IW.
C IW:    INTEGER array of length LIW. Intent(OUT). Used by the
C      routine as workspace.
C INFO:    INTEGER array of length 10.  Intent(OUT). On exit,
C      a negative value of INFO(1) is used to signal a fatal
C      error in the input data, a positive value of INFO(1)
C      indicates that a warning has been issued, and a
C      zero value is used to indicate a successful call.
C      In cases of error, further information is held in INFO(2).
C      For warnings, further information is
C      provided in INFO(3) to INFO(6).  INFO(7) to INFO(10) are not
C      currently used and are set to zero.
C      Possible nonzero values of INFO(1):
C       -1 -  The restriction ICNTL(2) = 0, 1, or 2 violated.
C         Value of ICNTL(2) is given by INFO(2).
C       -2 -  NC.LE.0. Value of NC is given by INFO(2).
C       -3 -  Error in NR. Value of NR is given by INFO(2).
C       -4 -  NE.LE.0. Value of NE is given by INFO(2).
C       -5 -  LJCN too small. Min. value of LJCN is given by INFO(2).
C       -6 -  LA too small. Min. value of LA is given by INFO(2).
C       -7 -  LIW too small. Value of LIW is given by INFO(2).
C       -8 -  LIP too small. Value of LIP is given by INFO(2).
C       -9 -  The entries of IP not monotonic increasing.
C      -10 -  For each I, IRN(I) or JCN(I) out-of-range.
C      -11 -  ICNTL(6) is out of range.
C       +1 -  One or more duplicated entries. One copy of
C         each such entry is kept and, if ICNTL(3) = 0 and
C         ICNTL(6).GE.0, the values of these entries are
C         added together. If  ICNTL(3) = 0 and ICNTL(6).LT.0,
C         the value of the first occurrence of the entry is used.
C         Initially INFO(3) is set to zero. If an entry appears
C         k times, INFO(3) is incremented by k-1 and INFO(6)
C         is set to the revised number of entries in the
C         matrix.
C       +2 - One or more of the entries in IRN out-of-range. These
C         entries are removed by the routine.`INFO(4) is set to
C         the number of entries which were out-of-range and
C         INFO(6) is set to the revised number of entries in the
C         matrix.
C       +4 - One or more of the entries in JCN out-of-range. These
C         entries are removed by the routine. INFO(5) is set to
C         the number of entries which were out-of-range and
C         INFO(6) is set to the revised number of entries in the
C         matrix. Positive values of INFO(1) are summed so that
C         the user can identify all warnings.
C
C   .. Scalar Arguments ..
    INTEGER LA,LIP,LIW,LJCN,NC,NE,NR
C   ..
C   .. Array Arguments ..
    DOUBLE PRECISION A(LA)
    INTEGER ICNTL(10),IP(LIP),INFO(10),IRN(NE),IW(LIW),JCN(LJCN)
C   ..
C   .. Local Scalars ..
    INTEGER I,ICNTL1,ICNTL2,ICNTL3,ICNTL6,LAA
    INTEGER IDUP,IOUT,IUP,JOUT,LP,MP,KNE,PART
    LOGICAL LCHECK
C   ..
C   .. External Subroutines ..
    EXTERNAL MC59BD,MC59CD,MC59DD,MC59ED,MC59FD
C   ..
C   .. Intrinsic Functions ..
    INTRINSIC MAX
C   ..
C   .. Executable Statements ..

C Initialise
    DO 10 I = 1,10
       INFO(I) = 0
   10 CONTINUE

    ICNTL1 = ICNTL(1)
    ICNTL2 = ICNTL(2)
    ICNTL3 = ICNTL(3)
    ICNTL6 = ICNTL(6)
    LCHECK = (ICNTL1.EQ.0)
C Streams for errors/warnings
    LP = ICNTL(4)
    MP = ICNTL(5)

C  Check the input data
    IF (ICNTL2.GT.2 .OR. ICNTL2.LT.0) THEN
       INFO(1) = -1
       INFO(2) = ICNTL2
       IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9010) ICNTL2
       END IF
       GO TO 70
    END IF

    IF (ICNTL6.GT.2 .OR. ICNTL6.LT.-2) THEN
       INFO(1) = -11
       INFO(2) = ICNTL6
       IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9150) ICNTL6
       END IF
       GO TO 70
    END IF
C For real matrices, symmetric = Hermitian so only
C have to distinguish between unsymmetric (ICNTL6 = 0) and
C symmetric (ICNTL6.ne.0)

    IF (NC.LT.1) THEN
      INFO(1) = -2
      INFO(2) = NC
      IF (LP.GT.0) THEN
      WRITE (LP,FMT=9000) INFO(1)
      WRITE (LP,FMT=9020) NC
      END IF
      GO TO 70
    END IF

    IF (NR.LT.1) THEN
      INFO(1) = -3
      INFO(2) = NR
      IF (LP.GT.0) THEN
      WRITE (LP,FMT=9000) INFO(1)
      WRITE (LP,FMT=9030) NR
      END IF
      GO TO 70
    END IF

    IF (ICNTL6.NE.0 .AND. NR.NE.NC) THEN
      INFO(1) = -3
      INFO(2) = NR
      IF (LP.GT.0) THEN
      WRITE (LP,FMT=9000) INFO(1)
      WRITE (LP,FMT=9035) NC,NR
      END IF
      GO TO 70
    END IF

    IF (NE.LT.1) THEN
      INFO(1) = -4
      INFO(2) = NE
      IF (LP.GT.0) THEN
      WRITE (LP,FMT=9000) INFO(1)
      WRITE (LP,FMT=9040) NE
      END IF
      GO TO 70
    END IF

    IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.1) THEN
      IF (LJCN.LT.NE) THEN
      INFO(1) = -5
      INFO(2) = NE
      END IF
    ELSE
      IF (LJCN.LT.1) THEN
      INFO(1) = -5
      INFO(2) = 1
      END IF
    END IF
    IF (INFO(1).EQ.-5) THEN
       IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9050) LJCN,INFO(2)
       END IF
       GO TO 70
    END IF

    IF (ICNTL3.EQ.0) THEN
      IF (LA.LT.NE) THEN
      INFO(1) = -6
      INFO(2) = NE
      END IF
    ELSE
      IF (LA.LT.1) THEN
      INFO(1) = -6
      INFO(2) = 1
      END IF
    END IF
    IF (INFO(1).EQ.-6) THEN
       IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9060) LA,INFO(2)
       END IF
       GO TO 70
    END IF

    IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.2) THEN
      IF (LIP.LT.NC+1) THEN
      INFO(1) = -7
      INFO(2) = NC+1
      END IF
    ELSE IF (LIP.LT.MAX(NR,NC)+1) THEN
      INFO(1) = -7
      INFO(2) = MAX(NR,NC)+1
    END IF
    IF (INFO(1).EQ.-7) THEN
      IF (LP.GT.0) THEN
      WRITE (LP,FMT=9000) INFO(1)
      WRITE (LP,FMT=9065) LIP,INFO(2)
      END IF
      GO TO 70
    END IF

C Check workspace sufficient
    IF (LIW.LT.MAX(NR,NC)+1) THEN
      INFO(1) = -8
      INFO(2) = MAX(NR,NC)+1
      IF (LP.GT.0) THEN
      WRITE (LP,FMT=9000) INFO(1)
      WRITE (LP,FMT=9070) LIW,INFO(2)
      END IF
      GO TO 70
    END IF

    LAA = NE
    IF (ICNTL3.NE.0) LAA = 1
C Initialise counts of number of out-of-range entries and duplicates
    IOUT = 0
    JOUT = 0
    IDUP = 0
    IUP = 0

C PART is used by MC59BD to indicate if upper or lower or
C all of matrix is required.
C PART =  0 : unsymmetric case, whole matrix wanted
C PART =  1 : symmetric case, lower triangular part of matrix wanted
C PART = -1 : symmetric case, upper triangular part of matrix wanted
    PART = 0
    IF (ICNTL6.NE.0) PART = 1

    IF (ICNTL2.EQ.0) THEN

C Order directly by columns
C On exit from MC59BD, KNE holds number of entries in matrix
C after removal of out-of-range entries. If no data checking, KNE = NE.
      CALL MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LAA,A,IP,IW,
   +          IOUT,JOUT,KNE)
C Return if ALL entries out-of-range.
      IF (KNE.EQ.0) GO TO 50

C Check for duplicates
      IF (LCHECK) CALL MC59ED(NC,NR,NE,IRN,LIP,IP,LAA,A,IW,IDUP,
   &                KNE,ICNTL6)

    ELSE IF (ICNTL2.EQ.1) THEN

C First order by rows.
C Interchanged roles of IRN and JCN, so set PART = -1
C if matrix is symmetric case
      IF (ICNTL6.NE.0) PART = -1
      CALL MC59BD(LCHECK,PART,NR,NC,NE,JCN,IRN,LAA,A,IW,IP,
   +          JOUT,IOUT,KNE)
C Return if ALL entries out-of-range.
      IF (KNE.EQ.0) GO TO 50

C At this point, JCN and IW hold column indices and row pointers
C Optionally, check for duplicates.
      IF (LCHECK) CALL MC59ED(NR,NC,NE,JCN,NR+1,IW,LAA,A,IP,
   &                IDUP,KNE,ICNTL6)

C Now order by columns and by rows within each column
      CALL MC59CD(NC,NR,KNE,IRN,JCN,LAA,A,IP,IW)

    ELSE IF (ICNTL2.EQ.2) THEN
C Input is using IP, IRN.
C Optionally check for duplicates and remove out-of-range entries
      IF (LCHECK) THEN
      CALL MC59FD(NC,NR,NE,IRN,NC+1,IP,LAA,A,LIW,IW,IDUP,
   +          IOUT,IUP,KNE,ICNTL6,INFO)
C Return if IP not monotonic.
      IF (INFO(1).EQ.-9) GO TO 40
C Return if ALL entries out-of-range.
      IF (KNE.EQ.0) GO TO 50
      ELSE
       KNE = NE
      END IF

C  Order by rows within each column
      CALL MC59DD(NC,KNE,IRN,IP,LAA,A)

    END IF

    INFO(3) = IDUP
    INFO(4) = IOUT
    INFO(5) = JOUT
    INFO(6) = KNE
    INFO(7) = IUP

C Set warning flag if out-of-range /duplicates found
    IF (IDUP.GT.0) INFO(1) = INFO(1) + 1
    IF (IOUT.GT.0) INFO(1) = INFO(1) + 2
    IF (JOUT.GT.0) INFO(1) = INFO(1) + 4
    IF (INFO(1).GT.0 .AND. MP.GT.0) THEN
      WRITE (MP,FMT=9080) INFO(1)
      IF (IOUT.GT.0) WRITE (MP,FMT=9090) IOUT
      IF (JOUT.GT.0) WRITE (MP,FMT=9110) JOUT
      IF (IDUP.GT.0) WRITE (MP,FMT=9100) IDUP
      IF (IUP.GT.0)  WRITE (MP,FMT=9130) IUP
    END IF
    GO TO 70

   40 INFO(3) = IDUP
    INFO(4) = IOUT
    INFO(7) = IUP
    IF (LP.GT.0) THEN
      WRITE (LP,FMT=9000) INFO(1)
      WRITE (LP,FMT=9140)
    END IF
    GO TO 70

   50 INFO(1) = -10
    INFO(4) = IOUT
    INFO(5) = JOUT
    INFO(2) = IOUT + JOUT
    IF (LP.GT.0) THEN
      WRITE (LP,FMT=9000) INFO(1)
      WRITE (LP,FMT=9120)
    END IF
   70 RETURN

 9000 FORMAT (/,' *** Error return from MC59AD *** INFO(1) = ',I3)
 9010 FORMAT (1X,'ICNTL(2) = ',I2,' is out of range')
 9020 FORMAT (1X,'NC = ',I6,' is out of range')
 9030 FORMAT (1X,'NR = ',I6,' is out of range')
 9035 FORMAT (1X,'Symmetric case. NC = ',I6,' but NR = ',I6)
 9040 FORMAT (1X,'NE = ',I10,' is out of range')
 9050 FORMAT (1X,'Increase LJCN from ',I10,' to at least ',I10)
 9060 FORMAT (1X,'Increase LA from ',I10,' to at least ',I10)
 9065 FORMAT (1X,'Increase LIP from ',I8,' to at least ',I10)
 9070 FORMAT (1X,'Increase LIW from ',I8,' to at least ',I10)
 9080 FORMAT (/,' *** Warning message from MC59AD *** INFO(1) = ',I3)
 9090 FORMAT (1X,I8,' entries in IRN supplied by the user were ',
   +     /,'     out of range and were ignored by the routine')
 9100 FORMAT (1X,I8,' duplicate entries were supplied by the user')
 9110 FORMAT (1X,I8,' entries in JCN supplied by the user were ',
   +     /,'     out of range and were ignored by the routine')
 9120 FORMAT (1X,'All entries out of range')
 9130 FORMAT (1X,I8,' of these entries were in the upper triangular ',
   +     /,'     part of matrix')
 9140 FORMAT (1X,'Entries in IP are not monotonic increasing')
 9150 FORMAT (1X,'ICNTL(6) = ',I2,' is out of range')
    END
C***********************************************************************
    SUBROUTINE MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LA,A,IP,IW,IOUT,
   +            JOUT,KNE)
C
C   To sort a sparse matrix from arbitrary order to
C   column order, unordered within each column. Optionally
C   checks for out-of-range entries in IRN,JCN.
C
C LCHECK - logical variable. Intent(IN). If true, check
C      for out-of-range indices.
C PART -   integer variable. Intent(IN)
C PART =  0 : unsymmetric case, whole matrix wanted
C PART =  1 : symmetric case, lower triangular part of matrix wanted
C         (ie IRN(K) .ge. JCN(K) on exit)
C PART = -1 : symmetric case, upper triangular part of matrix wanted
C         (ie IRN(K) .le. JCN(K) on exit)
C   NC - integer variable. Intent(IN)
C    - on entry must be set to the number of columns in the matrix
C   NR - integer variable. Intent(IN)
C    - on entry must be set to the number of rows in the matrix
C   NE - integer variable. Intent(IN)
C    - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(INOUT)
C    - on entry set to contain the row indices of the nonzeros
C      in arbitrary order.
C    - on exit, the entries in IRN are reordered so that the row
C      indices for column 1 precede those for column 2 and so on,
C      but the order within columns is arbitrary.
C  JCN - integer array of length NE. Intent(INOUT)
C    - on entry set to contain the column indices of the nonzeros
C    - JCN(K) must be the column index of
C      the entry in IRN(K)
C    - on exit, JCN(K) is the column index for the entry with
C      row index IRN(K) (K=1,...,NE).
C  LA  - integer variable which defines the length of the array A.
C      Intent(IN)
C   A  - real (double precision/complex/complex*16) array of length LA
C      Intent(INOUT)
C    - if LA > 1, the array must be of length NE, and A(K)
C      must be set to the value of the entry in (IRN(K), JCN(K));
C      on exit A is reordered in the same way as IRN
C    - if LA = 1, the array is not accessed
C  IP  - integer array of length NC+1. Intent(INOUT)
C    - not set on entry
C    - on exit, IP(J) contains the position in IRN (and A) of the
C      first entry in column J (J=1,...,NC)
C    - IP(NC+1) is set to NE+1
C  IW  - integer array of length NC+1.  Intent(INOUT)
C    - the array is used as workspace
C    - on exit IW(I) = IP(I) (so IW(I) points to the beginning
C      of column I).
C IOUT - integer variable. Intent(OUT). On exit, holds number
C      of entries in IRN found to be out-of-range
C JOUT - integer variable. Intent(OUT). On exit, holds number
C      of entries in JCN found to be out-of-range
C  KNE - integer variable. Intent(OUT). On exit, holds number
C      of entries in matrix after removal of out-of-range entries.
C      If no data checking, KNE = NE.

C   .. Scalar Arguments ..
    INTEGER LA,NC,NE,NR,IOUT,JOUT,KNE,PART
    LOGICAL LCHECK
C   ..
C   .. Array Arguments ..
    DOUBLE PRECISION A(LA)
    INTEGER IP(NC+1),IRN(NE),IW(NC+1),JCN(NE)
C   ..
C   .. Local Scalars ..
    DOUBLE PRECISION ACE,ACEP
    INTEGER I,ICE,ICEP,J,JCE,JCEP,K,L,LOC
C   ..
C   .. Executable Statements ..

C Initialise IW
    DO 10 J = 1,NC + 1
      IW(J) = 0
   10 CONTINUE

    KNE = 0
    IOUT = 0
    JOUT = 0
C Count the number of entries in each column and store in IW.
C We also allow checks for out-of-range indices
    IF (LCHECK) THEN
C Check data.
C Treat case of pattern only separately.
      IF (LA.GT.1) THEN
      IF (PART.EQ.0) THEN
C Unsymmetric
        DO 20 K = 1,NE
          I = IRN(K)
          J = JCN(K)
          IF (I.GT.NR .OR. I.LT.1) THEN
          IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
          IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
          ELSE IF (J.GT.NC .OR. J.LT.1) THEN
          JOUT = JOUT + 1
          ELSE
          KNE = KNE + 1
          IRN(KNE) = I
          JCN(KNE) = J
          A(KNE) = A(K)
          IW(J) = IW(J) + 1
          END IF
   20     CONTINUE
      ELSE IF (PART.EQ.1) THEN
C Symmetric, lower triangle
        DO 21 K = 1,NE
          I = IRN(K)
          J = JCN(K)
          IF (I.GT.NR .OR. I.LT.1) THEN
          IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
          IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
          ELSE IF (J.GT.NC .OR. J.LT.1) THEN
          JOUT = JOUT + 1
          ELSE
          KNE = KNE + 1
C Lower triangle ... swap if necessary
          IF (I.LT.J) THEN
            IRN(KNE) = J
            JCN(KNE) = I
            IW(I) = IW(I) + 1
          ELSE
            IRN(KNE) = I
            JCN(KNE) = J
            IW(J) = IW(J) + 1
          END IF
          A(KNE) = A(K)
          END IF
   21     CONTINUE
      ELSE IF (PART.EQ.-1) THEN
C Symmetric, upper triangle
        DO 22 K = 1,NE
          I = IRN(K)
          J = JCN(K)
          IF (I.GT.NR .OR. I.LT.1) THEN
          IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
          IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
          ELSE IF (J.GT.NC .OR. J.LT.1) THEN
          JOUT = JOUT + 1
          ELSE
          KNE = KNE + 1
C Upper triangle ... swap if necessary
          IF (I.GT.J) THEN
            IRN(KNE) = J
            JCN(KNE) = I
            IW(I) = IW(I) + 1
          ELSE
            IRN(KNE) = I
            JCN(KNE) = J
            IW(J) = IW(J) + 1
          END IF
          A(KNE) = A(K)
          END IF
   22     CONTINUE
      END IF
      ELSE
C Pattern only
      IF (PART.EQ.0) THEN
        DO 25 K = 1,NE
          I = IRN(K)
          J = JCN(K)
          IF (I.GT.NR .OR. I.LT.1) THEN
          IOUT = IOUT + 1
          IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
          ELSE IF (J.GT.NC .OR. J.LT.1) THEN
          JOUT = JOUT + 1
          ELSE
          KNE = KNE + 1
          IRN(KNE) = I
          JCN(KNE) = J
          IW(J) = IW(J) + 1
          END IF
   25     CONTINUE
      ELSE IF (PART.EQ.1) THEN
        DO 26 K = 1,NE
          I = IRN(K)
          J = JCN(K)
          IF (I.GT.NR .OR. I.LT.1) THEN
          IOUT = IOUT + 1
          IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
          ELSE IF (J.GT.NC .OR. J.LT.1) THEN
          JOUT = JOUT + 1
          ELSE
          KNE = KNE + 1
C Lower triangle ... swap if necessary
          IF (I.LT.J) THEN
            IRN(KNE) = J
            JCN(KNE) = I
            IW(I) = IW(I) + 1
          ELSE
            IRN(KNE) = I
            JCN(KNE) = J
            IW(J) = IW(J) + 1
          END IF
          END IF
   26     CONTINUE
      ELSE IF (PART.EQ.-1) THEN
        DO 27 K = 1,NE
          I = IRN(K)
          J = JCN(K)
          IF (I.GT.NR .OR. I.LT.1) THEN
          IOUT = IOUT + 1
          IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
          ELSE IF (J.GT.NC .OR. J.LT.1) THEN
          JOUT = JOUT + 1
          ELSE
          KNE = KNE + 1
C Upper triangle ... swap if necessary
          IF (I.GT.J) THEN
            IRN(KNE) = J
            JCN(KNE) = I
            IW(I) = IW(I) + 1
          ELSE
            IRN(KNE) = I
            JCN(KNE) = J
            IW(J) = IW(J) + 1
          END IF
          END IF
   27     CONTINUE
      END IF
      END IF
C Return if ALL entries out-of-range.
      IF (KNE.EQ.0) GO TO 130

    ELSE

C No checks
      KNE = NE
      IF (PART.EQ.0) THEN
      DO 30 K = 1,NE
        J = JCN(K)
        IW(J) = IW(J) + 1
   30   CONTINUE
      ELSE IF (PART.EQ.1) THEN
      DO 35 K = 1,NE
        I = IRN(K)
        J = JCN(K)
C Lower triangle ... swap if necessary
        IF (I.LT.J) THEN
         IRN(K) = J
         JCN(K) = I
         IW(I) = IW(I) + 1
        ELSE
          IW(J) = IW(J) + 1
        END IF
   35   CONTINUE
      ELSE IF (PART.EQ.-1) THEN
      DO 36 K = 1,NE
        I = IRN(K)
        J = JCN(K)
C Upper triangle ... swap if necessary
        IF (I.GT.J) THEN
         IRN(K) = J
         JCN(K) = I
         IW(I) = IW(I) + 1
        ELSE
          IW(J) = IW(J) + 1
        END IF
   36   CONTINUE
      END IF
    END IF

C KNE is now the number of nonzero entries in matrix.

C Put into IP and IW the positions where each column
C would begin in a compressed collection with the columns
C in natural order.

    IP(1) = 1
    DO 37 J = 2,NC + 1
      IP(J) = IW(J-1) + IP(J-1)
      IW(J-1) = IP(J-1)
   37 CONTINUE

C Reorder the elements into column order.
C Fill in each column from the front, and as a new entry is placed
C in column K increase the pointer IW(K) by one.

    IF (LA.EQ.1) THEN
C Pattern only
      DO 70 L = 1,NC
      DO 60 K = IW(L),IP(L+1) - 1
        ICE = IRN(K)
        JCE = JCN(K)
        DO 40 J = 1,NE
          IF (JCE.EQ.L) GO TO 50
          LOC = IW(JCE)
          JCEP = JCN(LOC)
          ICEP = IRN(LOC)
          IW(JCE) = LOC + 1
          JCN(LOC) = JCE
          IRN(LOC) = ICE
          JCE = JCEP
          ICE = ICEP
   40     CONTINUE
   50     JCN(K) = JCE
        IRN(K) = ICE
   60   CONTINUE
   70   CONTINUE
    ELSE

      DO 120 L = 1,NC
      DO 110 K = IW(L),IP(L+1) - 1
        ICE = IRN(K)
        JCE = JCN(K)
        ACE = A(K)
        DO 90 J = 1,NE
          IF (JCE.EQ.L) GO TO 100
          LOC = IW(JCE)
          JCEP = JCN(LOC)
          ICEP = IRN(LOC)
          IW(JCE) = LOC + 1
          JCN(LOC) = JCE
          IRN(LOC) = ICE
          JCE = JCEP
          ICE = ICEP
          ACEP = A(LOC)
          A(LOC) = ACE
          ACE = ACEP
   90     CONTINUE
  100     JCN(K) = JCE
        IRN(K) = ICE
        A(K) = ACE
  110   CONTINUE
  120   CONTINUE
    END IF

  130 CONTINUE

    RETURN
    END
C
C**********************************************************
    SUBROUTINE MC59CD(NC,NR,NE,IRN,JCN,LA,A,IP,IW)
C
C   To sort a sparse matrix stored by rows,
C   unordered within each row, to ordering by columns, with
C   ordering by rows within each column.
C
C   NC - integer variable. Intent(IN)
C    - on entry must be set to the number of columns in the matrix
C   NR - integer variable. Intent(IN)
C    - on entry must be set to the number of rows in the matrix
C  NE - integer variable. Intent(IN)
C    - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(OUT).
C    - not set on entry.
C    - on exit,  IRN holds row indices with the row
C      indices for column 1 preceding those for column 2 and so on,
C      with ordering by rows within each column.
C  JCN - integer array of length NE. Intent(INOUT)
C    - on entry set to contain the column indices of the nonzeros
C      with indices for column 1 preceding those for column 2
C      and so on, with the order within columns is arbitrary.
C    - on exit, contents destroyed.
C  LA  - integer variable which defines the length of the array A.
C      Intent(IN)
C   A  - real (double precision/complex/complex*16) array of length LA
C      Intent(INOUT)
C    - if LA > 1, the array must be of length NE, and A(K)
C      must be set to the value of the entry in JCN(K);
C      on exit A, A(K) holds the value of the entry in IRN(K).
C    - if LA = 1, the array is not accessed
C  IP  - integer array of length NC+1. Intent(INOUT)
C    - not set on entry
C    - on exit, IP(J) contains the position in IRN (and A) of the
C      first entry in column J (J=1,...,NC)
C    - IP(NC+1) is set to NE+1
C  IW  - integer array of length NR+1.  Intent(IN)
C    - on entry, must be set on entry so that IW(J) points to the
C      position in JCN of the first entry in row J, J=1,...,NR, and
C      IW(NR+1) must be set to NE+1
C
C   .. Scalar Arguments ..
    INTEGER LA,NC,NE,NR
C   ..
C   .. Array Arguments ..
    DOUBLE PRECISION A(LA)
    INTEGER IP(NC+1),IRN(NE),IW(NR+1),JCN(NE)
C   ..
C   .. Local Scalars ..
    DOUBLE PRECISION ACE,ACEP
    INTEGER I,ICE,ICEP,J,J1,J2,K,L,LOC,LOCP
C   ..
C   .. Executable Statements ..

C  Count the number of entries in each column

    DO 10 J = 1,NC
      IP(J) = 0
   10 CONTINUE

    IF (LA.GT.1) THEN

      DO 20 K = 1,NE
      I = JCN(K)
      IP(I) = IP(I) + 1
      IRN(K) = JCN(K)
   20   CONTINUE
      IP(NC+1) = NE + 1

C  Set IP so that IP(I) points to the first entry in column I+1

      IP(1) = IP(1) + 1
      DO 30 J = 2,NC
      IP(J) = IP(J) + IP(J-1)
   30   CONTINUE

      DO 50 I = NR,1,-1
      J1 = IW(I)
      J2 = IW(I+1) - 1
      DO 40 J = J1,J2
        K = IRN(J)
        L = IP(K) - 1
        JCN(J) = L
        IRN(J) = I
        IP(K) = L
   40   CONTINUE
   50   CONTINUE
      IP(NC+1) = NE + 1
      DO 70 J = 1,NE
      LOC = JCN(J)
      IF (LOC.EQ.0) GO TO 70
      ICE = IRN(J)
      ACE = A(J)
      JCN(J) = 0
      DO 60 K = 1,NE
        LOCP = JCN(LOC)
        ICEP = IRN(LOC)
        ACEP = A(LOC)
        JCN(LOC) = 0
        IRN(LOC) = ICE
        A(LOC) = ACE
        IF (LOCP.EQ.0) GO TO 70
        ICE = ICEP
        ACE = ACEP
        LOC = LOCP
   60   CONTINUE
   70   CONTINUE
    ELSE

C Pattern only

C  Count the number of entries in each column

      DO 90 K = 1,NE
      I = JCN(K)
      IP(I) = IP(I) + 1
   90   CONTINUE
      IP(NC+1) = NE + 1

C  Set IP so that IP(I) points to the first entry in column I+1

      IP(1) = IP(1) + 1
      DO 100 J = 2,NC
      IP(J) = IP(J) + IP(J-1)
  100   CONTINUE

      DO 120 I = NR,1,-1
      J1 = IW(I)
      J2 = IW(I+1) - 1
      DO 110 J = J1,J2
        K = JCN(J)
        L = IP(K) - 1
        IRN(L) = I
        IP(K) = L
  110   CONTINUE
  120   CONTINUE

    END IF

    RETURN
    END

C**********************************************************

    SUBROUTINE MC59DD(NC,NE,IRN,IP,LA,A)
C
C To sort from arbitrary order within each column to order
C by increasing row index. Note: this is taken from MC20B/BD.
C
C   NC - integer variable. Intent(IN)
C    - on entry must be set to the number of columns in the matrix
C   NE - integer variable. Intent(IN)
C    - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(INOUT)
C    - on entry set to contain the row indices of the nonzeros
C      ordered so that the row
C      indices for column 1 precede those for column 2 and so on,
C      but the order within columns is arbitrary.
C      On exit, the order within each column is by increasing
C      row indices.
C   LA - integer variable which defines the length of the array A.
C      Intent(IN)
C    A - real (double precision/complex/complex*16) array of length LA
C      Intent(INOUT)
C    - if LA > 1, the array must be of length NE, and A(K)
C      must be set to the value of the entry in IRN(K);
C      on exit A is reordered in the same way as IRN
C    - if LA = 1, the array is not accessed
C  IP  - integer array of length NC. Intent(IN)
C    - on entry, IP(J) contains the position in IRN (and A) of the
C      first entry in column J (J=1,...,NC)
C   . .
C   .. Scalar Arguments ..
    INTEGER LA,NC,NE
C   ..
C   .. Array Arguments ..
    DOUBLE PRECISION A(LA)
    INTEGER IRN(NE),IP(NC)
C   ..
C   .. Local Scalars ..
    DOUBLE PRECISION ACE
    INTEGER ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
C   ..
C   .. Intrinsic Functions ..
    INTRINSIC ABS
C   ..
C   .. Executable Statements ..

C Jump if pattern only.
    IF (LA.GT.1) THEN
      KMAX = NE
      DO 50 JJ = 1,NC
      J = NC + 1 - JJ
      KLO = IP(J) + 1
      IF (KLO.GT.KMAX) GO TO 40
      KOR = KMAX
      DO 30 KDUMMY = KLO,KMAX
C Items KOR, KOR+1, .... ,KMAX are in order
        ACE = A(KOR-1)
        ICE = IRN(KOR-1)
        DO 10 K = KOR,KMAX
          IK = IRN(K)
          IF (ABS(ICE).LE.ABS(IK)) GO TO 20
          IRN(K-1) = IK
          A(K-1) = A(K)
   10     CONTINUE
        K = KMAX + 1
   20     IRN(K-1) = ICE
        A(K-1) = ACE
        KOR = KOR - 1
   30   CONTINUE
C Next column
   40   KMAX = KLO - 2
   50   CONTINUE
    ELSE

C Pattern only.
      KMAX = NE
      DO 150 JJ = 1,NC
      J = NC + 1 - JJ
      KLO = IP(J) + 1
      IF (KLO.GT.KMAX) GO TO 140
      KOR = KMAX
      DO 130 KDUMMY = KLO,KMAX
C Items KOR, KOR+1, .... ,KMAX are in order
        ICE = IRN(KOR-1)
        DO 110 K = KOR,KMAX
          IK = IRN(K)
          IF (ABS(ICE).LE.ABS(IK)) GO TO 120
          IRN(K-1) = IK
  110     CONTINUE
        K = KMAX + 1
  120     IRN(K-1) = ICE
        KOR = KOR - 1
  130   CONTINUE
C Next column
  140   KMAX = KLO - 2
  150   CONTINUE
    END IF
    END
C***********************************************************************

    SUBROUTINE MC59ED(NC,NR,NE,IRN,LIP,IP,LA,A,IW,IDUP,KNE,ICNTL6)

C Checks IRN for duplicate entries.
C On exit, IDUP holds number of duplicates found and KNE is number
C of entries in matrix after removal of duplicates
C   . .
C   .. Scalar Arguments ..
    INTEGER ICNTL6,IDUP,KNE,LIP,LA,NC,NR,NE
C   ..
C   .. Array Arguments ..
    DOUBLE PRECISION A(LA)
    INTEGER IRN(NE),IP(LIP),IW(NR)
C   ..
C   .. Local Scalars ..
    INTEGER I,J,K,KSTART,KSTOP,NZJ

    IDUP = 0
    KNE = 0
C Initialise IW
    DO 10 I = 1,NR
      IW(I) = 0
   10 CONTINUE

    KSTART = IP(1)
    IF (LA.GT.1) THEN
C Matrix entries considered
      NZJ = 0
      DO 30 J = 1,NC
      KSTOP = IP(J+1)
      IP(J+1) = IP(J)
      DO 20 K = KSTART,KSTOP - 1
        I = IRN(K)
        IF (IW(I).LE.NZJ) THEN
          KNE = KNE + 1
          IRN(KNE) = I
          A(KNE) = A(K)
          IP(J+1) = IP(J+1) + 1
          IW(I) = KNE
        ELSE
C We have a duplicate in column J
          IDUP = IDUP + 1
C If requested, sum duplicates
          IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
        END IF
   20   CONTINUE
      KSTART = KSTOP
      NZJ = KNE
   30   CONTINUE

    ELSE

C Pattern only
      DO 50 J = 1,NC
      KSTOP = IP(J+1)
      IP(J+1) = IP(J)
      DO 40 K = KSTART,KSTOP - 1
        I = IRN(K)
        IF (IW(I).LT.J) THEN
          KNE = KNE + 1
          IRN(KNE) = I
          IP(J+1) = IP(J+1) + 1
          IW(I) = J
        ELSE
C  We have a duplicate in column J
          IDUP = IDUP + 1
        END IF
   40   CONTINUE
      KSTART = KSTOP
   50   CONTINUE
    END IF

    RETURN
    END
C***********************************************************************

    SUBROUTINE MC59FD(NC,NR,NE,IRN,LIP,IP,LA,A,LIW,IW,IDUP,IOUT,
   +            IUP,KNE,ICNTL6,INFO)

C Checks IRN for duplicate and out-of-range entries.
C For symmetric matrix, also checks NO entries lie in upper triangle.
C Also checks IP is monotonic.
C On exit:
C IDUP holds number of duplicates found
C IOUT holds number of out-of-range entries
C For symmetric matrix, IUP holds number of entries in upper
C triangular part.
C KNE holds number of entries in matrix after removal of
C out-of-range and duplicate entries.
C Note: this is similar to MC59ED except it also checks IP is
C monotonic and removes out-of-range entries in IRN.
C   . .
C   .. Scalar Arguments ..
    INTEGER ICNTL6,IDUP,IOUT,IUP,KNE,LA,LIP,LIW,NC,NR,NE
C   ..
C   .. Array Arguments ..
    DOUBLE PRECISION A(LA)
    INTEGER IRN(NE),IP(LIP),IW(LIW),INFO(2)
C   ..
C   .. Local Scalars ..
    INTEGER I,J,K,KSTART,KSTOP,NZJ,LOWER

    IDUP = 0
    IOUT = 0
    IUP = 0
    KNE = 0
C Initialise IW
    DO 10 I = 1,NR
      IW(I) = 0
   10 CONTINUE

    KSTART = IP(1)
    LOWER = 1
    IF (LA.GT.1) THEN
      NZJ = 0
      DO 30 J = 1,NC
C In symmetric case, entries out-of-range if they lie
C in upper triangular part.
      IF (ICNTL6.NE.0) LOWER = J
      KSTOP = IP(J+1)
      IF (KSTART.GT.KSTOP) THEN
        INFO(1) = -9
        INFO(2) = J
        RETURN
      END IF
      IP(J+1) = IP(J)
      DO 20 K = KSTART,KSTOP - 1
        I = IRN(K)
C Check for out-of-range
        IF (I.GT.NR .OR. I.LT.LOWER) THEN
          IOUT = IOUT + 1
C In symmetric case, check if entry is out-of-range because
C it lies in upper triangular part.
          IF (ICNTL6.NE.0 .AND. I.LT.J) IUP = IUP + 1
        ELSE IF (IW(I).LE.NZJ) THEN
          KNE = KNE + 1
          IRN(KNE) = I
          A(KNE) = A(K)
          IP(J+1) = IP(J+1) + 1
          IW(I) = KNE
        ELSE
C  We have a duplicate in column J
          IDUP = IDUP + 1
C If requested, sum duplicates
          IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
        END IF
   20   CONTINUE
      KSTART = KSTOP
      NZJ = KNE
   30   CONTINUE

    ELSE

C Pattern only
      DO 50 J = 1,NC
C In symmetric case, entries out-of-range if lie
C in upper triangular part.
      IF (ICNTL6.NE.0) LOWER = J
      KSTOP = IP(J+1)
      IF (KSTART.GT.KSTOP) THEN
        INFO(1) = -9
        INFO(2) = J
        RETURN
      END IF
      IP(J+1) = IP(J)
      DO  40 K = KSTART,KSTOP - 1
        I = IRN(K)
C Check for out-of-range
        IF (I.GT.NR .OR. I.LT.LOWER) THEN
          IOUT = IOUT + 1
          IF (ICNTL6.NE.0 .AND. I.GT.1) IUP = IUP + 1
        ELSE IF (IW(I).LT.J) THEN
          KNE = KNE + 1
          IRN(KNE) = I
          IP(J+1) = IP(J+1) + 1
          IW(I) = J
        ELSE
C  We have a duplicate in column J
          IDUP = IDUP + 1
        END IF
   40   CONTINUE
      KSTART = KSTOP
   50   CONTINUE
    END IF

    RETURN
    END

* COPYRIGHT (c) 1977 AEA Technology
* Original date 8 Oct 1992
C######8/10/92 Toolpack tool decs employed.
C######8/10/92 D version created by name change only.
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

    SUBROUTINE MC21A(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
C   .. Scalar Arguments ..
    INTEGER LICN,N,NUMNZ
C   ..
C   .. Array Arguments ..
    INTEGER ICN(LICN),IP(N),IPERM(N),IW(N,4),LENR(N)
C   ..
C   .. External Subroutines ..
    EXTERNAL MC21B
C   ..
C   .. Executable Statements ..
    CALL MC21B(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),IW(1,3),
   +       IW(1,4))
    RETURN
C
    END
    SUBROUTINE MC21B(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
C   PR(I) IS THE PREVIOUS ROW TO I IN THE DEPTH FIRST SEARCH.
C IT IS USED AS A WORK ARRAY IN THE SORTING ALGORITHM.
C   ELEMENTS (IPERM(I),I) I=1, ... N  ARE NON-ZERO AT THE END OF THE
C ALGORITHM UNLESS N ASSIGNMENTS HAVE NOT BEEN MADE.  IN WHICH CASE
C (IPERM(I),I) WILL BE ZERO FOR N-NUMNZ ENTRIES.
C   CV(I) IS THE MOST RECENT ROW EXTENSION AT WHICH COLUMN I
C WAS VISITED.
C   ARP(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED WHEN LOOKING FOR A CHEAP ASSIGNMENT.
C   OUT(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED DURING ONE PASS THROUGH THE MAIN LOOP.
C
C   INITIALIZATION OF ARRAYS.
C   .. Scalar Arguments ..
    INTEGER LICN,N,NUMNZ
C   ..
C   .. Array Arguments ..
    INTEGER ARP(N),CV(N),ICN(LICN),IP(N),IPERM(N),LENR(N),OUT(N),PR(N)
C   ..
C   .. Local Scalars ..
    INTEGER I,II,IN1,IN2,IOUTK,J,J1,JORD,K,KK
C   ..
C   .. Executable Statements ..
    DO 10 I = 1,N
      ARP(I) = LENR(I) - 1
      CV(I) = 0
      IPERM(I) = 0
   10 CONTINUE
    NUMNZ = 0
C
C
C   MAIN LOOP.
C   EACH PASS ROUND THIS LOOP EITHER RESULTS IN A NEW ASSIGNMENT
C OR GIVES A ROW WITH NO ASSIGNMENT.
    DO 100 JORD = 1,N
      J = JORD
      PR(J) = -1
      DO 70 K = 1,JORD
C LOOK FOR A CHEAP ASSIGNMENT
      IN1 = ARP(J)
      IF (IN1.LT.0) GO TO 30
      IN2 = IP(J) + LENR(J) - 1
      IN1 = IN2 - IN1
      DO 20 II = IN1,IN2
        I = ICN(II)
        IF (IPERM(I).EQ.0) GO TO 80
   20   CONTINUE
C   NO CHEAP ASSIGNMENT IN ROW.
      ARP(J) = -1
C   BEGIN LOOKING FOR ASSIGNMENT CHAIN STARTING WITH ROW J.
   30   CONTINUE
      OUT(J) = LENR(J) - 1
C INNER LOOP.  EXTENDS CHAIN BY ONE OR BACKTRACKS.
      DO 60 KK = 1,JORD
        IN1 = OUT(J)
        IF (IN1.LT.0) GO TO 50
        IN2 = IP(J) + LENR(J) - 1
        IN1 = IN2 - IN1
C FORWARD SCAN.
        DO 40 II = IN1,IN2
          I = ICN(II)
          IF (CV(I).EQ.JORD) GO TO 40
C   COLUMN I HAS NOT YET BEEN ACCESSED DURING THIS PASS.
          J1 = J
          J = IPERM(I)
          CV(I) = JORD
          PR(J) = J1
          OUT(J1) = IN2 - II - 1
          GO TO 70
C
   40     CONTINUE
C
C   BACKTRACKING STEP.
   50     CONTINUE
        J = PR(J)
        IF (J.EQ.-1) GO TO 100
   60   CONTINUE
C
   70   CONTINUE
C
C   NEW ASSIGNMENT IS MADE.
   80   CONTINUE
      IPERM(I) = J
      ARP(J) = IN2 - II - 1
      NUMNZ = NUMNZ + 1
      DO 90 K = 1,JORD
      J = PR(J)
      IF (J.EQ.-1) GO TO 100
      II = IP(J) + LENR(J) - OUT(J) - 2
      I = ICN(II)
      IPERM(I) = J
   90   CONTINUE
C
  100 CONTINUE
C
C   IF MATRIX IS STRUCTURALLY SINGULAR, WE NOW COMPLETE THE
C PERMUTATION IPERM.
    IF (NUMNZ.EQ.N) RETURN
    DO 110 I = 1,N
      ARP(I) = 0
  110 CONTINUE
    K = 0
    DO 130 I = 1,N
      IF (IPERM(I).NE.0) GO TO 120
      K = K + 1
      OUT(K) = I
      GO TO 130
C
  120   CONTINUE
      J = IPERM(I)
      ARP(J) = I
  130 CONTINUE
    K = 0
    DO 140 I = 1,N
      IF (ARP(I).NE.0) GO TO 140
      K = K + 1
      IOUTK = OUT(K)
      IPERM(IOUTK) = I
  140 CONTINUE
    RETURN
C
    END
* COPYRIGHT (c) 1976 AEA Technology
* Original date 21 Jan 1993
C 8 August 2000: CONTINUEs given to DOs.
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C
C 12th July 2004 Version 1.0.0. Version numbering added.

    SUBROUTINE MC22A(N,ICN,A,NZ,LENROW,IP,IQ,IW,IW1)
C   .. Scalar Arguments ..
    INTEGER N,NZ
C   ..
C   .. Array Arguments ..
    REAL A(NZ)
    INTEGER ICN(NZ),IP(N),IQ(N),IW(N,2),IW1(NZ),LENROW(N)
C   ..
C   .. Local Scalars ..
    REAL AVAL
    INTEGER I,ICHAIN,IOLD,IPOS,J,J2,JJ,JNUM,JVAL,LENGTH,NEWPOS
C   ..
C   .. Intrinsic Functions ..
    INTRINSIC IABS
C   ..
C   .. Executable Statements ..
    IF (NZ.LE.0) GO TO 1000
    IF (N.LE.0) GO TO 1000
C SET START OF ROW I IN IW(I,1) AND LENROW(I) IN IW(I,2)
    IW(1,1) = 1
    IW(1,2) = LENROW(1)
    DO 10 I = 2,N
      IW(I,1) = IW(I-1,1) + LENROW(I-1)
      IW(I,2) = LENROW(I)
   10 CONTINUE
C PERMUTE LENROW ACCORDING TO IP.  SET OFF-SETS FOR NEW POSITION
C   OF ROW IOLD IN IW(IOLD,1) AND PUT OLD ROW INDICES IN IW1 IN
C   POSITIONS CORRESPONDING TO THE NEW POSITION OF THIS ROW IN A/ICN.
    JJ = 1
    DO 20 I = 1,N
      IOLD = IP(I)
      IOLD = IABS(IOLD)
      LENGTH = IW(IOLD,2)
      LENROW(I) = LENGTH
      IF (LENGTH.EQ.0) GO TO 20
      IW(IOLD,1) = IW(IOLD,1) - JJ
      J2 = JJ + LENGTH - 1
      DO 15 J = JJ,J2
      IW1(J) = IOLD
   15   CONTINUE
      JJ = J2 + 1
   20 CONTINUE
C SET INVERSE PERMUTATION TO IQ IN IW(.,2).
    DO 30 I = 1,N
      IOLD = IQ(I)
      IOLD = IABS(IOLD)
      IW(IOLD,2) = I
   30 CONTINUE
C PERMUTE A AND ICN IN PLACE, CHANGING TO NEW COLUMN NUMBERS.
C
C ***   MAIN LOOP   ***
C EACH PASS THROUGH THIS LOOP PLACES A CLOSED CHAIN OF COLUMN INDICES
C   IN THEIR NEW (AND FINAL) POSITIONS ... THIS IS RECORDED BY
C   SETTING THE IW1 ENTRY TO ZERO SO THAT ANY WHICH ARE SUBSEQUENTLY
C   ENCOUNTERED DURING THIS MAJOR SCAN CAN BE BYPASSED.
    DO 200 I = 1,NZ
      IOLD = IW1(I)
      IF (IOLD.EQ.0) GO TO 200
      IPOS = I
      JVAL = ICN(I)
C IF ROW IOLD IS IN SAME POSITIONS AFTER PERMUTATION GO TO 150.
      IF (IW(IOLD,1).EQ.0) GO TO 150
      AVAL = A(I)
C **  CHAIN LOOP  **
C EACH PASS THROUGH THIS LOOP PLACES ONE (PERMUTED) COLUMN INDEX
C   IN ITS FINAL POSITION  .. VIZ. IPOS.
      DO 100 ICHAIN = 1,NZ
C NEWPOS IS THE ORIGINAL POSITION IN A/ICN OF THE ELEMENT TO BE PLACED
C IN POSITION IPOS.  IT IS ALSO THE POSITION OF THE NEXT ELEMENT IN
C   THE CHAIN.
      NEWPOS = IPOS + IW(IOLD,1)
C IS CHAIN COMPLETE ?
      IF (NEWPOS.EQ.I) GO TO 130
      A(IPOS) = A(NEWPOS)
      JNUM = ICN(NEWPOS)
      ICN(IPOS) = IW(JNUM,2)
      IPOS = NEWPOS
      IOLD = IW1(IPOS)
      IW1(IPOS) = 0
C **  END OF CHAIN LOOP  **
  100   CONTINUE
  130   A(IPOS) = AVAL
  150   ICN(IPOS) = IW(JVAL,2)
C ***   END OF MAIN LOOP   ***
  200 CONTINUE
C
 1000 RETURN

    END
* COPYRIGHT (c) 1987 AEA Technology
* Original date 10 Feb 1993
C     Toolpack tool decs employed.
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.

    SUBROUTINE MC34A(N,IRN,JCOLST,YESA,A,IW)
C THIS SUBROUTINE ACCEPTS AS INPUT THE STANDARD DATA STRUCTURE FOR
C   A SYMMETRIC MATRIX STORED AS A LOWER TRIANGLE AND PRODUCES
C   AS OUTPUT THE SYMMETRIC MATRIX HELD IN THE SAME DATA
C   STRUCTURE AS A GENERAL MATRIX.
C N IS AN INTEGER VARIABLE THAT MUST BE SET BY THE USER TO THE
C   ORDER OF THE MATRIX. NOT ALTERED BY THE ROUTINE
C   RESTRICTION (IBM VERSION ONLY): N LE 32767.
C IRN IS AN INTEGER (INTEGER*2 IN IBM VERSION) ARRAY THAT
C   MUST BE SET BY THE USER TO HOLD THE ROW INDICES OF THE LOWER
C   TRIANGULAR PART OF THE SYMMETRIC MATRIX.  THE ENTRIES OF A
C   SINGLE COLUMN ARE CONTIGUOUS. THE ENTRIES OF COLUMN J
C   PRECEDE THOSE OF COLUMN J+1 (J_=_1, ..., N-1), AND THERE IS
C   NO WASTED SPACE BETWEEN COLUMNS. ROW INDICES WITHIN A COLUMN
C   MAY BE IN ANY ORDER.  ON EXIT IT WILL HAVE THE SAME MEANING
C   BUT WILL BE CHANGED TO HOLD THE ROW INDICES OF ENTRIES IN
C   THE EXPANDED STRUCTURE.  DIAGONAL ENTRIES NEED NOT BE
C   PRESENT. THE NEW ROW INDICES ADDED IN THE UPPER TRIANGULAR
C   PART WILL BE IN ORDER FOR EACH COLUMN AND WILL PRECEDE THE
C   ROW INDICES FOR THE LOWER TRIANGULAR PART WHICH WILL REMAIN
C   IN THE INPUT ORDER.
C JCOLST IS AN INTEGER ARRAY OF LENGTH N+1 THAT MUST BE SET BY
C   THE USER SO THAT JCOLST(J) IS THE POSITION IN ARRAYS IRN AND
C   A OF THE FIRST ENTRY IN COLUMN J (J_=_1, ..., N).
C   JCOLST(N+1) MUST BE SET TO ONE MORE THAN THE TOTAL NUMBER OF
C   ENTRIES.  ON EXIT, JCOLST(J) WILL HAVE THE SAME MEANING BUT
C   WILL BE CHANGED TO POINT TO THE POSITION OF THE FIRST ENTRY
C   OF COLUMN J IN THE EXPANDED STRUCTURE. THE NEW VALUE OF
C   JCOLST(N+1) WILL BE ONE GREATER THAN THE NUMBER OF ENTRIES
C   IN THE EXPANDED STRUCTURE.
C YESA IS A LOGICAL VARIABLE THAT MUST BE SET TO .TRUE. IF THE
C   USER DESIRES TO GENERATE THE EXPANDED FORM FOR THE VALUES ALSO.
C   IF YESA IS .FALSE., THE ARRAY A WILL NOT BE REFERENCED.  IT IS
C   NOT ALTERED BY THE ROUTINE.
C A IS A REAL (DOUBLE PRECISION IN THE D VERSION) ARRAY THAT
C   CAN BE SET BY THE USER SO THAT A(K) HOLDS THE VALUE OF THE
C   ENTRY IN POSITION K OF IRN, {K = 1, _..._ JCOLST(N+1)-1}.
C   ON EXIT, IF YESA IS .TRUE., THE ARRAY WILL HOLD THE VALUES
C   OF THE ENTRIES IN THE EXPANDED STRUCTURE CORRESPONDING TO
C   THE OUTPUT VALUES OF IRN.   IF YESA IS .FALSE., THE ARRAY IS
C   NOT ACCESSED BY THE SUBROUTINE.
C IW IS AN INTEGER (INTEGER*2 IN IBM VERSION) ARRAY OF LENGTH
C   N THAT WILL BE USED AS WORKSPACE.
C
C CKP1 IS A LOCAL VARIABLE USED AS A RUNNING POINTER.
C OLDTAU IS NUMBER OF ENTRIES IN SYMMETRIC STORAGE.
C   .. Scalar Arguments ..
    INTEGER N
    LOGICAL YESA
C   ..
C   .. Array Arguments ..
    REAL A(*)
    INTEGER IRN(*),IW(*),JCOLST(*)
C   ..
C   .. Local Scalars ..
    INTEGER CKP1,I,I1,I2,II,IPKP1,IPOS,J,JSTART,LENK,NDIAG,NEWTAU,
   +      OLDTAU
C   ..
C   .. Executable Statements ..
C
    OLDTAU = JCOLST(N+1) - 1
C INITIALIZE WORK ARRAY
    DO 5 I = 1,N
      IW(I) = 0
    5 CONTINUE
C
C IW(J) IS SET EQUAL TO THE TOTAL NUMBER OF ENTRIES IN COLUMN J
C   OF THE EXPANDED SYMMETRIC MATRIX.
C NDIAG COUNTS NUMBER OF DIAGONAL ENTRIES PRESENT
    NDIAG = 0
    DO 20 J = 1,N
      I1 = JCOLST(J)
      I2 = JCOLST(J+1) - 1
      IW(J) = IW(J) + I2 - I1 + 1
      DO 10 II = I1,I2
      I = IRN(II)
      IF (I.NE.J) THEN
        IW(I) = IW(I) + 1

      ELSE
        NDIAG = NDIAG + 1
      END IF

   10   CONTINUE
   20 CONTINUE
C
C NEWTAU IS NUMBER OF ENTRIES IN EXPANDED STORAGE.
    NEWTAU = 2*OLDTAU - NDIAG
C IPKP1 POINTS TO POSITION AFTER END OF COLUMN BEING CURRENTLY
C   PROCESSED
    IPKP1 = OLDTAU + 1
C CKP1 POINTS TO POSITION AFTER END OF SAME COLUMN IN EXPANDED
C   STRUCTURE
    CKP1 = NEWTAU + 1
C GO THROUGH THE ARRAY IN THE REVERSE ORDER PLACING LOWER TRIANGULAR
C   ELEMENTS IN THE APPROPRIATE SLOTS.
    DO 40 J = N,1,-1
      I1 = JCOLST(J)
      I2 = IPKP1
C LENK IS NUMBER OF ENTRIES IN COLUMN J OF ORIGINAL STRUCTURE
      LENK = I2 - I1
C JSTART IS RUNNING POINTER TO POSITION IN NEW STRUCTURE
      JSTART = CKP1
C SET IKP1 FOR NEXT COLUMN
      IPKP1 = I1
      I2 = I2 - 1
C RUN THROUGH COLUMNS IN REVERSE ORDER
C LOWER TRIANGULAR PART OF COLUMN MOVED TO END OF SAME COLUMN IN
C   EXPANDED FORM
      DO 30 II = I2,I1,-1
      JSTART = JSTART - 1
      IF (YESA) A(JSTART) = A(II)
      IRN(JSTART) = IRN(II)
   30   CONTINUE
C JCOLST IS SET TO POSITION OF FIRST ENTRY IN LOWER TRIANGULAR PART OF
C   COLUMN J IN EXPANDED FORM
      JCOLST(J) = JSTART
C SET CKP1 FOR NEXT COLUMN
      CKP1 = CKP1 - IW(J)
C RESET IW(J) TO NUMBER OF ENTRIES IN LOWER TRIANGLE OF COLUMN.
      IW(J) = LENK
   40 CONTINUE
C
C AGAIN SWEEP THROUGH THE COLUMNS IN THE REVERSE ORDER, THIS
C   TIME WHEN ONE IS HANDLING COLUMN J THE UPPER TRIANGULAR
C   ELEMENTS A(J,I) ARE PUT IN POSITION.
    DO 80 J = N,1,-1
      I1 = JCOLST(J)
      I2 = JCOLST(J) + IW(J) - 1
C RUN DOWN COLUMN IN ORDER
C NOTE THAT I IS ALWAYS GREATER THAN OR EQUAL TO J
      DO 60 II = I1,I2
      I = IRN(II)
      IF (I.EQ.J) GO TO 60
      JCOLST(I) = JCOLST(I) - 1
      IPOS = JCOLST(I)
      IF (YESA) A(IPOS) = A(II)
      IRN(IPOS) = J
   60   CONTINUE
   80 CONTINUE
    JCOLST(N+1) = NEWTAU + 1
    RETURN

    END
* COPYRIGHT (c) 1988 AEA Technology
* Original date 17 Feb 2005

C 17th February 2005 Version 1.0.0. Replacement for FD05.

    REAL FUNCTION FD15A(T)
C----------------------------------------------------------------
C  Fortran 77 implementation of the Fortran 90 intrinsic
C    functions: EPSILON, TINY, HUGE and RADIX.  Note that
C    the RADIX result is returned as REAL.
C
C  The CHARACTER argument specifies the type of result:
C     
C   'E'  smallest positive real number: 1.0 + FD15A > 1.0, i.e.
C      EPSILON(REAL)
C   'T'  smallest full precision positive real number, i.e.
C      TINY(REAL)
C   'H'  largest finite positive real number, i.e.
C      HUGE(REAL)
C   'R'  the base of the floating point arithematic, i.e.
C      RADIX(REAL)
C
C    any other value gives a result of zero.
C----------------------------------------------------------------
    CHARACTER T

    IF ( T.EQ.'E' ) THEN
       FD15A = EPSILON(1.0)
    ELSE IF ( T.EQ.'T' ) THEN
       FD15A = TINY(1.0)
    ELSE IF ( T.EQ.'H' ) THEN
       FD15A = HUGE(1.0)
    ELSE IF ( T.EQ.'R' ) THEN
       FD15A = REAL(RADIX(1.0))
    ELSE
       FD15A = 0.0
    ENDIF
    RETURN
    END
C COPYRIGHT (c) 1995 Timothy A. Davis, Patrick Amestoy and
C       Council for the Central Laboratory of the Research Councils
C Original date 30 November 1995
C  April 2001: call to MC49 changed to MC59 to make routine threadsafe
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.
C 23 May 2007 Version 1.1.0. Absolute value of hash taken to cover the
C        case of integer overflow.
C        Comments with character in column 2 corrected.
C 2 August 2007 Version 2.0.0 Dense row handling added, error & warning
C        messages added, iovflo added, interface changed. MC47I/ID
C        added.
C 31 October 2007 Version 2.1.0 Corrected tree formation when handling
C        full variables
C

    SUBROUTINE MC47I(ICNTL)
    INTEGER ICNTL(10)
C ICNTL is an INTEGER array of length 10 that contains control
C   parameters and must be set by the user. Default values are set
C   by MA57ID.
C
C   ICNTL(1) is the stream number for error messages. Printing is
C   suppressed if ICNTL(1)<0. The default is 6.
C
C   ICNTL(2) is the stream number for warning messages. Printing is
C   suppressed if ICNTL(2)<0. The default is 6.
C
C   ICNTL(3) is the stream number for printing matrix data on entry
C   to and exit from MC47A/AD. Printing is suppressed if ICNTL(3)<0.
C   The default is -1.
C
C   ICNTL(4) controls the choice of AMD algorithm
C   =-1 Classical MC47B (AMD) algorithm (no dense row detection)
C   = 0 Only exactly dense rows in the reduced matrix are selected.
C   = 1 Corresponds to automatic setting of the minimum density
C   requirement.
C   The default value is 1.
C
C   ICNTL(5) defines the largest positive
C   integer that your computer can represent (-iovflo should also
C   be representable). HUGE(1) in Fortran 95. The default value is
C   2139062143
C
C   ICNTL(6) to ICNTL(10) are set to zero by MA57ID but are not
C   currently used by MC47.
C
C Local variables
    INTEGER I

    ICNTL(1) = 6
    ICNTL(2) = 6
    ICNTL(3) = -1
    ICNTL(4) = 1
    ICNTL(5) = 2139062143

    DO 100 I=6,10
      ICNTL(I) = 0
 100  CONTINUE
    RETURN
    END


    SUBROUTINE MC47A(N, NE, PE, IW, IWLEN,
   *    ICNTL,INFO, RINFO)
    INTEGER N, NE, PE(N+1), IWLEN, IW(IWLEN), INFO(10)
    INTEGER ICNTL(10)
    REAL  RINFO(10)
C N is an INTEGER variable that must be set by the user to the
C   order of the matrix A.  It is not altered by the subroutine.
C   Restriction: N >= 1.
C NE is an INTEGER variable that must be set by the user to the
C   number of entries in the matrix A.
C   It is not altered by the subroutine.
C PE is an INTEGER  array of size N+1 that must be set by the user.
C   If the user is supplying the entries by columns, then PE(i) must
C   hold the index in IW of the start of column i, i=1, ..., N and
C   PE(N+1) must be equal to NE+1.  If the user is supplying row and
C   column indices for the matrix, PE(1) must be negative.
C   On exit, PE will hold information on the matrix factors.
C IW is an INTEGER  array of length IWLEN that must be set by the user
C   to hold the pattern of the matrix A.  If PE(1) is positive,
C   IW(PE(J)), ..., IW(PE(J+1)-1) must hold the row indices of entries
C   in column J, J = 1, ..., N.
C   The entries within a column need not be in order.
C   If PE(1) is negative, then (IW(k), IW(NE+k)), k = 1, ..., NE, must
C   hold the row and column index of an entry.
C   Duplicates, out-of-range entries, diagonal entries, and entries
C   in upper triangle are ignored.
C   IW is used as workspace by the subroutine.  On exit, the
C   permutation generated by MC47 is held in IW(IWLEN-I+1), I=1, ... N
C   and is such that the kth column of the permuted matrix is column
C   IW(IWLEN-N+k) of the original matrix.
C   The inverse permutation is held in positions IWLEN-2N+1 to
C   IWLEN-N of IW, preceded by further information on the structure
C   of the factors
C IWLEN is an INTEGER variable. It must be set by the user to the length
C   of array IW and is not altered by the subroutine.
C   We recommend IWLEN > 2NE + 9N.   Restriction: IWLEN >= 2NE + 8N.
C ICNTL is an INTEGER array of length 10 that contains control
C   parameters and must be set by the user. Default values for the
C   components may be set by a call to MC47I.
C INFO is an INTEGER array of length 8 that need not be set by the user.
C   On return from MC47A/AD, a value of zero for INFO(1) indicates
C   that the subroutine has performed successfully.  Negative values
C   for INFO(1) signify a fatal error.  Possible values are:
C   -1    N < 1
C   -2    IWLEN < 2NE + 8N.
C   -3    Error in PE when input is by columns.
C   -4    Matrix is null (usually because all entries in upper
C       triangle.
C   There is one warning indicated by a positive value for INFO(1)
C   +1    Out-of-range index, duplicate, diagonal or entry in upper
C       triangle in input. Action taken is to ignore these entries.
C   The other entries of INFO give information to the user.
C   INFO(2) gives the number of compresses performed on the array IW.
C         A large value for this indicates that the ordering could
C         be found more quickly if IWLEN were increased.
C   INFO(3) gives the minimum necessary value for IWLEN for a
C         successful run of MC47A/AD on the same matrix as has just
C         been analysed, without the need for any compresses of IW.
C   INFO(4) gives the number of entries with row or column indices
C         that are out of range.  Any such entry is ignored.
C   INFO(5) gives the number of duplicate entries.  Any such entry is
C         ignored.
C   INFO(6) gives the number of entries in upper triangle.  Any such
C         entry is ignored.
C   INFO(7) gives the number of diagonal entries.  Any such entry is
C         ignored.
C   INFO(8) gives the number of restarts performed.
C RINFO is an INTEGER array of length 10 that need not be set by the
C   user. The other entries of RINFO give information to the user.
C   RINFO(1) gives forecast number of reals to hold the factorization
C   RINFO(2) gives the forecast number of flops required by the
C   factorization if no pivoting is performed.

C Local variables
    INTEGER DEGREE
    REAL  DUMMY(1)
    INTEGER ELEN,HEAD,I,II,I1,I2,J,LAST,LEN,LENIW,LP,MP,
   *      NEXT,NV,PFREE,W,WP
    INTEGER ICT59(10),INFO59(10),IOUT,JOUT,IDUP,JNFO(10)
C DEGREE is used to subdivide array IW
C DUMMY  is dummy real for call to MC34A/AD
C ELEN   is used to subdivide array IW
C HEAD   is used to subdivide array IW
C I    is DO loop variable
C IFLAG  is eror return from MC59A/AD
C II   is running index for entries of IW in column
C I1   is first location for entries of IW in column
C I2   is flast location for entries of IW in column
C J    is column index and DO loop variable
C LAST   is used to subdivide array IW
C LEN    is used to subdivide array IW
C LENIW  is space left in IW after allocation of work vectors
C LP   is a local copy of ICNTL(1)
C MP   is a local copy of ICNTL(3)
C NEXT   is used to subdivide array IW
C NV   is used to subdivide array IW
C PFREE  marks active length of IW for call to MC47B/BD
C W    is used to subdivide array IW
C WP   is a local copy of ICNTL(2)
C IOUT   number of indices out-of-range
C JOUT   number of column indices out-of-range (as detected by MC59A/AD)
C IDUP   number of duplicates


C Subroutines called
    EXTERNAL MC59A,MC34A,MC47B

C Initialize info
    DO 5 J = 1,10
       INFO(J) = 0
 5    CONTINUE

C Set LP,WP,MP
    LP = ICNTL(1)
    WP = ICNTL(2)
    MP = ICNTL(3)

C Check value of N
    IF (N.LT.1) THEN
      INFO(1) = -1
      IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10)')
   +     '**** Error return from MC47A **** INFO(1) =',INFO(1),
   +     'N has value ',N
      GO TO 1000
    ENDIF

C Error check to see if enough space in IW to get started
    IF (PE(1).LT.1) THEN
      IF (2*NE+N.GT.IWLEN) THEN
      INFO(1) = -2
       IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10)')
   +      '**** Error return from MC47A **** INFO(1) =',INFO(1),
   +      'IWLEN has value ',IWLEN
      GO TO 1000
      ENDIF
    ELSE
      IF (NE+N.GT.IWLEN) THEN
      INFO(1) = -2
       IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10)')
   +      '**** Error return from MC47A **** INFO(1) =',INFO(1),
   +      'IWLEN has value ',IWLEN
      GO TO 1000
      ENDIF
    ENDIF

C Diagnostic print
    IF (MP.GE.0) THEN
      WRITE(MP,'(/A)') 'Entry to MC47A/AD'
      WRITE(MP,'(A,I10,A,I10,A)') 'Matrix of order',N,' with',NE,
   *                  ' entries'
      IF (PE(1).LT.0)  THEN
      WRITE(MP,'(A)') 'Matrix input in coordinate form'
      WRITE(MP,'(A/(4(I8,I8)))') 'Row and column indices',
   *      (IW(I),IW(NE+I),I=1,NE)
      ELSE
      WRITE(MP,'(A)') 'Matrix input by columns'
      DO 10 J=1,N
        WRITE(MP,'(A,I4/(10I8))') 'Column',J,
   *                    (IW(I),I=PE(J),PE(J+1)-1)
   10   CONTINUE
      ENDIF
    ENDIF

C Divide workspace
    LAST   = IWLEN  - N + 1
    ELEN   = LAST   - N
    NV   = ELEN   - N
    W    = NV   - N
    DEGREE = W    - N
    HEAD   = DEGREE - N
    NEXT   = HEAD   - N
    LEN    = NEXT   - N
    LENIW = LEN-1

C Set counters for number of upper triangular entries and diagonals
C   present in input matrix.
C These will be removed for later processing.
    INFO(6) = 0
    INFO(7) = 0

    IF (PE(1).LT.0) THEN
C First remove diagonals (if present) and all entries in upper triangle
C Note that this may give out-of-range entries from MC59.
      DO 20 I=1,NE
      IF (IW(I).LE.IW(NE+I)) THEN
        IF (IW(I).EQ.IW(NE+I) .AND. IW(I).NE.0) THEN
          INFO(7) = INFO(7) + 1
        ELSE
          IF (IW(I).GT.0) INFO(6) = INFO(6) + 1
        ENDIF
        IW(I)=0
      ENDIF
   20   CONTINUE

C Call sort routine
      ICT59(1) = 0
      ICT59(2) = 1
      ICT59(3) = 1
      ICT59(4) = LP
      ICT59(5) = -1
      ICT59(6) = 0
      CALL MC59A(ICT59,N,N,NE,IW,NE,IW(NE+1),1,DUMMY,
   *          N+1,PE,N+1,IW(2*NE+1),INFO59)
C      IFLAG = INFO59(1)
      IDUP  = INFO59(3)
      IOUT  = INFO59(4)
      JOUT  = INFO59(5)
    ELSE

C Matrix already sorted by columns.
C Remove duplicates, out-of-range indices, and entries in upper
C     triangle.  First initialize counts.
      IDUP = 0
      IOUT = 0
      JOUT = 0

C Set array used to find duplicates
      DO 30 I = 1,N
      IW(NE+I) = 0
   30   CONTINUE

      DO 50 J=1,N
      I1 = PE(J)
      PE(J) = I1-(IOUT+IDUP)
      I2 = PE(J+1)-1
      IF (I2.LT.I1-1) THEN
        INFO(1) = -3
        GO TO 1000
      ENDIF
      DO 40 II = I1,I2
        I = IW(II)
        IF (I.LE.J .OR. I.GT.N) THEN
          IF (I.EQ.J) INFO(7) = INFO(7) + 1
          IF (I.GT.0 .AND. I.LT.J) INFO(6) = INFO(6) + 1
          IOUT = IOUT + 1
        ELSE
          IF (IW(NE+I).EQ.J) THEN
C Duplicate found
          IDUP = IDUP + 1
          ELSE
          IW(NE+I)=J
          IW(II-(IOUT+IDUP)) = I
          ENDIF
        ENDIF
   40   CONTINUE
   50   CONTINUE
      PE(N+1) = NE - (IOUT+IDUP) + 1
    ENDIF

C Set flags for duplicates or out-of-range entries
C Check if there were duplicates
    IF (IDUP.GT.0) THEN
      INFO(1) = 1
      INFO(4) = IDUP
      IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
   +     '**** Warning from MC47A **** INFO(1) =',INFO(1),
   +     'Number of duplicates found: ',INFO(4)
    ELSE
      INFO(4) = 0
    ENDIF
C Check for out of range entries
    IF (IOUT+ JOUT - INFO(7) .GT.0 ) THEN
      INFO(1) = 1
      INFO(5) = IOUT + JOUT - INFO(7)
      IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
   +     '**** Warning from MC47A **** INFO(1) =',INFO(1),
   +     'Number of out of range entries found and ignored: ',
   +     INFO(5)
    ELSE
      INFO(5) = 0
    ENDIF

C Check for entries in upper triangle
    IF (INFO(6).GT.0) THEN
       INFO(1) = 1
      IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
   +     '**** Warning from MC47A **** INFO(1) =',INFO(1),
   +     'Number of entries in upper triangle found and ignored: ',
   +      INFO(6)
    ENDIF

C Check for entries in diagonals
    IF (INFO(7).GT.0) THEN
       INFO(1) = 1
      IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
   +     '**** Warning from MC47A **** INFO(1) =',INFO(1),
   +     'Number of entries in diagonals found and ignored: ',
   +      INFO(7)
    ENDIF

C Check for null matrix
C Usually happens if wrong triangle is input
    IF (NE-(IOUT+IDUP).EQ.0) THEN
      INFO(1) = -4
      IF (LP.GE.0) WRITE(LP,'(/A,I3/A)')
   +     '**** Error return from MC47A **** INFO(1) =',INFO(1),
   +     'Matrix is null'
      GO TO 1000
    ENDIF

C Generate matrix in expanded form
C First check there is sufficient space in IW
    IF (LENIW.LT.2*(PE(N+1)-1)) THEN
      INFO(1) = -2
       IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10/A,I10)')
   +      '**** Error return from MC47A **** INFO(1) =',INFO(1),
   +      'IWLEN has value ',IWLEN,
   +      'Should be at least', 2*(PE(N+1)-1)+8*N
      GO TO 1000
    ENDIF

    CALL MC34A(N,IW,PE,.FALSE.,DUMMY,IW(W))
    PFREE = PE(N+1)

C Set length array for MC47B/BD
    DO 60 I=1,N
      IW(LEN+I-1) = PE(I+1) - PE(I)
   60 CONTINUE

C Call to approximate minimum degree subroutine.
    CALL MC47B(N,LENIW,PE,PFREE,IW(LEN),IW,IW(NV),
   *        IW(ELEN),IW(LAST),IW(DEGREE),
   *        IW(HEAD),IW(NEXT),IW(W), ICNTL,JNFO, RINFO)

    INFO(2) = JNFO(1)
    INFO(3) = PFREE+8*N
    INFO(8) = JNFO(2)

C Print diagnostics
    IF (MP.GE.0) THEN
      WRITE(MP,'(/A)') 'Exit from MC47A/AD'
      WRITE(MP,'(A/(7I10))') 'INFO(1-10):',(INFO(I),I=1,10)
      WRITE(MP,'(A/(8I10))') 'Parent array',(PE(I),I=1,N)
      WRITE(MP,'(A/(8I10))') 'Permutation',(IW(ELEN+I-1),I=1,N)
      WRITE(MP,'(A/(8I10))') 'Inverse permutation',
   *               (IW(LAST+I-1),I=1,N)
      WRITE(MP,'(A/(8I10))') 'Degree array',(IW(NV+I-1),I=1,N)
    ENDIF

 1000 RETURN
    END

C ====================================================================
C ====================================================================
C ====================================================================

    SUBROUTINE MC47B (N, IWLEN, PE, PFREE, LEN, IW, NV,
   $             ELEN, LAST, DEGREE,
   $             HEAD, DENXT, W, ICNTL, JNFO, RJNFO)

    INTEGER N, IWLEN, PE(N), PFREE, LEN(N), IW(IWLEN), NV(N),
   $      ELEN(N), LAST(N),  DEGREE(N),
   $       HEAD(N), DENXT(N), W(N), ICNTL(10), JNFO(10)

    REAL  RJNFO(10)

C -------------------------------------------------------------------
C AMDD is a modified version of
C   MC47B:  Approximate Minimum (UMFPACK/MA38-style, external) Degree
C      ordering algorithm, with aggresive absorption
C designed to automatically detect and exploit dense
C rows in the reduced matrix at any step of the minimum degree.
C
C We use the term Le to denote the set of all supervariables in element
C E.
C **** Reword below***
C A row is declared as full if none of its entries can be guaranteed
C   to be zero.
C A row is quasi dense if at most N-THRESM-1 of its entries can be
C   guaranteed to be zero at it has not been recognized as being full
C   in the calculation so far.
C A row is dense if it is either full or quasi dense.
C A row is sparse if it is not dense.
C -------------------------------------------------------------------
C
C
C N must be set to the matrix order. It is not altered.
C   Restriction:  N .ge. 1
C
C IWLEN must be set to the length of IW. It is not altered. On input,
C   the matrix is stored in IW (1..PFREE-1).
C   *** We do not recommend running this algorithm with ***
C   ***    IWLEN .LT. PFREE + N.              ***
C   *** Better performance will be obtained if      ***
C   ***    IWLEN .GE. PFREE + N               ***
C   *** or better yet                     ***
C   ***    IWLEN .GT. 1.2 * PFREE             ***
C   Restriction: IWLEN .GE. PFREE-1
C
C PE(i) must be set to the the index in IW of the start of row I, or be
C   zero if row I has no off-diagonal entries. During execution,
C   it is used for both supervariables and elements:
C     * Principal supervariable I:  index into IW of the
C         list of supervariable I.  A supervariable
C         represents one or more rows of the matrix
C         with identical pattern.
C     * Non-principal supervariable I:  if I has been absorbed
C         into another supervariable J, then PE(I) = -J.
C         That is, J has the same pattern as I.
C         Note that J might later be absorbed into another
C         supervariable J2, in which case PE(I) is still -J,
C         and PE(J) = -J2.
C     * Unabsorbed element E:  the index into IW of the list
C         of element E.  Element E is created when
C         the supervariable of the same name is selected as
C         the pivot.
C     * Absorbed element E:  if element E is absorbed into element
C         E2, then PE(E) = -E2.  This occurs when one of its
C         variables is eliminated and when the pattern of
C         E (that is, Le) is found to be a subset of the pattern
C         of E2 (that is, Le2).  If element E is "null" (it has
C         no entries outside its pivot block), then PE(E) = 0.
C
C   On output, PE holds the assembly tree/forest, which implicitly
C   represents a pivot order with identical fill-in as the actual
C   order (via a depth-first search of the tree). If NV(I) .GT. 0,
C   then I represents a node in the assembly tree, and the parent of
C   I is -PE(I), or zero if I is a root. If NV(I)=0, then (I,-PE(I))
C   represents an edge in a subtree, the root of which is a node in
C   the assembly tree.
C
C PFREE must be set to the position in IW of the first free variable.
C   During execution, additional data is placed in IW, and PFREE is
C   modified so that components  of IW from PFREE are free.
C   On output, PFREE is set equal to the size of IW that would have
C   caused no compressions to occur.  If NCMPA is zero, then
C   PFREE (on output) is less than or equal to IWLEN, and the space
C   IW(PFREE+1 ... IWLEN) was not used. Otherwise, PFREE (on output)
C   is greater than IWLEN, and all the memory in IW was used.
C
C LEN(I) must be set to hold the number of entries in row I of the
C   matrix, excluding the diagonal.  The contents of LEN(1..N) are
C   undefined on output.
C
C IW(1..PFREE-1) must be set to  hold the patterns of the rows of
C   the matrix.  The matrix must be symmetric, and both upper and
C   lower triangular parts must be present.  The diagonal must not be
C   present.  Row I is held as follows:
C         IW(PE(I)...PE(I) + LEN(I) - 1) must hold the list of
C         column indices for entries in row I (simple
C         supervariables), excluding the diagonal.  All
C         supervariables start with one row/column each
C         (supervariable I is just row I). If LEN(I) is zero on
C         input, then PE(I) is ignored on input. Note that the
C         rows need not be in any particular order, and there may
C         be empty space between the rows.
C   During execution, the supervariable I experiences fill-in. This
C   is represented by constructing a list of the elements that cause
C   fill-in in supervariable I:
C         IE(PE(i)...PE(I) + ELEN(I) - 1) is the list of elements
C         that contain I. This list is kept short by removing
C         absorbed elements. IW(PE(I)+ELEN(I)...PE(I)+LEN(I)-1)
C         is the list of supervariables in I. This list is kept
C         short by removing nonprincipal variables, and any entry
C         J that is also contained in at least one of the
C         elements in the list for I.
C   When supervariable I is selected as pivot, we create an element E
C   of the same name (E=I):
C         IE(PE(E)..PE(E)+LEN(E)-1) is the list of supervariables
C          in element E.
C   An element represents the fill-in that occurs when supervariable
C   I is selected as pivot.
C   CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
C   The contents of IW are undefined on output.
C
C NV(I) need not be set. During execution, ABS(NV(I)) is equal to the
C   number of rows represented by the principal supervariable I. If I
C   is a nonprincipal variable, then NV(I) = 0. Initially, NV(I) = 1
C   for all I.  NV(I) .LT. 0 signifies that I is a principal variable
C   in the pattern Lme of the current pivot element ME. On output,
C   NV(E) holds the true degree of element E at the time it was
C   created (including the diagonal part).
C
C ELEN(I) need not be set. See the description of IW above. At the
C   start of execution, ELEN(I) is set to zero. For a supervariable,
C   ELEN(I) is the number of elements in the list for supervariable
C   I. For an element, ELEN(E) is the negation of the position in the
C   pivot sequence of the supervariable that generated it. ELEN(I)=0
C   if I is nonprincipal.
C   On output ELEN(1..N) holds the inverse permutation (the same
C   as the 'INVP' argument in Sparspak). That is, if K = ELEN(I),
C   then row I is the Kth pivot row.  Row I of A appears as the
C   (ELEN(I))-th row in the permuted matrix, PAP^T.
C
C LAST(I) need not be set on input. In a degree list, LAST(I) is the
C   supervariable preceding I, or zero if I is the head of the list.
C   In a hash bucket, LAST(I) is the hash key for I. LAST(HEAD(HASH))
C   is also used as the head of a hash bucket if HEAD(HASH) contains
C   a degree list (see HEAD, below).
C   On output, LAST(1..N) holds the permutation (the same as the
C   'PERM' argument in Sparspak). That is, if I = LAST(K), then row I
C   is the Kth pivot row.  Row LAST(K) of A is the K-th row in the
C   permuted matrix, PAP^T.
C
C
C DEGREE need not be set on input. If I is a supervariable and sparse,
C   then DEGREE(I) holds the current approximation of the external
C   degree of row I (an upper bound). The external degree is the
C   number of entries in row I, minus ABS(NV(I)) (the diagonal
C   part). The bound is equal to the external degree if ELEN(I) is
C   less than or equal to two. We also use the term "external degree"
C   for elements E to refer to |Le \ Lme|. If I is full in the reduced
C   matrix, then DEGREE(I)=N+1. If I is dense in the reduced matrix,
C   then DEGREE(I)=N+1+last_approximate_external_deg of I.
C   All dense rows are stored in the list pointed by HEAD(N).
C   Quasi dense rows are stored first, and are followed by full rows
C   in the reduced matrix. LASTD holds the last row in
C   this list of dense rows or is zero if the list is empty.
C
C HEAD(DEG) need not be set on input. HEAD is used for degree lists.
C   HEAD(DEG) is the first supervariable in a degree list (all
C   supervariables I in a degree list DEG have the same approximate
C   degree, namely, DEG = DEGREE(I)). If the list DEG is empty then
C   HEAD(DEG) = 0.
C   During supervariable detection HEAD(HASH) also serves as a
C   pointer to a hash bucket.
C   If HEAD(HASH) .GT. 0, there is a degree list of degree HASH. The
C   hash bucket head pointer is LAST(HEAD(HASH)).
C   If HEAD(HASH) = 0, then the degree list and hash bucket are
C   both empty.
C   If HEAD(HASH) .LT. 0, then the degree list is empty, and
C   -HEAD(HASH) is the head of the hash bucket.
C   After supervariable detection is complete, all hash buckets are
C   empty, and the (LAST(HEAD(HASH)) = 0) condition is restored for
C   the non-empty degree lists.
C
C DENXT(I) need not be set on input. For supervariable I, DENXT(I) is
C   the supervariable following I in a link list, or zero if I is
C   the last in the list. Used for two kinds of lists: degree lists
C   and hash buckets (a supervariable can be in only one kind of
C   list at a time). For element E, DENXT(E) is the number of
C   variables with dense or full rows in the element E.
C
C W(I) need not be set on input. The flag array W determines the status
C   of elements and variables, and the external degree of elements.
C   For elements:
C      if W(E) = 0, then the element E is absorbed.
C      if W(E) .GE. WFLG, then W(E)-WFLG is the size of the set
C         |Le \ Lme|, in terms of nonzeros (the sum of ABS(NV(I))
C         for each principal variable I that is both in the
C         pattern of element E and NOT in the pattern of the
C         current pivot element, ME).
C      if WFLG .GT. WE(E) .GT. 0, then E is not absorbed and has
C         not yet been seen in the scan of the element lists in
C         the computation of |Le\Lme| in loop 150 below.
C         ***SD: change comment to remove reference to label***
C   For variables:
C      during supervariable detection, if W(J) .NE. WFLG then J is
C      not in the pattern of variable I.
C   The W array is initialized by setting W(I) = 1 for all I, and by
C   setting WFLG = 2. It is reinitialized if WFLG becomes too large
C   (to ensure that WFLG+N does not cause integer overflow).
C
C ICNTL is an INTEGER array of length 10 that contains control
C   parameters and must be set by the user. Default values for the
C   components may be set by a call to MC47I.
C
C RJNFO is an REAL (REAL  in  D version) array of length 7
C   that need not be set by the user. This array supplies information
C   on the execution of MC47B.
C   RJNFO(1) gives forecast number of reals to hold the factorization
C   RJNFO(2) gives the forecast number of flops required by the
C   factorization if no pivoting is performed.
C
C Local variables:
C ---------------
    INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, HASH, HMOD, I,
   $   IDUMMY, ILAST, INEXT, IOVFLO,J, JDUMMY, JLAST, JNEXT, K,
   $   KNT1, KNT2, KNT3, LASTD,  LENJ, LN, MAXMEM, ME,
   $   MEM, MINDEG, NBD, NCMPA, NDME, NEL, NELME, NEWMEM,
   $   NFULL, NLEFT, NRLADU, NVI, NVJ, NVPIV, P, P1, P2, P3, PDST,
   $   PEE, PEE1, PEND, PJ, PME, PME1, PME2, PN, PSRC, RSTRT,
   $   SLENME, THRESH, THRESM, WE, WFLG, WNVI,X
   $
    REAL  RELDEN, SM, STD, OPS
    LOGICAL IDENSE
C
C DEG:      the degree of a variable or element
C DEGME:    size (no. of variables), |Lme|, of the current element,
C         ME (= DEGREE(ME))
C DEXT:     external degree, |Le \ Lme|, of some element E
C DMAX:     largest |Le| seen so far
C E:      an element
C ELENME:   the length, ELEN(ME), of element list of pivotal var.
C ELN:      the length, ELEN(...), of an element list
C HASH:     the computed value of the hash function
C HMOD:     the hash function is computed modulo HMOD = MAX(1,N-1)
C I:      a supervariable
C IDUMMY:   loop counter
C ILAST:    the entry in a link list preceding I
C INEXT:    the entry in a link list following I
C IOVFLO:   local copy of ICNTL(5)
C J:      a supervariable
C JDUMMY:   loop counter
C JLAST:    the entry in a link list preceding J
C JNEXT:    the entry in a link list, or path, following J
C K:      the pivot order of an element or variable
C KNT1:     loop counter used during element construction
C KNT2:     loop counter used during element construction
C KNT3:     loop counter used during element construction
C LASTD:    index of the last row in the list of dense rows
C LENJ:     LEN(J)
C LN:       length of a supervariable list
C MAXMEM:   amount of memory needed for no compressions
C ME:       current supervariable being eliminated, and the
C             current element created by eliminating that
C             supervariable
C MEM:      memory in use assuming no compressions have occurred
C MINDEG:   current approximate minimum degree
C NBD:      total number of dense rows selected
C NCMPA:    counter for the number of times IW was compressed
C NDME  :   number of dense rows adjacent to me
C NEL:      number of pivots selected so far
C NELME:    number of pivots selected when reaching the root
C NEWMEM:   amount of new memory needed for current pivot element
C NFULL:    total number of full rows detected.
C NLEFT:    N-NEL, the number of nonpivotal rows/columns remaining
C NRLADU:   counter for the forecast number of reals in matrix factor
C NVI:      the number of variables in a supervariable I (= NV(I))
C NVJ:      the number of variables in a supervariable J (= NV(J))
C NVPIV:    number of pivots in current element
C P:      pointer into lots of things
C P1:       pe (i) for some variable i (start of element list)
C P2:       pe (i) + elen (i) -  1 for some var. i (end of el. list)
C P3:       index of first supervariable in clean list
C PJ:       pointer into an element or variable
C PDST:     destination pointer, for compression
C PEE:      pointer into element E
C PEE1:     pointer into element E
C PEND:     end of memory to compress
C PME:      pointer into the current element (PME1...PME2)
C PME1:     the current element, ME, is stored in IW(PME1...PME2)
C PME2:     the end of the current element
C PN:       pointer into a "clean" variable, also used to compress
C PSRC:     source pointer, for compression
C RSTRT:    counter for the number of restarts carried out
C SLENME:   number of variables in variable list of pivotal variable
C THRESH:   local copy of ICNTL(4)
C THRESM :    local integer holding the threshold used to detect quasi
C         dense rows. When quasi dense rows are reintegrated in the
C         graph to be processed then THRESM is modified.
C WE:       W(E)
C WFLG:     used for flagging the W array.  See description of W.
C WNVI:     WFLG-NV(I)
C X:      either a supervariable or an element
C
C OPS:      counter for forecast number of flops
C RELDEN :    holds average density to set THRESM automatically
C SM:       counter used for forming standard deviation
C STD:      standard deviation
C
C IDENSE is true if supervariable I is dense
C
C -------------------------------------------------------------------
C  FUNCTIONS CALLED:
C -------------------------------------------------------------------
    INTRINSIC MAX, MIN, MOD
C ====================================================================
C  INITIALIZATIONS
C ====================================================================

    DO 2 I = 1,10
       RJNFO(I) = 0.0
       JNFO(I) = 0
 2    CONTINUE
    DMAX = 0
    HMOD = MAX (1, N-1)
    IOVFLO = ICNTL(5)
    LASTD = 0
    MEM = PFREE - 1
    MAXMEM = MEM
    MINDEG = 1
    NBD   = 0
    NCMPA = 0
    NEL = 0
    NFULL  = 0
    NRLADU = 0
    RSTRT = 0
    OPS = 0.00
    THRESH = ICNTL(4)
    WFLG = 2

C   ------------------------------------------------------
C   Experiments with automatic setting of parameter THRESH.
C   ------------------------------------------------------
    IF (THRESH.GT.0) THEN
       THRESM  = 0
       RELDEN = 0.0
       SM = 0
C      ----------------------------------------------------------
C      initialize arrays and eliminate rows with no off-diag. nz.
C      ----------------------------------------------------------
       DO 5 I=1,N
        THRESM = MAX(THRESM, LEN(I))
        IF (LEN(I).GT.0) THEN
         RELDEN = RELDEN + LEN(I)
         SM = SM + (LEN(I) * LEN(I))
        END IF
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        DEGREE (I) = LEN (I)
        IF (DEGREE(I) .EQ. 0) THEN
         NEL = NEL + 1
         ELEN (I) = -NEL
         PE (I) = 0
         W (I) = 0
         NRLADU = NRLADU + 1
         OPS = OPS + 1
        ELSE
         W (I) = 1
         ELEN (I) = 0
        ENDIF
 5     CONTINUE
       IF (N .EQ. NEL) GOTO 265

       RELDEN = RELDEN/(N-NEL)
C      RELDEN holds average row length
       SM = SM/(N-NEL-NFULL) - RELDEN*RELDEN
       STD = SQRT(ABS(SM))
C      STD holds standard deviation of the row lengths
       IF (STD .LE. RELDEN) THEN
        THRESM = -1
       ELSE
        THRESM = INT(9*RELDEN + 0.5*STD*((STD/(RELDEN+0.01))**1.5)+
   *       2*RELDEN*RELDEN/(STD+0.01) +1)
       END IF
C   ------------------------------------------------------
C   end automatic setting of THRESM
C   ------------------------------------------------------

    ELSE
       THRESM = THRESH
       DO 10 I = 1, N
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        DEGREE (I) = LEN (I)
        IF (DEGREE(I) .EQ. 0) THEN
         NEL = NEL + 1
         ELEN (I) = -NEL
         PE (I) = 0
         W (I) = 0
         NRLADU = NRLADU + 1
         OPS = OPS + 1
        ELSE
         W (I) = 1
         ELEN (I) = 0
        ENDIF
 10    CONTINUE
    ENDIF
    IF (THRESM.GE.0) THEN
       IF (THRESM.GE.N) THEN
C       full rows only
        THRESM = -1
       ELSE IF (THRESM.EQ.0) THEN
        THRESM = N
       ENDIF
    ENDIF

C   ----------------------------------------------------------------
C   initialize degree lists
C   ----------------------------------------------------------------
    DO 20 I = 1, N
       DEG = DEGREE (I)
       IF (DEG .GT. 0) THEN
C       ----------------------------------------------------------
C       place i in the degree list corresponding to its degree
C       or in the dense row list if i is dense
C       ----------------------------------------------------------
C       test for row density
        IF ( (THRESM.GE.0) .AND.
   &       (DEG+1.GE.THRESM.OR.DEG+1.GE.N-NEL )) THEN
C          I is dense and will be inserted in the degree
C          list of N
         NBD = NBD+1
         IF (DEG+1.NE.N-NEL) THEN
C           I is quasi dense
            DEGREE(I) = DEGREE(I)+N+1
C           insert I at the beginning of degree list of n
            DEG = N
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            DENXT (I) = INEXT
            HEAD (DEG) = I
            LAST(I)  = 0
            IF (LASTD.EQ.0) THEN
             LASTD=I
            END IF
         ELSE
C           I is full
            NFULL = NFULL+1
            DEGREE(I) = N+1
C           insert I at the end of degree list of n
            DEG = N
            IF (LASTD.EQ.0) THEN
C            degree list is empty
             LASTD   = I
             HEAD(DEG) = I
             DENXT(I)   = 0
             LAST(I)   = 0
            ELSE
C             IF (NFULL.EQ.1) THEN
C               First full row encountered
                DENXT(LASTD) = I
                LAST(I)   = LASTD
                LASTD     = I
                DENXT(I)   = 0
C             ELSE
C               Absorb I into LASTD (first full row found)
C                PE(I) = - LASTD
C                NV(LASTD) = NV(LASTD) + NV(I)
C                NV(I) = 0
C                ELEN(I) = 0
C             END IF
            ENDIF
         ENDIF
        ELSE
C          place i in the degree list corresponding to its degree
         INEXT = HEAD (DEG)
         IF (INEXT .NE. 0) LAST (INEXT) = I
         DENXT (I) = INEXT
         HEAD (DEG) = I
        ENDIF
       ENDIF
 20   CONTINUE

C   We suppress dense row selection if none of them was found in A
C   in the 1st pass
    IF (NBD.EQ.0 .AND. THRESH.GT.0) THEN
       THRESM = -1
    END IF
C
C ====================================================================
C  WHILE (selecting pivots) DO
C ====================================================================

 30   IF (NEL .LT. N) THEN

C ==================================================================
C  GET PIVOT OF MINIMUM APPROXIMATE DEGREE
C ==================================================================
C     -------------------------------------------------------------
C     find next supervariable for elimination
C     -------------------------------------------------------------
       DO 40 DEG = MINDEG, N
        ME = HEAD (DEG)
        IF (ME .GT. 0) GO TO 50
 40    CONTINUE
 50    MINDEG = DEG
       IF (DEG.LT.N)  THEN
C     -------------------------------------------------------------
C     remove chosen variable from linked list
C     -------------------------------------------------------------
        INEXT = DENXT (ME)
        IF (INEXT .NE. 0) LAST (INEXT) = 0
        HEAD (DEG) = INEXT
       ELSE
        IF (DEGREE(ME).EQ.N+1) GO TO 263
C DEGREE(ME).GT.N+1 so ME is quasi dense
C RESTARTING STRATEGY
C       FOR EACH  quasi dense row d
C       1/ insert d in the degree list according to the
C        value degree(d)-(N+1) (updating MINDEG)
C       2/ Build the adjacency list of d in the quotient graph
C          update DENXT(e_me)= DENXT(e_me)-NV(ME)
C       4/ get back to min degree process
C
C       THRESM > 0 because quasi dense rows were selected
C       While loop: ME is the current dense row
C       make sure that WFLG is not too large
        RSTRT = RSTRT + 1
        RELDEN = 0.0
        SM = 0
        IF (WFLG .GT. IOVFLO-NBD-1) THEN
         DO  51 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
 51        CONTINUE
         WFLG = 2
        END IF
        WFLG = WFLG + 1
        DO 57 IDUMMY = 1,N

C       ---------------------------------------------------------
C       remove chosen variable from link list
C       ---------------------------------------------------------
         INEXT = DENXT (ME)
         IF (INEXT .NE. 0) THEN
            LAST (INEXT) = 0
         ELSE
            LASTD = 0
         ENDIF
C       ----------------------------------------------------------
c       build adjacency list of ME in quotient graph
C       and calculate its external degree in ndense(me)
C       ----------------------------------------------------------
         DENXT(ME) = 0
C          Flag ME as having been considered in this calculation
         W(ME)    = WFLG
         P1 = PE(ME)
         P2 = P1 + LEN(ME) -1
C       LN-1 holds the pointer in IW to last elt/var in adj list
C          of ME.  LEN(ME) will then be set to LN-P1
C       ELN-1 hold the pointer in IW to  last elt in in adj list
C          of ME.  ELEN(ME) will then be set to ELN-P1
C       element adjacent to ME
         LN     = P1
         ELN    = P1
         DO 55 P=P1,P2
            E= IW(P)
            IF (W(E).EQ.WFLG) GO TO 55
            W(E) = WFLG
C         -------------------------------------------
C         Ensure that E is an unabsorbed element or a quasi dense
C           row and flag it
C         -------------------------------------------
            DO 52 JDUMMY = 1,N
             IF ( PE(E) .GE. 0 ) GOTO 53
             E = -PE(E)
             IF (W(E) .EQ.WFLG) GOTO 55
             W(E) = WFLG
 52         CONTINUE
 53         IF (ELEN(E).LT.0) THEN
C            E is a new element in adj(ME)
             DENXT(E) = DENXT(E) - NV(ME)
C            Move first entry in ME's list of adjacent variables
C            to the end
             IW(LN) = IW(ELN)
C            Place E at end of ME's list of adjacent elements
             IW(ELN) = E
             LN  = LN+1
             ELN = ELN + 1
C            update ndense of ME with all unflagged dense
C            rows in E
             PEE1 = PE(E)
             DO 54 PEE = PEE1, PEE1+LEN(E)-1
                X = IW(PEE)
                IF ((ELEN(X).GE.0).AND.(W(X).NE.WFLG)) THEN
C                X is a dense row
                 DENXT(ME) = DENXT(ME) + NV(X)
                 W(X) = WFLG
                ENDIF
 54            CONTINUE
            ELSE
C            E is a dense row
             DENXT(ME) = DENXT(ME) + NV(E)
C            Place E at end of ME's list of adjacent variables
             IW(LN)=E
             LN = LN+1
            ENDIF
 55        CONTINUE

C       ----------------------------------------------
C       DEGREE(ME)-(N+1) holds last external degree computed
C       when ME was detected as dense
C       DENXT(ME) is the exact external degree of ME
C       ----------------------------------------------
         WFLG   = WFLG + 1
         LEN(ME)  = LN-P1
         ELEN(ME) = ELN- P1
         NDME = DENXT(ME)+NV(ME)
C        If we want to select ME as full (NDME.EQ.NBD)
C        or quasi dense (NDME.GE.THRESM) then
C        denxt(of elements adjacent to ME) should be updated
         IF (DENXT(ME).EQ.0) DENXT(ME) =1
C       ---------------------------------------------------------
C       place ME in the degree list of DENXT(ME), update DEGREE
C       ---------------------------------------------------------
C         IF (DEGREE(ME)+NV(ME) .LT. NBD  ) THEN
         IF (NDME .LT. NBD) THEN
C           ME is not full
            RELDEN = RELDEN + NV(ME)*NDME
            SM = SM + NV(ME)*NDME*NDME
            DEGREE(ME) = DENXT(ME)
            DEG = DEGREE(ME)
            MINDEG = MIN(DEG,MINDEG)
            JNEXT = HEAD(DEG)
            IF (JNEXT.NE. 0) LAST (JNEXT) = ME
            DENXT(ME) = JNEXT
            HEAD(DEG) = ME
         ELSE
C           ME is full
            DEGREE(ME) = N+1
            DEG = DENXT(ME)
            MINDEG = MIN(DEG,MINDEG)
            DEG = N

C           Update DENXT of all elements in the list of elements
C           adjacent to ME
            P1 = PE(ME)
            P2 = P1 + ELEN(ME) - 1
            DO 56 PJ=P1,P2
             E= IW(PJ)
             DENXT (E) = DENXT(E) + NV(ME)
 56         CONTINUE
C            insert ME in the list of dense rows
            DEG = N
C           ME at the end of the list
            NFULL = NFULL +NV(ME)
            IF (LASTD.EQ.0) THEN
C                degree list is empty
             LASTD   = ME
             HEAD(N) = ME
             DENXT(ME)   = 0
             LAST(ME)   = 0
             IF (INEXT.EQ.0) INEXT = LASTD
            ELSE
C             IF (NFULL.EQ.NV(ME)) THEN
C                 First full row encountered
                DENXT(LASTD) = ME
                LAST(ME)   = LASTD
                LASTD      = ME
                DENXT(ME)   = 0
                IF (INEXT.EQ.0) INEXT = LASTD
C             ELSE
C             Absorb ME into LASTD (first full row found)
C                PE(ME) = - LASTD
C                NV(LASTD) = NV(LASTD) + NV(ME)
C                NV(ME) = 0
C                ELEN(ME) = 0
C             END IF
            ENDIF
         END IF

C       ------------------------------
C       process next quasi dense row
C       ------------------------------
         ME    = INEXT
         IF (ME.EQ.0) GO TO 58
         IF (DEGREE(ME).LE.(N+1) ) GOTO 58
 57       CONTINUE
 58       HEAD (N) = ME
C       ---------------------------------------
C       update dense row selection strategy
C       -------------------------------------
        IF (NBD.EQ.NFULL) THEN
         RELDEN = 0
         SM = 0
        ELSE
         RELDEN = (RELDEN + NFULL*NBD)/(NBD)
         SM = (SM + NFULL*NBD*NBD)/(NBD) - RELDEN*RELDEN
        END IF
        STD = SQRT(ABS(SM))
        THRESM = INT(9*RELDEN+0.5*STD*((STD/(RELDEN + 0.01))**1.5)
   *       + 2*RELDEN*RELDEN/(STD+0.01) +1)
        THRESM = MIN(THRESM,NBD)
        IF (THRESM.GE.NBD) THEN
         THRESM = N
        END IF
        NBD = NFULL
C       get back to min degree elimination loop
        GOTO 30
C       -------------------------------------------------------------
C       -------------------------------------------------------------
       ENDIF
C     -------------------------------------------------------------
C     me represents the elimination of pivots nel+1 to nel+nv(me).
C     place me itself as the first in this set.  It will be moved
C     to the nel+nv(me) position when the permutation vectors are
C     computed.
C     -------------------------------------------------------------
       ELENME = ELEN (ME)
       ELEN (ME) = - (NEL + 1)
       NVPIV = NV (ME)
       NEL = NEL + NVPIV
       DENXT(ME) = 0

C ====================================================================
C  CONSTRUCT NEW ELEMENT
C ====================================================================
C
C     -------------------------------------------------------------
C     At this point, me is the pivotal supervariable.  It will be
C     converted into the current element.  Scan list of the
C     pivotal supervariable, me, setting tree pointers and
C     constructing new list of supervariables for the new element,
C     me.  p is a pointer to the current position in the old list.
C     -------------------------------------------------------------
C
C     flag the variable "me" as being in the front by negating nv(me)
       NV (ME) = -NVPIV
       DEGME = 0
       IF (ELENME .EQ. 0) THEN
C       ----------------------------------------------------------
C       There are no elements involved.
C       Construct the new element in place.
C       ----------------------------------------------------------
        PME1 = PE (ME)
        PME2 = PME1 - 1
        DO 60 P = PME1, PME1 + LEN (ME) - 1
         I = IW (P)
         NVI = NV (I)
         IF (NVI .GT. 0) THEN
C         ----------------------------------------------------
C         i is a principal variable not yet placed in the
C         generated element. Store i in new list
C         ----------------------------------------------------
            DEGME = DEGME + NVI
C           flag i as being in Lme by negating nv (i)
            NV (I) = -NVI
            PME2 = PME2 + 1
            IW (PME2) = I

C         ----------------------------------------------------
C         remove variable i from degree list.
C         ----------------------------------------------------
C         only done for sparse rows
            IF (DEGREE(I).LE.N) THEN
             ILAST = LAST (I)
             INEXT = DENXT (I)
             IF (INEXT .NE. 0) LAST (INEXT) = ILAST
             IF (ILAST .NE. 0) THEN
                DENXT (ILAST) = INEXT
             ELSE
C               i is at the head of the degree list
                HEAD (DEGREE (I)) = INEXT
             ENDIF
            ELSE
C            Dense rows remain dense so do not remove from list
             DENXT(ME) = DENXT(ME) + NVI
            ENDIF
         ENDIF
 60       CONTINUE
C       this element takes no new memory in iw:
        NEWMEM = 0
       ELSE
C       ----------------------------------------------------------
C       construct the new element in empty space, iw (pfree ...)
C       ----------------------------------------------------------
        P  = PE (ME)
        PME1 = PFREE
        SLENME = LEN (ME) - ELENME
        DO 120 KNT1 = 1, ELENME
C          search the elements in me.
         E = IW (P)
         P = P + 1
         PJ = PE (E)
         LN = LEN (E)
C       -------------------------------------------------------
C       search for different supervariables and add them to the
C       new list, compressing when necessary.
C       -------------------------------------------------------
         DO 110 KNT2 = 1, LN
            I = IW (PJ)
            PJ = PJ + 1
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
C         -------------------------------------------------
C         compress iw, if necessary
C         -------------------------------------------------
             IF (PFREE .GT. IWLEN) THEN
C            prepare for compressing iw by adjusting
C            pointers and lengths so that the lists being
C            searched in the inner and outer loops contain
C            only the remaining entries.
C            ***** SD: Seperate compression subroutine tried
C              but found to be inefficient in comparison ****
                PE (ME) = P
                LEN (ME) = LEN (ME) - KNT1
C               Check if anything left in supervariable ME
                IF (LEN (ME) .EQ. 0) PE (ME) = 0
                PE (E) = PJ
                LEN (E) = LN - KNT2
C               Check if anything left in element E
                IF (LEN (E) .EQ. 0) PE (E) = 0
                NCMPA = NCMPA + 1
C               store first item in pe
C               set first entry to -item
                DO 70 J = 1, N
                 PN = PE (J)
                 IF (PN .GT. 0) THEN
                  PE (J) = IW (PN)
                  IW (PN) = -J
                 ENDIF
 70             CONTINUE

C               psrc/pdst point to source/destination
                PDST = 1
                PSRC = 1
                PEND = PME1 - 1

C               while loop:
                DO 91 IDUMMY = 1, IWLEN
                 IF (PSRC .GT. PEND) GO TO 95
C                search for next negative entry
                 J = -IW (PSRC)
                 PSRC = PSRC + 1
                 IF (J .GT. 0) THEN
                  IW (PDST) = PE (J)
                  PE (J) = PDST
                  PDST = PDST + 1
C             copy from source to destination
                  LENJ = LEN (J)
                  DO 90 KNT3 = 0, LENJ - 2
                     IW (PDST + KNT3) = IW (PSRC + KNT3)
 90                 CONTINUE
                  PDST = PDST + LENJ - 1
                  PSRC = PSRC + LENJ - 1
                 ENDIF
 91             END DO

C               move the new partially-constructed element
 95             P1 = PDST
                DO 100 PSRC = PME1, PFREE - 1
                 IW (PDST) = IW (PSRC)
                 PDST = PDST + 1
 100            CONTINUE
                PME1 = P1
                PFREE = PDST
                PJ = PE (E)
                P = PE (ME)
             ENDIF

C         -------------------------------------------------
C         i is a principal variable not yet placed in Lme
C         store i in new list
C         -------------------------------------------------
             DEGME = DEGME + NVI
C            flag i as being in Lme by negating nv (i)
             NV (I) = -NVI
             IW (PFREE) = I
             PFREE = PFREE + 1

C         -------------------------------------------------
C         remove variable i from degree link list
C         -------------------------------------------------
C            only done for sparse rows
             IF (DEGREE(I).LE.N) THEN
                ILAST = LAST (I)
                INEXT = DENXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                 DENXT (ILAST) = INEXT
                ELSE
C               i is at the head of the degree list
                 HEAD (DEGREE (I)) = INEXT
                ENDIF
             ELSE
C            Dense rows remain dense so do not remove from list
                DENXT(ME) = DENXT(ME) + NVI
             ENDIF
            ENDIF
 110       CONTINUE

C          set tree pointer and flag to indicate element e is
C          absorbed into new element me (the parent of e is me)
            PE (E) = -ME
            W (E) = 0
 120      CONTINUE

C       search the supervariables in me.
        KNT1 = ELENME + 1
        E = ME
        PJ = P
        LN = SLENME

C       -------------------------------------------------------
C       search for different supervariables and add them to the
C       new list, compressing when necessary.
C       -------------------------------------------------------
        DO 126 KNT2 = 1, LN
         I = IW (PJ)
         PJ = PJ + 1
         NVI = NV (I)
         IF (NVI .GT. 0) THEN
C         -------------------------------------------------
C         compress iw, if necessary
C         -------------------------------------------------
            IF (PFREE .GT. IWLEN) THEN
C           prepare for compressing iw by adjusting
C           pointers and lengths so that the lists being
C           searched in the inner and outer loops contain
C           only the remaining entries.
             PE (ME) = P
             LEN (ME) = LEN (ME) - KNT1
C            Check if anything left in supervariable ME
             IF (LEN (ME) .EQ. 0) PE (ME) = 0
             PE (E) = PJ
             LEN (E) = LN - KNT2
C            Check if anything left in element E
             IF (LEN (E) .EQ. 0) PE (E) = 0
             NCMPA = NCMPA + 1
C            store first item in pe
C            set first entry to -item
             DO 121 J = 1, N
                PN = PE (J)
                IF (PN .GT. 0) THEN
                 PE (J) = IW (PN)
                 IW (PN) = -J
                ENDIF
 121           CONTINUE

C            psrc/pdst point to source/destination
             PDST = 1
             PSRC = 1
             PEND = PME1 - 1

C           while loop:
C 122          CONTINUE
             DO 123 IDUMMY = 1,IWLEN
                IF (PSRC .GT. PEND) GO TO 124
C               search for next negative entry
                J = -IW (PSRC)
                PSRC = PSRC + 1
                IF (J .GT. 0) THEN
                 IW (PDST) = PE (J)
                 PE (J) = PDST
                 PDST = PDST + 1
C                copy from source to destination
                 LENJ = LEN (J)
                 DO 122 KNT3 = 0, LENJ - 2
                  IW (PDST + KNT3) = IW (PSRC + KNT3)
 122               CONTINUE
                 PDST = PDST + LENJ - 1
                 PSRC = PSRC + LENJ - 1
                ENDIF
 123           END DO

C           move the new partially-constructed element
 124           P1 = PDST
             DO 125 PSRC = PME1, PFREE - 1
                IW (PDST) = IW (PSRC)
                PDST = PDST + 1
 125           CONTINUE
             PME1 = P1
             PFREE = PDST
             PJ = PE (E)
             P = PE (ME)
            END IF

C         -------------------------------------------------
C         i is a principal variable not yet placed in Lme
C         store i in new list
C         -------------------------------------------------
            DEGME = DEGME + NVI
C           flag i as being in Lme by negating nv (i)
            NV (I) = -NVI
            IW (PFREE) = I
            PFREE = PFREE + 1

C         -------------------------------------------------
C         remove variable i from degree link list
C         -------------------------------------------------
C           only done for sparse rows
            IF (DEGREE(I).LE.N) THEN
             ILAST = LAST (I)
             INEXT = DENXT (I)
             IF (INEXT .NE. 0) LAST (INEXT) = ILAST
             IF (ILAST .NE. 0) THEN
                DENXT (ILAST) = INEXT
             ELSE
C           i is at the head of the degree list
                HEAD (DEGREE (I)) = INEXT
             ENDIF
            ELSE
C            Dense rows remain dense so do not remove from list
             DENXT(ME) = DENXT(ME) + NVI
            ENDIF
         ENDIF
 126       CONTINUE

        PME2 = PFREE - 1
C       this element takes newmem new memory in iw (possibly zero)
        NEWMEM = PFREE - PME1
        MEM = MEM + NEWMEM
        MAXMEM = MAX (MAXMEM, MEM)
       ENDIF

C     -------------------------------------------------------------
C     me has now been converted into an element in iw (pme1..pme2)
C     -------------------------------------------------------------
C     degme holds the external degree of new element
       DEGREE (ME) = DEGME
       PE (ME) = PME1
       LEN (ME) = PME2 - PME1 + 1

C     -------------------------------------------------------------
C     make sure that wflg is not too large.  With the current
C     value of wflg, wflg+n must not cause integer overflow
C     -------------------------------------------------------------
       IF (WFLG .GT. IOVFLO-N) THEN
        DO 130 X = 1, N
         IF (W (X) .NE. 0) W (X) = 1
 130      CONTINUE
        WFLG = 2
       ENDIF

C ====================================================================
C   COMPUTE (w(e) - wflg) = |Le(G')\Lme(G')| FOR ALL ELEMENTS
C   where G' is the subgraph of G containing just the sparse rows)
C ====================================================================
C     -------------------------------------------------------------
C     Scan 1:  compute the external degrees of elements touched
C     with respect to the current element.  That is:
C        (w (e) - wflg) = |Le \ Lme|
C     for each element e involving a supervariable in Lme.
C     The notation Le refers to the pattern (list of
C     supervariables) of a previous element e, where e is not yet
C     absorbed, stored in iw (pe (e) + 1 ... pe (e) + iw (pe (e))).
C     The notation Lme refers to the pattern of the current element
C     (stored in iw (pme1..pme2)).
C     aggressive absorption is possible only if DENXT(ME) = NBD
C     which is true when only full rows have been selected.
C     -------------------------------------------------------------
       IF (NBD.GT.0) THEN
C       Dense rows have been found
        DO 150 PME = PME1, PME2
         I = IW (PME)
c          skip dense rows
         IF (DEGREE(I).GT.N) GOTO 150
         ELN = ELEN (I)
         IF (ELN .GT. 0) THEN
C          note that nv (i) has been negated to denote i in Lme:
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
             E = IW (P)
             WE = W (E)
             IF (WE .GE. WFLG) THEN
C            unabsorbed element e has been seen in this loop
                WE = WE - NVI
             ELSE IF (WE .NE. 0) THEN
C            e is an unabsorbed element - this is
C            the first we have seen e in all of Scan 1
                WE = DEGREE (E) + WNVI - DENXT(E)
             ENDIF
             W (E) = WE
 140          CONTINUE
         ENDIF
 150      CONTINUE
       ELSE
C       No dense rows have been found
        DO 152 PME = PME1, PME2
         I = IW (PME)
         ELN = ELEN (I)
         IF (ELN .GT. 0) THEN
C       note that nv (i) has been negated to denote i in Lme:
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 151 P = PE (I), PE (I) + ELN - 1
             E = IW (P)
             WE = W (E)
             IF (WE .GE. WFLG) THEN
C            unabsorbed element e has been seen in this loop
                WE = WE - NVI
             ELSE IF (WE .NE. 0) THEN
C            e is an unabsorbed element - this is
C            the first we have seen e in all of Scan 1
                WE = DEGREE (E) + WNVI
             ENDIF
             W (E) = WE
 151          CONTINUE
         ENDIF
 152      CONTINUE
       END IF

C ====================================================================
C  DEGREE UPDATE AND ELEMENT ABSORPTION
C ====================================================================

C     -------------------------------------------------------------
C     Scan 2:  for each sparse i in Lme, sum up the external degrees
C     of each Le for the elements e appearing within i, plus the
C     supervariables in i.  Place i in hash list.
C     -------------------------------------------------------------
       IF (NBD.GT.0) THEN
C       Dense rows have been found
        DO 180 PME = PME1, PME2
         I = IW (PME)
C          skip dense rows
         IF (DEGREE(I).GT.N) GOTO 180
C          remove absorbed elements from the list for i
         P1 = PE (I)
         P2 = P1 + ELEN (I) - 1
         PN = P1
         HASH = 0
         DEG = 0

C       ----------------------------------------------------------
C       scan the element list associated with supervariable i
C       ----------------------------------------------------------
         DO 160 P = P1, P2
            E = IW (P)
C           dext = | Le | - | (Le \cap Lme)\D | - DENXT(e)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
             DEG = DEG + DEXT
             IW (PN) = E
             PN = PN + 1
             HASH = HASH+E
            ELSE IF ((DEXT .EQ. 0) .AND.
   &            (DENXT(ME).EQ.NBD)) THEN
C         aggressive absorption: e is not adjacent to me, but
C         |Le(G') \ Lme(G')| is 0 and all dense rows
C         are in me, so absorb it into me
             PE (E) = -ME
             W (E)  = 0
            ELSE IF (DEXT.EQ.0) THEN
             IW(PN) = E
             PN   = PN+1
             HASH = HASH + E
            ENDIF
 160       CONTINUE

C       count the number of elements in i (including me):
         ELEN (I) = PN - P1 + 1

C       ----------------------------------------------------------
C       scan the supervariables in the list associated with i
C       ----------------------------------------------------------
         P3 = PN
         DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
C           j is unabsorbed, and not in Lme.
C           add to degree and add to new list
C           add degree only of sparse rows.
             IF (DEGREE(J).LE.N) DEG=DEG+NVJ
             IW (PN) = J
             PN = PN + 1
             HASH = HASH + J
            ENDIF
 170       CONTINUE

C       ----------------------------------------------------------
C       update the degree and check for mass elimination
C       ----------------------------------------------------------
         IF ((DEG .EQ. 0).AND.(DENXT(ME).EQ.NBD)) THEN
C          mass elimination only possible when all dense rows
C          are in ME
C       -------------------------------------------------------
C       mass elimination - supervariable i can be eliminated
C       -------------------------------------------------------
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
         ELSE
C       -------------------------------------------------------
C       update the upper-bound degree of i
C       A bound for the new external degree is the old bound plus
C       the size of the generated element
C       -------------------------------------------------------
C       the following degree does not yet include the size
C       of the current element, which is added later:
            DEGREE(I) = MIN (DEG+NBD-DENXT(ME), DEGREE(I))

C       -------------------------------------------------------
C       add me to the list for i
C       -------------------------------------------------------
C          move first supervariable to end of list
            IW (PN) = IW (P3)
C          move first element to end of element part of list
            IW (P3) = IW (P1)
C          add new element to front of list.
            IW (P1) = ME
C          store the new length of the list in len (i)
            LEN (I) = PN - P1 + 1

C       -------------------------------------------------------
C       place in hash bucket.  Save hash key of i in last (i).
C       -------------------------------------------------------
            HASH = ABS(MOD (HASH, HMOD)) + 1
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
C          the degree list is empty, hash head is -j
             DENXT (I) = -J
             HEAD (HASH) = -I
            ELSE
C          degree list is not empty - has j as its head
C          last is hash head
             DENXT (I) = LAST (J)
             LAST (J) = I
            ENDIF
            LAST (I) = HASH
         ENDIF
 180      CONTINUE
       ELSE
C       No dense rows have been found
        DO 183 PME = PME1, PME2
         I = IW (PME)
C          remove absorbed elements from the list for i
         P1 = PE (I)
         P2 = P1 + ELEN (I) - 1
         PN = P1
         HASH = 0
         DEG = 0

C          -------------------------------------------------------
C          scan the element list associated with supervariable i
C          -------------------------------------------------------
         DO 181 P = P1, P2
            E = IW (P)
C           dext = | Le | - | (Le \cap Lme)\D | - DENXT(e)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
             DEG = DEG + DEXT
             IW (PN) = E
             PN = PN + 1
             HASH = HASH + E
            ELSE IF (DEXT .EQ. 0) THEN
C            aggressive absorption: e is not adjacent to me, but
C            |Le(G') \ Lme(G')| is 0, so absorb it into me
             PE (E) = -ME
             W (E)  = 0
            ENDIF
 181       CONTINUE

C       count the number of elements in i (including me):
         ELEN (I) = PN - P1 + 1

C       ----------------------------------------------------------
C       scan the supervariables in the list associated with i
C       ----------------------------------------------------------
         P3 = PN
         DO 182 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
C           j is unabsorbed, and not in Lme.
C           add to degree and add to new list
             DEG=DEG+NVJ
             IW (PN) = J
             PN = PN + 1
             HASH = HASH + J
            ENDIF
 182       CONTINUE

C       ----------------------------------------------------------
C       update the degree and check for mass elimination
C       ----------------------------------------------------------
         IF (DEG .EQ. 0) THEN
C       -------------------------------------------------------
C       mass elimination - supervariable i can be eliminated
C       -------------------------------------------------------
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
         ELSE
C       -------------------------------------------------------
C       update the upper-bound degree of i
C       A bound for the new external degree is the old bound plus
C       the size of the generated element
C       -------------------------------------------------------
C
C       the following degree does not yet include the size
C       of the current element, which is added later:
            DEGREE(I) = MIN (DEG,  DEGREE(I))

C       -------------------------------------------------------
C       add me to the list for i
C       -------------------------------------------------------
C          move first supervariable to end of list
            IW (PN) = IW (P3)
C          move first element to end of element part of list
            IW (P3) = IW (P1)
C          add new element to front of list.
            IW (P1) = ME
C          store the new length of the list in len (i)
            LEN (I) = PN - P1 + 1

C       -------------------------------------------------------
C       place in hash bucket.  Save hash key of i in last (i).
C       -------------------------------------------------------
            HASH = ABS(MOD (HASH, HMOD)) + 1
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
C          the degree list is empty, hash head is -j
             DENXT (I) = -J
             HEAD (HASH) = -I
            ELSE
C          degree list is not empty - has j as its head
C          last is hash head
             DENXT (I) = LAST (J)
             LAST (J) = I
            ENDIF
            LAST (I) = HASH
         ENDIF
 183      CONTINUE
       END IF
       DEGREE (ME) = DEGME

C     -------------------------------------------------------------
C     Clear the counter array, w (...), by incrementing wflg.
C     -------------------------------------------------------------
       DMAX = MAX (DMAX, DEGME)
       WFLG = WFLG + DMAX

C      make sure that wflg+n does not cause integer overflow
       IF (WFLG .GE. IOVFLO - N) THEN
        DO 190 X = 1, N
         IF (W (X) .NE. 0) W (X) = 1
 190      CONTINUE
        WFLG = 2
       ENDIF
C      at this point, w (1..n) .lt. wflg holds

C ====================================================================
C  SUPERVARIABLE DETECTION
C ====================================================================
       DO 250 PME = PME1, PME2
        I = IW (PME)
        IF ( (NV(I).GE.0) .OR. (DEGREE(I).GT.N) ) GO TO 250
C       only done for sparse rows
C       replace i by head of its hash bucket, and set the hash
C       bucket header to zero

C       -------------------------------------------------------
C       examine all hash buckets with 2 or more variables.  We
C       do this by examing all unique hash keys for super-
C       variables in the pattern Lme of the current element, me
C       -------------------------------------------------------
        HASH = LAST (I)
C       let i = head of hash bucket, and empty the hash bucket
        J = HEAD (HASH)
        IF (J .EQ. 0) GO TO 250
        IF (J .LT. 0) THEN
C         degree list is empty
         I = -J
         HEAD (HASH) = 0
        ELSE
C         degree list is not empty, restore last () of head
         I = LAST (J)
         LAST (J) = 0
        ENDIF
        IF (I .EQ. 0) GO TO 250

C       while loop:
        DO 247 JDUMMY = 1,N
         IF (DENXT (I) .EQ. 0) GO TO 250
C         ----------------------------------------------------
C         this bucket has one or more variables following i.
C         scan all of them to see if i can absorb any entries
C         that follow i in hash bucket.  Scatter i into w.
C         ----------------------------------------------------
         LN = LEN (I)
         ELN = ELEN (I)
C          do not flag the first element in the list (me)
         DO 210 P = PE (I) + 1, PE (I) + LN - 1
            W (IW (P)) = WFLG
 210       CONTINUE

C         ----------------------------------------------------
C         scan every other entry j following i in bucket
C         ----------------------------------------------------
         JLAST = I
         J = DENXT (I)

C         while loop:
         DO 245 IDUMMY=1,N
            IF (J .EQ. 0) GO TO 246

C         -------------------------------------------------
C         check if j and i have identical nonzero pattern
C         -------------------------------------------------
C         jump if i and j do not have same size data structure
            IF (LEN (J) .NE. LN) GO TO 240
C         jump if i and j do not have same number adj elts
            IF (ELEN (J) .NE. ELN) GO TO 240
C         do not flag the first element in the list (me)

            DO 230 P = PE (J) + 1, PE (J) + LN - 1
C           jump if an entry (iw(p)) is in j but not in i
             IF (W (IW (P)) .NE. WFLG) GO TO 240
 230          CONTINUE

C           -------------------------------------------------
C           found it!  j can be absorbed into i
C           -------------------------------------------------
            PE (J) = -I
C           both nv (i) and nv (j) are negated since they
C           are in Lme, and the absolute values of each
C           are the number of variables in i and j:
            NV (I) = NV (I) + NV (J)
            NV (J) = 0
            ELEN (J) = 0
C           delete j from hash bucket
            J = DENXT (J)
            DENXT (JLAST) = J
            GO TO 245

C           -------------------------------------------------
 240          CONTINUE
C         j cannot be absorbed into i
C         -------------------------------------------------
            JLAST = J
            J = DENXT (J)
 245       CONTINUE

C         ----------------------------------------------------
C         no more variables can be absorbed into i
C         go to next i in bucket and clear flag array
C         ----------------------------------------------------
 246       WFLG = WFLG + 1
         I = DENXT (I)
         IF (I .EQ. 0) GO TO 250
 247       CONTINUE
 250    CONTINUE

C ====================================================================
C  RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMENT
C  Squeeze out absorbed variables
C ====================================================================
      P = PME1
      NLEFT = N - NEL
      DO 260 PME = PME1, PME2
         I = IW (PME)
         NVI = -NV (I)
         IF (NVI .LE. 0) GO TO 260
C        i is a principal variable in Lme
C        restore nv (i) to signify that i is principal
         NV (I) = NVI
         IF (DEGREE(I).GT.N) GO TO 258
C       -------------------------------------------------------
C       compute the external degree (add size of current elem)
C       -------------------------------------------------------
         DEG = MIN (DEGREE (I)+ DEGME - NVI, NLEFT - NVI)
         DEGREE (I) = DEG
         IDENSE = .FALSE.
C       -------------------
C       Dense row detection
C       -------------------
         IF (THRESM.GE.0) THEN
C       DEGME is exact external degree of pivot ME |Le\Ve|,
C       DEG is is approx external degree of I
          IF ((DEG+NVI .GE. THRESM).OR.
   &         (DEG+NVI .GE. NLEFT)) THEN
             IF (THRESM.EQ.N) THEN
C            We must be sure that I is full in reduced matrix
              IF ((ELEN(I).LE.2) .AND.((DEG+NVI).EQ.NLEFT)
   &             .AND. NBD.EQ.NFULL ) THEN
C                DEG approximation is exact and I is dense
               DEGREE(I) = N+1
               IDENSE = .TRUE.
              ENDIF
             ELSE
C             relaxed dense row detection
              IDENSE = .TRUE.
              IF ((ELEN(I).LE.2).AND. ((DEG+NVI).EQ.NLEFT)
   &             .AND. NBD.EQ.NFULL ) THEN
               DEGREE(I) = N+1
              ELSE
               DEGREE(I) = N+1+DEGREE(I)
              ENDIF
             ENDIF
          ENDIF
          IF (IDENSE) THEN
C            update DENXT of all elements in the list of element
C            adjacent to I (including ME).
             P1 = PE(I)
             P2 = P1 + ELEN(I) - 1
             DO 255 PJ=P1,P2
              E= IW(PJ)
              DENXT (E) = DENXT(E) + NVI
 255         CONTINUE
C            insert I in the list of dense rows
             NBD = NBD+NVI
             DEG = N
             IF (DEGREE(I).EQ.N+1) THEN
c             insert I at the end of the list
              NFULL = NFULL +NVI
              IF (LASTD.EQ.0) THEN
C                degree list is empty
               LASTD   = I
               HEAD(DEG) = I
               DENXT(I)   = 0
               LAST(I)   = 0
              ELSE
C               IF (NFULL.EQ.NVI) THEN
C                 First full row encountered
                  DENXT(LASTD) = I
                  LAST(I)   = LASTD
                  LASTD     = I
                  DENXT(I)   = 0
C               ELSE
C                 Absorb I into LASTD (first full row found)
C                  PE(I) = - LASTD
C                  NV(LASTD) = NV(LASTD) + NV(I)
C                  NV(I) = 0
C                  ELEN(I) = 0
C               END IF
              ENDIF
             ELSE
C             insert I at the beginning of the list
              INEXT = HEAD(DEG)
              IF (INEXT .NE. 0) LAST (INEXT) = I
              DENXT (I) = INEXT
              HEAD (DEG) = I
              LAST(I)    = 0
              IF (LASTD.EQ.0) LASTD=I
             ENDIF
C         end of IDENSE=true
          ENDIF
C        end of THRESM>0
         ENDIF

         IF (.NOT.IDENSE) THEN
C       -------------------------------------------------------
C       place the supervariable at the head of the degree list
C       -------------------------------------------------------
          INEXT = HEAD (DEG)
          IF (INEXT .NE. 0) LAST (INEXT) = I
          DENXT (I) = INEXT
          LAST (I) = 0
          HEAD (DEG) = I
         ENDIF
C       -------------------------------------------------------
C       save the new degree, and find the minimum degree
C       -------------------------------------------------------
         MINDEG = MIN (MINDEG, DEG)
 258       CONTINUE
C       -------------------------------------------------------
C       place the supervariable in the element pattern
C       -------------------------------------------------------
         IW (P) = I
         P = P + 1
 260    CONTINUE

C =====================================================================
C  FINALIZE THE NEW ELEMENT
C =====================================================================
      OPS = OPS + DEGME*NVPIV + DEGME * NVPIV*NVPIV +
   *       DEGME*DEGME*NVPIV + NVPIV*NVPIV*NVPIV/3 +
   *       NVPIV*NVPIV/2 + NVPIV/6 + NVPIV
      NRLADU = NRLADU + (NVPIV*(NVPIV+1))/2 + (DEGME*NVPIV)
      NV (ME) = NVPIV + DEGME
C       nv (me) is now the degree of pivot (including diagonal part)
C       save the length of the list for the new element me
      LEN (ME) = P - PME1
      IF (LEN (ME) .EQ. 0) THEN
C        there is nothing left of the current pivot element
         PE (ME) = 0
         W (ME) = 0
      ENDIF
      IF (NEWMEM .NE. 0) THEN
C        element was not constructed in place: deallocate part
C        of it (final size is less than or equal to newmem,
C        since newly nonprincipal variables have been removed).
         PFREE = P
         MEM = MEM - NEWMEM + LEN (ME)
      ENDIF

C =====================================================================
C     END WHILE (selecting pivots)
      GO TO 30
     ENDIF
C =====================================================================
     GO TO 265

C    We have only full rows that we amalgamate at the root
C    node and ME = LASTD
C    Perform mass elimination of all full rows
 263   NELME    = -(NEL+1)
     DO 264 X=1,N
      IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
C        X is an unabsorbed element
         PE(X) = -ME
      ELSEIF (DEGREE(X).EQ.N+1) THEN
C        X is a dense row, absorb it in ME (mass elimination)
         NEL   = NEL + NV(X)
         PE(X) = -ME
         ELEN(X) = 0
         NV(X) = 0
      ENDIF
 264   CONTINUE
C    ME is the root node
     ELEN(ME) = NELME
     NV(ME)   = NBD
     NRLADU = NRLADU + (NBD*(NBD+1))/2
     OPS = OPS + NBD*NBD*NBD/3 + NBD*NBD/2 + NBD/6 + NBD
     PE(ME)   = 0

 265   CONTINUE
C ===================================================================
C  COMPUTE THE PERMUTATION VECTORS
C ===================================================================

C   ----------------------------------------------------------------
C   The time taken by the following code is O(n).  At this
C   point, elen (e) = -k has been done for all elements e,
C   and elen (i) = 0 has been done for all nonprincipal
C   variables i.  At this point, there are no principal
C   supervariables left, and all elements are absorbed.
C   ----------------------------------------------------------------
C
C   ----------------------------------------------------------------
C   compute the ordering of unordered nonprincipal variables
C   ----------------------------------------------------------------

     DO 290 I = 1, N
      IF (ELEN (I) .EQ. 0) THEN
C       ----------------------------------------------------------
C       i is an un-ordered row.  Traverse the tree from i until
C       reaching an element, e.  The element, e, was the
C       principal supervariable of i and all nodes in the path
C       from i to when e was selected as pivot.
C       ----------------------------------------------------------
         J = -PE (I)
C       while (j is a variable) do:
         DO 270 JDUMMY = 1,N
          IF (ELEN (J) .LT. 0) GO TO 275
          J = -PE (J)
 270       CONTINUE
 275       E = J
C       ----------------------------------------------------------
C       get the current pivot ordering of e
C       ----------------------------------------------------------
         K = -ELEN (E)

C       ----------------------------------------------------------
C       traverse the path again from i to e, and compress the
C       path (all nodes point to e).  Path compression allows
C       this code to compute in O(n) time.  Order the unordered
C       nodes in the path, and place the element e at the end.
C       ----------------------------------------------------------
         J = I
C        while (j is a variable) do:
         DO 280 IDUMMY = 1,N
          IF (ELEN (J) .LT. 0) GO TO 285
          JNEXT = -PE (J)
          PE (J) = -E
          IF (ELEN (J) .EQ. 0) THEN
C            j is an unordered row
             ELEN (J) = K
             K = K + 1
          ENDIF
          J = JNEXT
 280       CONTINUE
C        leave elen (e) negative, so we know it is an element
 285       ELEN (E) = -K
      ENDIF
 290   CONTINUE

C   ----------------------------------------------------------------
C   reset the inverse permutation (elen (1..n)) to be positive,
C   and compute the permutation (last (1..n)).
C   ----------------------------------------------------------------
     DO 300 I = 1, N
      K = ABS (ELEN (I))
      LAST (K) = I
      ELEN (I) = K
 300   CONTINUE

C ====================================================================
C  RETURN THE MEMORY USAGE IN IW AND SET INFORMATION ARRAYS
C ====================================================================
C   If maxmem is less than or equal to iwlen, then no compressions
C   occurred, and iw (maxmem+1 ... iwlen) was unused.  Otherwise
C   compressions did occur, and iwlen would have had to have been
C   greater than or equal to maxmem for no compressions to occur.
C   Return the value of maxmem in the pfree argument.

     RJNFO(1) = OPS
     RJNFO(2) = NRLADU
     JNFO(1) = NCMPA
     JNFO(2) = RSTRT
     PFREE = MAXMEM

     RETURN
     END
C ====================================================================
C ====================================================================
C ====================================================================
! COPYRIGHT (c) 1999 Council for the Central Laboratory
!            of the Research Councils
! Original date July 1999
! AUTHORS Iain Duff (i.duff@rl.ac.uk) and
!       Jacko Koster (jacko.koster@uninett.no)
!
! Version 1.6.0
! See ChangeLog for version history.
!


    SUBROUTINE MC64I(ICNTL)
    IMPLICIT NONE

C  Purpose
C  =======
C
C  The components of the array ICNTL control the action of MC64A/AD.
C  Default values for these are set in this subroutine.
C
C  Parameters
C  ==========
C
    INTEGER ICNTL(10)
C
C  Local variables
    INTEGER I
C
C    ICNTL(1) has default value 6.
C   It is the output stream for error messages. If it
C   is negative, these messages will be suppressed.
C
C    ICNTL(2) has default value 6.
C   It is the output stream for warning messages.
C   If it is negative, these messages are suppressed.
C
C    ICNTL(3) has default value -1.
C   It is the output stream for monitoring printing.
C   If it is negative, these messages are suppressed.
C
C    ICNTL(4) has default value 0.
C   If left at the defaut value, the incoming data is checked for
C   out-of-range indices and duplicates.  Setting ICNTL(4) to any
C   other will avoid the checks but is likely to cause problems
C   later if out-of-range indices or duplicates are present.
C   The user should only set ICNTL(4) non-zero, if the data is
C   known to avoid these problems.
C
C    ICNTL(5) to ICNTL(10) are not used by MC64A/AD but are set to
C   zero in this routine.

C Initialization of the ICNTL array.
    ICNTL(1) = 6
    ICNTL(2) = 6
    ICNTL(3) = -1
    DO 10 I = 4,10
      ICNTL(I) = 0
   10 CONTINUE

    RETURN
    END

C**********************************************************************
    SUBROUTINE MC64A(JOB,N,NE,IP,IRN,A,NUM,CPERM,LIW,IW,LDW,DW,
   &       ICNTL,INFO)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
C  Purpose
C  =======
C
C This subroutine attempts to find a column permutation for an NxN
C sparse matrix A = {a_ij} that makes the permuted matrix have N
C entries on its diagonal.
C If the matrix is structurally nonsingular, the subroutine optionally
C returns a column permutation that maximizes the smallest element
C on the diagonal, maximizes the sum of the diagonal entries, or
C maximizes the product of the diagonal entries of the permuted matrix.
C For the latter option, the subroutine also finds scaling factors
C that may be used to scale the matrix so that the nonzero diagonal
C entries of the permuted matrix are one in absolute value and all the
C off-diagonal entries are less than or equal to one in absolute value.
C The natural logarithms of the scaling factors u(i), i=1..N, for the
C rows and v(j), j=1..N, for the columns are returned so that the
C scaled matrix B = {b_ij} has entries b_ij = a_ij * EXP(u_i + v_j).
C
C  Parameters
C  ==========
C
    INTEGER JOB,N,NE,NUM,LIW,LDW
    INTEGER IP(N+1),IRN(NE),CPERM(N),IW(LIW),ICNTL(10),INFO(10)
    REAL A(NE),DW(LDW)
C
C JOB is an INTEGER variable which must be set by the user to
C control the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   1 Compute a column permutation of the matrix so that the
C   permuted matrix has as many entries on its diagonal as possible.
C   The values on the diagonal are of arbitrary size. HSL subroutine
C   MC21A/AD is used for this. See [1].
C   2 Compute a column permutation of the matrix so that the smallest
C   value on the diagonal of the permuted matrix is maximized.
C   See [3].
C   3 Compute a column permutation of the matrix so that the smallest
C   value on the diagonal of the permuted matrix is maximized.
C   The algorithm differs from the one used for JOB = 2 and may
C   have quite a different performance. See [2].
C   4 Compute a column permutation of the matrix so that the sum
C   of the diagonal entries of the permuted matrix is maximized.
C   See [3].
C   5 Compute a column permutation of the matrix so that the product
C   of the diagonal entries of the permuted matrix is maximized
C   and vectors to scale the matrix so that the nonzero diagonal
C   entries of the permuted matrix are one in absolute value and
C   all the off-diagonal entries are less than or equal to one in
C   absolute value. See [3].
C  Restriction: 1 <= JOB <= 5.
C
C N is an INTEGER variable which must be set by the user to the
C   order of the matrix A. It is not altered by the subroutine.
C   Restriction: N >= 1.
C
C NE is an INTEGER variable which must be set by the user to the
C   number of entries in the matrix. It is not altered by the
C   subroutine.
C   Restriction: NE >= 1.
C
C IP is an INTEGER array of length N+1.
C   IP(J), J=1..N, must be set by the user to the position in array IRN
C   of the first row index of an entry in column J. IP(N+1) must be set
C   to NE+1. It is not altered by the subroutine.
C
C IRN is an INTEGER array of length NE.
C   IRN(K), K=1..NE, must be set by the user to hold the row indices of
C   the entries of the matrix. Those belonging to column J must be
C   stored contiguously in the positions IP(J)..IP(J+1)-1. The ordering
C   of the row indices within each column is unimportant. Repeated
C   entries are not allowed. The array IRN is not altered by the
C   subroutine.
C
C A is a REAL array of length NE.
C   The user must set A(K), K=1..NE, to the numerical value of the
C   entry that corresponds to IRN(K).
C   It is not used by the subroutine when JOB = 1.
C   It is not altered by the subroutine.
C
C NUM is an INTEGER variable that need not be set by the user.
C   On successful exit, NUM will be the number of entries on the
C   diagonal of the permuted matrix.
C   If NUM < N, the matrix is structurally singular.
C
C CPERM is an INTEGER array of length N that need not be set by the
C   user. On successful exit, CPERM contains the column permutation.
C   Column ABS(CPERM(J)) of the original matrix is column J in the
C   permuted matrix, J=1..N. For the N-NUM  entries of CPERM that are
C   not matched the permutation array is set negative so that a full
C   permutation of the matrix is obtained even in the structurally
C   singular case.
C
C LIW is an INTEGER variable that must be set by the user to
C   the dimension of array IW. It is not altered by the subroutine.
C   Restriction:
C   JOB = 1 :  LIW >= 5N
C   JOB = 2 :  LIW >= 4N
C   JOB = 3 :  LIW >= 10N + NE
C   JOB = 4 :  LIW >= 5N
C   JOB = 5 :  LIW >= 5N
C
C IW is an INTEGER array of length LIW that is used for workspace.
C
C LDW is an INTEGER variable that must be set by the user to the
C   dimension of array DW. It is not altered by the subroutine.
C   Restriction:
C   JOB = 1 :  LDW is not used
C   JOB = 2 :  LDW >= N
C   JOB = 3 :  LDW >= NE
C   JOB = 4 :  LDW >= 2N + NE
C   JOB = 5 :  LDW >= 3N + NE
C
C DW is a REAL array of length LDW
C   that is used for workspace. If JOB = 5, on return,
C   DW(i) contains u_i, i=1..N, and DW(N+j) contains v_j, j=1..N.
C
C ICNTL is an INTEGER array of length 10. Its components control the
C   output of MC64A/AD and must be set by the user before calling
C   MC64A/AD. They are not altered by the subroutine.
C
C   ICNTL(1) must be set to specify the output stream for
C   error messages. If ICNTL(1) < 0, messages are suppressed.
C   The default value set by MC46I/ID is 6.
C
C   ICNTL(2) must be set by the user to specify the output stream for
C   warning messages. If ICNTL(2) < 0, messages are suppressed.
C   The default value set by MC46I/ID is 6.
C
C   ICNTL(3) must be set by the user to specify the output stream for
C   diagnostic messages. If ICNTL(3) < 0, messages are suppressed.
C   The default value set by MC46I/ID is -1.
C
C   ICNTL(4) must be set by the user to a value other than 0 to avoid
C   checking of the input data.
C   The default value set by MC46I/ID is 0.
C
C INFO is an INTEGER array of length 10 which need not be set by the
C   user. INFO(1) is set non-negative to indicate success. A negative
C   value is returned if an error occurred, a positive value if a
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : successful entry (for structurally singular matrix).
C   +2 : the returned scaling factors are large and may cause
C      overflow when used to scale the matrix.
C      (For JOB = 5 entry only.)
C   -1 : JOB < 1 or JOB > 5.  Value of JOB held in INFO(2).
C   -2 : N < 1.  Value of N held in INFO(2).
C   -3 : NE < 1. Value of NE held in INFO(2).
C   -4 : the defined length LIW violates the restriction on LIW.
C      Value of LIW required given by INFO(2).
C   -5 : the defined length LDW violates the restriction on LDW.
C      Value of LDW required given by INFO(2).
C   -6 : entries are found whose row indices are out of range. INFO(2)
C      contains the index of a column in which such an entry is found.
C   -7 : repeated entries are found. INFO(2) contains the index of a
C      column in which such entries are found.
C  INFO(3) to INFO(10) are not currently used and are set to zero by
C      the routine.
C
C References:
C  [1]  I. S. Duff, (1981),
C     "Algorithm 575. Permutations for a zero-free diagonal",
C     ACM Trans. Math. Software 7(3), 387-390.
C  [2]  I. S. Duff and J. Koster, (1998),
C     "The design and use of algorithms for permuting large
C     entries to the diagonal of sparse matrices",
C     SIAM J. Matrix Anal. Appl., vol. 20, no. 4, pp. 889-901.
C  [3]  I. S. Duff and J. Koster, (1999),
C     "On algorithms for permuting large entries to the diagonal
C     of sparse matrices",
C     Technical Report RAL-TR-1999-030, RAL, Oxfordshire, England.

C Local variables and parameters
    INTEGER I,J,K
    REAL FACT,ZERO,RINF
    PARAMETER (ZERO=0.0D+00)
C External routines and functions
    EXTERNAL MC21A,MC64B,MC64R,MC64S,MC64W
C Intrinsic functions
    INTRINSIC ABS,LOG

C Set RINF to largest positive real number (infinity)
    RINF = HUGE(RINF)

C Check value of JOB
    IF (JOB.LT.1 .OR. JOB.GT.5) THEN
      INFO(1) = -1
      INFO(2) = JOB
      IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
      GO TO 99
    ENDIF
C Check value of N
    IF (N.LT.1) THEN
      INFO(1) = -2
      INFO(2) = N
      IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
      GO TO 99
    ENDIF
C Check value of NE
    IF (NE.LT.1) THEN
      INFO(1) = -3
      INFO(2) = NE
      IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NE',NE
      GO TO 99
    ENDIF
C Check LIW
    IF (JOB.EQ.1) K = 5*N
    IF (JOB.EQ.2) K = 4*N
    IF (JOB.EQ.3) K = 10*N + NE
    IF (JOB.EQ.4) K = 5*N
    IF (JOB.EQ.5) K = 5*N
    IF (LIW.LT.K) THEN
      INFO(1) = -4
      INFO(2) = K
      IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
      GO TO 99
    ENDIF
C Check LDW
C If JOB = 1, do not check
    IF (JOB.GT.1) THEN
      IF (JOB.EQ.2) K = N
      IF (JOB.EQ.3) K = NE
      IF (JOB.EQ.4) K = 2*N + NE
      IF (JOB.EQ.5) K = 3*N + NE
      IF (LDW.LT.K) THEN
      INFO(1) = -5
      INFO(2) = K
      IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
      GO TO 99
      ENDIF
    ENDIF
    IF (ICNTL(4).EQ.0) THEN
C Check row indices. Use IW(1:N) as workspace
      DO 3 I = 1,N
      IW(I) = 0
    3   CONTINUE
      DO 6 J = 1,N
      DO 4 K = IP(J),IP(J+1)-1
        I = IRN(K)
C Check for row indices that are out of range
        IF (I.LT.1 .OR. I.GT.N) THEN
          INFO(1) = -6
          INFO(2) = J
          IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
          GO TO 99
        ENDIF
C Check for repeated row indices within a column
        IF (IW(I).EQ.J) THEN
          INFO(1) = -7
          INFO(2) = J
          IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I
          GO TO 99
        ELSE
          IW(I) = J
        ENDIF
    4   CONTINUE
    6   CONTINUE
    ENDIF

C Print diagnostics on input
    IF (ICNTL(3).GE.0) THEN
      WRITE(ICNTL(3),9020) JOB,N,NE
      WRITE(ICNTL(3),9021) (IP(J),J=1,N+1)
      WRITE(ICNTL(3),9022) (IRN(J),J=1,NE)
      IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(J),J=1,NE)
    ENDIF

C Set components of INFO to zero
    DO 8 I=1,10
      INFO(I) = 0
    8 CONTINUE

C Compute maximum matching with MC21A/AD
    IF (JOB.EQ.1) THEN
C Put length of column J in IW(J)
      DO 10 J = 1,N
      IW(J) = IP(J+1) - IP(J)
   10   CONTINUE
C IW(N+1:5N) is workspace
      CALL MC21A(N,IRN,NE,IP,IW(1),CPERM,NUM,IW(N+1))
      GO TO 90
    ENDIF

C Compute bottleneck matching
    IF (JOB.EQ.2) THEN
C IW(1:5N), DW(1:N) are workspaces
      CALL MC64B(N,NE,IP,IRN,A,CPERM,NUM,
   &   IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),DW)
      GO TO 90
    ENDIF

C Compute bottleneck matching
    IF (JOB.EQ.3) THEN
C Copy IRN(K) into IW(K), ABS(A(K)) into DW(K), K=1..NE
      DO 20 K = 1,NE
      IW(K) = IRN(K)
      DW(K) = ABS(A(K))
   20   CONTINUE
C Sort entries in each column by decreasing value.
      CALL MC64R(N,NE,IP,IW,DW)
C IW(NE+1:NE+10N) is workspace
      CALL MC64S(N,NE,IP,IW(1),DW,CPERM,NUM,IW(NE+1),
   &   IW(NE+N+1),IW(NE+2*N+1),IW(NE+3*N+1),IW(NE+4*N+1),
   &   IW(NE+5*N+1),IW(NE+6*N+1))
      GO TO 90
    ENDIF

    IF (JOB.EQ.4) THEN
      DO 50 J = 1,N
      FACT = ZERO
      DO 30 K = IP(J),IP(J+1)-1
        IF (ABS(A(K)).GT.FACT) FACT = ABS(A(K))
   30   CONTINUE
      DO 40 K = IP(J),IP(J+1)-1
        DW(2*N+K) = FACT - ABS(A(K))
   40   CONTINUE
   50   CONTINUE
C B = DW(2N+1:2N+NE); IW(1:5N) and DW(1:2N) are workspaces
      CALL MC64W(N,NE,IP,IRN,DW(2*N+1),CPERM,NUM,
   &   IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
   &   DW(1),DW(N+1))
      GO TO 90
    ENDIF

    IF (JOB.EQ.5) THEN
      DO 75 J = 1,N
      FACT = ZERO
      DO 60 K = IP(J),IP(J+1)-1
        DW(3*N+K) = ABS(A(K))
        IF (DW(3*N+K).GT.FACT) FACT = DW(3*N+K)
   60   CONTINUE
      DW(2*N+J) = FACT
      IF (FACT.NE.ZERO) THEN
        FACT = LOG(FACT)
      ELSE
        FACT = RINF/N
      ENDIF
      DO 70 K = IP(J),IP(J+1)-1
        IF (DW(3*N+K).NE.ZERO) THEN
          DW(3*N+K) = FACT - LOG(DW(3*N+K))
        ELSE
          DW(3*N+K) = RINF/N
        ENDIF
   70   CONTINUE
   75   CONTINUE
C B = DW(3N+1:3N+NE); IW(1:5N) and DW(1:2N) are workspaces
      CALL MC64W(N,NE,IP,IRN,DW(3*N+1),CPERM,NUM,
   &   IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
   &   DW(1),DW(N+1))
      IF (NUM.EQ.N) THEN
      DO 80 J = 1,N
        IF (DW(2*N+J).NE.ZERO) THEN
          DW(N+J) = DW(N+J) - LOG(DW(2*N+J))
        ELSE
          DW(N+J) = ZERO
        ENDIF
   80   CONTINUE
      ENDIF
C Check size of scaling factors
      FACT = 0.5*LOG(RINF)
      DO 86 J = 1,N
      IF (DW(J).LT.FACT .AND. DW(N+J).LT.FACT) GO TO 86
      INFO(1) = 2
      GO TO 90
   86   CONTINUE
C     GO TO 90
    ENDIF

   90 IF (INFO(1).EQ.0 .AND. NUM.LT.N) THEN
C Matrix is structurally singular, return with warning
      INFO(1) = 1
      IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9011) INFO(1)
    ENDIF
    IF (INFO(1).EQ.2) THEN
C Scaling factors are large, return with warning
      IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9012) INFO(1)
    ENDIF

C Print diagnostics on output
    IF (ICNTL(3).GE.0) THEN
      WRITE(ICNTL(3),9030) (INFO(J),J=1,2)
      WRITE(ICNTL(3),9031) NUM
      WRITE(ICNTL(3),9032) (CPERM(J),J=1,N)
      IF (JOB.EQ.5) THEN
      WRITE(ICNTL(3),9033) (DW(J),J=1,N)
      WRITE(ICNTL(3),9034) (DW(N+J),J=1,N)
      ENDIF
    ENDIF

C Return from subroutine.
   99 RETURN

 9001 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2,
   &      ' because ',(A),' = ',I10)
 9004 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
   &      '      LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
   &      '      LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
   &      '      Column ',I8,
   &      ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
   &      '      Column ',I8,
   &      ' contains two or more entries with row index ',I8)
 9011 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
   &      '      The matrix is structurally singular.')
 9012 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
   &      '      Some scaling factors may be too large.')
 9020 FORMAT (' ****** Input parameters for MC64A/AD:'/
   &      ' JOB = ',I8/' N   = ',I8/' NE  = ',I8)
 9021 FORMAT (' IP(1:N+1)  = ',8I8/(14X,8I8))
 9022 FORMAT (' IRN(1:NE)  = ',8I8/(14X,8I8))
 9023 FORMAT (' A(1:NE)    = ',4(1PE14.4)/(14X,4(1PE14.4)))
 9030 FORMAT (' ****** Output parameters for MC64A/AD:'/
   &      ' INFO(1:2)  = ',2I8)
 9031 FORMAT (' NUM      = ',I8)
 9032 FORMAT (' CPERM(1:N) = ',8I8/(14X,8I8))
 9033 FORMAT (' DW(1:N)    = ',5(F11.3)/(14X,5(F11.3)))
 9034 FORMAT (' DW(N+1:2N) = ',5(F11.3)/(14X,5(F11.3)))
    END

C**********************************************************************
    SUBROUTINE MC64B(N,NE,IP,IRN,A,IPERM,NUM,JPERM,PR,Q,L,D)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER N,NE,NUM
    INTEGER IP(N+1),IRN(NE),IPERM(N),JPERM(N),PR(N),Q(N),L(N)
    REAL A(NE),D(N)

C N, NE, IP, IRN are described in MC64A/AD.
C A is a REAL array of length
C   NE. A(K), K=1..NE, must be set to the value of the entry
C   that corresponds to IRN(K). It is not altered.
C IPERM is an INTEGER array of length N. On exit, it contains the
C    matching: IPERM(I) = 0 or row I is matched to column IPERM(I).
C NUM is INTEGER variable. On exit, it contains the cardinality of the
C    matching stored in IPERM.
C IW is an INTEGER work array of length 4N.
C DW is a REAL array of length N.

C Local variables
    INTEGER I,II,J,JJ,JORD,Q0,QLEN,IDUM,JDUM,ISP,JSP,
   &      K,KK,KK1,KK2,I0,UP,LOW,LPOS
    REAL CSP,DI,DNEW,DQ0,AI,A0,BV
C Local parameters
    REAL RINF,ZERO,MINONE
    PARAMETER (ZERO=0.0D+0,MINONE=-1.0D+0)
C Intrinsic functions
    INTRINSIC ABS,MIN
C External subroutines and/or functions
    EXTERNAL MC64D,MC64E,MC64F


C Set RINF to largest positive real number
    RINF = HUGE(RINF)

C Initialization
    NUM = 0
    BV = RINF
    DO 10 K = 1,N
      IPERM(K) = 0
      JPERM(K) = 0
      PR(K) = IP(K)
      D(K) = ZERO
   10 CONTINUE
C Scan columns of matrix;
    DO 20 J = 1,N
      A0 = MINONE
      DO 30 K = IP(J),IP(J+1)-1
      I = IRN(K)
      AI = ABS(A(K))
      IF (AI.GT.D(I)) D(I) = AI
      IF (JPERM(J).NE.0) GO TO 30
      IF (AI.GE.BV) THEN
        A0 = BV
        IF (IPERM(I).NE.0) GO TO 30
        JPERM(J) = I
        IPERM(I) = J
        NUM = NUM + 1
      ELSE
        IF (AI.LE.A0) GO TO 30
        A0 = AI
        I0 = I
      ENDIF
   30   CONTINUE
      IF (A0.NE.MINONE .AND. A0.LT.BV) THEN
      BV = A0
      IF (IPERM(I0).NE.0) GO TO 20
      IPERM(I0) = J
      JPERM(J) = I0
      NUM = NUM + 1
      ENDIF
   20 CONTINUE
C Update BV with smallest of all the largest maximum absolute values
C of the rows.
    DO 25 I = 1,N
      BV = MIN(BV,D(I))
   25 CONTINUE
    IF (NUM.EQ.N) GO TO 1000
C Rescan unassigned columns; improve initial assignment
    DO 95 J = 1,N
      IF (JPERM(J).NE.0) GO TO 95
      DO 50 K = IP(J),IP(J+1)-1
      I = IRN(K)
      AI = ABS(A(K))
      IF (AI.LT.BV) GO TO 50
      IF (IPERM(I).EQ.0) GO TO 90
      JJ = IPERM(I)
      KK1 = PR(JJ)
      KK2 = IP(JJ+1) - 1
      IF (KK1.GT.KK2) GO TO 50
      DO 70 KK = KK1,KK2
        II = IRN(KK)
        IF (IPERM(II).NE.0) GO TO 70
        IF (ABS(A(KK)).GE.BV) GO TO 80
   70   CONTINUE
      PR(JJ) = KK2 + 1
   50   CONTINUE
      GO TO 95
   80   JPERM(JJ) = II
      IPERM(II) = JJ
      PR(JJ) = KK + 1
   90   NUM = NUM + 1
      JPERM(J) = I
      IPERM(I) = J
      PR(J) = K + 1
   95 CONTINUE
    IF (NUM.EQ.N) GO TO 1000

C Prepare for main loop
    DO 99 I = 1,N
      D(I) = MINONE
      L(I) = 0
   99 CONTINUE

C Main loop ... each pass round this loop is similar to Dijkstra's
C algorithm for solving the single source shortest path problem

    DO 100 JORD = 1,N

      IF (JPERM(JORD).NE.0) GO TO 100
      QLEN = 0
      LOW = N + 1
      UP = N + 1
C CSP is cost of shortest path to any unassigned row
C ISP is matrix position of unassigned row element in shortest path
C JSP is column index of unassigned row element in shortest path
      CSP = MINONE
C Build shortest path tree starting from unassigned column JORD
      J = JORD
      PR(J) = -1

C Scan column J
      DO 115 K = IP(J),IP(J+1)-1
      I = IRN(K)
      DNEW = ABS(A(K))
      IF (CSP.GE.DNEW) GO TO 115
      IF (IPERM(I).EQ.0) THEN
C Row I is unassigned; update shortest path info
        CSP = DNEW
        ISP = I
        JSP = J
        IF (CSP.GE.BV) GO TO 160
      ELSE
        D(I) = DNEW
        IF (DNEW.GE.BV) THEN
C Add row I to Q2
          LOW = LOW - 1
          Q(LOW) = I
        ELSE
C Add row I to Q, and push it
          QLEN = QLEN + 1
          L(I) = QLEN
          CALL MC64D(I,N,Q,D,L,1)
        ENDIF
        JJ = IPERM(I)
        PR(JJ) = J
      ENDIF
  115   CONTINUE

      DO 150 JDUM = 1,NUM
C If Q2 is empty, extract new rows from Q
      IF (LOW.EQ.UP) THEN
        IF (QLEN.EQ.0) GO TO 160
        I = Q(1)
        IF (CSP.GE.D(I)) GO TO 160
        BV = D(I)
        DO 152 IDUM = 1,N
          CALL MC64E(QLEN,N,Q,D,L,1)
          L(I) = 0
          LOW = LOW - 1
          Q(LOW) = I
          IF (QLEN.EQ.0) GO TO 153
          I = Q(1)
          IF (D(I).NE.BV) GO TO 153
  152     CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C Move row Q0
  153   UP = UP - 1
      Q0 = Q(UP)
      DQ0 = D(Q0)
      L(Q0) = UP
C Scan column that matches with row Q0
      J = IPERM(Q0)
      DO 155 K = IP(J),IP(J+1)-1
        I = IRN(K)
C Update D(I)
        IF (L(I).GE.UP) GO TO 155
        DNEW = MIN(DQ0,ABS(A(K)))
        IF (CSP.GE.DNEW) GO TO 155
        IF (IPERM(I).EQ.0) THEN
C Row I is unassigned; update shortest path info
          CSP = DNEW
          ISP = I
          JSP = J
          IF (CSP.GE.BV) GO TO 160
        ELSE
          DI = D(I)
          IF (DI.GE.BV .OR. DI.GE.DNEW) GO TO 155
          D(I) = DNEW
          IF (DNEW.GE.BV) THEN
C Delete row I from Q (if necessary); add row I to Q2
          IF (DI.NE.MINONE) THEN
            LPOS = L(I)
            CALL MC64F(LPOS,QLEN,N,Q,D,L,1)
          ENDIF
          L(I) = 0
          LOW = LOW - 1
          Q(LOW) = I
          ELSE
C Add row I to Q (if necessary); push row I up Q
          IF (DI.EQ.MINONE) THEN
            QLEN = QLEN + 1
            L(I) = QLEN
          ENDIF
          CALL MC64D(I,N,Q,D,L,1)
          ENDIF
C Update tree
          JJ = IPERM(I)
          PR(JJ) = J
        ENDIF
  155   CONTINUE
  150   CONTINUE

C If CSP = MINONE, no augmenting path is found
  160   IF (CSP.EQ.MINONE) GO TO 190
C Update bottleneck value
      BV = MIN(BV,CSP)
C Find augmenting path by tracing backward in PR; update IPERM,JPERM
      NUM = NUM + 1
      I = ISP
      J = JSP
      DO 170 JDUM = 1,NUM+1
      I0 = JPERM(J)
      JPERM(J) = I
      IPERM(I) = J
      J = PR(J)
      IF (J.EQ.-1) GO TO 190
      I = I0
  170   CONTINUE
C End of dummy loop; this point is never reached
  190   DO 191 KK = UP,N
      I = Q(KK)
      D(I) = MINONE
      L(I) = 0
  191   CONTINUE
      DO 192 KK = LOW,UP-1
      I = Q(KK)
      D(I) = MINONE
  192   CONTINUE
      DO 193 KK = 1,QLEN
      I = Q(KK)
      D(I) = MINONE
      L(I) = 0
  193   CONTINUE

  100 CONTINUE
C End of main loop

C BV is bottleneck value of final matching
    IF (NUM.EQ.N) GO TO 1000

C Matrix is structurally singular, complete IPERM.
C JPERM, PR are work arrays
    DO 300 J = 1,N
      JPERM(J) = 0
  300 CONTINUE
    K = 0
    DO 310 I = 1,N
      IF (IPERM(I).EQ.0) THEN
      K = K + 1
      PR(K) = I
      ELSE
      J = IPERM(I)
      JPERM(J) = I
      ENDIF
  310 CONTINUE
    K = 0
    DO 320 I = 1,N
      IF (JPERM(I).NE.0) GO TO 320
      K = K + 1
      JDUM = PR(K)
      IPERM(JDUM) = - I
  320 CONTINUE

 1000 RETURN
    END

C**********************************************************************
    SUBROUTINE MC64D(I,N,Q,D,L,IWAY)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER I,N,IWAY
    INTEGER Q(N),L(N)
    REAL D(N)

C Variables N,Q,D,L are described in MC64B/BD
C IF IWAY is equal to 1, then
C node I is pushed from its current position upwards
C IF IWAY is not equal to 1, then
C node I is pushed from its current position downwards

C Local variables and parameters
    INTEGER IDUM,K,POS,POSK,QK
    PARAMETER (K=2)
    REAL DI

    POS = L(I)
    IF (POS.LE.1) GO TO 20
    DI = D(I)
C POS is index of current position of I in the tree
    IF (IWAY.EQ.1) THEN
      DO 10 IDUM = 1,N
      POSK = POS/K
      QK = Q(POSK)
      IF (DI.LE.D(QK)) GO TO 20
      Q(POS) = QK
      L(QK) = POS
      POS = POSK
      IF (POS.LE.1) GO TO 20
   10   CONTINUE
C End of dummy loop; this point is never reached
    ELSE
      DO 15 IDUM = 1,N
      POSK = POS/K
      QK = Q(POSK)
      IF (DI.GE.D(QK)) GO TO 20
      Q(POS) = QK
      L(QK) = POS
      POS = POSK
      IF (POS.LE.1) GO TO 20
   15   CONTINUE
C End of dummy loop; this point is never reached
    ENDIF
C End of dummy if; this point is never reached
   20 Q(POS) = I
    L(I) = POS

    RETURN
    END

C**********************************************************************
    SUBROUTINE MC64E(QLEN,N,Q,D,L,IWAY)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER QLEN,N,IWAY
    INTEGER Q(N),L(N)
    REAL D(N)

C Variables QLEN,N,Q,D,L are described in MC64B/BD (IWAY = 1) or
C   MC64W/WD (IWAY = 2)
C The root node is deleted from the binary heap.

C Local variables and parameters
    INTEGER I,IDUM,K,POS,POSK
    PARAMETER (K=2)
    REAL DK,DR,DI

C Move last element to begin of Q
    I = Q(QLEN)
    DI = D(I)
    QLEN = QLEN - 1
    POS = 1
    IF (IWAY.EQ.1) THEN
      DO 10 IDUM = 1,N
      POSK = K*POS
      IF (POSK.GT.QLEN) GO TO 20
      DK = D(Q(POSK))
      IF (POSK.LT.QLEN) THEN
        DR = D(Q(POSK+1))
        IF (DK.LT.DR) THEN
          POSK = POSK + 1
          DK = DR
        ENDIF
      ENDIF
      IF (DI.GE.DK) GO TO 20
C Exchange old last element with larger priority child
      Q(POS) = Q(POSK)
      L(Q(POS)) = POS
      POS = POSK
   10   CONTINUE
C End of dummy loop; this point is never reached
    ELSE
      DO 15 IDUM = 1,N
      POSK = K*POS
      IF (POSK.GT.QLEN) GO TO 20
      DK = D(Q(POSK))
      IF (POSK.LT.QLEN) THEN
        DR = D(Q(POSK+1))
        IF (DK.GT.DR) THEN
          POSK = POSK + 1
          DK = DR
        ENDIF
      ENDIF
      IF (DI.LE.DK) GO TO 20
C Exchange old last element with smaller child
      Q(POS) = Q(POSK)
      L(Q(POS)) = POS
      POS = POSK
   15   CONTINUE
C End of dummy loop; this point is never reached
    ENDIF
C End of dummy if; this point is never reached
   20 Q(POS) = I
    L(I) = POS

    RETURN
    END

C**********************************************************************
    SUBROUTINE MC64F(POS0,QLEN,N,Q,D,L,IWAY)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER POS0,QLEN,N,IWAY
    INTEGER Q(N),L(N)
    REAL D(N)

C Variables QLEN,N,Q,D,L are described in MC64B/BD (IWAY = 1) or
C   MC64W (IWAY = 2).
C Move last element in the heap

    INTEGER I,IDUM,K,POS,POSK,QK
    PARAMETER (K=2)
    REAL DK,DR,DI

C Quick return, if possible
    IF (QLEN.EQ.POS0) THEN
      QLEN = QLEN - 1
      RETURN
    ENDIF

C Move last element from queue Q to position POS0
C POS is current position of node I in the tree
    I = Q(QLEN)
    DI = D(I)
    QLEN = QLEN - 1
    POS = POS0
    IF (IWAY.EQ.1) THEN
      IF (POS.LE.1) GO TO 20
      DO 10 IDUM = 1,N
      POSK = POS/K
      QK = Q(POSK)
      IF (DI.LE.D(QK)) GO TO 20
      Q(POS) = QK
      L(QK) = POS
      POS = POSK
      IF (POS.LE.1) GO TO 20
   10   CONTINUE
C End of dummy loop; this point is never reached
   20   Q(POS) = I
      L(I) = POS
      IF (POS.NE.POS0) RETURN
      DO 30 IDUM = 1,N
      POSK = K*POS
      IF (POSK.GT.QLEN) GO TO 40
      DK = D(Q(POSK))
      IF (POSK.LT.QLEN) THEN
        DR = D(Q(POSK+1))
        IF (DK.LT.DR) THEN
          POSK = POSK + 1
          DK = DR
        ENDIF
      ENDIF
      IF (DI.GE.DK) GO TO 40
      QK = Q(POSK)
      Q(POS) = QK
      L(QK) = POS
      POS = POSK
   30   CONTINUE
C End of dummy loop; this point is never reached
    ELSE
      IF (POS.LE.1) GO TO 34
      DO 32 IDUM = 1,N
      POSK = POS/K
      QK = Q(POSK)
      IF (DI.GE.D(QK)) GO TO 34
      Q(POS) = QK
      L(QK) = POS
      POS = POSK
      IF (POS.LE.1) GO TO 34
   32   CONTINUE
C End of dummy loop; this point is never reached
   34   Q(POS) = I
      L(I) = POS
      IF (POS.NE.POS0) RETURN
      DO 36 IDUM = 1,N
      POSK = K*POS
      IF (POSK.GT.QLEN) GO TO 40
      DK = D(Q(POSK))
      IF (POSK.LT.QLEN) THEN
        DR = D(Q(POSK+1))
        IF (DK.GT.DR) THEN
          POSK = POSK + 1
          DK = DR
        ENDIF
      ENDIF
      IF (DI.LE.DK) GO TO 40
      QK = Q(POSK)
      Q(POS) = QK
      L(QK) = POS
      POS = POSK
   36   CONTINUE
C End of dummy loop; this point is never reached
    ENDIF
C End of dummy if; this point is never reached
   40 Q(POS) = I
    L(I) = POS

    RETURN
    END

C**********************************************************************
    SUBROUTINE MC64R(N,NE,IP,IRN,A)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER N,NE
    INTEGER IP(N+1),IRN(NE)
    REAL A(NE)

C This subroutine sorts the entries in each column of the
C sparse matrix (defined by N,NE,IP,IRN,A) by decreasing
C numerical value.

C Local constants
    INTEGER THRESH,TDLEN
    PARAMETER (THRESH=15,TDLEN=50)
C Local variables
    INTEGER J,IPJ,K,LEN,R,S,HI,FIRST,MID,LAST,TD
    REAL HA,KEY
C Local arrays
    INTEGER TODO(TDLEN)

    DO 100 J = 1,N
      LEN = IP(J+1) - IP(J)
      IF (LEN.LE.1) GO TO 100
      IPJ = IP(J)

C Sort array roughly with partial quicksort
      IF (LEN.LT.THRESH) GO TO 400
      TODO(1) = IPJ
      TODO(2) = IPJ + LEN
      TD = 2
  500   CONTINUE
      FIRST = TODO(TD-1)
      LAST = TODO(TD)
C KEY is the smallest of two values present in interval [FIRST,LAST)
      KEY = A((FIRST+LAST)/2)
      DO 475 K = FIRST,LAST-1
      HA = A(K)
      IF (HA.EQ.KEY) GO TO 475
      IF (HA.GT.KEY) GO TO 470
      KEY = HA
      GO TO 470
  475   CONTINUE
C Only one value found in interval, so it is already sorted
      TD = TD - 2
      GO TO 425

C Reorder interval [FIRST,LAST) such that entries before MID are gt KEY
  470   MID = FIRST
      DO 450 K = FIRST,LAST-1
      IF (A(K).LE.KEY) GO TO 450
      HA = A(MID)
      A(MID) = A(K)
      A(K) = HA
      HI = IRN(MID)
      IRN(MID) = IRN(K)
      IRN(K) = HI
      MID = MID + 1
  450   CONTINUE
C Both subintervals [FIRST,MID), [MID,LAST) are nonempty
C Stack the longest of the two subintervals first
      IF (MID-FIRST.GE.LAST-MID) THEN
      TODO(TD+2) = LAST
      TODO(TD+1) = MID
      TODO(TD) = MID
C      TODO(TD-1) = FIRST
      ELSE
      TODO(TD+2) = MID
      TODO(TD+1) = FIRST
      TODO(TD) = LAST
      TODO(TD-1) = MID
      ENDIF
      TD = TD + 2

  425   CONTINUE
      IF (TD.EQ.0) GO TO 400
C There is still work to be done
      IF (TODO(TD)-TODO(TD-1).GE.THRESH) GO TO 500
C Next interval is already short enough for straightforward insertion
      TD = TD - 2
      GO TO 425

C Complete sorting with straightforward insertion
  400   DO 200 R = IPJ+1,IPJ+LEN-1
      IF (A(R-1) .LT. A(R)) THEN
        HA = A(R)
        HI = IRN(R)
        A(R) = A(R-1)
        IRN(R) = IRN(R-1)
        DO 300 S = R-1,IPJ+1,-1
          IF (A(S-1) .LT. HA) THEN
          A(S) = A(S-1)
          IRN(S) = IRN(S-1)
          ELSE
          A(S) = HA
          IRN(S) = HI
          GO TO 200
          END IF
  300     CONTINUE
        A(IPJ) = HA
        IRN(IPJ) = HI
      END IF
  200   CONTINUE

  100 CONTINUE

    RETURN
    END

C**********************************************************************
    SUBROUTINE MC64S(N,NE,IP,IRN,A,IPERM,NUMX,
   &       W,LEN,LENL,LENH,FC,IW,IW4)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER N,NE,NUMX
    INTEGER IP(N+1),IRN(NE),IPERM(N),
   &      W(N),LEN(N),LENL(N),LENH(N),FC(N),IW(N),IW4(4*N)
    REAL A(NE)

C N, NE, IP, IRN, are described in MC64A/AD.
C A is a REAL array of length NE.
C   A(K), K=1..NE, must be set to the value of the entry that
C   corresponds to IRN(k). The entries in each column must be
C   non-negative and ordered by decreasing value.
C IPERM is an INTEGER array of length N. On exit, it contains the
C   bottleneck matching: IPERM(I) - 0 or row I is matched to column
C   IPERM(I).
C NUMX is an INTEGER variable. On exit, it contains the cardinality
C   of the matching stored in IPERM.
C IW is an INTEGER work array of length 10N.

C FC is an integer array of length N that contains the list of
C   unmatched columns.
C LEN(J), LENL(J), LENH(J) are integer arrays of length N that point
C   to entries in matrix column J.
C   In the matrix defined by the column parts IP(J)+LENL(J) we know
C   a matching does not exist; in the matrix defined by the column
C   parts IP(J)+LENH(J) we know one exists.
C   LEN(J) lies between LENL(J) and LENH(J) and determines the matrix
C   that is tested for a maximum matching.
C W is an integer array of length N and contains the indices of the
C   columns for which LENL ne LENH.
C WLEN is number of indices stored in array W.
C IW is integer work array of length N.
C IW4 is integer work array of length 4N used by MC64U/UD.

    INTEGER NUM,NVAL,WLEN,II,I,J,K,L,CNT,MOD,IDUM1,IDUM2,IDUM3
    REAL BVAL,BMIN,BMAX,RINF
    EXTERNAL MC64Q,MC64U

C BMIN and BMAX are such that a maximum matching exists for the input
C   matrix in which all entries smaller than BMIN are dropped.
C   For BMAX, a maximum matching does not exist.
C BVAL is a value between BMIN and BMAX.
C CNT is the number of calls made to MC64U/UD so far.
C NUM is the cardinality of last matching found.

C Set RINF to largest positive real number
    RINF = HUGE(RINF)

C Compute a first maximum matching from scratch on whole matrix.
    DO 20 J = 1,N
      FC(J) = J
      IW(J) = 0
      LEN(J) = IP(J+1) - IP(J)
   20 CONTINUE
C The first call to MC64U/UD
    CNT = 1
    MOD = 1
    NUMX = 0
    CALL MC64U(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUMX,N,
   &        IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))

C IW contains a maximum matching of length NUMX.
    NUM = NUMX

    IF (NUM.NE.N) THEN
C Matrix is structurally singular
      BMAX = RINF
    ELSE
C Matrix is structurally nonsingular, NUM=NUMX=N;
C Set BMAX just above the smallest of all the maximum absolute
C values of the columns
      BMAX = RINF
      DO 30 J = 1,N
      BVAL = 0.0
      DO 25 K = IP(J),IP(J+1)-1
        IF (A(K).GT.BVAL) BVAL = A(K)
   25   CONTINUE
      IF (BVAL.LT.BMAX) BMAX = BVAL
   30   CONTINUE
      BMAX = 1.001 * BMAX
    ENDIF

C Initialize BVAL,BMIN
    BVAL = 0.0
    BMIN = 0.0
C Initialize LENL,LEN,LENH,W,WLEN according to BMAX.
C Set LEN(J), LENH(J) just after last entry in column J.
C Set LENL(J) just after last entry in column J with value ge BMAX.
    WLEN = 0
    DO 48 J = 1,N
      L = IP(J+1) - IP(J)
      LENH(J) = L
      LEN(J) = L
      DO 45 K = IP(J),IP(J+1)-1
      IF (A(K).LT.BMAX) GO TO 46
   45   CONTINUE
C Column J is empty or all entries are ge BMAX
      K = IP(J+1)
   46   LENL(J) = K - IP(J)
C Add J to W if LENL(J) ne LENH(J)
      IF (LENL(J).EQ.L) GO TO 48
      WLEN = WLEN + 1
      W(WLEN) = J
   48 CONTINUE

C Main loop
    DO 90 IDUM1 = 1,NE
      IF (NUM.EQ.NUMX) THEN
C We have a maximum matching in IW; store IW in IPERM
      DO 50 I = 1,N
        IPERM(I) = IW(I)
   50   CONTINUE
C Keep going round this loop until matching IW is no longer maximum.
      DO 80 IDUM2 = 1,NE
        BMIN = BVAL
        IF (BMAX .EQ. BMIN) GO TO 99
C Find splitting value BVAL
        CALL MC64Q(IP,LENL,LEN,W,WLEN,A,NVAL,BVAL)
        IF (NVAL.LE.1) GO TO 99
C Set LEN such that all matrix entries with value lt BVAL are
C discarded. Store old LEN in LENH. Do this for all columns W(K).
C Each step, either K is incremented or WLEN is decremented.
        K = 1
        DO 70 IDUM3 = 1,N
          IF (K.GT.WLEN) GO TO 71
          J = W(K)
          DO 55 II = IP(J)+LEN(J)-1,IP(J)+LENL(J),-1
          IF (A(II).GE.BVAL) GO TO 60
          I = IRN(II)
          IF (IW(I).NE.J) GO TO 55
C Remove entry from matching
          IW(I) = 0
          NUM = NUM - 1
          FC(N-NUM) = J
   55       CONTINUE
   60       LENH(J) = LEN(J)
C IP(J)+LEN(J)-1 is last entry in column ge BVAL
          LEN(J) = II - IP(J) + 1
C If LENH(J) = LENL(J), remove J from W
          IF (LENL(J).EQ.LENH(J)) THEN
          W(K) = W(WLEN)
          WLEN = WLEN - 1
          ELSE
          K = K + 1
          ENDIF
   70     CONTINUE
   71     IF (NUM.LT.NUMX) GO TO 81
   80   CONTINUE
C End of dummy loop; this point is never reached
C Set mode for next call to MC64U/UD
   81   MOD = 1
      ELSE
C We do not have a maximum matching in IW.
      BMAX = BVAL
C BMIN is the bottleneck value of a maximum matching;
C for BMAX the matching is not maximum, so BMAX>BMIN
C      IF (BMAX .EQ. BMIN) GO TO 99
C Find splitting value BVAL
      CALL MC64Q(IP,LEN,LENH,W,WLEN,A,NVAL,BVAL)
      IF (NVAL.EQ.0. OR. BVAL.EQ.BMIN) GO TO 99
C Set LEN such that all matrix entries with value ge BVAL are
C inside matrix. Store old LEN in LENL. Do this for all columns W(K).
C Each step, either K is incremented or WLEN is decremented.
      K = 1
      DO 87 IDUM3 = 1,N
        IF (K.GT.WLEN) GO TO 88
        J = W(K)
        DO 85 II = IP(J)+LEN(J),IP(J)+LENH(J)-1
          IF (A(II).LT.BVAL) GO TO 86
   85     CONTINUE
   86     LENL(J) = LEN(J)
        LEN(J) = II - IP(J)
        IF (LENL(J).EQ.LENH(J)) THEN
          W(K) = W(WLEN)
          WLEN = WLEN - 1
        ELSE
          K = K + 1
        ENDIF
   87   CONTINUE
C End of dummy loop; this point is never reached
C Set mode for next call to MC64U/UD
   88   MOD = 0
      ENDIF
      CNT = CNT + 1
      CALL MC64U(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUM,NUMX,
   &        IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))

C IW contains maximum matching of length NUM
   90 CONTINUE
C End of dummy loop; this point is never reached

C BMIN is bottleneck value of final matching
   99 IF (NUMX.EQ.N) GO TO 1000
C The matrix is structurally singular, complete IPERM
C W, IW are work arrays
    DO 300 J = 1,N
      W(J) = 0
  300 CONTINUE
    K = 0
    DO 310 I = 1,N
      IF (IPERM(I).EQ.0) THEN
      K = K + 1
      IW(K) = I
      ELSE
      J = IPERM(I)
      W(J) = I
      ENDIF
  310 CONTINUE
    K = 0
    DO 320 J = 1,N
      IF (W(J).NE.0) GO TO 320
      K = K + 1
      IDUM1 = IW(K)
      IPERM(IDUM1) =  - J
  320 CONTINUE

 1000 RETURN
    END

C**********************************************************************
    SUBROUTINE MC64Q(IP,LENL,LENH,W,WLEN,A,NVAL,VAL)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER WLEN,NVAL
    INTEGER IP(*),LENL(*),LENH(*),W(*)
    REAL A(*),VAL

C This routine searches for at most XX different numerical values
C in the columns W(1:WLEN). XX>=2.
C Each column J is scanned between IP(J)+LENL(J) and IP(J)+LENH(J)-1
C until XX values are found or all columns have been considered.
C On output, NVAL is the number of different values that is found
C and SPLIT(1:NVAL) contains the values in decreasing order.
C If NVAL > 0, the routine returns VAL = SPLIT((NVAL+1)/2).
C
    INTEGER XX,J,K,II,S,POS
    PARAMETER (XX=10)
    REAL SPLIT(XX),HA

C Scan columns in W(1:WLEN). For each encountered value, if value not
C already present in SPLIT(1:NVAL), insert value such that SPLIT
C remains sorted by decreasing value.
C The sorting is done by straightforward insertion; therefore the use
C of this routine should be avoided for large XX (XX < 20).
    NVAL = 0
    DO 10 K = 1,WLEN
      J = W(K)
      DO 15 II = IP(J)+LENL(J),IP(J)+LENH(J)-1
      HA = A(II)
      IF (NVAL.EQ.0) THEN
        SPLIT(1) = HA
        NVAL = 1
      ELSE
C Check presence of HA in SPLIT
        DO 20 S = NVAL,1,-1
          IF (SPLIT(S).EQ.HA) GO TO 15
          IF (SPLIT(S).GT.HA) THEN
          POS = S + 1
          GO TO 21
          ENDIF
  20      CONTINUE
        POS = 1
C The insertion
  21      DO 22 S = NVAL,POS,-1
          SPLIT(S+1) = SPLIT(S)
  22      CONTINUE
        SPLIT(POS) = HA
        NVAL = NVAL + 1
      ENDIF
C Exit loop if XX values are found
      IF (NVAL.EQ.XX) GO TO 11
  15    CONTINUE
  10  CONTINUE
C Determine VAL
  11  IF (NVAL.GT.0) VAL = SPLIT((NVAL+1)/2)

    RETURN
    END

C**********************************************************************
    SUBROUTINE MC64U(ID,MOD,N,IRN,LIRN,IP,LENC,FC,IPERM,NUM,NUMX,
   &       PR,ARP,CV,OUT)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER ID,MOD,N,LIRN,NUM,NUMX
    INTEGER ARP(N),CV(N),IRN(LIRN),IP(N),
   &      FC(N),IPERM(N),LENC(N),OUT(N),PR(N)

C PR(J) is the previous column to J in the depth first search.
C   Array PR is used as workspace in the sorting algorithm.
C Elements (I,IPERM(I)) I=1,..,N are entries at the end of the
C   algorithm unless N assignments have not been made in which case
C   N-NUM pairs (I,IPERM(I)) will not be entries in the matrix.
C CV(I) is the most recent loop number (ID+JORD) at which row I
C   was visited.
C ARP(J) is the number of entries in column J which have been scanned
C   when looking for a cheap assignment.
C OUT(J) is one less than the number of entries in column J which have
C   not been scanned during one pass through the main loop.
C NUMX is maximum possible size of matching.

    INTEGER I,II,IN1,IN2,J,J1,JORD,K,KK,LAST,NFC,
   &      NUM0,NUM1,NUM2,ID0,ID1

    IF (ID.EQ.1) THEN
C The first call to MC64U/UD.
C Initialize CV and ARP; parameters MOD, NUMX are not accessed
      DO 5 I = 1,N
      CV(I) = 0
      ARP(I) = 0
    5   CONTINUE
      NUM1 = N
      NUM2 = N
    ELSE
C Not the first call to MC64U/UD.
C Re-initialize ARP if entries were deleted since last call to MC64U/UD
      IF (MOD.EQ.1) THEN
      DO 8 I = 1,N
        ARP(I) = 0
    8   CONTINUE
      ENDIF
      NUM1 = NUMX
      NUM2 = N - NUMX
    ENDIF
    NUM0 = NUM

C NUM0 is size of input matching
C NUM1 is maximum possible size of matching
C NUM2 is maximum allowed number of unassigned rows/columns
C NUM is size of current matching

C Quick return if possible
C    IF (NUM.EQ.N) GO TO 199
C NFC is number of rows/columns that could not be assigned
    NFC = 0
C Integers ID0+1 to ID0+N are unique numbers for call ID to MC64U/UD,
C so 1st call uses 1..N, 2nd call uses N+1..2N, etc
    ID0 = (ID-1)*N

C Main loop. Each pass round this loop either results in a new
C assignment or gives a column with no assignment

    DO 100 JORD = NUM0+1,N

C Each pass uses unique number ID1
      ID1 = ID0 + JORD
C J is unmatched column
      J = FC(JORD-NUM0)
      PR(J) = -1
      DO 70 K = 1,JORD
C Look for a cheap assignment
      IF (ARP(J).GE.LENC(J)) GO TO 30
      IN1 = IP(J) + ARP(J)
      IN2 = IP(J) + LENC(J) - 1
      DO 20 II = IN1,IN2
        I = IRN(II)
        IF (IPERM(I).EQ.0) GO TO 80
   20   CONTINUE
C No cheap assignment in row
      ARP(J) = LENC(J)
C Begin looking for assignment chain starting with row J
   30   OUT(J) = LENC(J) - 1
C Inner loop.  Extends chain by one or backtracks
      DO 60 KK = 1,JORD
        IN1 = OUT(J)
        IF (IN1.LT.0) GO TO 50
        IN2 = IP(J) + LENC(J) - 1
        IN1 = IN2 - IN1
C Forward scan
        DO 40 II = IN1,IN2
          I = IRN(II)
          IF (CV(I).EQ.ID1) GO TO 40
C Column J has not yet been accessed during this pass
          J1 = J
          J = IPERM(I)
          CV(I) = ID1
          PR(J) = J1
          OUT(J1) = IN2 - II - 1
          GO TO 70
   40     CONTINUE
C Backtracking step.
   50     J1 = PR(J)
        IF (J1.EQ.-1) THEN
C No augmenting path exists for column J.
          NFC = NFC + 1
          FC(NFC) = J
          IF (NFC.GT.NUM2) THEN
C A matching of maximum size NUM1 is not possible
          LAST = JORD
          GO TO 101
          ENDIF
          GO TO 100
        ENDIF
        J = J1
   60   CONTINUE
C End of dummy loop; this point is never reached
   70   CONTINUE
C End of dummy loop; this point is never reached

C New assignment is made.
   80   IPERM(I) = J
      ARP(J) = II - IP(J) + 1
      NUM = NUM + 1
      DO 90 K = 1,JORD
      J = PR(J)
      IF (J.EQ.-1) GO TO 95
      II = IP(J) + LENC(J) - OUT(J) - 2
      I = IRN(II)
      IPERM(I) = J
   90   CONTINUE
C End of dummy loop; this point is never reached

   95   IF (NUM.EQ.NUM1) THEN
C A matching of maximum size NUM1 is found
      LAST = JORD
      GO TO 101
      ENDIF
C
  100 CONTINUE

C All unassigned columns have been considered
    LAST = N

C Now, a transversal is computed or is not possible.
C Complete FC before returning.
  101 DO 110 JORD = LAST+1,N
      NFC = NFC + 1
      FC(NFC) = FC(JORD-NUM0)
  110 CONTINUE

C  199 RETURN
    RETURN
    END

C**********************************************************************
    SUBROUTINE MC64W(N,NE,IP,IRN,A,IPERM,NUM,
   &       JPERM,OUT,PR,Q,L,U,D)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER N,NE,NUM
    INTEGER IP(N+1),IRN(NE),IPERM(N),
   &      JPERM(N),OUT(N),PR(N),Q(N),L(N)
    REAL A(NE),U(N),D(N)

C N, NE, IP, IRN are described in MC64A/AD.
C A is a REAL array of length NE.
C   A(K), K=1..NE, must be set to the value of the entry that
C   corresponds to IRN(K). It is not altered.
C   All values A(K) must be non-negative.
C IPERM is an INTEGER array of length N. On exit, it contains the
C   weighted matching: IPERM(I) = 0 or row I is matched to column
C   IPERM(I).
C NUM is an INTEGER variable. On exit, it contains the cardinality of
C   the matching stored in IPERM.
C IW is an INTEGER work array of length 5N.
C DW is a REAL array of length 2N.
C   On exit, U = D(1:N) contains the dual row variable and
C   V = D(N+1:2N) contains the dual column variable. If the matrix
C   is structurally nonsingular (NUM = N), the following holds:
C    U(I)+V(J) <= A(I,J)  if IPERM(I) |= J
C    U(I)+V(J)  = A(I,J)  if IPERM(I)  = J
C    U(I) = 0  if IPERM(I) = 0
C    V(J) = 0  if there is no I for which IPERM(I) = J

C Local variables
    INTEGER I,I0,II,J,JJ,JORD,Q0,QLEN,JDUM,ISP,JSP,
   &      K,K0,K1,K2,KK,KK1,KK2,UP,LOW,LPOS
    REAL CSP,DI,DMIN,DNEW,DQ0,VJ
C Local parameters
    REAL RINF,ZERO
    PARAMETER (ZERO=0.0D+0)
C External subroutines and/or functions
    EXTERNAL MC64D,MC64E,MC64F


C Set RINF to largest positive real number
    RINF = HUGE(RINF)

C Initialization
    NUM = 0
    DO 10 K = 1,N
      U(K) = RINF
      D(K) = ZERO
      IPERM(K) = 0
      JPERM(K) = 0
      PR(K) = IP(K)
      L(K) = 0
   10 CONTINUE
C Initialize U(I)
    DO 30 J = 1,N
      DO 20 K = IP(J),IP(J+1)-1
      I = IRN(K)
      IF (A(K).GT.U(I)) GO TO 20
      U(I) = A(K)
      IPERM(I) = J
      L(I) = K
   20   CONTINUE
   30 CONTINUE
    DO 40 I = 1,N
      J = IPERM(I)
      IF (J.EQ.0) GO TO 40
C Row I is not empty
      IPERM(I) = 0
      IF (JPERM(J).NE.0) GO TO 40
C Don't choose cheap assignment from dense columns
      IF (IP(J+1)-IP(J) .GT. N/10 .AND. N.GT.50) GO TO 40
C Assignment of column J to row I
      NUM = NUM + 1
      IPERM(I) = J
      JPERM(J) = L(I)
   40 CONTINUE
    IF (NUM.EQ.N) GO TO 1000
C Scan unassigned columns; improve assignment
    DO 95 J = 1,N
C JPERM(J) ne 0 iff column J is already assigned
      IF (JPERM(J).NE.0) GO TO 95
      K1 = IP(J)
      K2 = IP(J+1) - 1
C Continue only if column J is not empty
      IF (K1.GT.K2) GO TO 95
C     VJ = RINF
C Changes made to allow for NaNs
      I0 = IRN(K1)
      VJ = A(K1) - U(I0)
      K0 = K1
      DO 50 K = K1+1,K2
      I = IRN(K)
      DI = A(K) - U(I)
      IF (DI.GT.VJ) GO TO 50
      IF (DI.LT.VJ .OR. DI.EQ.RINF) GO TO 55
      IF (IPERM(I).NE.0 .OR. IPERM(I0).EQ.0) GO TO 50
   55   VJ = DI
      I0 = I
      K0 = K
   50   CONTINUE
      D(J) = VJ
      K = K0
      I = I0
      IF (IPERM(I).EQ.0) GO TO 90
      DO 60 K = K0,K2
      I = IRN(K)
      IF (A(K)-U(I).GT.VJ) GO TO 60
      JJ = IPERM(I)
C Scan remaining part of assigned column JJ
      KK1 = PR(JJ)
      KK2 = IP(JJ+1) - 1
      IF (KK1.GT.KK2) GO TO 60
      DO 70 KK = KK1,KK2
        II = IRN(KK)
        IF (IPERM(II).GT.0) GO TO 70
        IF (A(KK)-U(II).LE.D(JJ)) GO TO 80
   70   CONTINUE
      PR(JJ) = KK2 + 1
   60   CONTINUE
      GO TO 95
   80   JPERM(JJ) = KK
      IPERM(II) = JJ
      PR(JJ) = KK + 1
   90   NUM = NUM + 1
      JPERM(J) = K
      IPERM(I) = J
      PR(J) = K + 1
   95 CONTINUE
    IF (NUM.EQ.N) GO TO 1000

C Prepare for main loop
    DO 99 I = 1,N
      D(I) = RINF
      L(I) = 0
   99 CONTINUE

C Main loop ... each pass round this loop is similar to Dijkstra's
C algorithm for solving the single source shortest path problem

    DO 100 JORD = 1,N

      IF (JPERM(JORD).NE.0) GO TO 100
C JORD is next unmatched column
C DMIN is the length of shortest path in the tree
      DMIN = RINF
      QLEN = 0
      LOW = N + 1
      UP = N + 1
C CSP is the cost of the shortest augmenting path to unassigned row
C IRN(ISP). The corresponding column index is JSP.
      CSP = RINF
C Build shortest path tree starting from unassigned column (root) JORD
      J = JORD
      PR(J) = -1

C Scan column J
      DO 115 K = IP(J),IP(J+1)-1
      I = IRN(K)
      DNEW = A(K) - U(I)
      IF (DNEW.GE.CSP) GO TO 115
      IF (IPERM(I).EQ.0) THEN
        CSP = DNEW
        ISP = K
        JSP = J
      ELSE
        IF (DNEW.LT.DMIN) DMIN = DNEW
        D(I) = DNEW
        QLEN = QLEN + 1
        Q(QLEN) = K
      ENDIF
  115   CONTINUE
C Initialize heap Q and Q2 with rows held in Q(1:QLEN)
      Q0 = QLEN
      QLEN = 0
      DO 120 KK = 1,Q0
      K = Q(KK)
      I = IRN(K)
      IF (CSP.LE.D(I)) THEN
        D(I) = RINF
        GO TO 120
      ENDIF
      IF (D(I).LE.DMIN) THEN
        LOW = LOW - 1
        Q(LOW) = I
        L(I) = LOW
      ELSE
        QLEN = QLEN + 1
        L(I) = QLEN
        CALL MC64D(I,N,Q,D,L,2)
      ENDIF
C Update tree
      JJ = IPERM(I)
      OUT(JJ) = K
      PR(JJ) = J
  120   CONTINUE

      DO 150 JDUM = 1,NUM

C If Q2 is empty, extract rows from Q
      IF (LOW.EQ.UP) THEN
        IF (QLEN.EQ.0) GO TO 160
        I = Q(1)
        IF (D(I).GE.CSP) GO TO 160
        DMIN = D(I)
  152     CALL MC64E(QLEN,N,Q,D,L,2)
        LOW = LOW - 1
        Q(LOW) = I
        L(I) = LOW
        IF (QLEN.EQ.0) GO TO 153
        I = Q(1)
        IF (D(I).GT.DMIN) GO TO 153
        GO TO 152
      ENDIF
C Q0 is row whose distance D(Q0) to the root is smallest
  153   Q0 = Q(UP-1)
      DQ0 = D(Q0)
C Exit loop if path to Q0 is longer than the shortest augmenting path
      IF (DQ0.GE.CSP) GO TO 160
      UP = UP - 1

C Scan column that matches with row Q0
      J = IPERM(Q0)
      VJ = DQ0 - A(JPERM(J)) + U(Q0)
      DO 155 K = IP(J),IP(J+1)-1
        I = IRN(K)
        IF (L(I).GE.UP) GO TO 155
C DNEW is new cost
        DNEW = VJ + A(K)-U(I)
C Do not update D(I) if DNEW ge cost of shortest path
        IF (DNEW.GE.CSP) GO TO 155
        IF (IPERM(I).EQ.0) THEN
C Row I is unmatched; update shortest path info
          CSP = DNEW
          ISP = K
          JSP = J
        ELSE
C Row I is matched; do not update D(I) if DNEW is larger
          DI = D(I)
          IF (DI.LE.DNEW) GO TO 155
          IF (L(I).GE.LOW) GO TO 155
          D(I) = DNEW
          IF (DNEW.LE.DMIN) THEN
          LPOS = L(I)
          IF (LPOS.NE.0)
   *        CALL MC64F(LPOS,QLEN,N,Q,D,L,2)
          LOW = LOW - 1
          Q(LOW) = I
          L(I) = LOW
          ELSE
          IF (L(I).EQ.0) THEN
            QLEN = QLEN + 1
            L(I) = QLEN
          ENDIF
          CALL MC64D(I,N,Q,D,L,2)
          ENDIF
C Update tree
          JJ = IPERM(I)
          OUT(JJ) = K
          PR(JJ) = J
        ENDIF
  155   CONTINUE
  150   CONTINUE

C If CSP = RINF, no augmenting path is found
  160   IF (CSP.EQ.RINF) GO TO 190
C Find augmenting path by tracing backward in PR; update IPERM,JPERM
      NUM = NUM + 1
      I = IRN(ISP)
      IPERM(I) = JSP
      JPERM(JSP) = ISP
      J = JSP
      DO 170 JDUM = 1,NUM
      JJ = PR(J)
      IF (JJ.EQ.-1) GO TO 180
      K = OUT(J)
      I = IRN(K)
      IPERM(I) = JJ
      JPERM(JJ) = K
      J = JJ
  170   CONTINUE
C End of dummy loop; this point is never reached

C Update U for rows in Q(UP:N)
  180   DO 185 KK = UP,N
      I = Q(KK)
      U(I) = U(I) + D(I) - CSP
  185   CONTINUE
  190   DO 191 KK = LOW,N
      I = Q(KK)
      D(I) = RINF
      L(I) = 0
  191   CONTINUE
      DO 193 KK = 1,QLEN
      I = Q(KK)
      D(I) = RINF
      L(I) = 0
  193   CONTINUE

  100 CONTINUE
C End of main loop


C Set dual column variable in D(1:N)
 1000 DO 200 J = 1,N
      K = JPERM(J)
      IF (K.NE.0) THEN
      D(J) = A(K) - U(IRN(K))
      ELSE
      D(J) = ZERO
      ENDIF
      IF (IPERM(J).EQ.0) U(J) = ZERO
  200 CONTINUE

    IF (NUM.EQ.N) GO TO 1100

C The matrix is structurally singular, complete IPERM.
C JPERM, OUT are work arrays
    DO 300 J = 1,N
      JPERM(J) = 0
  300 CONTINUE
    K = 0
    DO 310 I = 1,N
      IF (IPERM(I).EQ.0) THEN
      K = K + 1
      OUT(K) = I
      ELSE
      J = IPERM(I)
      JPERM(J) = I
      ENDIF
  310 CONTINUE
    K = 0
    DO 320 J = 1,N
      IF (JPERM(J).NE.0) GO TO 320
      K = K + 1
      JDUM = OUT(K)
      IPERM(JDUM) = - J
  320 CONTINUE
 1100 RETURN
    END


* COPYRIGHT (c) 1993 AEA Technology and
* Council for the Central Laboratory of the Research Councils
C Original date 18 June 2001
C  June 2001: threadsafe version of MF41
C
C 12th July 2004 Version 1.0.0. Version numbering added.
C 3 December 2009. Version 1.1.0. CMPLX replaced by DCMPLX in D version.

    SUBROUTINE MF71A(N,KASE,X,EST,W,KEEP)
C
C    MF71AD ESTIMATES THE 1-NORM OF A SQUARE, DOUBLE PRECISION
C    COMPLEX MATRIX A.
C    REVERSE COMMUNICATION IS USED FOR EVALUATING
C    MATRIX-VECTOR PRODUCTS.
C
C
C       N     INTEGER
C           THE ORDER OF THE MATRIX.  N .GE. 1.
C
C       KASE    INTEGER
C           SET INITIALLY TO ZERO . IF N .LE. 0 SET TO -1
C           ON INTERMEDIATE RETURN
C           = 1 OR 2.
C          ON FINAL RETURN
C           =  0  ,IF SUCCESS
C           = -1  ,IF N .LE.0
C
C       X     COMPLEX*16 ARRAY OF DIMENSION (N)
C           IF 1-NORM IS REQUIRED
C           MUST BE OVERWRITTEN BY
C
C              A*X,         IF KASE=1,
C              CTRANS(A)*X,   IF KASE=2,
C
C           WHERE  CTRANS(A)  IS THE CONJUGATE TRANSPOSE,
C           AND MF71 MUST BE RE-CALLED, WITH ALL THE OTHER
C           PARAMETERS UNCHANGED.
C           IF INFINITY-NORM IS REQUIRED
C           MUST BE OVERWRITTEN BY
C
C              CTRANS(A)*X,   IF KASE=1,
C              A*X,         IF KASE=2,
C
C           AND MF71 MUST BE RE-CALLED, WITH ALL THE OTHER
C           PARAMETERS UNCHANGED.
C
C       EST   DOUBLE PRECISION
C           CONTAINS AN ESTIMATE (A LOWER BOUND) FOR NORM(A).
C
C       W     COMPLEX*16 ARRAY OF DIMENSION (N)
C           = A*V,   WHERE  EST = NORM(W)/NORM(V)
C                (V  IS NOT RETURNED).
C
C       KEEP    INTEGER ARRAY LENGTH 5 USED TO PRESERVE PRIVATE
C           DATA, JUMP, ITER, J AND JLAST BETWEEN CALLS,
C           KEEP(5) IS SPARE.
C
C    REFERENCE
C    N.J. HIGHAM (1987) FORTRAN CODES FOR ESTIMATING
C    THE ONE-NORM OF A
C    REAL OR COMPLEX MATRIX, WITH APPLICATIONS
C    TO CONDITION  ESTIMATION, ACM TRANS MATH SOFTWARE,VOL 14,
C    1988, PP381-396.
C
C   .. Parameters ..
    INTEGER ITMAX
    PARAMETER (ITMAX=5)
    REAL ZERO,ONE,TWO
    PARAMETER (ZERO=0.0,ONE=1.0,TWO=2.0)
    COMPLEX CZERO,CONE
    PARAMETER (CZERO= (0.0,0.0),CONE= (1.0,0.0))
C   ..
C   .. Scalar Arguments ..
    REAL EST
    INTEGER KASE,N
C   ..
C   .. Array Arguments ..
    COMPLEX W(N),X(N)
    INTEGER KEEP(5)
C   ..
C   .. Local Scalars ..
    REAL ALTSGN,TEMP,XMAX
    INTEGER I,ITER,J,JLAST,JUMP
C   ..
C   .. Intrinsic Functions ..
    INTRINSIC ABS,CMPLX,REAL
C   ..
C
    IF (N.LE.0) THEN
       KASE = -1
       RETURN

    END IF

    IF (KASE.EQ.0) THEN
       DO 10 I = 1,N
        X(I) = CMPLX(ONE/REAL(N))
   10    CONTINUE
       KASE = 1
       JUMP = 1
       KEEP(1) = JUMP
       KEEP(2) = 0
       KEEP(3) = 0
       KEEP(4) = 0
       RETURN

    END IF
C
    JUMP  = KEEP(1)
    ITER  = KEEP(2)
    J   = KEEP(3)
    JLAST = KEEP(4)
C
    GO TO (20,40,80,110,150) JUMP
C
C    ................ ENTRY   (JUMP = 1)
C
   20 CONTINUE
    IF (N.EQ.1) THEN
       W(1) = X(1)
       EST = ABS(W(1))
C       ... QUIT
       GO TO 180

    END IF
C
    DO 30 I = 1,N
       IF (ABS(X(I)).EQ.ZERO) THEN
        X(I) = CONE
       ELSE
        X(I) = X(I)/CMPLX(ABS(X(I)))
       END IF
   30 CONTINUE
    KASE = 2
    JUMP = 2
    GO TO 1010
C
C    ................ ENTRY   (JUMP = 2)
C   FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
C
   40 CONTINUE
C FIND MAX ENTRY IN X
    J = 1
    XMAX = ABS(X(J))
    DO 50 I = 2,N
       IF (ABS(X(I)).LE.XMAX) GO TO 50
       J = I
       XMAX = ABS(X(J))
   50 CONTINUE
    ITER = 2
C
C    MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
C
   60 CONTINUE
    DO 70 I = 1,N
       X(I) = CZERO
   70 CONTINUE
    X(J) = CONE
    KASE = 1
    JUMP = 3
    GO TO 1010
C
C    ................ ENTRY   (JUMP = 3)
C   X HAS BEEN OVERWRITTEN BY A*X.
C
   80 CONTINUE
C
C   COPY X INTO W
    DO 90 I = 1,N
       W(I) = X(I)
   90 CONTINUE
C
    DO 100 I = 1,N
       IF (ABS(X(I)).EQ.ZERO) THEN
        X(I) = CONE
       ELSE
        X(I) = X(I)/CMPLX(ABS(X(I)))
       END IF
  100 CONTINUE
    KASE = 2
    JUMP = 4
    GO TO 1010
C
C    ................ ENTRY   (JUMP = 4)
C
  110 CONTINUE
    JLAST = J
C FIND MAX ENTRY IN X
    J = 1
    XMAX = ABS(X(J))
    DO 120 I = 2,N
       IF (ABS(X(I)).LE.XMAX) GO TO 120
       J = I
       XMAX = ABS(X(J))
  120 CONTINUE
    IF (ABS(X(JLAST)).NE.ABS(X(J)) .AND.
   +    (ITER.LT.ITMAX)) THEN
       ITER = ITER + 1
       GO TO 60

    END IF
C
C    ITERATION COMPLETE.  FINAL STAGE.
C
    EST = ZERO
    DO 130 I = 1,N
       EST = EST + ABS(W(I))
  130 CONTINUE
C
    ALTSGN = ONE
    DO 140 I = 1,N
       X(I) = CMPLX(ALTSGN* (ONE+REAL(I-1)/REAL(N-1)))
       ALTSGN = -ALTSGN
  140 CONTINUE
    KASE = 1
    JUMP = 5
    GO TO 1010
C
C    ................ ENTRY   (JUMP = 5)
C
  150 CONTINUE
    TEMP = ZERO
    DO 160 I = 1,N
       TEMP = TEMP + ABS(X(I))
  160 CONTINUE
    TEMP = TWO*TEMP/REAL(3*N)
    IF (TEMP.GT.EST) THEN
C
C    COPY X INTO W
C
       DO 170 I = 1,N
        W(I) = X(I)
  170    CONTINUE
       EST = TEMP
    END IF
C
  180 KASE = 0
C
 1010 CONTINUE
    KEEP(1) = JUMP
    KEEP(2) = ITER
    KEEP(3) = J
    KEEP(4) = JLAST
    RETURN
C
    END
* COPYRIGHT (c) 1977 AEA Technology
* Original date 8 Oct 1992
C######8/10/92 Toolpack tool decs employed.
C######8/10/92 D version created by name change only.
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

    SUBROUTINE MC21AD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
C   .. Scalar Arguments ..
    INTEGER LICN,N,NUMNZ
C   ..
C   .. Array Arguments ..
    INTEGER ICN(LICN),IP(N),IPERM(N),IW(N,4),LENR(N)
C   ..
C   .. External Subroutines ..
    EXTERNAL MC21BD
C   ..
C   .. Executable Statements ..
    CALL MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),
   +        IW(1,3),IW(1,4))
    RETURN
C
    END
    SUBROUTINE MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
C   PR(I) IS THE PREVIOUS ROW TO I IN THE DEPTH FIRST SEARCH.
C IT IS USED AS A WORK ARRAY IN THE SORTING ALGORITHM.
C   ELEMENTS (IPERM(I),I) I=1, ... N  ARE NON-ZERO AT THE END OF THE
C ALGORITHM UNLESS N ASSIGNMENTS HAVE NOT BEEN MADE.  IN WHICH CASE
C (IPERM(I),I) WILL BE ZERO FOR N-NUMNZ ENTRIES.
C   CV(I) IS THE MOST RECENT ROW EXTENSION AT WHICH COLUMN I
C WAS VISITED.
C   ARP(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED WHEN LOOKING FOR A CHEAP ASSIGNMENT.
C   OUT(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED DURING ONE PASS THROUGH THE MAIN LOOP.
C
C   INITIALIZATION OF ARRAYS.
C   .. Scalar Arguments ..
    INTEGER LICN,N,NUMNZ
C   ..
C   .. Array Arguments ..
    INTEGER ARP(N),CV(N),ICN(LICN),IP(N),IPERM(N),LENR(N),OUT(N),PR(N)
C   ..
C   .. Local Scalars ..
    INTEGER I,II,IN1,IN2,IOUTK,J,J1,JORD,K,KK
C   ..
C   .. Executable Statements ..
    DO 10 I = 1,N
      ARP(I) = LENR(I) - 1
      CV(I) = 0
      IPERM(I) = 0
   10 CONTINUE
    NUMNZ = 0
C
C
C   MAIN LOOP.
C   EACH PASS ROUND THIS LOOP EITHER RESULTS IN A NEW ASSIGNMENT
C OR GIVES A ROW WITH NO ASSIGNMENT.
    DO 100 JORD = 1,N
      J = JORD
      PR(J) = -1
      DO 70 K = 1,JORD
C LOOK FOR A CHEAP ASSIGNMENT
      IN1 = ARP(J)
      IF (IN1.LT.0) GO TO 30
      IN2 = IP(J) + LENR(J) - 1
      IN1 = IN2 - IN1
      DO 20 II = IN1,IN2
        I = ICN(II)
        IF (IPERM(I).EQ.0) GO TO 80
   20   CONTINUE
C   NO CHEAP ASSIGNMENT IN ROW.
      ARP(J) = -1
C   BEGIN LOOKING FOR ASSIGNMENT CHAIN STARTING WITH ROW J.
   30   CONTINUE
      OUT(J) = LENR(J) - 1
C INNER LOOP.  EXTENDS CHAIN BY ONE OR BACKTRACKS.
      DO 60 KK = 1,JORD
        IN1 = OUT(J)
        IF (IN1.LT.0) GO TO 50
        IN2 = IP(J) + LENR(J) - 1
        IN1 = IN2 - IN1
C FORWARD SCAN.
        DO 40 II = IN1,IN2
          I = ICN(II)
          IF (CV(I).EQ.JORD) GO TO 40
C   COLUMN I HAS NOT YET BEEN ACCESSED DURING THIS PASS.
          J1 = J
          J = IPERM(I)
          CV(I) = JORD
          PR(J) = J1
          OUT(J1) = IN2 - II - 1
          GO TO 70
C
   40     CONTINUE
C
C   BACKTRACKING STEP.
   50     CONTINUE
        J = PR(J)
        IF (J.EQ.-1) GO TO 100
   60   CONTINUE
C
   70   CONTINUE
C
C   NEW ASSIGNMENT IS MADE.
   80   CONTINUE
      IPERM(I) = J
      ARP(J) = IN2 - II - 1
      NUMNZ = NUMNZ + 1
      DO 90 K = 1,JORD
      J = PR(J)
      IF (J.EQ.-1) GO TO 100
      II = IP(J) + LENR(J) - OUT(J) - 2
      I = ICN(II)
      IPERM(I) = J
   90   CONTINUE
C
  100 CONTINUE
C
C   IF MATRIX IS STRUCTURALLY SINGULAR, WE NOW COMPLETE THE
C PERMUTATION IPERM.
    IF (NUMNZ.EQ.N) RETURN
    DO 110 I = 1,N
      ARP(I) = 0
  110 CONTINUE
    K = 0
    DO 130 I = 1,N
      IF (IPERM(I).NE.0) GO TO 120
      K = K + 1
      OUT(K) = I
      GO TO 130
C
  120   CONTINUE
      J = IPERM(I)
      ARP(J) = I
  130 CONTINUE
    K = 0
    DO 140 I = 1,N
      IF (ARP(I).NE.0) GO TO 140
      K = K + 1
      IOUTK = OUT(K)
      IPERM(IOUTK) = I
  140 CONTINUE
    RETURN
C
    END
* COPYRIGHT (c) 1976 AEA Technology
* Original date 21 Jan 1993
C 8 August 2000: CONTINUEs given to DOs.
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C
C 12th July 2004 Version 1.0.0. Version numbering added.

    SUBROUTINE MC22AD(N,ICN,A,NZ,LENROW,IP,IQ,IW,IW1)
C   .. Scalar Arguments ..
    INTEGER N,NZ
C   ..
C   .. Array Arguments ..
    DOUBLE PRECISION A(NZ)
    INTEGER ICN(NZ),IP(N),IQ(N),IW(N,2),IW1(NZ),LENROW(N)
C   ..
C   .. Local Scalars ..
    DOUBLE PRECISION AVAL
    INTEGER I,ICHAIN,IOLD,IPOS,J,J2,JJ,JNUM,JVAL,LENGTH,NEWPOS
C   ..
C   .. Intrinsic Functions ..
    INTRINSIC IABS
C   ..
C   .. Executable Statements ..
    IF (NZ.LE.0) GO TO 1000
    IF (N.LE.0) GO TO 1000
C SET START OF ROW I IN IW(I,1) AND LENROW(I) IN IW(I,2)
    IW(1,1) = 1
    IW(1,2) = LENROW(1)
    DO 10 I = 2,N
      IW(I,1) = IW(I-1,1) + LENROW(I-1)
      IW(I,2) = LENROW(I)
   10 CONTINUE
C PERMUTE LENROW ACCORDING TO IP.  SET OFF-SETS FOR NEW POSITION
C   OF ROW IOLD IN IW(IOLD,1) AND PUT OLD ROW INDICES IN IW1 IN
C   POSITIONS CORRESPONDING TO THE NEW POSITION OF THIS ROW IN A/ICN.
    JJ = 1
    DO 20 I = 1,N
      IOLD = IP(I)
      IOLD = IABS(IOLD)
      LENGTH = IW(IOLD,2)
      LENROW(I) = LENGTH
      IF (LENGTH.EQ.0) GO TO 20
      IW(IOLD,1) = IW(IOLD,1) - JJ
      J2 = JJ + LENGTH - 1
      DO 15 J = JJ,J2
      IW1(J) = IOLD
   15   CONTINUE
      JJ = J2 + 1
   20 CONTINUE
C SET INVERSE PERMUTATION TO IQ IN IW(.,2).
    DO 30 I = 1,N
      IOLD = IQ(I)
      IOLD = IABS(IOLD)
      IW(IOLD,2) = I
   30 CONTINUE
C PERMUTE A AND ICN IN PLACE, CHANGING TO NEW COLUMN NUMBERS.
C
C ***   MAIN LOOP   ***
C EACH PASS THROUGH THIS LOOP PLACES A CLOSED CHAIN OF COLUMN INDICES
C   IN THEIR NEW (AND FINAL) POSITIONS ... THIS IS RECORDED BY
C   SETTING THE IW1 ENTRY TO ZERO SO THAT ANY WHICH ARE SUBSEQUENTLY
C   ENCOUNTERED DURING THIS MAJOR SCAN CAN BE BYPASSED.
    DO 200 I = 1,NZ
      IOLD = IW1(I)
      IF (IOLD.EQ.0) GO TO 200
      IPOS = I
      JVAL = ICN(I)
C IF ROW IOLD IS IN SAME POSITIONS AFTER PERMUTATION GO TO 150.
      IF (IW(IOLD,1).EQ.0) GO TO 150
      AVAL = A(I)
C **  CHAIN LOOP  **
C EACH PASS THROUGH THIS LOOP PLACES ONE (PERMUTED) COLUMN INDEX
C   IN ITS FINAL POSITION  .. VIZ. IPOS.
      DO 100 ICHAIN = 1,NZ
C NEWPOS IS THE ORIGINAL POSITION IN A/ICN OF THE ELEMENT TO BE PLACED
C IN POSITION IPOS.  IT IS ALSO THE POSITION OF THE NEXT ELEMENT IN
C   THE CHAIN.
      NEWPOS = IPOS + IW(IOLD,1)
C IS CHAIN COMPLETE ?
      IF (NEWPOS.EQ.I) GO TO 130
      A(IPOS) = A(NEWPOS)
      JNUM = ICN(NEWPOS)
      ICN(IPOS) = IW(JNUM,2)
      IPOS = NEWPOS
      IOLD = IW1(IPOS)
      IW1(IPOS) = 0
C **  END OF CHAIN LOOP  **
  100   CONTINUE
  130   A(IPOS) = AVAL
  150   ICN(IPOS) = IW(JVAL,2)
C ***   END OF MAIN LOOP   ***
  200 CONTINUE
C
 1000 RETURN

    END
* COPYRIGHT (c) 1987 AEA Technology
* Original date 10 Feb 1993
C     Toolpack tool decs employed.
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.

    SUBROUTINE MC34AD(N,IRN,JCOLST,YESA,A,IW)
C THIS SUBROUTINE ACCEPTS AS INPUT THE STANDARD DATA STRUCTURE FOR
C   A SYMMETRIC MATRIX STORED AS A LOWER TRIANGLE AND PRODUCES
C   AS OUTPUT THE SYMMETRIC MATRIX HELD IN THE SAME DATA
C   STRUCTURE AS A GENERAL MATRIX.
C N IS AN INTEGER VARIABLE THAT MUST BE SET BY THE USER TO THE
C   ORDER OF THE MATRIX. NOT ALTERED BY THE ROUTINE
C   RESTRICTION (IBM VERSION ONLY): N LE 32767.
C IRN IS AN INTEGER (INTEGER*2 IN IBM VERSION) ARRAY THAT
C   MUST BE SET BY THE USER TO HOLD THE ROW INDICES OF THE LOWER
C   TRIANGULAR PART OF THE SYMMETRIC MATRIX.  THE ENTRIES OF A
C   SINGLE COLUMN ARE CONTIGUOUS. THE ENTRIES OF COLUMN J
C   PRECEDE THOSE OF COLUMN J+1 (J_=_1, ..., N-1), AND THERE IS
C   NO WASTED SPACE BETWEEN COLUMNS. ROW INDICES WITHIN A COLUMN
C   MAY BE IN ANY ORDER.  ON EXIT IT WILL HAVE THE SAME MEANING
C   BUT WILL BE CHANGED TO HOLD THE ROW INDICES OF ENTRIES IN
C   THE EXPANDED STRUCTURE.  DIAGONAL ENTRIES NEED NOT BE
C   PRESENT. THE NEW ROW INDICES ADDED IN THE UPPER TRIANGULAR
C   PART WILL BE IN ORDER FOR EACH COLUMN AND WILL PRECEDE THE
C   ROW INDICES FOR THE LOWER TRIANGULAR PART WHICH WILL REMAIN
C   IN THE INPUT ORDER.
C JCOLST IS AN INTEGER ARRAY OF LENGTH N+1 THAT MUST BE SET BY
C   THE USER SO THAT JCOLST(J) IS THE POSITION IN ARRAYS IRN AND
C   A OF THE FIRST ENTRY IN COLUMN J (J_=_1, ..., N).
C   JCOLST(N+1) MUST BE SET TO ONE MORE THAN THE TOTAL NUMBER OF
C   ENTRIES.  ON EXIT, JCOLST(J) WILL HAVE THE SAME MEANING BUT
C   WILL BE CHANGED TO POINT TO THE POSITION OF THE FIRST ENTRY
C   OF COLUMN J IN THE EXPANDED STRUCTURE. THE NEW VALUE OF
C   JCOLST(N+1) WILL BE ONE GREATER THAN THE NUMBER OF ENTRIES
C   IN THE EXPANDED STRUCTURE.
C YESA IS A LOGICAL VARIABLE THAT MUST BE SET TO .TRUE. IF THE
C   USER DESIRES TO GENERATE THE EXPANDED FORM FOR THE VALUES ALSO.
C   IF YESA IS .FALSE., THE ARRAY A WILL NOT BE REFERENCED.  IT IS
C   NOT ALTERED BY THE ROUTINE.
C A IS A REAL (DOUBLE PRECISION IN THE D VERSION) ARRAY THAT
C   CAN BE SET BY THE USER SO THAT A(K) HOLDS THE VALUE OF THE
C   ENTRY IN POSITION K OF IRN, {K = 1, _..._ JCOLST(N+1)-1}.
C   ON EXIT, IF YESA IS .TRUE., THE ARRAY WILL HOLD THE VALUES
C   OF THE ENTRIES IN THE EXPANDED STRUCTURE CORRESPONDING TO
C   THE OUTPUT VALUES OF IRN.   IF YESA IS .FALSE., THE ARRAY IS
C   NOT ACCESSED BY THE SUBROUTINE.
C IW IS AN INTEGER (INTEGER*2 IN IBM VERSION) ARRAY OF LENGTH
C   N THAT WILL BE USED AS WORKSPACE.
C
C CKP1 IS A LOCAL VARIABLE USED AS A RUNNING POINTER.
C OLDTAU IS NUMBER OF ENTRIES IN SYMMETRIC STORAGE.
C   .. Scalar Arguments ..
    INTEGER N
    LOGICAL YESA
C   ..
C   .. Array Arguments ..
    DOUBLE PRECISION A(*)
    INTEGER IRN(*),IW(*),JCOLST(*)
C   ..
C   .. Local Scalars ..
    INTEGER CKP1,I,I1,I2,II,IPKP1,IPOS,J,JSTART,LENK,NDIAG,NEWTAU,
   +      OLDTAU
C   ..
C   .. Executable Statements ..
C
    OLDTAU = JCOLST(N+1) - 1
C INITIALIZE WORK ARRAY
    DO 5 I = 1,N
      IW(I) = 0
    5 CONTINUE
C
C IW(J) IS SET EQUAL TO THE TOTAL NUMBER OF ENTRIES IN COLUMN J
C   OF THE EXPANDED SYMMETRIC MATRIX.
C NDIAG COUNTS NUMBER OF DIAGONAL ENTRIES PRESENT
    NDIAG = 0
    DO 20 J = 1,N
      I1 = JCOLST(J)
      I2 = JCOLST(J+1) - 1
      IW(J) = IW(J) + I2 - I1 + 1
      DO 10 II = I1,I2
      I = IRN(II)
      IF (I.NE.J) THEN
        IW(I) = IW(I) + 1

      ELSE
        NDIAG = NDIAG + 1
      END IF

   10   CONTINUE
   20 CONTINUE
C
C NEWTAU IS NUMBER OF ENTRIES IN EXPANDED STORAGE.
    NEWTAU = 2*OLDTAU - NDIAG
C IPKP1 POINTS TO POSITION AFTER END OF COLUMN BEING CURRENTLY
C   PROCESSED
    IPKP1 = OLDTAU + 1
C CKP1 POINTS TO POSITION AFTER END OF SAME COLUMN IN EXPANDED
C   STRUCTURE
    CKP1 = NEWTAU + 1
C GO THROUGH THE ARRAY IN THE REVERSE ORDER PLACING LOWER TRIANGULAR
C   ELEMENTS IN THE APPROPRIATE SLOTS.
    DO 40 J = N,1,-1
      I1 = JCOLST(J)
      I2 = IPKP1
C LENK IS NUMBER OF ENTRIES IN COLUMN J OF ORIGINAL STRUCTURE
      LENK = I2 - I1
C JSTART IS RUNNING POINTER TO POSITION IN NEW STRUCTURE
      JSTART = CKP1
C SET IKP1 FOR NEXT COLUMN
      IPKP1 = I1
      I2 = I2 - 1
C RUN THROUGH COLUMNS IN REVERSE ORDER
C LOWER TRIANGULAR PART OF COLUMN MOVED TO END OF SAME COLUMN IN
C   EXPANDED FORM
      DO 30 II = I2,I1,-1
      JSTART = JSTART - 1
      IF (YESA) A(JSTART) = A(II)
      IRN(JSTART) = IRN(II)
   30   CONTINUE
C JCOLST IS SET TO POSITION OF FIRST ENTRY IN LOWER TRIANGULAR PART OF
C   COLUMN J IN EXPANDED FORM
      JCOLST(J) = JSTART
C SET CKP1 FOR NEXT COLUMN
      CKP1 = CKP1 - IW(J)
C RESET IW(J) TO NUMBER OF ENTRIES IN LOWER TRIANGLE OF COLUMN.
      IW(J) = LENK
   40 CONTINUE
C
C AGAIN SWEEP THROUGH THE COLUMNS IN THE REVERSE ORDER, THIS
C   TIME WHEN ONE IS HANDLING COLUMN J THE UPPER TRIANGULAR
C   ELEMENTS A(J,I) ARE PUT IN POSITION.
    DO 80 J = N,1,-1
      I1 = JCOLST(J)
      I2 = JCOLST(J) + IW(J) - 1
C RUN DOWN COLUMN IN ORDER
C NOTE THAT I IS ALWAYS GREATER THAN OR EQUAL TO J
      DO 60 II = I1,I2
      I = IRN(II)
      IF (I.EQ.J) GO TO 60
      JCOLST(I) = JCOLST(I) - 1
      IPOS = JCOLST(I)
      IF (YESA) A(IPOS) = A(II)
      IRN(IPOS) = J
   60   CONTINUE
   80 CONTINUE
    JCOLST(N+1) = NEWTAU + 1
    RETURN

    END
* COPYRIGHT (c) 1988 AEA Technology
* Original date 17 Feb 2005

C 17th February 2005 Version 1.0.0. Replacement for FD05.

    DOUBLE PRECISION FUNCTION FD15AD(T)
C----------------------------------------------------------------
C  Fortran 77 implementation of the Fortran 90 intrinsic
C    functions: EPSILON, TINY, HUGE and RADIX.  Note that
C    the RADIX result is returned as DOUBLE PRECISION.
C
C  The CHARACTER argument specifies the type of result:
C     
C   'E'  smallest positive real number: 1.0 + DC(1) > 1.0, i.e.
C      EPSILON(DOUBLE PRECISION)
C   'T'  smallest full precision positive real number, i.e.
C      TINY(DOUBLE PRECISION)
C   'H'  largest finite positive real number, i.e.
C      HUGE(DOUBLE PRECISION)
C   'R'  the base of the floating point arithematic, i.e.
C      RADIX(DOUBLE PRECISION)
C
C    any other value gives a result of zero.
C----------------------------------------------------------------
    CHARACTER T

    IF ( T.EQ.'E' ) THEN
       FD15AD = EPSILON(1.0D0)
    ELSE IF ( T.EQ.'T' ) THEN
       FD15AD = TINY(1.0D0)
    ELSE IF ( T.EQ.'H' ) THEN
       FD15AD = HUGE(1.0D0)
    ELSE IF ( T.EQ.'R' ) THEN
       FD15AD = DBLE(RADIX(1.0D0))
    ELSE
       FD15AD = 0.0D0
    ENDIF
    RETURN
    END
C COPYRIGHT (c) 1995 Timothy A. Davis, Patrick Amestoy and
C       Council for the Central Laboratory of the Research Councils
C Original date 30 November 1995
C  April 2001: call to MC49 changed to MC59 to make routine threadsafe
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.
C 23 May 2007 Version 1.1.0. Absolute value of hash taken to cover the
C        case of integer overflow.
C        Comments with character in column 2 corrected.
C 2 August 2007 Version 2.0.0 Dense row handling added, error & warning
C        messages added, iovflo added, interface changed. MC47I/ID
C        added.
C 31 October 2007 Version 2.1.0 Corrected tree formation when handling
C        full variables
C

    SUBROUTINE MC47ID(ICNTL)
    INTEGER ICNTL(10)
C ICNTL is an INTEGER array of length 10 that contains control
C   parameters and must be set by the user. Default values are set
C   by MA57ID.
C
C   ICNTL(1) is the stream number for error messages. Printing is
C   suppressed if ICNTL(1)<0. The default is 6.
C
C   ICNTL(2) is the stream number for warning messages. Printing is
C   suppressed if ICNTL(2)<0. The default is 6.
C
C   ICNTL(3) is the stream number for printing matrix data on entry
C   to and exit from MC47A/AD. Printing is suppressed if ICNTL(3)<0.
C   The default is -1.
C
C   ICNTL(4) controls the choice of AMD algorithm
C   =-1 Classical MC47B (AMD) algorithm (no dense row detection)
C   = 0 Only exactly dense rows in the reduced matrix are selected.
C   = 1 Corresponds to automatic setting of the minimum density
C   requirement.
C   The default value is 1.
C
C   ICNTL(5) defines the largest positive
C   integer that your computer can represent (-iovflo should also
C   be representable). HUGE(1) in Fortran 95. The default value is
C   2139062143
C
C   ICNTL(6) to ICNTL(10) are set to zero by MA57ID but are not
C   currently used by MC47.
C
C Local variables
    INTEGER I

    ICNTL(1) = 6
    ICNTL(2) = 6
    ICNTL(3) = -1
    ICNTL(4) = 1
    ICNTL(5) = 2139062143

    DO 100 I=6,10
      ICNTL(I) = 0
 100  CONTINUE
    RETURN
    END


    SUBROUTINE MC47AD(N, NE, PE, IW, IWLEN,
   *    ICNTL,INFO, RINFO)
    INTEGER N, NE, PE(N+1), IWLEN, IW(IWLEN), INFO(10)
    INTEGER ICNTL(10)
    DOUBLE PRECISION RINFO(10)
C N is an INTEGER variable that must be set by the user to the
C   order of the matrix A.  It is not altered by the subroutine.
C   Restriction: N >= 1.
C NE is an INTEGER variable that must be set by the user to the
C   number of entries in the matrix A.
C   It is not altered by the subroutine.
C PE is an INTEGER  array of size N+1 that must be set by the user.
C   If the user is supplying the entries by columns, then PE(i) must
C   hold the index in IW of the start of column i, i=1, ..., N and
C   PE(N+1) must be equal to NE+1.  If the user is supplying row and
C   column indices for the matrix, PE(1) must be negative.
C   On exit, PE will hold information on the matrix factors.
C IW is an INTEGER  array of length IWLEN that must be set by the user
C   to hold the pattern of the matrix A.  If PE(1) is positive,
C   IW(PE(J)), ..., IW(PE(J+1)-1) must hold the row indices of entries
C   in column J, J = 1, ..., N.
C   The entries within a column need not be in order.
C   If PE(1) is negative, then (IW(k), IW(NE+k)), k = 1, ..., NE, must
C   hold the row and column index of an entry.
C   Duplicates, out-of-range entries, diagonal entries, and entries
C   in upper triangle are ignored.
C   IW is used as workspace by the subroutine.  On exit, the
C   permutation generated by MC47 is held in IW(IWLEN-I+1), I=1, ... N
C   and is such that the kth column of the permuted matrix is column
C   IW(IWLEN-N+k) of the original matrix.
C   The inverse permutation is held in positions IWLEN-2N+1 to
C   IWLEN-N of IW, preceded by further information on the structure
C   of the factors
C IWLEN is an INTEGER variable. It must be set by the user to the length
C   of array IW and is not altered by the subroutine.
C   We recommend IWLEN > 2NE + 9N.   Restriction: IWLEN >= 2NE + 8N.
C ICNTL is an INTEGER array of length 10 that contains control
C   parameters and must be set by the user. Default values for the
C   components may be set by a call to MC47ID.
C INFO is an INTEGER array of length 8 that need not be set by the user.
C   On return from MC47A/AD, a value of zero for INFO(1) indicates
C   that the subroutine has performed successfully.  Negative values
C   for INFO(1) signify a fatal error.  Possible values are:
C   -1    N < 1
C   -2    IWLEN < 2NE + 8N.
C   -3    Error in PE when input is by columns.
C   -4    Matrix is null (usually because all entries in upper
C       triangle.
C   There is one warning indicated by a positive value for INFO(1)
C   +1    Out-of-range index, duplicate, diagonal or entry in upper
C       triangle in input. Action taken is to ignore these entries.
C   The other entries of INFO give information to the user.
C   INFO(2) gives the number of compresses performed on the array IW.
C         A large value for this indicates that the ordering could
C         be found more quickly if IWLEN were increased.
C   INFO(3) gives the minimum necessary value for IWLEN for a
C         successful run of MC47A/AD on the same matrix as has just
C         been analysed, without the need for any compresses of IW.
C   INFO(4) gives the number of entries with row or column indices
C         that are out of range.  Any such entry is ignored.
C   INFO(5) gives the number of duplicate entries.  Any such entry is
C         ignored.
C   INFO(6) gives the number of entries in upper triangle.  Any such
C         entry is ignored.
C   INFO(7) gives the number of diagonal entries.  Any such entry is
C         ignored.
C   INFO(8) gives the number of restarts performed.
C RINFO is an INTEGER array of length 10 that need not be set by the
C   user. The other entries of RINFO give information to the user.
C   RINFO(1) gives forecast number of reals to hold the factorization
C   RINFO(2) gives the forecast number of flops required by the
C   factorization if no pivoting is performed.

C Local variables
    INTEGER DEGREE
    DOUBLE PRECISION DUMMY(1)
    INTEGER ELEN,HEAD,I,II,I1,I2,J,LAST,LEN,LENIW,LP,MP,
   *      NEXT,NV,PFREE,W,WP
    INTEGER ICT59(10),INFO59(10),IOUT,JOUT,IDUP,JNFO(10)
C DEGREE is used to subdivide array IW
C DUMMY  is dummy real for call to MC34A/AD
C ELEN   is used to subdivide array IW
C HEAD   is used to subdivide array IW
C I    is DO loop variable
C IFLAG  is eror return from MC59A/AD
C II   is running index for entries of IW in column
C I1   is first location for entries of IW in column
C I2   is flast location for entries of IW in column
C J    is column index and DO loop variable
C LAST   is used to subdivide array IW
C LEN    is used to subdivide array IW
C LENIW  is space left in IW after allocation of work vectors
C LP   is a local copy of ICNTL(1)
C MP   is a local copy of ICNTL(3)
C NEXT   is used to subdivide array IW
C NV   is used to subdivide array IW
C PFREE  marks active length of IW for call to MC47B/BD
C W    is used to subdivide array IW
C WP   is a local copy of ICNTL(2)
C IOUT   number of indices out-of-range
C JOUT   number of column indices out-of-range (as detected by MC59A/AD)
C IDUP   number of duplicates


C Subroutines called
    EXTERNAL MC59AD,MC34AD,MC47BD

C Initialize info
    DO 5 J = 1,10
       INFO(J) = 0
 5    CONTINUE

C Set LP,WP,MP
    LP = ICNTL(1)
    WP = ICNTL(2)
    MP = ICNTL(3)

C Check value of N
    IF (N.LT.1) THEN
      INFO(1) = -1
      IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10)')
   +     '**** Error return from MC47AD **** INFO(1) =',INFO(1),
   +     'N has value ',N
      GO TO 1000
    ENDIF

C Error check to see if enough space in IW to get started
    IF (PE(1).LT.1) THEN
      IF (2*NE+N.GT.IWLEN) THEN
      INFO(1) = -2
       IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10)')
   +      '**** Error return from MC47AD **** INFO(1) =',INFO(1),
   +      'IWLEN has value ',IWLEN
      GO TO 1000
      ENDIF
    ELSE
      IF (NE+N.GT.IWLEN) THEN
      INFO(1) = -2
       IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10)')
   +      '**** Error return from MC47AD **** INFO(1) =',INFO(1),
   +      'IWLEN has value ',IWLEN
      GO TO 1000
      ENDIF
    ENDIF

C Diagnostic print
    IF (MP.GE.0) THEN
      WRITE(MP,'(/A)') 'Entry to MC47A/AD'
      WRITE(MP,'(A,I10,A,I10,A)') 'Matrix of order',N,' with',NE,
   *                  ' entries'
      IF (PE(1).LT.0)  THEN
      WRITE(MP,'(A)') 'Matrix input in coordinate form'
      WRITE(MP,'(A/(4(I8,I8)))') 'Row and column indices',
   *      (IW(I),IW(NE+I),I=1,NE)
      ELSE
      WRITE(MP,'(A)') 'Matrix input by columns'
      DO 10 J=1,N
        WRITE(MP,'(A,I4/(10I8))') 'Column',J,
   *                    (IW(I),I=PE(J),PE(J+1)-1)
   10   CONTINUE
      ENDIF
    ENDIF

C Divide workspace
    LAST   = IWLEN  - N + 1
    ELEN   = LAST   - N
    NV   = ELEN   - N
    W    = NV   - N
    DEGREE = W    - N
    HEAD   = DEGREE - N
    NEXT   = HEAD   - N
    LEN    = NEXT   - N
    LENIW = LEN-1

C Set counters for number of upper triangular entries and diagonals
C   present in input matrix.
C These will be removed for later processing.
    INFO(6) = 0
    INFO(7) = 0

    IF (PE(1).LT.0) THEN
C First remove diagonals (if present) and all entries in upper triangle
C Note that this may give out-of-range entries from MC59.
      DO 20 I=1,NE
      IF (IW(I).LE.IW(NE+I)) THEN
        IF (IW(I).EQ.IW(NE+I) .AND. IW(I).NE.0) THEN
          INFO(7) = INFO(7) + 1
        ELSE
          IF (IW(I).GT.0) INFO(6) = INFO(6) + 1
        ENDIF
        IW(I)=0
      ENDIF
   20   CONTINUE

C Call sort routine
      ICT59(1) = 0
      ICT59(2) = 1
      ICT59(3) = 1
      ICT59(4) = LP
      ICT59(5) = -1
      ICT59(6) = 0
      CALL MC59AD(ICT59,N,N,NE,IW,NE,IW(NE+1),1,DUMMY,
   *          N+1,PE,N+1,IW(2*NE+1),INFO59)
C      IFLAG = INFO59(1)
      IDUP  = INFO59(3)
      IOUT  = INFO59(4)
      JOUT  = INFO59(5)
    ELSE

C Matrix already sorted by columns.
C Remove duplicates, out-of-range indices, and entries in upper
C     triangle.  First initialize counts.
      IDUP = 0
      IOUT = 0
      JOUT = 0

C Set array used to find duplicates
      DO 30 I = 1,N
      IW(NE+I) = 0
   30   CONTINUE

      DO 50 J=1,N
      I1 = PE(J)
      PE(J) = I1-(IOUT+IDUP)
      I2 = PE(J+1)-1
      IF (I2.LT.I1-1) THEN
        INFO(1) = -3
        GO TO 1000
      ENDIF
      DO 40 II = I1,I2
        I = IW(II)
        IF (I.LE.J .OR. I.GT.N) THEN
          IF (I.EQ.J) INFO(7) = INFO(7) + 1
          IF (I.GT.0 .AND. I.LT.J) INFO(6) = INFO(6) + 1
          IOUT = IOUT + 1
        ELSE
          IF (IW(NE+I).EQ.J) THEN
C Duplicate found
          IDUP = IDUP + 1
          ELSE
          IW(NE+I)=J
          IW(II-(IOUT+IDUP)) = I
          ENDIF
        ENDIF
   40   CONTINUE
   50   CONTINUE
      PE(N+1) = NE - (IOUT+IDUP) + 1
    ENDIF

C Set flags for duplicates or out-of-range entries
C Check if there were duplicates
    IF (IDUP.GT.0) THEN
      INFO(1) = 1
      INFO(4) = IDUP
      IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
   +     '**** Warning from MC47AD **** INFO(1) =',INFO(1),
   +     'Number of duplicates found: ',INFO(4)
    ELSE
      INFO(4) = 0
    ENDIF
C Check for out of range entries
    IF (IOUT+ JOUT - INFO(7) .GT.0 ) THEN
      INFO(1) = 1
      INFO(5) = IOUT + JOUT - INFO(7)
      IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
   +     '**** Warning from MC47AD **** INFO(1) =',INFO(1),
   +     'Number of out of range entries found and ignored: ',
   +     INFO(5)
    ELSE
      INFO(5) = 0
    ENDIF

C Check for entries in upper triangle
    IF (INFO(6).GT.0) THEN
       INFO(1) = 1
      IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
   +     '**** Warning from MC47AD **** INFO(1) =',INFO(1),
   +     'Number of entries in upper triangle found and ignored: ',
   +      INFO(6)
    ENDIF

C Check for entries in diagonals
    IF (INFO(7).GT.0) THEN
       INFO(1) = 1
      IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
   +     '**** Warning from MC47AD **** INFO(1) =',INFO(1),
   +     'Number of entries in diagonals found and ignored: ',
   +      INFO(7)
    ENDIF

C Check for null matrix
C Usually happens if wrong triangle is input
    IF (NE-(IOUT+IDUP).EQ.0) THEN
      INFO(1) = -4
      IF (LP.GE.0) WRITE(LP,'(/A,I3/A)')
   +     '**** Error return from MC47AD **** INFO(1) =',INFO(1),
   +     'Matrix is null'
      GO TO 1000
    ENDIF

C Generate matrix in expanded form
C First check there is sufficient space in IW
    IF (LENIW.LT.2*(PE(N+1)-1)) THEN
      INFO(1) = -2
       IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10/A,I10)')
   +      '**** Error return from MC47AD **** INFO(1) =',INFO(1),
   +      'IWLEN has value ',IWLEN,
   +      'Should be at least', 2*(PE(N+1)-1)+8*N
      GO TO 1000
    ENDIF

    CALL MC34AD(N,IW,PE,.FALSE.,DUMMY,IW(W))
    PFREE = PE(N+1)

C Set length array for MC47B/BD
    DO 60 I=1,N
      IW(LEN+I-1) = PE(I+1) - PE(I)
   60 CONTINUE

C Call to approximate minimum degree subroutine.
    CALL MC47BD(N,LENIW,PE,PFREE,IW(LEN),IW,IW(NV),
   *        IW(ELEN),IW(LAST),IW(DEGREE),
   *        IW(HEAD),IW(NEXT),IW(W), ICNTL,JNFO, RINFO)

    INFO(2) = JNFO(1)
    INFO(3) = PFREE+8*N
    INFO(8) = JNFO(2)

C Print diagnostics
    IF (MP.GE.0) THEN
      WRITE(MP,'(/A)') 'Exit from MC47A/AD'
      WRITE(MP,'(A/(7I10))') 'INFO(1-10):',(INFO(I),I=1,10)
      WRITE(MP,'(A/(8I10))') 'Parent array',(PE(I),I=1,N)
      WRITE(MP,'(A/(8I10))') 'Permutation',(IW(ELEN+I-1),I=1,N)
      WRITE(MP,'(A/(8I10))') 'Inverse permutation',
   *               (IW(LAST+I-1),I=1,N)
      WRITE(MP,'(A/(8I10))') 'Degree array',(IW(NV+I-1),I=1,N)
    ENDIF

 1000 RETURN
    END

C ====================================================================
C ====================================================================
C ====================================================================

    SUBROUTINE MC47BD (N, IWLEN, PE, PFREE, LEN, IW, NV,
   $             ELEN, LAST, DEGREE,
   $             HEAD, DENXT, W, ICNTL, JNFO, RJNFO)

    INTEGER N, IWLEN, PE(N), PFREE, LEN(N), IW(IWLEN), NV(N),
   $      ELEN(N), LAST(N),  DEGREE(N),
   $       HEAD(N), DENXT(N), W(N), ICNTL(10), JNFO(10)

    DOUBLE PRECISION RJNFO(10)

C -------------------------------------------------------------------
C AMDD is a modified version of
C   MC47B:  Approximate Minimum (UMFPACK/MA38-style, external) Degree
C      ordering algorithm, with aggresive absorption
C designed to automatically detect and exploit dense
C rows in the reduced matrix at any step of the minimum degree.
C
C We use the term Le to denote the set of all supervariables in element
C E.
C **** Reword below***
C A row is declared as full if none of its entries can be guaranteed
C   to be zero.
C A row is quasi dense if at most N-THRESM-1 of its entries can be
C   guaranteed to be zero at it has not been recognized as being full
C   in the calculation so far.
C A row is dense if it is either full or quasi dense.
C A row is sparse if it is not dense.
C -------------------------------------------------------------------
C
C
C N must be set to the matrix order. It is not altered.
C   Restriction:  N .ge. 1
C
C IWLEN must be set to the length of IW. It is not altered. On input,
C   the matrix is stored in IW (1..PFREE-1).
C   *** We do not recommend running this algorithm with ***
C   ***    IWLEN .LT. PFREE + N.              ***
C   *** Better performance will be obtained if      ***
C   ***    IWLEN .GE. PFREE + N               ***
C   *** or better yet                     ***
C   ***    IWLEN .GT. 1.2 * PFREE             ***
C   Restriction: IWLEN .GE. PFREE-1
C
C PE(i) must be set to the the index in IW of the start of row I, or be
C   zero if row I has no off-diagonal entries. During execution,
C   it is used for both supervariables and elements:
C     * Principal supervariable I:  index into IW of the
C         list of supervariable I.  A supervariable
C         represents one or more rows of the matrix
C         with identical pattern.
C     * Non-principal supervariable I:  if I has been absorbed
C         into another supervariable J, then PE(I) = -J.
C         That is, J has the same pattern as I.
C         Note that J might later be absorbed into another
C         supervariable J2, in which case PE(I) is still -J,
C         and PE(J) = -J2.
C     * Unabsorbed element E:  the index into IW of the list
C         of element E.  Element E is created when
C         the supervariable of the same name is selected as
C         the pivot.
C     * Absorbed element E:  if element E is absorbed into element
C         E2, then PE(E) = -E2.  This occurs when one of its
C         variables is eliminated and when the pattern of
C         E (that is, Le) is found to be a subset of the pattern
C         of E2 (that is, Le2).  If element E is "null" (it has
C         no entries outside its pivot block), then PE(E) = 0.
C
C   On output, PE holds the assembly tree/forest, which implicitly
C   represents a pivot order with identical fill-in as the actual
C   order (via a depth-first search of the tree). If NV(I) .GT. 0,
C   then I represents a node in the assembly tree, and the parent of
C   I is -PE(I), or zero if I is a root. If NV(I)=0, then (I,-PE(I))
C   represents an edge in a subtree, the root of which is a node in
C   the assembly tree.
C
C PFREE must be set to the position in IW of the first free variable.
C   During execution, additional data is placed in IW, and PFREE is
C   modified so that components  of IW from PFREE are free.
C   On output, PFREE is set equal to the size of IW that would have
C   caused no compressions to occur.  If NCMPA is zero, then
C   PFREE (on output) is less than or equal to IWLEN, and the space
C   IW(PFREE+1 ... IWLEN) was not used. Otherwise, PFREE (on output)
C   is greater than IWLEN, and all the memory in IW was used.
C
C LEN(I) must be set to hold the number of entries in row I of the
C   matrix, excluding the diagonal.  The contents of LEN(1..N) are
C   undefined on output.
C
C IW(1..PFREE-1) must be set to  hold the patterns of the rows of
C   the matrix.  The matrix must be symmetric, and both upper and
C   lower triangular parts must be present.  The diagonal must not be
C   present.  Row I is held as follows:
C         IW(PE(I)...PE(I) + LEN(I) - 1) must hold the list of
C         column indices for entries in row I (simple
C         supervariables), excluding the diagonal.  All
C         supervariables start with one row/column each
C         (supervariable I is just row I). If LEN(I) is zero on
C         input, then PE(I) is ignored on input. Note that the
C         rows need not be in any particular order, and there may
C         be empty space between the rows.
C   During execution, the supervariable I experiences fill-in. This
C   is represented by constructing a list of the elements that cause
C   fill-in in supervariable I:
C         IE(PE(i)...PE(I) + ELEN(I) - 1) is the list of elements
C         that contain I. This list is kept short by removing
C         absorbed elements. IW(PE(I)+ELEN(I)...PE(I)+LEN(I)-1)
C         is the list of supervariables in I. This list is kept
C         short by removing nonprincipal variables, and any entry
C         J that is also contained in at least one of the
C         elements in the list for I.
C   When supervariable I is selected as pivot, we create an element E
C   of the same name (E=I):
C         IE(PE(E)..PE(E)+LEN(E)-1) is the list of supervariables
C          in element E.
C   An element represents the fill-in that occurs when supervariable
C   I is selected as pivot.
C   CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
C   The contents of IW are undefined on output.
C
C NV(I) need not be set. During execution, ABS(NV(I)) is equal to the
C   number of rows represented by the principal supervariable I. If I
C   is a nonprincipal variable, then NV(I) = 0. Initially, NV(I) = 1
C   for all I.  NV(I) .LT. 0 signifies that I is a principal variable
C   in the pattern Lme of the current pivot element ME. On output,
C   NV(E) holds the true degree of element E at the time it was
C   created (including the diagonal part).
C
C ELEN(I) need not be set. See the description of IW above. At the
C   start of execution, ELEN(I) is set to zero. For a supervariable,
C   ELEN(I) is the number of elements in the list for supervariable
C   I. For an element, ELEN(E) is the negation of the position in the
C   pivot sequence of the supervariable that generated it. ELEN(I)=0
C   if I is nonprincipal.
C   On output ELEN(1..N) holds the inverse permutation (the same
C   as the 'INVP' argument in Sparspak). That is, if K = ELEN(I),
C   then row I is the Kth pivot row.  Row I of A appears as the
C   (ELEN(I))-th row in the permuted matrix, PAP^T.
C
C LAST(I) need not be set on input. In a degree list, LAST(I) is the
C   supervariable preceding I, or zero if I is the head of the list.
C   In a hash bucket, LAST(I) is the hash key for I. LAST(HEAD(HASH))
C   is also used as the head of a hash bucket if HEAD(HASH) contains
C   a degree list (see HEAD, below).
C   On output, LAST(1..N) holds the permutation (the same as the
C   'PERM' argument in Sparspak). That is, if I = LAST(K), then row I
C   is the Kth pivot row.  Row LAST(K) of A is the K-th row in the
C   permuted matrix, PAP^T.
C
C
C DEGREE need not be set on input. If I is a supervariable and sparse,
C   then DEGREE(I) holds the current approximation of the external
C   degree of row I (an upper bound). The external degree is the
C   number of entries in row I, minus ABS(NV(I)) (the diagonal
C   part). The bound is equal to the external degree if ELEN(I) is
C   less than or equal to two. We also use the term "external degree"
C   for elements E to refer to |Le \ Lme|. If I is full in the reduced
C   matrix, then DEGREE(I)=N+1. If I is dense in the reduced matrix,
C   then DEGREE(I)=N+1+last_approximate_external_deg of I.
C   All dense rows are stored in the list pointed by HEAD(N).
C   Quasi dense rows are stored first, and are followed by full rows
C   in the reduced matrix. LASTD holds the last row in
C   this list of dense rows or is zero if the list is empty.
C
C HEAD(DEG) need not be set on input. HEAD is used for degree lists.
C   HEAD(DEG) is the first supervariable in a degree list (all
C   supervariables I in a degree list DEG have the same approximate
C   degree, namely, DEG = DEGREE(I)). If the list DEG is empty then
C   HEAD(DEG) = 0.
C   During supervariable detection HEAD(HASH) also serves as a
C   pointer to a hash bucket.
C   If HEAD(HASH) .GT. 0, there is a degree list of degree HASH. The
C   hash bucket head pointer is LAST(HEAD(HASH)).
C   If HEAD(HASH) = 0, then the degree list and hash bucket are
C   both empty.
C   If HEAD(HASH) .LT. 0, then the degree list is empty, and
C   -HEAD(HASH) is the head of the hash bucket.
C   After supervariable detection is complete, all hash buckets are
C   empty, and the (LAST(HEAD(HASH)) = 0) condition is restored for
C   the non-empty degree lists.
C
C DENXT(I) need not be set on input. For supervariable I, DENXT(I) is
C   the supervariable following I in a link list, or zero if I is
C   the last in the list. Used for two kinds of lists: degree lists
C   and hash buckets (a supervariable can be in only one kind of
C   list at a time). For element E, DENXT(E) is the number of
C   variables with dense or full rows in the element E.
C
C W(I) need not be set on input. The flag array W determines the status
C   of elements and variables, and the external degree of elements.
C   For elements:
C      if W(E) = 0, then the element E is absorbed.
C      if W(E) .GE. WFLG, then W(E)-WFLG is the size of the set
C         |Le \ Lme|, in terms of nonzeros (the sum of ABS(NV(I))
C         for each principal variable I that is both in the
C         pattern of element E and NOT in the pattern of the
C         current pivot element, ME).
C      if WFLG .GT. WE(E) .GT. 0, then E is not absorbed and has
C         not yet been seen in the scan of the element lists in
C         the computation of |Le\Lme| in loop 150 below.
C         ***SD: change comment to remove reference to label***
C   For variables:
C      during supervariable detection, if W(J) .NE. WFLG then J is
C      not in the pattern of variable I.
C   The W array is initialized by setting W(I) = 1 for all I, and by
C   setting WFLG = 2. It is reinitialized if WFLG becomes too large
C   (to ensure that WFLG+N does not cause integer overflow).
C
C ICNTL is an INTEGER array of length 10 that contains control
C   parameters and must be set by the user. Default values for the
C   components may be set by a call to MC47ID.
C
C RJNFO is an REAL (DOUBLE PRECISION in  D version) array of length 7
C   that need not be set by the user. This array supplies information
C   on the execution of MC47BD.
C   RJNFO(1) gives forecast number of reals to hold the factorization
C   RJNFO(2) gives the forecast number of flops required by the
C   factorization if no pivoting is performed.
C
C Local variables:
C ---------------
    INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, HASH, HMOD, I,
   $   IDUMMY, ILAST, INEXT, IOVFLO,J, JDUMMY, JLAST, JNEXT, K,
   $   KNT1, KNT2, KNT3, LASTD,  LENJ, LN, MAXMEM, ME,
   $   MEM, MINDEG, NBD, NCMPA, NDME, NEL, NELME, NEWMEM,
   $   NFULL, NLEFT, NRLADU, NVI, NVJ, NVPIV, P, P1, P2, P3, PDST,
   $   PEE, PEE1, PEND, PJ, PME, PME1, PME2, PN, PSRC, RSTRT,
   $   SLENME, THRESH, THRESM, WE, WFLG, WNVI,X
   $
    DOUBLE PRECISION RELDEN, SM, STD, OPS
    LOGICAL IDENSE
C
C DEG:      the degree of a variable or element
C DEGME:    size (no. of variables), |Lme|, of the current element,
C         ME (= DEGREE(ME))
C DEXT:     external degree, |Le \ Lme|, of some element E
C DMAX:     largest |Le| seen so far
C E:      an element
C ELENME:   the length, ELEN(ME), of element list of pivotal var.
C ELN:      the length, ELEN(...), of an element list
C HASH:     the computed value of the hash function
C HMOD:     the hash function is computed modulo HMOD = MAX(1,N-1)
C I:      a supervariable
C IDUMMY:   loop counter
C ILAST:    the entry in a link list preceding I
C INEXT:    the entry in a link list following I
C IOVFLO:   local copy of ICNTL(5)
C J:      a supervariable
C JDUMMY:   loop counter
C JLAST:    the entry in a link list preceding J
C JNEXT:    the entry in a link list, or path, following J
C K:      the pivot order of an element or variable
C KNT1:     loop counter used during element construction
C KNT2:     loop counter used during element construction
C KNT3:     loop counter used during element construction
C LASTD:    index of the last row in the list of dense rows
C LENJ:     LEN(J)
C LN:       length of a supervariable list
C MAXMEM:   amount of memory needed for no compressions
C ME:       current supervariable being eliminated, and the
C             current element created by eliminating that
C             supervariable
C MEM:      memory in use assuming no compressions have occurred
C MINDEG:   current approximate minimum degree
C NBD:      total number of dense rows selected
C NCMPA:    counter for the number of times IW was compressed
C NDME  :   number of dense rows adjacent to me
C NEL:      number of pivots selected so far
C NELME:    number of pivots selected when reaching the root
C NEWMEM:   amount of new memory needed for current pivot element
C NFULL:    total number of full rows detected.
C NLEFT:    N-NEL, the number of nonpivotal rows/columns remaining
C NRLADU:   counter for the forecast number of reals in matrix factor
C NVI:      the number of variables in a supervariable I (= NV(I))
C NVJ:      the number of variables in a supervariable J (= NV(J))
C NVPIV:    number of pivots in current element
C P:      pointer into lots of things
C P1:       pe (i) for some variable i (start of element list)
C P2:       pe (i) + elen (i) -  1 for some var. i (end of el. list)
C P3:       index of first supervariable in clean list
C PJ:       pointer into an element or variable
C PDST:     destination pointer, for compression
C PEE:      pointer into element E
C PEE1:     pointer into element E
C PEND:     end of memory to compress
C PME:      pointer into the current element (PME1...PME2)
C PME1:     the current element, ME, is stored in IW(PME1...PME2)
C PME2:     the end of the current element
C PN:       pointer into a "clean" variable, also used to compress
C PSRC:     source pointer, for compression
C RSTRT:    counter for the number of restarts carried out
C SLENME:   number of variables in variable list of pivotal variable
C THRESH:   local copy of ICNTL(4)
C THRESM :    local integer holding the threshold used to detect quasi
C         dense rows. When quasi dense rows are reintegrated in the
C         graph to be processed then THRESM is modified.
C WE:       W(E)
C WFLG:     used for flagging the W array.  See description of W.
C WNVI:     WFLG-NV(I)
C X:      either a supervariable or an element
C
C OPS:      counter for forecast number of flops
C RELDEN :    holds average density to set THRESM automatically
C SM:       counter used for forming standard deviation
C STD:      standard deviation
C
C IDENSE is true if supervariable I is dense
C
C -------------------------------------------------------------------
C  FUNCTIONS CALLED:
C -------------------------------------------------------------------
    INTRINSIC MAX, MIN, MOD
C ====================================================================
C  INITIALIZATIONS
C ====================================================================

    DO 2 I = 1,10
       RJNFO(I) = 0.0
       JNFO(I) = 0
 2    CONTINUE
    DMAX = 0
    HMOD = MAX (1, N-1)
    IOVFLO = ICNTL(5)
    LASTD = 0
    MEM = PFREE - 1
    MAXMEM = MEM
    MINDEG = 1
    NBD   = 0
    NCMPA = 0
    NEL = 0
    NFULL  = 0
    NRLADU = 0
    RSTRT = 0
    OPS = 0.00
    THRESH = ICNTL(4)
    WFLG = 2

C   ------------------------------------------------------
C   Experiments with automatic setting of parameter THRESH.
C   ------------------------------------------------------
    IF (THRESH.GT.0) THEN
       THRESM  = 0
       RELDEN = 0.0
       SM = 0
C      ----------------------------------------------------------
C      initialize arrays and eliminate rows with no off-diag. nz.
C      ----------------------------------------------------------
       DO 5 I=1,N
        THRESM = MAX(THRESM, LEN(I))
        IF (LEN(I).GT.0) THEN
         RELDEN = RELDEN + LEN(I)
         SM = SM + (LEN(I) * LEN(I))
        END IF
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        DEGREE (I) = LEN (I)
        IF (DEGREE(I) .EQ. 0) THEN
         NEL = NEL + 1
         ELEN (I) = -NEL
         PE (I) = 0
         W (I) = 0
         NRLADU = NRLADU + 1
         OPS = OPS + 1
        ELSE
         W (I) = 1
         ELEN (I) = 0
        ENDIF
 5     CONTINUE
       IF (N .EQ. NEL) GOTO 265

       RELDEN = RELDEN/(N-NEL)
C      RELDEN holds average row length
       SM = SM/(N-NEL-NFULL) - RELDEN*RELDEN
       STD = SQRT(ABS(SM))
C      STD holds standard deviation of the row lengths
       IF (STD .LE. RELDEN) THEN
        THRESM = -1
       ELSE
        THRESM = INT(9*RELDEN + 0.5*STD*((STD/(RELDEN+0.01))**1.5)+
   *       2*RELDEN*RELDEN/(STD+0.01) +1)
       END IF
C   ------------------------------------------------------
C   end automatic setting of THRESM
C   ------------------------------------------------------

    ELSE
       THRESM = THRESH
       DO 10 I = 1, N
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        DEGREE (I) = LEN (I)
        IF (DEGREE(I) .EQ. 0) THEN
         NEL = NEL + 1
         ELEN (I) = -NEL
         PE (I) = 0
         W (I) = 0
         NRLADU = NRLADU + 1
         OPS = OPS + 1
        ELSE
         W (I) = 1
         ELEN (I) = 0
        ENDIF
 10    CONTINUE
    ENDIF
    IF (THRESM.GE.0) THEN
       IF (THRESM.GE.N) THEN
C       full rows only
        THRESM = -1
       ELSE IF (THRESM.EQ.0) THEN
        THRESM = N
       ENDIF
    ENDIF

C   ----------------------------------------------------------------
C   initialize degree lists
C   ----------------------------------------------------------------
    DO 20 I = 1, N
       DEG = DEGREE (I)
       IF (DEG .GT. 0) THEN
C       ----------------------------------------------------------
C       place i in the degree list corresponding to its degree
C       or in the dense row list if i is dense
C       ----------------------------------------------------------
C       test for row density
        IF ( (THRESM.GE.0) .AND.
   &       (DEG+1.GE.THRESM.OR.DEG+1.GE.N-NEL )) THEN
C          I is dense and will be inserted in the degree
C          list of N
         NBD = NBD+1
         IF (DEG+1.NE.N-NEL) THEN
C           I is quasi dense
            DEGREE(I) = DEGREE(I)+N+1
C           insert I at the beginning of degree list of n
            DEG = N
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            DENXT (I) = INEXT
            HEAD (DEG) = I
            LAST(I)  = 0
            IF (LASTD.EQ.0) THEN
             LASTD=I
            END IF
         ELSE
C           I is full
            NFULL = NFULL+1
            DEGREE(I) = N+1
C           insert I at the end of degree list of n
            DEG = N
            IF (LASTD.EQ.0) THEN
C            degree list is empty
             LASTD   = I
             HEAD(DEG) = I
             DENXT(I)   = 0
             LAST(I)   = 0
            ELSE
C             IF (NFULL.EQ.1) THEN
C               First full row encountered
                DENXT(LASTD) = I
                LAST(I)   = LASTD
                LASTD     = I
                DENXT(I)   = 0
C             ELSE
C               Absorb I into LASTD (first full row found)
C                PE(I) = - LASTD
C                NV(LASTD) = NV(LASTD) + NV(I)
C                NV(I) = 0
C                ELEN(I) = 0
C             END IF
            ENDIF
         ENDIF
        ELSE
C          place i in the degree list corresponding to its degree
         INEXT = HEAD (DEG)
         IF (INEXT .NE. 0) LAST (INEXT) = I
         DENXT (I) = INEXT
         HEAD (DEG) = I
        ENDIF
       ENDIF
 20   CONTINUE

C   We suppress dense row selection if none of them was found in A
C   in the 1st pass
    IF (NBD.EQ.0 .AND. THRESH.GT.0) THEN
       THRESM = -1
    END IF
C
C ====================================================================
C  WHILE (selecting pivots) DO
C ====================================================================

 30   IF (NEL .LT. N) THEN

C ==================================================================
C  GET PIVOT OF MINIMUM APPROXIMATE DEGREE
C ==================================================================
C     -------------------------------------------------------------
C     find next supervariable for elimination
C     -------------------------------------------------------------
       DO 40 DEG = MINDEG, N
        ME = HEAD (DEG)
        IF (ME .GT. 0) GO TO 50
 40    CONTINUE
 50    MINDEG = DEG
       IF (DEG.LT.N)  THEN
C     -------------------------------------------------------------
C     remove chosen variable from linked list
C     -------------------------------------------------------------
        INEXT = DENXT (ME)
        IF (INEXT .NE. 0) LAST (INEXT) = 0
        HEAD (DEG) = INEXT
       ELSE
        IF (DEGREE(ME).EQ.N+1) GO TO 263
C DEGREE(ME).GT.N+1 so ME is quasi dense
C RESTARTING STRATEGY
C       FOR EACH  quasi dense row d
C       1/ insert d in the degree list according to the
C        value degree(d)-(N+1) (updating MINDEG)
C       2/ Build the adjacency list of d in the quotient graph
C          update DENXT(e_me)= DENXT(e_me)-NV(ME)
C       4/ get back to min degree process
C
C       THRESM > 0 because quasi dense rows were selected
C       While loop: ME is the current dense row
C       make sure that WFLG is not too large
        RSTRT = RSTRT + 1
        RELDEN = 0.0
        SM = 0
        IF (WFLG .GT. IOVFLO-NBD-1) THEN
         DO  51 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
 51        CONTINUE
         WFLG = 2
        END IF
        WFLG = WFLG + 1
        DO 57 IDUMMY = 1,N

C       ---------------------------------------------------------
C       remove chosen variable from link list
C       ---------------------------------------------------------
         INEXT = DENXT (ME)
         IF (INEXT .NE. 0) THEN
            LAST (INEXT) = 0
         ELSE
            LASTD = 0
         ENDIF
C       ----------------------------------------------------------
c       build adjacency list of ME in quotient graph
C       and calculate its external degree in ndense(me)
C       ----------------------------------------------------------
         DENXT(ME) = 0
C          Flag ME as having been considered in this calculation
         W(ME)    = WFLG
         P1 = PE(ME)
         P2 = P1 + LEN(ME) -1
C       LN-1 holds the pointer in IW to last elt/var in adj list
C          of ME.  LEN(ME) will then be set to LN-P1
C       ELN-1 hold the pointer in IW to  last elt in in adj list
C          of ME.  ELEN(ME) will then be set to ELN-P1
C       element adjacent to ME
         LN     = P1
         ELN    = P1
         DO 55 P=P1,P2
            E= IW(P)
            IF (W(E).EQ.WFLG) GO TO 55
            W(E) = WFLG
C         -------------------------------------------
C         Ensure that E is an unabsorbed element or a quasi dense
C           row and flag it
C         -------------------------------------------
            DO 52 JDUMMY = 1,N
             IF ( PE(E) .GE. 0 ) GOTO 53
             E = -PE(E)
             IF (W(E) .EQ.WFLG) GOTO 55
             W(E) = WFLG
 52         CONTINUE
 53         IF (ELEN(E).LT.0) THEN
C            E is a new element in adj(ME)
             DENXT(E) = DENXT(E) - NV(ME)
C            Move first entry in ME's list of adjacent variables
C            to the end
             IW(LN) = IW(ELN)
C            Place E at end of ME's list of adjacent elements
             IW(ELN) = E
             LN  = LN+1
             ELN = ELN + 1
C            update ndense of ME with all unflagged dense
C            rows in E
             PEE1 = PE(E)
             DO 54 PEE = PEE1, PEE1+LEN(E)-1
                X = IW(PEE)
                IF ((ELEN(X).GE.0).AND.(W(X).NE.WFLG)) THEN
C                X is a dense row
                 DENXT(ME) = DENXT(ME) + NV(X)
                 W(X) = WFLG
                ENDIF
 54            CONTINUE
            ELSE
C            E is a dense row
             DENXT(ME) = DENXT(ME) + NV(E)
C            Place E at end of ME's list of adjacent variables
             IW(LN)=E
             LN = LN+1
            ENDIF
 55        CONTINUE

C       ----------------------------------------------
C       DEGREE(ME)-(N+1) holds last external degree computed
C       when ME was detected as dense
C       DENXT(ME) is the exact external degree of ME
C       ----------------------------------------------
         WFLG   = WFLG + 1
         LEN(ME)  = LN-P1
         ELEN(ME) = ELN- P1
         NDME = DENXT(ME)+NV(ME)
C        If we want to select ME as full (NDME.EQ.NBD)
C        or quasi dense (NDME.GE.THRESM) then
C        denxt(of elements adjacent to ME) should be updated
         IF (DENXT(ME).EQ.0) DENXT(ME) =1
C       ---------------------------------------------------------
C       place ME in the degree list of DENXT(ME), update DEGREE
C       ---------------------------------------------------------
C         IF (DEGREE(ME)+NV(ME) .LT. NBD  ) THEN
         IF (NDME .LT. NBD) THEN
C           ME is not full
            RELDEN = RELDEN + NV(ME)*NDME
            SM = SM + NV(ME)*NDME*NDME
            DEGREE(ME) = DENXT(ME)
            DEG = DEGREE(ME)
            MINDEG = MIN(DEG,MINDEG)
            JNEXT = HEAD(DEG)
            IF (JNEXT.NE. 0) LAST (JNEXT) = ME
            DENXT(ME) = JNEXT
            HEAD(DEG) = ME
         ELSE
C           ME is full
            DEGREE(ME) = N+1
            DEG = DENXT(ME)
            MINDEG = MIN(DEG,MINDEG)
            DEG = N

C           Update DENXT of all elements in the list of elements
C           adjacent to ME
            P1 = PE(ME)
            P2 = P1 + ELEN(ME) - 1
            DO 56 PJ=P1,P2
             E= IW(PJ)
             DENXT (E) = DENXT(E) + NV(ME)
 56         CONTINUE
C            insert ME in the list of dense rows
            DEG = N
C           ME at the end of the list
            NFULL = NFULL +NV(ME)
            IF (LASTD.EQ.0) THEN
C                degree list is empty
             LASTD   = ME
             HEAD(N) = ME
             DENXT(ME)   = 0
             LAST(ME)   = 0
             IF (INEXT.EQ.0) INEXT = LASTD
            ELSE
C             IF (NFULL.EQ.NV(ME)) THEN
C                 First full row encountered
                DENXT(LASTD) = ME
                LAST(ME)   = LASTD
                LASTD      = ME
                DENXT(ME)   = 0
                IF (INEXT.EQ.0) INEXT = LASTD
C             ELSE
C             Absorb ME into LASTD (first full row found)
C                PE(ME) = - LASTD
C                NV(LASTD) = NV(LASTD) + NV(ME)
C                NV(ME) = 0
C                ELEN(ME) = 0
C             END IF
            ENDIF
         END IF

C       ------------------------------
C       process next quasi dense row
C       ------------------------------
         ME    = INEXT
         IF (ME.EQ.0) GO TO 58
         IF (DEGREE(ME).LE.(N+1) ) GOTO 58
 57       CONTINUE
 58       HEAD (N) = ME
C       ---------------------------------------
C       update dense row selection strategy
C       -------------------------------------
        IF (NBD.EQ.NFULL) THEN
         RELDEN = 0
         SM = 0
        ELSE
         RELDEN = (RELDEN + NFULL*NBD)/(NBD)
         SM = (SM + NFULL*NBD*NBD)/(NBD) - RELDEN*RELDEN
        END IF
        STD = SQRT(ABS(SM))
        THRESM = INT(9*RELDEN+0.5*STD*((STD/(RELDEN + 0.01))**1.5)
   *       + 2*RELDEN*RELDEN/(STD+0.01) +1)
        THRESM = MIN(THRESM,NBD)
        IF (THRESM.GE.NBD) THEN
         THRESM = N
        END IF
        NBD = NFULL
C       get back to min degree elimination loop
        GOTO 30
C       -------------------------------------------------------------
C       -------------------------------------------------------------
       ENDIF
C     -------------------------------------------------------------
C     me represents the elimination of pivots nel+1 to nel+nv(me).
C     place me itself as the first in this set.  It will be moved
C     to the nel+nv(me) position when the permutation vectors are
C     computed.
C     -------------------------------------------------------------
       ELENME = ELEN (ME)
       ELEN (ME) = - (NEL + 1)
       NVPIV = NV (ME)
       NEL = NEL + NVPIV
       DENXT(ME) = 0

C ====================================================================
C  CONSTRUCT NEW ELEMENT
C ====================================================================
C
C     -------------------------------------------------------------
C     At this point, me is the pivotal supervariable.  It will be
C     converted into the current element.  Scan list of the
C     pivotal supervariable, me, setting tree pointers and
C     constructing new list of supervariables for the new element,
C     me.  p is a pointer to the current position in the old list.
C     -------------------------------------------------------------
C
C     flag the variable "me" as being in the front by negating nv(me)
       NV (ME) = -NVPIV
       DEGME = 0
       IF (ELENME .EQ. 0) THEN
C       ----------------------------------------------------------
C       There are no elements involved.
C       Construct the new element in place.
C       ----------------------------------------------------------
        PME1 = PE (ME)
        PME2 = PME1 - 1
        DO 60 P = PME1, PME1 + LEN (ME) - 1
         I = IW (P)
         NVI = NV (I)
         IF (NVI .GT. 0) THEN
C         ----------------------------------------------------
C         i is a principal variable not yet placed in the
C         generated element. Store i in new list
C         ----------------------------------------------------
            DEGME = DEGME + NVI
C           flag i as being in Lme by negating nv (i)
            NV (I) = -NVI
            PME2 = PME2 + 1
            IW (PME2) = I

C         ----------------------------------------------------
C         remove variable i from degree list.
C         ----------------------------------------------------
C         only done for sparse rows
            IF (DEGREE(I).LE.N) THEN
             ILAST = LAST (I)
             INEXT = DENXT (I)
             IF (INEXT .NE. 0) LAST (INEXT) = ILAST
             IF (ILAST .NE. 0) THEN
                DENXT (ILAST) = INEXT
             ELSE
C               i is at the head of the degree list
                HEAD (DEGREE (I)) = INEXT
             ENDIF
            ELSE
C            Dense rows remain dense so do not remove from list
             DENXT(ME) = DENXT(ME) + NVI
            ENDIF
         ENDIF
 60       CONTINUE
C       this element takes no new memory in iw:
        NEWMEM = 0
       ELSE
C       ----------------------------------------------------------
C       construct the new element in empty space, iw (pfree ...)
C       ----------------------------------------------------------
        P  = PE (ME)
        PME1 = PFREE
        SLENME = LEN (ME) - ELENME
        DO 120 KNT1 = 1, ELENME
C          search the elements in me.
         E = IW (P)
         P = P + 1
         PJ = PE (E)
         LN = LEN (E)
C       -------------------------------------------------------
C       search for different supervariables and add them to the
C       new list, compressing when necessary.
C       -------------------------------------------------------
         DO 110 KNT2 = 1, LN
            I = IW (PJ)
            PJ = PJ + 1
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
C         -------------------------------------------------
C         compress iw, if necessary
C         -------------------------------------------------
             IF (PFREE .GT. IWLEN) THEN
C            prepare for compressing iw by adjusting
C            pointers and lengths so that the lists being
C            searched in the inner and outer loops contain
C            only the remaining entries.
C            ***** SD: Seperate compression subroutine tried
C              but found to be inefficient in comparison ****
                PE (ME) = P
                LEN (ME) = LEN (ME) - KNT1
C               Check if anything left in supervariable ME
                IF (LEN (ME) .EQ. 0) PE (ME) = 0
                PE (E) = PJ
                LEN (E) = LN - KNT2
C               Check if anything left in element E
                IF (LEN (E) .EQ. 0) PE (E) = 0
                NCMPA = NCMPA + 1
C               store first item in pe
C               set first entry to -item
                DO 70 J = 1, N
                 PN = PE (J)
                 IF (PN .GT. 0) THEN
                  PE (J) = IW (PN)
                  IW (PN) = -J
                 ENDIF
 70             CONTINUE

C               psrc/pdst point to source/destination
                PDST = 1
                PSRC = 1
                PEND = PME1 - 1

C               while loop:
                DO 91 IDUMMY = 1, IWLEN
                 IF (PSRC .GT. PEND) GO TO 95
C                search for next negative entry
                 J = -IW (PSRC)
                 PSRC = PSRC + 1
                 IF (J .GT. 0) THEN
                  IW (PDST) = PE (J)
                  PE (J) = PDST
                  PDST = PDST + 1
C             copy from source to destination
                  LENJ = LEN (J)
                  DO 90 KNT3 = 0, LENJ - 2
                     IW (PDST + KNT3) = IW (PSRC + KNT3)
 90                 CONTINUE
                  PDST = PDST + LENJ - 1
                  PSRC = PSRC + LENJ - 1
                 ENDIF
 91             END DO

C               move the new partially-constructed element
 95             P1 = PDST
                DO 100 PSRC = PME1, PFREE - 1
                 IW (PDST) = IW (PSRC)
                 PDST = PDST + 1
 100            CONTINUE
                PME1 = P1
                PFREE = PDST
                PJ = PE (E)
                P = PE (ME)
             ENDIF

C         -------------------------------------------------
C         i is a principal variable not yet placed in Lme
C         store i in new list
C         -------------------------------------------------
             DEGME = DEGME + NVI
C            flag i as being in Lme by negating nv (i)
             NV (I) = -NVI
             IW (PFREE) = I
             PFREE = PFREE + 1

C         -------------------------------------------------
C         remove variable i from degree link list
C         -------------------------------------------------
C            only done for sparse rows
             IF (DEGREE(I).LE.N) THEN
                ILAST = LAST (I)
                INEXT = DENXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                 DENXT (ILAST) = INEXT
                ELSE
C               i is at the head of the degree list
                 HEAD (DEGREE (I)) = INEXT
                ENDIF
             ELSE
C            Dense rows remain dense so do not remove from list
                DENXT(ME) = DENXT(ME) + NVI
             ENDIF
            ENDIF
 110       CONTINUE

C          set tree pointer and flag to indicate element e is
C          absorbed into new element me (the parent of e is me)
            PE (E) = -ME
            W (E) = 0
 120      CONTINUE

C       search the supervariables in me.
        KNT1 = ELENME + 1
        E = ME
        PJ = P
        LN = SLENME

C       -------------------------------------------------------
C       search for different supervariables and add them to the
C       new list, compressing when necessary.
C       -------------------------------------------------------
        DO 126 KNT2 = 1, LN
         I = IW (PJ)
         PJ = PJ + 1
         NVI = NV (I)
         IF (NVI .GT. 0) THEN
C         -------------------------------------------------
C         compress iw, if necessary
C         -------------------------------------------------
            IF (PFREE .GT. IWLEN) THEN
C           prepare for compressing iw by adjusting
C           pointers and lengths so that the lists being
C           searched in the inner and outer loops contain
C           only the remaining entries.
             PE (ME) = P
             LEN (ME) = LEN (ME) - KNT1
C            Check if anything left in supervariable ME
             IF (LEN (ME) .EQ. 0) PE (ME) = 0
             PE (E) = PJ
             LEN (E) = LN - KNT2
C            Check if anything left in element E
             IF (LEN (E) .EQ. 0) PE (E) = 0
             NCMPA = NCMPA + 1
C            store first item in pe
C            set first entry to -item
             DO 121 J = 1, N
                PN = PE (J)
                IF (PN .GT. 0) THEN
                 PE (J) = IW (PN)
                 IW (PN) = -J
                ENDIF
 121           CONTINUE

C            psrc/pdst point to source/destination
             PDST = 1
             PSRC = 1
             PEND = PME1 - 1

C           while loop:
C 122          CONTINUE
             DO 123 IDUMMY = 1,IWLEN
                IF (PSRC .GT. PEND) GO TO 124
C               search for next negative entry
                J = -IW (PSRC)
                PSRC = PSRC + 1
                IF (J .GT. 0) THEN
                 IW (PDST) = PE (J)
                 PE (J) = PDST
                 PDST = PDST + 1
C                copy from source to destination
                 LENJ = LEN (J)
                 DO 122 KNT3 = 0, LENJ - 2
                  IW (PDST + KNT3) = IW (PSRC + KNT3)
 122               CONTINUE
                 PDST = PDST + LENJ - 1
                 PSRC = PSRC + LENJ - 1
                ENDIF
 123           END DO

C           move the new partially-constructed element
 124           P1 = PDST
             DO 125 PSRC = PME1, PFREE - 1
                IW (PDST) = IW (PSRC)
                PDST = PDST + 1
 125           CONTINUE
             PME1 = P1
             PFREE = PDST
             PJ = PE (E)
             P = PE (ME)
            END IF

C         -------------------------------------------------
C         i is a principal variable not yet placed in Lme
C         store i in new list
C         -------------------------------------------------
            DEGME = DEGME + NVI
C           flag i as being in Lme by negating nv (i)
            NV (I) = -NVI
            IW (PFREE) = I
            PFREE = PFREE + 1

C         -------------------------------------------------
C         remove variable i from degree link list
C         -------------------------------------------------
C           only done for sparse rows
            IF (DEGREE(I).LE.N) THEN
             ILAST = LAST (I)
             INEXT = DENXT (I)
             IF (INEXT .NE. 0) LAST (INEXT) = ILAST
             IF (ILAST .NE. 0) THEN
                DENXT (ILAST) = INEXT
             ELSE
C           i is at the head of the degree list
                HEAD (DEGREE (I)) = INEXT
             ENDIF
            ELSE
C            Dense rows remain dense so do not remove from list
             DENXT(ME) = DENXT(ME) + NVI
            ENDIF
         ENDIF
 126       CONTINUE

        PME2 = PFREE - 1
C       this element takes newmem new memory in iw (possibly zero)
        NEWMEM = PFREE - PME1
        MEM = MEM + NEWMEM
        MAXMEM = MAX (MAXMEM, MEM)
       ENDIF

C     -------------------------------------------------------------
C     me has now been converted into an element in iw (pme1..pme2)
C     -------------------------------------------------------------
C     degme holds the external degree of new element
       DEGREE (ME) = DEGME
       PE (ME) = PME1
       LEN (ME) = PME2 - PME1 + 1

C     -------------------------------------------------------------
C     make sure that wflg is not too large.  With the current
C     value of wflg, wflg+n must not cause integer overflow
C     -------------------------------------------------------------
       IF (WFLG .GT. IOVFLO-N) THEN
        DO 130 X = 1, N
         IF (W (X) .NE. 0) W (X) = 1
 130      CONTINUE
        WFLG = 2
       ENDIF

C ====================================================================
C   COMPUTE (w(e) - wflg) = |Le(G')\Lme(G')| FOR ALL ELEMENTS
C   where G' is the subgraph of G containing just the sparse rows)
C ====================================================================
C     -------------------------------------------------------------
C     Scan 1:  compute the external degrees of elements touched
C     with respect to the current element.  That is:
C        (w (e) - wflg) = |Le \ Lme|
C     for each element e involving a supervariable in Lme.
C     The notation Le refers to the pattern (list of
C     supervariables) of a previous element e, where e is not yet
C     absorbed, stored in iw (pe (e) + 1 ... pe (e) + iw (pe (e))).
C     The notation Lme refers to the pattern of the current element
C     (stored in iw (pme1..pme2)).
C     aggressive absorption is possible only if DENXT(ME) = NBD
C     which is true when only full rows have been selected.
C     -------------------------------------------------------------
       IF (NBD.GT.0) THEN
C       Dense rows have been found
        DO 150 PME = PME1, PME2
         I = IW (PME)
c          skip dense rows
         IF (DEGREE(I).GT.N) GOTO 150
         ELN = ELEN (I)
         IF (ELN .GT. 0) THEN
C          note that nv (i) has been negated to denote i in Lme:
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
             E = IW (P)
             WE = W (E)
             IF (WE .GE. WFLG) THEN
C            unabsorbed element e has been seen in this loop
                WE = WE - NVI
             ELSE IF (WE .NE. 0) THEN
C            e is an unabsorbed element - this is
C            the first we have seen e in all of Scan 1
                WE = DEGREE (E) + WNVI - DENXT(E)
             ENDIF
             W (E) = WE
 140          CONTINUE
         ENDIF
 150      CONTINUE
       ELSE
C       No dense rows have been found
        DO 152 PME = PME1, PME2
         I = IW (PME)
         ELN = ELEN (I)
         IF (ELN .GT. 0) THEN
C       note that nv (i) has been negated to denote i in Lme:
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 151 P = PE (I), PE (I) + ELN - 1
             E = IW (P)
             WE = W (E)
             IF (WE .GE. WFLG) THEN
C            unabsorbed element e has been seen in this loop
                WE = WE - NVI
             ELSE IF (WE .NE. 0) THEN
C            e is an unabsorbed element - this is
C            the first we have seen e in all of Scan 1
                WE = DEGREE (E) + WNVI
             ENDIF
             W (E) = WE
 151          CONTINUE
         ENDIF
 152      CONTINUE
       END IF

C ====================================================================
C  DEGREE UPDATE AND ELEMENT ABSORPTION
C ====================================================================

C     -------------------------------------------------------------
C     Scan 2:  for each sparse i in Lme, sum up the external degrees
C     of each Le for the elements e appearing within i, plus the
C     supervariables in i.  Place i in hash list.
C     -------------------------------------------------------------
       IF (NBD.GT.0) THEN
C       Dense rows have been found
        DO 180 PME = PME1, PME2
         I = IW (PME)
C          skip dense rows
         IF (DEGREE(I).GT.N) GOTO 180
C          remove absorbed elements from the list for i
         P1 = PE (I)
         P2 = P1 + ELEN (I) - 1
         PN = P1
         HASH = 0
         DEG = 0

C       ----------------------------------------------------------
C       scan the element list associated with supervariable i
C       ----------------------------------------------------------
         DO 160 P = P1, P2
            E = IW (P)
C           dext = | Le | - | (Le \cap Lme)\D | - DENXT(e)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
             DEG = DEG + DEXT
             IW (PN) = E
             PN = PN + 1
             HASH = HASH+E
            ELSE IF ((DEXT .EQ. 0) .AND.
   &            (DENXT(ME).EQ.NBD)) THEN
C         aggressive absorption: e is not adjacent to me, but
C         |Le(G') \ Lme(G')| is 0 and all dense rows
C         are in me, so absorb it into me
             PE (E) = -ME
             W (E)  = 0
            ELSE IF (DEXT.EQ.0) THEN
             IW(PN) = E
             PN   = PN+1
             HASH = HASH + E
            ENDIF
 160       CONTINUE

C       count the number of elements in i (including me):
         ELEN (I) = PN - P1 + 1

C       ----------------------------------------------------------
C       scan the supervariables in the list associated with i
C       ----------------------------------------------------------
         P3 = PN
         DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
C           j is unabsorbed, and not in Lme.
C           add to degree and add to new list
C           add degree only of sparse rows.
             IF (DEGREE(J).LE.N) DEG=DEG+NVJ
             IW (PN) = J
             PN = PN + 1
             HASH = HASH + J
            ENDIF
 170       CONTINUE

C       ----------------------------------------------------------
C       update the degree and check for mass elimination
C       ----------------------------------------------------------
         IF ((DEG .EQ. 0).AND.(DENXT(ME).EQ.NBD)) THEN
C          mass elimination only possible when all dense rows
C          are in ME
C       -------------------------------------------------------
C       mass elimination - supervariable i can be eliminated
C       -------------------------------------------------------
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
         ELSE
C       -------------------------------------------------------
C       update the upper-bound degree of i
C       A bound for the new external degree is the old bound plus
C       the size of the generated element
C       -------------------------------------------------------
C       the following degree does not yet include the size
C       of the current element, which is added later:
            DEGREE(I) = MIN (DEG+NBD-DENXT(ME), DEGREE(I))

C       -------------------------------------------------------
C       add me to the list for i
C       -------------------------------------------------------
C          move first supervariable to end of list
            IW (PN) = IW (P3)
C          move first element to end of element part of list
            IW (P3) = IW (P1)
C          add new element to front of list.
            IW (P1) = ME
C          store the new length of the list in len (i)
            LEN (I) = PN - P1 + 1

C       -------------------------------------------------------
C       place in hash bucket.  Save hash key of i in last (i).
C       -------------------------------------------------------
            HASH = ABS(MOD (HASH, HMOD)) + 1
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
C          the degree list is empty, hash head is -j
             DENXT (I) = -J
             HEAD (HASH) = -I
            ELSE
C          degree list is not empty - has j as its head
C          last is hash head
             DENXT (I) = LAST (J)
             LAST (J) = I
            ENDIF
            LAST (I) = HASH
         ENDIF
 180      CONTINUE
       ELSE
C       No dense rows have been found
        DO 183 PME = PME1, PME2
         I = IW (PME)
C          remove absorbed elements from the list for i
         P1 = PE (I)
         P2 = P1 + ELEN (I) - 1
         PN = P1
         HASH = 0
         DEG = 0

C          -------------------------------------------------------
C          scan the element list associated with supervariable i
C          -------------------------------------------------------
         DO 181 P = P1, P2
            E = IW (P)
C           dext = | Le | - | (Le \cap Lme)\D | - DENXT(e)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
             DEG = DEG + DEXT
             IW (PN) = E
             PN = PN + 1
             HASH = HASH + E
            ELSE IF (DEXT .EQ. 0) THEN
C            aggressive absorption: e is not adjacent to me, but
C            |Le(G') \ Lme(G')| is 0, so absorb it into me
             PE (E) = -ME
             W (E)  = 0
            ENDIF
 181       CONTINUE

C       count the number of elements in i (including me):
         ELEN (I) = PN - P1 + 1

C       ----------------------------------------------------------
C       scan the supervariables in the list associated with i
C       ----------------------------------------------------------
         P3 = PN
         DO 182 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
C           j is unabsorbed, and not in Lme.
C           add to degree and add to new list
             DEG=DEG+NVJ
             IW (PN) = J
             PN = PN + 1
             HASH = HASH + J
            ENDIF
 182       CONTINUE

C       ----------------------------------------------------------
C       update the degree and check for mass elimination
C       ----------------------------------------------------------
         IF (DEG .EQ. 0) THEN
C       -------------------------------------------------------
C       mass elimination - supervariable i can be eliminated
C       -------------------------------------------------------
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
         ELSE
C       -------------------------------------------------------
C       update the upper-bound degree of i
C       A bound for the new external degree is the old bound plus
C       the size of the generated element
C       -------------------------------------------------------
C
C       the following degree does not yet include the size
C       of the current element, which is added later:
            DEGREE(I) = MIN (DEG,  DEGREE(I))

C       -------------------------------------------------------
C       add me to the list for i
C       -------------------------------------------------------
C          move first supervariable to end of list
            IW (PN) = IW (P3)
C          move first element to end of element part of list
            IW (P3) = IW (P1)
C          add new element to front of list.
            IW (P1) = ME
C          store the new length of the list in len (i)
            LEN (I) = PN - P1 + 1

C       -------------------------------------------------------
C       place in hash bucket.  Save hash key of i in last (i).
C       -------------------------------------------------------
            HASH = ABS(MOD (HASH, HMOD)) + 1
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
C          the degree list is empty, hash head is -j
             DENXT (I) = -J
             HEAD (HASH) = -I
            ELSE
C          degree list is not empty - has j as its head
C          last is hash head
             DENXT (I) = LAST (J)
             LAST (J) = I
            ENDIF
            LAST (I) = HASH
         ENDIF
 183      CONTINUE
       END IF
       DEGREE (ME) = DEGME

C     -------------------------------------------------------------
C     Clear the counter array, w (...), by incrementing wflg.
C     -------------------------------------------------------------
       DMAX = MAX (DMAX, DEGME)
       WFLG = WFLG + DMAX

C      make sure that wflg+n does not cause integer overflow
       IF (WFLG .GE. IOVFLO - N) THEN
        DO 190 X = 1, N
         IF (W (X) .NE. 0) W (X) = 1
 190      CONTINUE
        WFLG = 2
       ENDIF
C      at this point, w (1..n) .lt. wflg holds

C ====================================================================
C  SUPERVARIABLE DETECTION
C ====================================================================
       DO 250 PME = PME1, PME2
        I = IW (PME)
        IF ( (NV(I).GE.0) .OR. (DEGREE(I).GT.N) ) GO TO 250
C       only done for sparse rows
C       replace i by head of its hash bucket, and set the hash
C       bucket header to zero

C       -------------------------------------------------------
C       examine all hash buckets with 2 or more variables.  We
C       do this by examing all unique hash keys for super-
C       variables in the pattern Lme of the current element, me
C       -------------------------------------------------------
        HASH = LAST (I)
C       let i = head of hash bucket, and empty the hash bucket
        J = HEAD (HASH)
        IF (J .EQ. 0) GO TO 250
        IF (J .LT. 0) THEN
C         degree list is empty
         I = -J
         HEAD (HASH) = 0
        ELSE
C         degree list is not empty, restore last () of head
         I = LAST (J)
         LAST (J) = 0
        ENDIF
        IF (I .EQ. 0) GO TO 250

C       while loop:
        DO 247 JDUMMY = 1,N
         IF (DENXT (I) .EQ. 0) GO TO 250
C         ----------------------------------------------------
C         this bucket has one or more variables following i.
C         scan all of them to see if i can absorb any entries
C         that follow i in hash bucket.  Scatter i into w.
C         ----------------------------------------------------
         LN = LEN (I)
         ELN = ELEN (I)
C          do not flag the first element in the list (me)
         DO 210 P = PE (I) + 1, PE (I) + LN - 1
            W (IW (P)) = WFLG
 210       CONTINUE

C         ----------------------------------------------------
C         scan every other entry j following i in bucket
C         ----------------------------------------------------
         JLAST = I
         J = DENXT (I)

C         while loop:
         DO 245 IDUMMY=1,N
            IF (J .EQ. 0) GO TO 246

C         -------------------------------------------------
C         check if j and i have identical nonzero pattern
C         -------------------------------------------------
C         jump if i and j do not have same size data structure
            IF (LEN (J) .NE. LN) GO TO 240
C         jump if i and j do not have same number adj elts
            IF (ELEN (J) .NE. ELN) GO TO 240
C         do not flag the first element in the list (me)

            DO 230 P = PE (J) + 1, PE (J) + LN - 1
C           jump if an entry (iw(p)) is in j but not in i
             IF (W (IW (P)) .NE. WFLG) GO TO 240
 230          CONTINUE

C           -------------------------------------------------
C           found it!  j can be absorbed into i
C           -------------------------------------------------
            PE (J) = -I
C           both nv (i) and nv (j) are negated since they
C           are in Lme, and the absolute values of each
C           are the number of variables in i and j:
            NV (I) = NV (I) + NV (J)
            NV (J) = 0
            ELEN (J) = 0
C           delete j from hash bucket
            J = DENXT (J)
            DENXT (JLAST) = J
            GO TO 245

C           -------------------------------------------------
 240          CONTINUE
C         j cannot be absorbed into i
C         -------------------------------------------------
            JLAST = J
            J = DENXT (J)
 245       CONTINUE

C         ----------------------------------------------------
C         no more variables can be absorbed into i
C         go to next i in bucket and clear flag array
C         ----------------------------------------------------
 246       WFLG = WFLG + 1
         I = DENXT (I)
         IF (I .EQ. 0) GO TO 250
 247       CONTINUE
 250    CONTINUE

C ====================================================================
C  RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMENT
C  Squeeze out absorbed variables
C ====================================================================
      P = PME1
      NLEFT = N - NEL
      DO 260 PME = PME1, PME2
         I = IW (PME)
         NVI = -NV (I)
         IF (NVI .LE. 0) GO TO 260
C        i is a principal variable in Lme
C        restore nv (i) to signify that i is principal
         NV (I) = NVI
         IF (DEGREE(I).GT.N) GO TO 258
C       -------------------------------------------------------
C       compute the external degree (add size of current elem)
C       -------------------------------------------------------
         DEG = MIN (DEGREE (I)+ DEGME - NVI, NLEFT - NVI)
         DEGREE (I) = DEG
         IDENSE = .FALSE.
C       -------------------
C       Dense row detection
C       -------------------
         IF (THRESM.GE.0) THEN
C       DEGME is exact external degree of pivot ME |Le\Ve|,
C       DEG is is approx external degree of I
          IF ((DEG+NVI .GE. THRESM).OR.
   &         (DEG+NVI .GE. NLEFT)) THEN
             IF (THRESM.EQ.N) THEN
C            We must be sure that I is full in reduced matrix
              IF ((ELEN(I).LE.2) .AND.((DEG+NVI).EQ.NLEFT)
   &             .AND. NBD.EQ.NFULL ) THEN
C                DEG approximation is exact and I is dense
               DEGREE(I) = N+1
               IDENSE = .TRUE.
              ENDIF
             ELSE
C             relaxed dense row detection
              IDENSE = .TRUE.
              IF ((ELEN(I).LE.2).AND. ((DEG+NVI).EQ.NLEFT)
   &             .AND. NBD.EQ.NFULL ) THEN
               DEGREE(I) = N+1
              ELSE
               DEGREE(I) = N+1+DEGREE(I)
              ENDIF
             ENDIF
          ENDIF
          IF (IDENSE) THEN
C            update DENXT of all elements in the list of element
C            adjacent to I (including ME).
             P1 = PE(I)
             P2 = P1 + ELEN(I) - 1
             DO 255 PJ=P1,P2
              E= IW(PJ)
              DENXT (E) = DENXT(E) + NVI
 255         CONTINUE
C            insert I in the list of dense rows
             NBD = NBD+NVI
             DEG = N
             IF (DEGREE(I).EQ.N+1) THEN
c             insert I at the end of the list
              NFULL = NFULL +NVI
              IF (LASTD.EQ.0) THEN
C                degree list is empty
               LASTD   = I
               HEAD(DEG) = I
               DENXT(I)   = 0
               LAST(I)   = 0
              ELSE
C               IF (NFULL.EQ.NVI) THEN
C                 First full row encountered
                  DENXT(LASTD) = I
                  LAST(I)   = LASTD
                  LASTD     = I
                  DENXT(I)   = 0
C               ELSE
C                 Absorb I into LASTD (first full row found)
C                  PE(I) = - LASTD
C                  NV(LASTD) = NV(LASTD) + NV(I)
C                  NV(I) = 0
C                  ELEN(I) = 0
C               END IF
              ENDIF
             ELSE
C             insert I at the beginning of the list
              INEXT = HEAD(DEG)
              IF (INEXT .NE. 0) LAST (INEXT) = I
              DENXT (I) = INEXT
              HEAD (DEG) = I
              LAST(I)    = 0
              IF (LASTD.EQ.0) LASTD=I
             ENDIF
C         end of IDENSE=true
          ENDIF
C        end of THRESM>0
         ENDIF

         IF (.NOT.IDENSE) THEN
C       -------------------------------------------------------
C       place the supervariable at the head of the degree list
C       -------------------------------------------------------
          INEXT = HEAD (DEG)
          IF (INEXT .NE. 0) LAST (INEXT) = I
          DENXT (I) = INEXT
          LAST (I) = 0
          HEAD (DEG) = I
         ENDIF
C       -------------------------------------------------------
C       save the new degree, and find the minimum degree
C       -------------------------------------------------------
         MINDEG = MIN (MINDEG, DEG)
 258       CONTINUE
C       -------------------------------------------------------
C       place the supervariable in the element pattern
C       -------------------------------------------------------
         IW (P) = I
         P = P + 1
 260    CONTINUE

C =====================================================================
C  FINALIZE THE NEW ELEMENT
C =====================================================================
      OPS = OPS + DEGME*NVPIV + DEGME * NVPIV*NVPIV +
   *       DEGME*DEGME*NVPIV + NVPIV*NVPIV*NVPIV/3 +
   *       NVPIV*NVPIV/2 + NVPIV/6 + NVPIV
      NRLADU = NRLADU + (NVPIV*(NVPIV+1))/2 + (DEGME*NVPIV)
      NV (ME) = NVPIV + DEGME
C       nv (me) is now the degree of pivot (including diagonal part)
C       save the length of the list for the new element me
      LEN (ME) = P - PME1
      IF (LEN (ME) .EQ. 0) THEN
C        there is nothing left of the current pivot element
         PE (ME) = 0
         W (ME) = 0
      ENDIF
      IF (NEWMEM .NE. 0) THEN
C        element was not constructed in place: deallocate part
C        of it (final size is less than or equal to newmem,
C        since newly nonprincipal variables have been removed).
         PFREE = P
         MEM = MEM - NEWMEM + LEN (ME)
      ENDIF

C =====================================================================
C     END WHILE (selecting pivots)
      GO TO 30
     ENDIF
C =====================================================================
     GO TO 265

C    We have only full rows that we amalgamate at the root
C    node and ME = LASTD
C    Perform mass elimination of all full rows
 263   NELME    = -(NEL+1)
     DO 264 X=1,N
      IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
C        X is an unabsorbed element
         PE(X) = -ME
      ELSEIF (DEGREE(X).EQ.N+1) THEN
C        X is a dense row, absorb it in ME (mass elimination)
         NEL   = NEL + NV(X)
         PE(X) = -ME
         ELEN(X) = 0
         NV(X) = 0
      ENDIF
 264   CONTINUE
C    ME is the root node
     ELEN(ME) = NELME
     NV(ME)   = NBD
     NRLADU = NRLADU + (NBD*(NBD+1))/2
     OPS = OPS + NBD*NBD*NBD/3 + NBD*NBD/2 + NBD/6 + NBD
     PE(ME)   = 0

 265   CONTINUE
C ===================================================================
C  COMPUTE THE PERMUTATION VECTORS
C ===================================================================

C   ----------------------------------------------------------------
C   The time taken by the following code is O(n).  At this
C   point, elen (e) = -k has been done for all elements e,
C   and elen (i) = 0 has been done for all nonprincipal
C   variables i.  At this point, there are no principal
C   supervariables left, and all elements are absorbed.
C   ----------------------------------------------------------------
C
C   ----------------------------------------------------------------
C   compute the ordering of unordered nonprincipal variables
C   ----------------------------------------------------------------

     DO 290 I = 1, N
      IF (ELEN (I) .EQ. 0) THEN
C       ----------------------------------------------------------
C       i is an un-ordered row.  Traverse the tree from i until
C       reaching an element, e.  The element, e, was the
C       principal supervariable of i and all nodes in the path
C       from i to when e was selected as pivot.
C       ----------------------------------------------------------
         J = -PE (I)
C       while (j is a variable) do:
         DO 270 JDUMMY = 1,N
          IF (ELEN (J) .LT. 0) GO TO 275
          J = -PE (J)
 270       CONTINUE
 275       E = J
C       ----------------------------------------------------------
C       get the current pivot ordering of e
C       ----------------------------------------------------------
         K = -ELEN (E)

C       ----------------------------------------------------------
C       traverse the path again from i to e, and compress the
C       path (all nodes point to e).  Path compression allows
C       this code to compute in O(n) time.  Order the unordered
C       nodes in the path, and place the element e at the end.
C       ----------------------------------------------------------
         J = I
C        while (j is a variable) do:
         DO 280 IDUMMY = 1,N
          IF (ELEN (J) .LT. 0) GO TO 285
          JNEXT = -PE (J)
          PE (J) = -E
          IF (ELEN (J) .EQ. 0) THEN
C            j is an unordered row
             ELEN (J) = K
             K = K + 1
          ENDIF
          J = JNEXT
 280       CONTINUE
C        leave elen (e) negative, so we know it is an element
 285       ELEN (E) = -K
      ENDIF
 290   CONTINUE

C   ----------------------------------------------------------------
C   reset the inverse permutation (elen (1..n)) to be positive,
C   and compute the permutation (last (1..n)).
C   ----------------------------------------------------------------
     DO 300 I = 1, N
      K = ABS (ELEN (I))
      LAST (K) = I
      ELEN (I) = K
 300   CONTINUE

C ====================================================================
C  RETURN THE MEMORY USAGE IN IW AND SET INFORMATION ARRAYS
C ====================================================================
C   If maxmem is less than or equal to iwlen, then no compressions
C   occurred, and iw (maxmem+1 ... iwlen) was unused.  Otherwise
C   compressions did occur, and iwlen would have had to have been
C   greater than or equal to maxmem for no compressions to occur.
C   Return the value of maxmem in the pfree argument.

     RJNFO(1) = OPS
     RJNFO(2) = NRLADU
     JNFO(1) = NCMPA
     JNFO(2) = RSTRT
     PFREE = MAXMEM

     RETURN
     END
C ====================================================================
C ====================================================================
C ====================================================================
! COPYRIGHT (c) 1999 Council for the Central Laboratory
!            of the Research Councils
! Original date July 1999
! AUTHORS Iain Duff (i.duff@rl.ac.uk) and
!       Jacko Koster (jacko.koster@uninett.no)
!
! Version 1.6.0
! See ChangeLog for version history.
!
    SUBROUTINE MC64ID(ICNTL)
    IMPLICIT NONE

C  Purpose
C  =======
C
C  The components of the array ICNTL control the action of MC64A/AD.
C  Default values for these are set in this subroutine.
C
C  Parameters
C  ==========
C
    INTEGER ICNTL(10)
C
C  Local variables
    INTEGER I
C
C    ICNTL(1) has default value 6.
C   It is the output stream for error messages. If it
C   is negative, these messages will be suppressed.
C
C    ICNTL(2) has default value 6.
C   It is the output stream for warning messages.
C   If it is negative, these messages are suppressed.
C
C    ICNTL(3) has default value -1.
C   It is the output stream for monitoring printing.
C   If it is negative, these messages are suppressed.
C
C    ICNTL(4) has default value 0.
C   If left at the defaut value, the incoming data is checked for
C   out-of-range indices and duplicates.  Setting ICNTL(4) to any
C   other will avoid the checks but is likely to cause problems
C   later if out-of-range indices or duplicates are present.
C   The user should only set ICNTL(4) non-zero, if the data is
C   known to avoid these problems.
C
C    ICNTL(5) to ICNTL(10) are not used by MC64A/AD but are set to
C   zero in this routine.

C Initialization of the ICNTL array.
    ICNTL(1) = 6
    ICNTL(2) = 6
    ICNTL(3) = -1
    DO 10 I = 4,10
      ICNTL(I) = 0
   10 CONTINUE

    RETURN
    END

C**********************************************************************
    SUBROUTINE MC64AD(JOB,N,NE,IP,IRN,A,NUM,CPERM,LIW,IW,LDW,DW,
   &       ICNTL,INFO)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
C  Purpose
C  =======
C
C This subroutine attempts to find a column permutation for an NxN
C sparse matrix A = {a_ij} that makes the permuted matrix have N
C entries on its diagonal.
C If the matrix is structurally nonsingular, the subroutine optionally
C returns a column permutation that maximizes the smallest element
C on the diagonal, maximizes the sum of the diagonal entries, or
C maximizes the product of the diagonal entries of the permuted matrix.
C For the latter option, the subroutine also finds scaling factors
C that may be used to scale the matrix so that the nonzero diagonal
C entries of the permuted matrix are one in absolute value and all the
C off-diagonal entries are less than or equal to one in absolute value.
C The natural logarithms of the scaling factors u(i), i=1..N, for the
C rows and v(j), j=1..N, for the columns are returned so that the
C scaled matrix B = {b_ij} has entries b_ij = a_ij * EXP(u_i + v_j).
C
C  Parameters
C  ==========
C
    INTEGER JOB,N,NE,NUM,LIW,LDW
    INTEGER IP(N+1),IRN(NE),CPERM(N),IW(LIW),ICNTL(10),INFO(10)
    DOUBLE PRECISION A(NE),DW(LDW)
C
C JOB is an INTEGER variable which must be set by the user to
C control the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   1 Compute a column permutation of the matrix so that the
C   permuted matrix has as many entries on its diagonal as possible.
C   The values on the diagonal are of arbitrary size. HSL subroutine
C   MC21A/AD is used for this. See [1].
C   2 Compute a column permutation of the matrix so that the smallest
C   value on the diagonal of the permuted matrix is maximized.
C   See [3].
C   3 Compute a column permutation of the matrix so that the smallest
C   value on the diagonal of the permuted matrix is maximized.
C   The algorithm differs from the one used for JOB = 2 and may
C   have quite a different performance. See [2].
C   4 Compute a column permutation of the matrix so that the sum
C   of the diagonal entries of the permuted matrix is maximized.
C   See [3].
C   5 Compute a column permutation of the matrix so that the product
C   of the diagonal entries of the permuted matrix is maximized
C   and vectors to scale the matrix so that the nonzero diagonal
C   entries of the permuted matrix are one in absolute value and
C   all the off-diagonal entries are less than or equal to one in
C   absolute value. See [3].
C  Restriction: 1 <= JOB <= 5.
C
C N is an INTEGER variable which must be set by the user to the
C   order of the matrix A. It is not altered by the subroutine.
C   Restriction: N >= 1.
C
C NE is an INTEGER variable which must be set by the user to the
C   number of entries in the matrix. It is not altered by the
C   subroutine.
C   Restriction: NE >= 1.
C
C IP is an INTEGER array of length N+1.
C   IP(J), J=1..N, must be set by the user to the position in array IRN
C   of the first row index of an entry in column J. IP(N+1) must be set
C   to NE+1. It is not altered by the subroutine.
C
C IRN is an INTEGER array of length NE.
C   IRN(K), K=1..NE, must be set by the user to hold the row indices of
C   the entries of the matrix. Those belonging to column J must be
C   stored contiguously in the positions IP(J)..IP(J+1)-1. The ordering
C   of the row indices within each column is unimportant. Repeated
C   entries are not allowed. The array IRN is not altered by the
C   subroutine.
C
C A is a DOUBLE PRECISION array of length NE.
C   The user must set A(K), K=1..NE, to the numerical value of the
C   entry that corresponds to IRN(K).
C   It is not used by the subroutine when JOB = 1.
C   It is not altered by the subroutine.
C
C NUM is an INTEGER variable that need not be set by the user.
C   On successful exit, NUM will be the number of entries on the
C   diagonal of the permuted matrix.
C   If NUM < N, the matrix is structurally singular.
C
C CPERM is an INTEGER array of length N that need not be set by the
C   user. On successful exit, CPERM contains the column permutation.
C   Column ABS(CPERM(J)) of the original matrix is column J in the
C   permuted matrix, J=1..N. For the N-NUM  entries of CPERM that are
C   not matched the permutation array is set negative so that a full
C   permutation of the matrix is obtained even in the structurally
C   singular case.
C
C LIW is an INTEGER variable that must be set by the user to
C   the dimension of array IW. It is not altered by the subroutine.
C   Restriction:
C   JOB = 1 :  LIW >= 5N
C   JOB = 2 :  LIW >= 4N
C   JOB = 3 :  LIW >= 10N + NE
C   JOB = 4 :  LIW >= 5N
C   JOB = 5 :  LIW >= 5N
C
C IW is an INTEGER array of length LIW that is used for workspace.
C
C LDW is an INTEGER variable that must be set by the user to the
C   dimension of array DW. It is not altered by the subroutine.
C   Restriction:
C   JOB = 1 :  LDW is not used
C   JOB = 2 :  LDW >= N
C   JOB = 3 :  LDW >= NE
C   JOB = 4 :  LDW >= 2N + NE
C   JOB = 5 :  LDW >= 3N + NE
C
C DW is a DOUBLE PRECISION array of length LDW
C   that is used for workspace. If JOB = 5, on return,
C   DW(i) contains u_i, i=1..N, and DW(N+j) contains v_j, j=1..N.
C
C ICNTL is an INTEGER array of length 10. Its components control the
C   output of MC64A/AD and must be set by the user before calling
C   MC64A/AD. They are not altered by the subroutine.
C
C   ICNTL(1) must be set to specify the output stream for
C   error messages. If ICNTL(1) < 0, messages are suppressed.
C   The default value set by MC46I/ID is 6.
C
C   ICNTL(2) must be set by the user to specify the output stream for
C   warning messages. If ICNTL(2) < 0, messages are suppressed.
C   The default value set by MC46I/ID is 6.
C
C   ICNTL(3) must be set by the user to specify the output stream for
C   diagnostic messages. If ICNTL(3) < 0, messages are suppressed.
C   The default value set by MC46I/ID is -1.
C
C   ICNTL(4) must be set by the user to a value other than 0 to avoid
C   checking of the input data.
C   The default value set by MC46I/ID is 0.
C
C INFO is an INTEGER array of length 10 which need not be set by the
C   user. INFO(1) is set non-negative to indicate success. A negative
C   value is returned if an error occurred, a positive value if a
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : successful entry (for structurally singular matrix).
C   +2 : the returned scaling factors are large and may cause
C      overflow when used to scale the matrix.
C      (For JOB = 5 entry only.)
C   -1 : JOB < 1 or JOB > 5.  Value of JOB held in INFO(2).
C   -2 : N < 1.  Value of N held in INFO(2).
C   -3 : NE < 1. Value of NE held in INFO(2).
C   -4 : the defined length LIW violates the restriction on LIW.
C      Value of LIW required given by INFO(2).
C   -5 : the defined length LDW violates the restriction on LDW.
C      Value of LDW required given by INFO(2).
C   -6 : entries are found whose row indices are out of range. INFO(2)
C      contains the index of a column in which such an entry is found.
C   -7 : repeated entries are found. INFO(2) contains the index of a
C      column in which such entries are found.
C  INFO(3) to INFO(10) are not currently used and are set to zero by
C      the routine.
C
C References:
C  [1]  I. S. Duff, (1981),
C     "Algorithm 575. Permutations for a zero-free diagonal",
C     ACM Trans. Math. Software 7(3), 387-390.
C  [2]  I. S. Duff and J. Koster, (1998),
C     "The design and use of algorithms for permuting large
C     entries to the diagonal of sparse matrices",
C     SIAM J. Matrix Anal. Appl., vol. 20, no. 4, pp. 889-901.
C  [3]  I. S. Duff and J. Koster, (1999),
C     "On algorithms for permuting large entries to the diagonal
C     of sparse matrices",
C     Technical Report RAL-TR-1999-030, RAL, Oxfordshire, England.

C Local variables and parameters
    INTEGER I,J,K
    DOUBLE PRECISION FACT,ZERO,RINF
    PARAMETER (ZERO=0.0D+00)
C External routines and functions
    EXTERNAL MC21AD,MC64BD,MC64RD,MC64SD,MC64WD
C Intrinsic functions
    INTRINSIC ABS,LOG

C Set RINF to largest positive real number (infinity)
    RINF = HUGE(RINF)

C Check value of JOB
    IF (JOB.LT.1 .OR. JOB.GT.5) THEN
      INFO(1) = -1
      INFO(2) = JOB
      IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
      GO TO 99
    ENDIF
C Check value of N
    IF (N.LT.1) THEN
      INFO(1) = -2
      INFO(2) = N
      IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
      GO TO 99
    ENDIF
C Check value of NE
    IF (NE.LT.1) THEN
      INFO(1) = -3
      INFO(2) = NE
      IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NE',NE
      GO TO 99
    ENDIF
C Check LIW
    IF (JOB.EQ.1) K = 5*N
    IF (JOB.EQ.2) K = 4*N
    IF (JOB.EQ.3) K = 10*N + NE
    IF (JOB.EQ.4) K = 5*N
    IF (JOB.EQ.5) K = 5*N
    IF (LIW.LT.K) THEN
      INFO(1) = -4
      INFO(2) = K
      IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
      GO TO 99
    ENDIF
C Check LDW
C If JOB = 1, do not check
    IF (JOB.GT.1) THEN
      IF (JOB.EQ.2) K = N
      IF (JOB.EQ.3) K = NE
      IF (JOB.EQ.4) K = 2*N + NE
      IF (JOB.EQ.5) K = 3*N + NE
      IF (LDW.LT.K) THEN
      INFO(1) = -5
      INFO(2) = K
      IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
      GO TO 99
      ENDIF
    ENDIF
    IF (ICNTL(4).EQ.0) THEN
C Check row indices. Use IW(1:N) as workspace
      DO 3 I = 1,N
      IW(I) = 0
    3   CONTINUE
      DO 6 J = 1,N
      DO 4 K = IP(J),IP(J+1)-1
        I = IRN(K)
C Check for row indices that are out of range
        IF (I.LT.1 .OR. I.GT.N) THEN
          INFO(1) = -6
          INFO(2) = J
          IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
          GO TO 99
        ENDIF
C Check for repeated row indices within a column
        IF (IW(I).EQ.J) THEN
          INFO(1) = -7
          INFO(2) = J
          IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I
          GO TO 99
        ELSE
          IW(I) = J
        ENDIF
    4   CONTINUE
    6   CONTINUE
    ENDIF

C Print diagnostics on input
    IF (ICNTL(3).GE.0) THEN
      WRITE(ICNTL(3),9020) JOB,N,NE
      WRITE(ICNTL(3),9021) (IP(J),J=1,N+1)
      WRITE(ICNTL(3),9022) (IRN(J),J=1,NE)
      IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(J),J=1,NE)
    ENDIF

C Set components of INFO to zero
    DO 8 I=1,10
      INFO(I) = 0
    8 CONTINUE

C Compute maximum matching with MC21A/AD
    IF (JOB.EQ.1) THEN
C Put length of column J in IW(J)
      DO 10 J = 1,N
      IW(J) = IP(J+1) - IP(J)
   10   CONTINUE
C IW(N+1:5N) is workspace
      CALL MC21AD(N,IRN,NE,IP,IW(1),CPERM,NUM,IW(N+1))
      GO TO 90
    ENDIF

C Compute bottleneck matching
    IF (JOB.EQ.2) THEN
C IW(1:5N), DW(1:N) are workspaces
      CALL MC64BD(N,NE,IP,IRN,A,CPERM,NUM,
   &   IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),DW)
      GO TO 90
    ENDIF

C Compute bottleneck matching
    IF (JOB.EQ.3) THEN
C Copy IRN(K) into IW(K), ABS(A(K)) into DW(K), K=1..NE
      DO 20 K = 1,NE
      IW(K) = IRN(K)
      DW(K) = ABS(A(K))
   20   CONTINUE
C Sort entries in each column by decreasing value.
      CALL MC64RD(N,NE,IP,IW,DW)
C IW(NE+1:NE+10N) is workspace
      CALL MC64SD(N,NE,IP,IW(1),DW,CPERM,NUM,IW(NE+1),
   &   IW(NE+N+1),IW(NE+2*N+1),IW(NE+3*N+1),IW(NE+4*N+1),
   &   IW(NE+5*N+1),IW(NE+6*N+1))
      GO TO 90
    ENDIF

    IF (JOB.EQ.4) THEN
      DO 50 J = 1,N
      FACT = ZERO
      DO 30 K = IP(J),IP(J+1)-1
        IF (ABS(A(K)).GT.FACT) FACT = ABS(A(K))
   30   CONTINUE
      DO 40 K = IP(J),IP(J+1)-1
        DW(2*N+K) = FACT - ABS(A(K))
   40   CONTINUE
   50   CONTINUE
C B = DW(2N+1:2N+NE); IW(1:5N) and DW(1:2N) are workspaces
      CALL MC64WD(N,NE,IP,IRN,DW(2*N+1),CPERM,NUM,
   &   IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
   &   DW(1),DW(N+1))
      GO TO 90
    ENDIF

    IF (JOB.EQ.5) THEN
      DO 75 J = 1,N
      FACT = ZERO
      DO 60 K = IP(J),IP(J+1)-1
        DW(3*N+K) = ABS(A(K))
        IF (DW(3*N+K).GT.FACT) FACT = DW(3*N+K)
   60   CONTINUE
      DW(2*N+J) = FACT
      IF (FACT.NE.ZERO) THEN
        FACT = LOG(FACT)
      ELSE
        FACT = RINF/N
      ENDIF
      DO 70 K = IP(J),IP(J+1)-1
        IF (DW(3*N+K).NE.ZERO) THEN
          DW(3*N+K) = FACT - LOG(DW(3*N+K))
        ELSE
          DW(3*N+K) = RINF/N
        ENDIF
   70   CONTINUE
   75   CONTINUE
C B = DW(3N+1:3N+NE); IW(1:5N) and DW(1:2N) are workspaces
      CALL MC64WD(N,NE,IP,IRN,DW(3*N+1),CPERM,NUM,
   &   IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
   &   DW(1),DW(N+1))
      IF (NUM.EQ.N) THEN
      DO 80 J = 1,N
        IF (DW(2*N+J).NE.ZERO) THEN
          DW(N+J) = DW(N+J) - LOG(DW(2*N+J))
        ELSE
          DW(N+J) = ZERO
        ENDIF
   80   CONTINUE
      ENDIF
C Check size of scaling factors
      FACT = 0.5*LOG(RINF)
      DO 86 J = 1,N
      IF (DW(J).LT.FACT .AND. DW(N+J).LT.FACT) GO TO 86
      INFO(1) = 2
      GO TO 90
   86   CONTINUE
C     GO TO 90
    ENDIF

   90 IF (INFO(1).EQ.0 .AND. NUM.LT.N) THEN
C Matrix is structurally singular, return with warning
      INFO(1) = 1
      IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9011) INFO(1)
    ENDIF
    IF (INFO(1).EQ.2) THEN
C Scaling factors are large, return with warning
      IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9012) INFO(1)
    ENDIF

C Print diagnostics on output
    IF (ICNTL(3).GE.0) THEN
      WRITE(ICNTL(3),9030) (INFO(J),J=1,2)
      WRITE(ICNTL(3),9031) NUM
      WRITE(ICNTL(3),9032) (CPERM(J),J=1,N)
      IF (JOB.EQ.5) THEN
      WRITE(ICNTL(3),9033) (DW(J),J=1,N)
      WRITE(ICNTL(3),9034) (DW(N+J),J=1,N)
      ENDIF
    ENDIF

C Return from subroutine.
   99 RETURN

 9001 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2,
   &      ' because ',(A),' = ',I10)
 9004 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
   &      '      LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
   &      '      LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
   &      '      Column ',I8,
   &      ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
   &      '      Column ',I8,
   &      ' contains two or more entries with row index ',I8)
 9011 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
   &      '      The matrix is structurally singular.')
 9012 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
   &      '      Some scaling factors may be too large.')
 9020 FORMAT (' ****** Input parameters for MC64A/AD:'/
   &      ' JOB = ',I8/' N   = ',I8/' NE  = ',I8)
 9021 FORMAT (' IP(1:N+1)  = ',8I8/(14X,8I8))
 9022 FORMAT (' IRN(1:NE)  = ',8I8/(14X,8I8))
 9023 FORMAT (' A(1:NE)    = ',4(1PD14.4)/(14X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC64A/AD:'/
   &      ' INFO(1:2)  = ',2I8)
 9031 FORMAT (' NUM      = ',I8)
 9032 FORMAT (' CPERM(1:N) = ',8I8/(14X,8I8))
 9033 FORMAT (' DW(1:N)    = ',5(F11.3)/(14X,5(F11.3)))
 9034 FORMAT (' DW(N+1:2N) = ',5(F11.3)/(14X,5(F11.3)))
    END

C**********************************************************************
    SUBROUTINE MC64BD(N,NE,IP,IRN,A,IPERM,NUM,JPERM,PR,Q,L,D)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER N,NE,NUM
    INTEGER IP(N+1),IRN(NE),IPERM(N),JPERM(N),PR(N),Q(N),L(N)
    DOUBLE PRECISION A(NE),D(N)

C N, NE, IP, IRN are described in MC64A/AD.
C A is a DOUBLE PRECISION array of length
C   NE. A(K), K=1..NE, must be set to the value of the entry
C   that corresponds to IRN(K). It is not altered.
C IPERM is an INTEGER array of length N. On exit, it contains the
C    matching: IPERM(I) = 0 or row I is matched to column IPERM(I).
C NUM is INTEGER variable. On exit, it contains the cardinality of the
C    matching stored in IPERM.
C IW is an INTEGER work array of length 4N.
C DW is a DOUBLE PRECISION array of length N.

C Local variables
    INTEGER I,II,J,JJ,JORD,Q0,QLEN,IDUM,JDUM,ISP,JSP,
   &      K,KK,KK1,KK2,I0,UP,LOW,LPOS
    DOUBLE PRECISION CSP,DI,DNEW,DQ0,AI,A0,BV
C Local parameters
    DOUBLE PRECISION RINF,ZERO,MINONE
    PARAMETER (ZERO=0.0D+0,MINONE=-1.0D+0)
C Intrinsic functions
    INTRINSIC ABS,MIN
C External subroutines and/or functions
    EXTERNAL MC64DD,MC64ED,MC64FD


C Set RINF to largest positive real number
    RINF = HUGE(RINF)

C Initialization
    NUM = 0
    BV = RINF
    DO 10 K = 1,N
      IPERM(K) = 0
      JPERM(K) = 0
      PR(K) = IP(K)
      D(K) = ZERO
   10 CONTINUE
C Scan columns of matrix;
    DO 20 J = 1,N
      A0 = MINONE
      DO 30 K = IP(J),IP(J+1)-1
      I = IRN(K)
      AI = ABS(A(K))
      IF (AI.GT.D(I)) D(I) = AI
      IF (JPERM(J).NE.0) GO TO 30
      IF (AI.GE.BV) THEN
        A0 = BV
        IF (IPERM(I).NE.0) GO TO 30
        JPERM(J) = I
        IPERM(I) = J
        NUM = NUM + 1
      ELSE
        IF (AI.LE.A0) GO TO 30
        A0 = AI
        I0 = I
      ENDIF
   30   CONTINUE
      IF (A0.NE.MINONE .AND. A0.LT.BV) THEN
      BV = A0
      IF (IPERM(I0).NE.0) GO TO 20
      IPERM(I0) = J
      JPERM(J) = I0
      NUM = NUM + 1
      ENDIF
   20 CONTINUE
C Update BV with smallest of all the largest maximum absolute values
C of the rows.
    DO 25 I = 1,N
      BV = MIN(BV,D(I))
   25 CONTINUE
    IF (NUM.EQ.N) GO TO 1000
C Rescan unassigned columns; improve initial assignment
    DO 95 J = 1,N
      IF (JPERM(J).NE.0) GO TO 95
      DO 50 K = IP(J),IP(J+1)-1
      I = IRN(K)
      AI = ABS(A(K))
      IF (AI.LT.BV) GO TO 50
      IF (IPERM(I).EQ.0) GO TO 90
      JJ = IPERM(I)
      KK1 = PR(JJ)
      KK2 = IP(JJ+1) - 1
      IF (KK1.GT.KK2) GO TO 50
      DO 70 KK = KK1,KK2
        II = IRN(KK)
        IF (IPERM(II).NE.0) GO TO 70
        IF (ABS(A(KK)).GE.BV) GO TO 80
   70   CONTINUE
      PR(JJ) = KK2 + 1
   50   CONTINUE
      GO TO 95
   80   JPERM(JJ) = II
      IPERM(II) = JJ
      PR(JJ) = KK + 1
   90   NUM = NUM + 1
      JPERM(J) = I
      IPERM(I) = J
      PR(J) = K + 1
   95 CONTINUE
    IF (NUM.EQ.N) GO TO 1000

C Prepare for main loop
    DO 99 I = 1,N
      D(I) = MINONE
      L(I) = 0
   99 CONTINUE

C Main loop ... each pass round this loop is similar to Dijkstra's
C algorithm for solving the single source shortest path problem

    DO 100 JORD = 1,N

      IF (JPERM(JORD).NE.0) GO TO 100
      QLEN = 0
      LOW = N + 1
      UP = N + 1
C CSP is cost of shortest path to any unassigned row
C ISP is matrix position of unassigned row element in shortest path
C JSP is column index of unassigned row element in shortest path
      CSP = MINONE
C Build shortest path tree starting from unassigned column JORD
      J = JORD
      PR(J) = -1

C Scan column J
      DO 115 K = IP(J),IP(J+1)-1
      I = IRN(K)
      DNEW = ABS(A(K))
      IF (CSP.GE.DNEW) GO TO 115
      IF (IPERM(I).EQ.0) THEN
C Row I is unassigned; update shortest path info
        CSP = DNEW
        ISP = I
        JSP = J
        IF (CSP.GE.BV) GO TO 160
      ELSE
        D(I) = DNEW
        IF (DNEW.GE.BV) THEN
C Add row I to Q2
          LOW = LOW - 1
          Q(LOW) = I
        ELSE
C Add row I to Q, and push it
          QLEN = QLEN + 1
          L(I) = QLEN
          CALL MC64DD(I,N,Q,D,L,1)
        ENDIF
        JJ = IPERM(I)
        PR(JJ) = J
      ENDIF
  115   CONTINUE

      DO 150 JDUM = 1,NUM
C If Q2 is empty, extract new rows from Q
      IF (LOW.EQ.UP) THEN
        IF (QLEN.EQ.0) GO TO 160
        I = Q(1)
        IF (CSP.GE.D(I)) GO TO 160
        BV = D(I)
        DO 152 IDUM = 1,N
          CALL MC64ED(QLEN,N,Q,D,L,1)
          L(I) = 0
          LOW = LOW - 1
          Q(LOW) = I
          IF (QLEN.EQ.0) GO TO 153
          I = Q(1)
          IF (D(I).NE.BV) GO TO 153
  152     CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C Move row Q0
  153   UP = UP - 1
      Q0 = Q(UP)
      DQ0 = D(Q0)
      L(Q0) = UP
C Scan column that matches with row Q0
      J = IPERM(Q0)
      DO 155 K = IP(J),IP(J+1)-1
        I = IRN(K)
C Update D(I)
        IF (L(I).GE.UP) GO TO 155
        DNEW = MIN(DQ0,ABS(A(K)))
        IF (CSP.GE.DNEW) GO TO 155
        IF (IPERM(I).EQ.0) THEN
C Row I is unassigned; update shortest path info
          CSP = DNEW
          ISP = I
          JSP = J
          IF (CSP.GE.BV) GO TO 160
        ELSE
          DI = D(I)
          IF (DI.GE.BV .OR. DI.GE.DNEW) GO TO 155
          D(I) = DNEW
          IF (DNEW.GE.BV) THEN
C Delete row I from Q (if necessary); add row I to Q2
          IF (DI.NE.MINONE) THEN
            LPOS = L(I)
            CALL MC64FD(LPOS,QLEN,N,Q,D,L,1)
          ENDIF
          L(I) = 0
          LOW = LOW - 1
          Q(LOW) = I
          ELSE
C Add row I to Q (if necessary); push row I up Q
          IF (DI.EQ.MINONE) THEN
            QLEN = QLEN + 1
            L(I) = QLEN
          ENDIF
          CALL MC64DD(I,N,Q,D,L,1)
          ENDIF
C Update tree
          JJ = IPERM(I)
          PR(JJ) = J
        ENDIF
  155   CONTINUE
  150   CONTINUE

C If CSP = MINONE, no augmenting path is found
  160   IF (CSP.EQ.MINONE) GO TO 190
C Update bottleneck value
      BV = MIN(BV,CSP)
C Find augmenting path by tracing backward in PR; update IPERM,JPERM
      NUM = NUM + 1
      I = ISP
      J = JSP
      DO 170 JDUM = 1,NUM+1
      I0 = JPERM(J)
      JPERM(J) = I
      IPERM(I) = J
      J = PR(J)
      IF (J.EQ.-1) GO TO 190
      I = I0
  170   CONTINUE
C End of dummy loop; this point is never reached
  190   DO 191 KK = UP,N
      I = Q(KK)
      D(I) = MINONE
      L(I) = 0
  191   CONTINUE
      DO 192 KK = LOW,UP-1
      I = Q(KK)
      D(I) = MINONE
  192   CONTINUE
      DO 193 KK = 1,QLEN
      I = Q(KK)
      D(I) = MINONE
      L(I) = 0
  193   CONTINUE

  100 CONTINUE
C End of main loop

C BV is bottleneck value of final matching
    IF (NUM.EQ.N) GO TO 1000

C Matrix is structurally singular, complete IPERM.
C JPERM, PR are work arrays
    DO 300 J = 1,N
      JPERM(J) = 0
  300 CONTINUE
    K = 0
    DO 310 I = 1,N
      IF (IPERM(I).EQ.0) THEN
      K = K + 1
      PR(K) = I
      ELSE
      J = IPERM(I)
      JPERM(J) = I
      ENDIF
  310 CONTINUE
    K = 0
    DO 320 I = 1,N
      IF (JPERM(I).NE.0) GO TO 320
      K = K + 1
      JDUM = PR(K)
      IPERM(JDUM) = - I
  320 CONTINUE

 1000 RETURN
    END

C**********************************************************************
    SUBROUTINE MC64DD(I,N,Q,D,L,IWAY)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER I,N,IWAY
    INTEGER Q(N),L(N)
    DOUBLE PRECISION D(N)

C Variables N,Q,D,L are described in MC64B/BD
C IF IWAY is equal to 1, then
C node I is pushed from its current position upwards
C IF IWAY is not equal to 1, then
C node I is pushed from its current position downwards

C Local variables and parameters
    INTEGER IDUM,K,POS,POSK,QK
    PARAMETER (K=2)
    DOUBLE PRECISION DI

    POS = L(I)
    IF (POS.LE.1) GO TO 20
    DI = D(I)
C POS is index of current position of I in the tree
    IF (IWAY.EQ.1) THEN
      DO 10 IDUM = 1,N
      POSK = POS/K
      QK = Q(POSK)
      IF (DI.LE.D(QK)) GO TO 20
      Q(POS) = QK
      L(QK) = POS
      POS = POSK
      IF (POS.LE.1) GO TO 20
   10   CONTINUE
C End of dummy loop; this point is never reached
    ELSE
      DO 15 IDUM = 1,N
      POSK = POS/K
      QK = Q(POSK)
      IF (DI.GE.D(QK)) GO TO 20
      Q(POS) = QK
      L(QK) = POS
      POS = POSK
      IF (POS.LE.1) GO TO 20
   15   CONTINUE
C End of dummy loop; this point is never reached
    ENDIF
C End of dummy if; this point is never reached
   20 Q(POS) = I
    L(I) = POS

    RETURN
    END

C**********************************************************************
    SUBROUTINE MC64ED(QLEN,N,Q,D,L,IWAY)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER QLEN,N,IWAY
    INTEGER Q(N),L(N)
    DOUBLE PRECISION D(N)

C Variables QLEN,N,Q,D,L are described in MC64B/BD (IWAY = 1) or
C   MC64W/WD (IWAY = 2)
C The root node is deleted from the binary heap.

C Local variables and parameters
    INTEGER I,IDUM,K,POS,POSK
    PARAMETER (K=2)
    DOUBLE PRECISION DK,DR,DI

C Move last element to begin of Q
    I = Q(QLEN)
    DI = D(I)
    QLEN = QLEN - 1
    POS = 1
    IF (IWAY.EQ.1) THEN
      DO 10 IDUM = 1,N
      POSK = K*POS
      IF (POSK.GT.QLEN) GO TO 20
      DK = D(Q(POSK))
      IF (POSK.LT.QLEN) THEN
        DR = D(Q(POSK+1))
        IF (DK.LT.DR) THEN
          POSK = POSK + 1
          DK = DR
        ENDIF
      ENDIF
      IF (DI.GE.DK) GO TO 20
C Exchange old last element with larger priority child
      Q(POS) = Q(POSK)
      L(Q(POS)) = POS
      POS = POSK
   10   CONTINUE
C End of dummy loop; this point is never reached
    ELSE
      DO 15 IDUM = 1,N
      POSK = K*POS
      IF (POSK.GT.QLEN) GO TO 20
      DK = D(Q(POSK))
      IF (POSK.LT.QLEN) THEN
        DR = D(Q(POSK+1))
        IF (DK.GT.DR) THEN
          POSK = POSK + 1
          DK = DR
        ENDIF
      ENDIF
      IF (DI.LE.DK) GO TO 20
C Exchange old last element with smaller child
      Q(POS) = Q(POSK)
      L(Q(POS)) = POS
      POS = POSK
   15   CONTINUE
C End of dummy loop; this point is never reached
    ENDIF
C End of dummy if; this point is never reached
   20 Q(POS) = I
    L(I) = POS

    RETURN
    END

C**********************************************************************
    SUBROUTINE MC64FD(POS0,QLEN,N,Q,D,L,IWAY)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER POS0,QLEN,N,IWAY
    INTEGER Q(N),L(N)
    DOUBLE PRECISION D(N)

C Variables QLEN,N,Q,D,L are described in MC64B/BD (IWAY = 1) or
C   MC64WD (IWAY = 2).
C Move last element in the heap

    INTEGER I,IDUM,K,POS,POSK,QK
    PARAMETER (K=2)
    DOUBLE PRECISION DK,DR,DI

C Quick return, if possible
    IF (QLEN.EQ.POS0) THEN
      QLEN = QLEN - 1
      RETURN
    ENDIF

C Move last element from queue Q to position POS0
C POS is current position of node I in the tree
    I = Q(QLEN)
    DI = D(I)
    QLEN = QLEN - 1
    POS = POS0
    IF (IWAY.EQ.1) THEN
      IF (POS.LE.1) GO TO 20
      DO 10 IDUM = 1,N
      POSK = POS/K
      QK = Q(POSK)
      IF (DI.LE.D(QK)) GO TO 20
      Q(POS) = QK
      L(QK) = POS
      POS = POSK
      IF (POS.LE.1) GO TO 20
   10   CONTINUE
C End of dummy loop; this point is never reached
   20   Q(POS) = I
      L(I) = POS
      IF (POS.NE.POS0) RETURN
      DO 30 IDUM = 1,N
      POSK = K*POS
      IF (POSK.GT.QLEN) GO TO 40
      DK = D(Q(POSK))
      IF (POSK.LT.QLEN) THEN
        DR = D(Q(POSK+1))
        IF (DK.LT.DR) THEN
          POSK = POSK + 1
          DK = DR
        ENDIF
      ENDIF
      IF (DI.GE.DK) GO TO 40
      QK = Q(POSK)
      Q(POS) = QK
      L(QK) = POS
      POS = POSK
   30   CONTINUE
C End of dummy loop; this point is never reached
    ELSE
      IF (POS.LE.1) GO TO 34
      DO 32 IDUM = 1,N
      POSK = POS/K
      QK = Q(POSK)
      IF (DI.GE.D(QK)) GO TO 34
      Q(POS) = QK
      L(QK) = POS
      POS = POSK
      IF (POS.LE.1) GO TO 34
   32   CONTINUE
C End of dummy loop; this point is never reached
   34   Q(POS) = I
      L(I) = POS
      IF (POS.NE.POS0) RETURN
      DO 36 IDUM = 1,N
      POSK = K*POS
      IF (POSK.GT.QLEN) GO TO 40
      DK = D(Q(POSK))
      IF (POSK.LT.QLEN) THEN
        DR = D(Q(POSK+1))
        IF (DK.GT.DR) THEN
          POSK = POSK + 1
          DK = DR
        ENDIF
      ENDIF
      IF (DI.LE.DK) GO TO 40
      QK = Q(POSK)
      Q(POS) = QK
      L(QK) = POS
      POS = POSK
   36   CONTINUE
C End of dummy loop; this point is never reached
    ENDIF
C End of dummy if; this point is never reached
   40 Q(POS) = I
    L(I) = POS

    RETURN
    END

C**********************************************************************
    SUBROUTINE MC64RD(N,NE,IP,IRN,A)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER N,NE
    INTEGER IP(N+1),IRN(NE)
    DOUBLE PRECISION A(NE)

C This subroutine sorts the entries in each column of the
C sparse matrix (defined by N,NE,IP,IRN,A) by decreasing
C numerical value.

C Local constants
    INTEGER THRESH,TDLEN
    PARAMETER (THRESH=15,TDLEN=50)
C Local variables
    INTEGER J,IPJ,K,LEN,R,S,HI,FIRST,MID,LAST,TD
    DOUBLE PRECISION HA,KEY
C Local arrays
    INTEGER TODO(TDLEN)

    DO 100 J = 1,N
      LEN = IP(J+1) - IP(J)
      IF (LEN.LE.1) GO TO 100
      IPJ = IP(J)

C Sort array roughly with partial quicksort
      IF (LEN.LT.THRESH) GO TO 400
      TODO(1) = IPJ
      TODO(2) = IPJ + LEN
      TD = 2
  500   CONTINUE
      FIRST = TODO(TD-1)
      LAST = TODO(TD)
C KEY is the smallest of two values present in interval [FIRST,LAST)
      KEY = A((FIRST+LAST)/2)
      DO 475 K = FIRST,LAST-1
      HA = A(K)
      IF (HA.EQ.KEY) GO TO 475
      IF (HA.GT.KEY) GO TO 470
      KEY = HA
      GO TO 470
  475   CONTINUE
C Only one value found in interval, so it is already sorted
      TD = TD - 2
      GO TO 425

C Reorder interval [FIRST,LAST) such that entries before MID are gt KEY
  470   MID = FIRST
      DO 450 K = FIRST,LAST-1
      IF (A(K).LE.KEY) GO TO 450
      HA = A(MID)
      A(MID) = A(K)
      A(K) = HA
      HI = IRN(MID)
      IRN(MID) = IRN(K)
      IRN(K) = HI
      MID = MID + 1
  450   CONTINUE
C Both subintervals [FIRST,MID), [MID,LAST) are nonempty
C Stack the longest of the two subintervals first
      IF (MID-FIRST.GE.LAST-MID) THEN
      TODO(TD+2) = LAST
      TODO(TD+1) = MID
      TODO(TD) = MID
C      TODO(TD-1) = FIRST
      ELSE
      TODO(TD+2) = MID
      TODO(TD+1) = FIRST
      TODO(TD) = LAST
      TODO(TD-1) = MID
      ENDIF
      TD = TD + 2

  425   CONTINUE
      IF (TD.EQ.0) GO TO 400
C There is still work to be done
      IF (TODO(TD)-TODO(TD-1).GE.THRESH) GO TO 500
C Next interval is already short enough for straightforward insertion
      TD = TD - 2
      GO TO 425

C Complete sorting with straightforward insertion
  400   DO 200 R = IPJ+1,IPJ+LEN-1
      IF (A(R-1) .LT. A(R)) THEN
        HA = A(R)
        HI = IRN(R)
        A(R) = A(R-1)
        IRN(R) = IRN(R-1)
        DO 300 S = R-1,IPJ+1,-1
          IF (A(S-1) .LT. HA) THEN
          A(S) = A(S-1)
          IRN(S) = IRN(S-1)
          ELSE
          A(S) = HA
          IRN(S) = HI
          GO TO 200
          END IF
  300     CONTINUE
        A(IPJ) = HA
        IRN(IPJ) = HI
      END IF
  200   CONTINUE

  100 CONTINUE

    RETURN
    END

C**********************************************************************
    SUBROUTINE MC64SD(N,NE,IP,IRN,A,IPERM,NUMX,
   &       W,LEN,LENL,LENH,FC,IW,IW4)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER N,NE,NUMX
    INTEGER IP(N+1),IRN(NE),IPERM(N),
   &      W(N),LEN(N),LENL(N),LENH(N),FC(N),IW(N),IW4(4*N)
    DOUBLE PRECISION A(NE)

C N, NE, IP, IRN, are described in MC64A/AD.
C A is a DOUBLE PRECISION array of length NE.
C   A(K), K=1..NE, must be set to the value of the entry that
C   corresponds to IRN(k). The entries in each column must be
C   non-negative and ordered by decreasing value.
C IPERM is an INTEGER array of length N. On exit, it contains the
C   bottleneck matching: IPERM(I) - 0 or row I is matched to column
C   IPERM(I).
C NUMX is an INTEGER variable. On exit, it contains the cardinality
C   of the matching stored in IPERM.
C IW is an INTEGER work array of length 10N.

C FC is an integer array of length N that contains the list of
C   unmatched columns.
C LEN(J), LENL(J), LENH(J) are integer arrays of length N that point
C   to entries in matrix column J.
C   In the matrix defined by the column parts IP(J)+LENL(J) we know
C   a matching does not exist; in the matrix defined by the column
C   parts IP(J)+LENH(J) we know one exists.
C   LEN(J) lies between LENL(J) and LENH(J) and determines the matrix
C   that is tested for a maximum matching.
C W is an integer array of length N and contains the indices of the
C   columns for which LENL ne LENH.
C WLEN is number of indices stored in array W.
C IW is integer work array of length N.
C IW4 is integer work array of length 4N used by MC64U/UD.

    INTEGER NUM,NVAL,WLEN,II,I,J,K,L,CNT,MOD,IDUM1,IDUM2,IDUM3
    DOUBLE PRECISION BVAL,BMIN,BMAX,RINF
    EXTERNAL MC64QD,MC64UD

C BMIN and BMAX are such that a maximum matching exists for the input
C   matrix in which all entries smaller than BMIN are dropped.
C   For BMAX, a maximum matching does not exist.
C BVAL is a value between BMIN and BMAX.
C CNT is the number of calls made to MC64U/UD so far.
C NUM is the cardinality of last matching found.

C Set RINF to largest positive real number
    RINF = HUGE(RINF)

C Compute a first maximum matching from scratch on whole matrix.
    DO 20 J = 1,N
      FC(J) = J
      IW(J) = 0
      LEN(J) = IP(J+1) - IP(J)
   20 CONTINUE
C The first call to MC64U/UD
    CNT = 1
    MOD = 1
    NUMX = 0
    CALL MC64UD(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUMX,N,
   &        IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))

C IW contains a maximum matching of length NUMX.
    NUM = NUMX

    IF (NUM.NE.N) THEN
C Matrix is structurally singular
      BMAX = RINF
    ELSE
C Matrix is structurally nonsingular, NUM=NUMX=N;
C Set BMAX just above the smallest of all the maximum absolute
C values of the columns
      BMAX = RINF
      DO 30 J = 1,N
      BVAL = 0.0
      DO 25 K = IP(J),IP(J+1)-1
        IF (A(K).GT.BVAL) BVAL = A(K)
   25   CONTINUE
      IF (BVAL.LT.BMAX) BMAX = BVAL
   30   CONTINUE
      BMAX = 1.001 * BMAX
    ENDIF

C Initialize BVAL,BMIN
    BVAL = 0.0
    BMIN = 0.0
C Initialize LENL,LEN,LENH,W,WLEN according to BMAX.
C Set LEN(J), LENH(J) just after last entry in column J.
C Set LENL(J) just after last entry in column J with value ge BMAX.
    WLEN = 0
    DO 48 J = 1,N
      L = IP(J+1) - IP(J)
      LENH(J) = L
      LEN(J) = L
      DO 45 K = IP(J),IP(J+1)-1
      IF (A(K).LT.BMAX) GO TO 46
   45   CONTINUE
C Column J is empty or all entries are ge BMAX
      K = IP(J+1)
   46   LENL(J) = K - IP(J)
C Add J to W if LENL(J) ne LENH(J)
      IF (LENL(J).EQ.L) GO TO 48
      WLEN = WLEN + 1
      W(WLEN) = J
   48 CONTINUE

C Main loop
    DO 90 IDUM1 = 1,NE
      IF (NUM.EQ.NUMX) THEN
C We have a maximum matching in IW; store IW in IPERM
      DO 50 I = 1,N
        IPERM(I) = IW(I)
   50   CONTINUE
C Keep going round this loop until matching IW is no longer maximum.
      DO 80 IDUM2 = 1,NE
        BMIN = BVAL
        IF (BMAX .EQ. BMIN) GO TO 99
C Find splitting value BVAL
        CALL MC64QD(IP,LENL,LEN,W,WLEN,A,NVAL,BVAL)
        IF (NVAL.LE.1) GO TO 99
C Set LEN such that all matrix entries with value lt BVAL are
C discarded. Store old LEN in LENH. Do this for all columns W(K).
C Each step, either K is incremented or WLEN is decremented.
        K = 1
        DO 70 IDUM3 = 1,N
          IF (K.GT.WLEN) GO TO 71
          J = W(K)
          DO 55 II = IP(J)+LEN(J)-1,IP(J)+LENL(J),-1
          IF (A(II).GE.BVAL) GO TO 60
          I = IRN(II)
          IF (IW(I).NE.J) GO TO 55
C Remove entry from matching
          IW(I) = 0
          NUM = NUM - 1
          FC(N-NUM) = J
   55       CONTINUE
   60       LENH(J) = LEN(J)
C IP(J)+LEN(J)-1 is last entry in column ge BVAL
          LEN(J) = II - IP(J) + 1
C If LENH(J) = LENL(J), remove J from W
          IF (LENL(J).EQ.LENH(J)) THEN
          W(K) = W(WLEN)
          WLEN = WLEN - 1
          ELSE
          K = K + 1
          ENDIF
   70     CONTINUE
   71     IF (NUM.LT.NUMX) GO TO 81
   80   CONTINUE
C End of dummy loop; this point is never reached
C Set mode for next call to MC64U/UD
   81   MOD = 1
      ELSE
C We do not have a maximum matching in IW.
      BMAX = BVAL
C BMIN is the bottleneck value of a maximum matching;
C for BMAX the matching is not maximum, so BMAX>BMIN
C      IF (BMAX .EQ. BMIN) GO TO 99
C Find splitting value BVAL
      CALL MC64QD(IP,LEN,LENH,W,WLEN,A,NVAL,BVAL)
      IF (NVAL.EQ.0. OR. BVAL.EQ.BMIN) GO TO 99
C Set LEN such that all matrix entries with value ge BVAL are
C inside matrix. Store old LEN in LENL. Do this for all columns W(K).
C Each step, either K is incremented or WLEN is decremented.
      K = 1
      DO 87 IDUM3 = 1,N
        IF (K.GT.WLEN) GO TO 88
        J = W(K)
        DO 85 II = IP(J)+LEN(J),IP(J)+LENH(J)-1
          IF (A(II).LT.BVAL) GO TO 86
   85     CONTINUE
   86     LENL(J) = LEN(J)
        LEN(J) = II - IP(J)
        IF (LENL(J).EQ.LENH(J)) THEN
          W(K) = W(WLEN)
          WLEN = WLEN - 1
        ELSE
          K = K + 1
        ENDIF
   87   CONTINUE
C End of dummy loop; this point is never reached
C Set mode for next call to MC64U/UD
   88   MOD = 0
      ENDIF
      CNT = CNT + 1
      CALL MC64UD(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUM,NUMX,
   &        IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))

C IW contains maximum matching of length NUM
   90 CONTINUE
C End of dummy loop; this point is never reached

C BMIN is bottleneck value of final matching
   99 IF (NUMX.EQ.N) GO TO 1000
C The matrix is structurally singular, complete IPERM
C W, IW are work arrays
    DO 300 J = 1,N
      W(J) = 0
  300 CONTINUE
    K = 0
    DO 310 I = 1,N
      IF (IPERM(I).EQ.0) THEN
      K = K + 1
      IW(K) = I
      ELSE
      J = IPERM(I)
      W(J) = I
      ENDIF
  310 CONTINUE
    K = 0
    DO 320 J = 1,N
      IF (W(J).NE.0) GO TO 320
      K = K + 1
      IDUM1 = IW(K)
      IPERM(IDUM1) =  - J
  320 CONTINUE

 1000 RETURN
    END

C**********************************************************************
    SUBROUTINE MC64QD(IP,LENL,LENH,W,WLEN,A,NVAL,VAL)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER WLEN,NVAL
    INTEGER IP(*),LENL(*),LENH(*),W(*)
    DOUBLE PRECISION A(*),VAL

C This routine searches for at most XX different numerical values
C in the columns W(1:WLEN). XX>=2.
C Each column J is scanned between IP(J)+LENL(J) and IP(J)+LENH(J)-1
C until XX values are found or all columns have been considered.
C On output, NVAL is the number of different values that is found
C and SPLIT(1:NVAL) contains the values in decreasing order.
C If NVAL > 0, the routine returns VAL = SPLIT((NVAL+1)/2).
C
    INTEGER XX,J,K,II,S,POS
    PARAMETER (XX=10)
    DOUBLE PRECISION SPLIT(XX),HA

C Scan columns in W(1:WLEN). For each encountered value, if value not
C already present in SPLIT(1:NVAL), insert value such that SPLIT
C remains sorted by decreasing value.
C The sorting is done by straightforward insertion; therefore the use
C of this routine should be avoided for large XX (XX < 20).
    NVAL = 0
    DO 10 K = 1,WLEN
      J = W(K)
      DO 15 II = IP(J)+LENL(J),IP(J)+LENH(J)-1
      HA = A(II)
      IF (NVAL.EQ.0) THEN
        SPLIT(1) = HA
        NVAL = 1
      ELSE
C Check presence of HA in SPLIT
        DO 20 S = NVAL,1,-1
          IF (SPLIT(S).EQ.HA) GO TO 15
          IF (SPLIT(S).GT.HA) THEN
          POS = S + 1
          GO TO 21
          ENDIF
  20      CONTINUE
        POS = 1
C The insertion
  21      DO 22 S = NVAL,POS,-1
          SPLIT(S+1) = SPLIT(S)
  22      CONTINUE
        SPLIT(POS) = HA
        NVAL = NVAL + 1
      ENDIF
C Exit loop if XX values are found
      IF (NVAL.EQ.XX) GO TO 11
  15    CONTINUE
  10  CONTINUE
C Determine VAL
  11  IF (NVAL.GT.0) VAL = SPLIT((NVAL+1)/2)

    RETURN
    END

C**********************************************************************
    SUBROUTINE MC64UD(ID,MOD,N,IRN,LIRN,IP,LENC,FC,IPERM,NUM,NUMX,
   &       PR,ARP,CV,OUT)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER ID,MOD,N,LIRN,NUM,NUMX
    INTEGER ARP(N),CV(N),IRN(LIRN),IP(N),
   &      FC(N),IPERM(N),LENC(N),OUT(N),PR(N)

C PR(J) is the previous column to J in the depth first search.
C   Array PR is used as workspace in the sorting algorithm.
C Elements (I,IPERM(I)) I=1,..,N are entries at the end of the
C   algorithm unless N assignments have not been made in which case
C   N-NUM pairs (I,IPERM(I)) will not be entries in the matrix.
C CV(I) is the most recent loop number (ID+JORD) at which row I
C   was visited.
C ARP(J) is the number of entries in column J which have been scanned
C   when looking for a cheap assignment.
C OUT(J) is one less than the number of entries in column J which have
C   not been scanned during one pass through the main loop.
C NUMX is maximum possible size of matching.

    INTEGER I,II,IN1,IN2,J,J1,JORD,K,KK,LAST,NFC,
   &      NUM0,NUM1,NUM2,ID0,ID1

    IF (ID.EQ.1) THEN
C The first call to MC64U/UD.
C Initialize CV and ARP; parameters MOD, NUMX are not accessed
      DO 5 I = 1,N
      CV(I) = 0
      ARP(I) = 0
    5   CONTINUE
      NUM1 = N
      NUM2 = N
    ELSE
C Not the first call to MC64U/UD.
C Re-initialize ARP if entries were deleted since last call to MC64U/UD
      IF (MOD.EQ.1) THEN
      DO 8 I = 1,N
        ARP(I) = 0
    8   CONTINUE
      ENDIF
      NUM1 = NUMX
      NUM2 = N - NUMX
    ENDIF
    NUM0 = NUM

C NUM0 is size of input matching
C NUM1 is maximum possible size of matching
C NUM2 is maximum allowed number of unassigned rows/columns
C NUM is size of current matching

C Quick return if possible
C    IF (NUM.EQ.N) GO TO 199
C NFC is number of rows/columns that could not be assigned
    NFC = 0
C Integers ID0+1 to ID0+N are unique numbers for call ID to MC64U/UD,
C so 1st call uses 1..N, 2nd call uses N+1..2N, etc
    ID0 = (ID-1)*N

C Main loop. Each pass round this loop either results in a new
C assignment or gives a column with no assignment

    DO 100 JORD = NUM0+1,N

C Each pass uses unique number ID1
      ID1 = ID0 + JORD
C J is unmatched column
      J = FC(JORD-NUM0)
      PR(J) = -1
      DO 70 K = 1,JORD
C Look for a cheap assignment
      IF (ARP(J).GE.LENC(J)) GO TO 30
      IN1 = IP(J) + ARP(J)
      IN2 = IP(J) + LENC(J) - 1
      DO 20 II = IN1,IN2
        I = IRN(II)
        IF (IPERM(I).EQ.0) GO TO 80
   20   CONTINUE
C No cheap assignment in row
      ARP(J) = LENC(J)
C Begin looking for assignment chain starting with row J
   30   OUT(J) = LENC(J) - 1
C Inner loop.  Extends chain by one or backtracks
      DO 60 KK = 1,JORD
        IN1 = OUT(J)
        IF (IN1.LT.0) GO TO 50
        IN2 = IP(J) + LENC(J) - 1
        IN1 = IN2 - IN1
C Forward scan
        DO 40 II = IN1,IN2
          I = IRN(II)
          IF (CV(I).EQ.ID1) GO TO 40
C Column J has not yet been accessed during this pass
          J1 = J
          J = IPERM(I)
          CV(I) = ID1
          PR(J) = J1
          OUT(J1) = IN2 - II - 1
          GO TO 70
   40     CONTINUE
C Backtracking step.
   50     J1 = PR(J)
        IF (J1.EQ.-1) THEN
C No augmenting path exists for column J.
          NFC = NFC + 1
          FC(NFC) = J
          IF (NFC.GT.NUM2) THEN
C A matching of maximum size NUM1 is not possible
          LAST = JORD
          GO TO 101
          ENDIF
          GO TO 100
        ENDIF
        J = J1
   60   CONTINUE
C End of dummy loop; this point is never reached
   70   CONTINUE
C End of dummy loop; this point is never reached

C New assignment is made.
   80   IPERM(I) = J
      ARP(J) = II - IP(J) + 1
      NUM = NUM + 1
      DO 90 K = 1,JORD
      J = PR(J)
      IF (J.EQ.-1) GO TO 95
      II = IP(J) + LENC(J) - OUT(J) - 2
      I = IRN(II)
      IPERM(I) = J
   90   CONTINUE
C End of dummy loop; this point is never reached

   95   IF (NUM.EQ.NUM1) THEN
C A matching of maximum size NUM1 is found
      LAST = JORD
      GO TO 101
      ENDIF
C
  100 CONTINUE

C All unassigned columns have been considered
    LAST = N

C Now, a transversal is computed or is not possible.
C Complete FC before returning.
  101 DO 110 JORD = LAST+1,N
      NFC = NFC + 1
      FC(NFC) = FC(JORD-NUM0)
  110 CONTINUE

C  199 RETURN
    RETURN
    END

C**********************************************************************
    SUBROUTINE MC64WD(N,NE,IP,IRN,A,IPERM,NUM,
   &       JPERM,OUT,PR,Q,L,U,D)
    IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C   Research Councils                           ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we       ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                         ***
C *** Any problems?   Contact ...                     ***
C   Iain Duff (I.Duff@rl.ac.uk) or                    ***
C   Jacko Koster (jacko.koster@uninett.no)                ***
C
    INTEGER N,NE,NUM
    INTEGER IP(N+1),IRN(NE),IPERM(N),
   &      JPERM(N),OUT(N),PR(N),Q(N),L(N)
    DOUBLE PRECISION A(NE),U(N),D(N)

C N, NE, IP, IRN are described in MC64A/AD.
C A is a DOUBLE PRECISION array of length NE.
C   A(K), K=1..NE, must be set to the value of the entry that
C   corresponds to IRN(K). It is not altered.
C   All values A(K) must be non-negative.
C IPERM is an INTEGER array of length N. On exit, it contains the
C   weighted matching: IPERM(I) = 0 or row I is matched to column
C   IPERM(I).
C NUM is an INTEGER variable. On exit, it contains the cardinality of
C   the matching stored in IPERM.
C IW is an INTEGER work array of length 5N.
C DW is a DOUBLE PRECISION array of length 2N.
C   On exit, U = D(1:N) contains the dual row variable and
C   V = D(N+1:2N) contains the dual column variable. If the matrix
C   is structurally nonsingular (NUM = N), the following holds:
C    U(I)+V(J) <= A(I,J)  if IPERM(I) |= J
C    U(I)+V(J)  = A(I,J)  if IPERM(I)  = J
C    U(I) = 0  if IPERM(I) = 0
C    V(J) = 0  if there is no I for which IPERM(I) = J

C Local variables
    INTEGER I,I0,II,J,JJ,JORD,Q0,QLEN,JDUM,ISP,JSP,
   &      K,K0,K1,K2,KK,KK1,KK2,UP,LOW,LPOS
    DOUBLE PRECISION CSP,DI,DMIN,DNEW,DQ0,VJ
C Local parameters
    DOUBLE PRECISION RINF,ZERO
    PARAMETER (ZERO=0.0D+0)
C External subroutines and/or functions
    EXTERNAL MC64DD,MC64ED,MC64FD


C Set RINF to largest positive real number
    RINF = HUGE(RINF)

C Initialization
    NUM = 0
    DO 10 K = 1,N
      U(K) = RINF
      D(K) = ZERO
      IPERM(K) = 0
      JPERM(K) = 0
      PR(K) = IP(K)
      L(K) = 0
   10 CONTINUE
C Initialize U(I)
    DO 30 J = 1,N
      DO 20 K = IP(J),IP(J+1)-1
      I = IRN(K)
      IF (A(K).GT.U(I)) GO TO 20
      U(I) = A(K)
      IPERM(I) = J
      L(I) = K
   20   CONTINUE
   30 CONTINUE
    DO 40 I = 1,N
      J = IPERM(I)
      IF (J.EQ.0) GO TO 40
C Row I is not empty
      IPERM(I) = 0
      IF (JPERM(J).NE.0) GO TO 40
C Don't choose cheap assignment from dense columns
      IF (IP(J+1)-IP(J) .GT. N/10 .AND. N.GT.50) GO TO 40
C Assignment of column J to row I
      NUM = NUM + 1
      IPERM(I) = J
      JPERM(J) = L(I)
   40 CONTINUE
    IF (NUM.EQ.N) GO TO 1000
C Scan unassigned columns; improve assignment
    DO 95 J = 1,N
C JPERM(J) ne 0 iff column J is already assigned
      IF (JPERM(J).NE.0) GO TO 95
      K1 = IP(J)
      K2 = IP(J+1) - 1
C Continue only if column J is not empty
      IF (K1.GT.K2) GO TO 95
C     VJ = RINF
C Changes made to allow for NaNs
      I0 = IRN(K1)
      VJ = A(K1) - U(I0)
      K0 = K1
      DO 50 K = K1+1,K2
      I = IRN(K)
      DI = A(K) - U(I)
      IF (DI.GT.VJ) GO TO 50
      IF (DI.LT.VJ .OR. DI.EQ.RINF) GO TO 55
      IF (IPERM(I).NE.0 .OR. IPERM(I0).EQ.0) GO TO 50
   55   VJ = DI
      I0 = I
      K0 = K
   50   CONTINUE
      D(J) = VJ
      K = K0
      I = I0
      IF (IPERM(I).EQ.0) GO TO 90
      DO 60 K = K0,K2
      I = IRN(K)
      IF (A(K)-U(I).GT.VJ) GO TO 60
      JJ = IPERM(I)
C Scan remaining part of assigned column JJ
      KK1 = PR(JJ)
      KK2 = IP(JJ+1) - 1
      IF (KK1.GT.KK2) GO TO 60
      DO 70 KK = KK1,KK2
        II = IRN(KK)
        IF (IPERM(II).GT.0) GO TO 70
        IF (A(KK)-U(II).LE.D(JJ)) GO TO 80
   70   CONTINUE
      PR(JJ) = KK2 + 1
   60   CONTINUE
      GO TO 95
   80   JPERM(JJ) = KK
      IPERM(II) = JJ
      PR(JJ) = KK + 1
   90   NUM = NUM + 1
      JPERM(J) = K
      IPERM(I) = J
      PR(J) = K + 1
   95 CONTINUE
    IF (NUM.EQ.N) GO TO 1000

C Prepare for main loop
    DO 99 I = 1,N
      D(I) = RINF
      L(I) = 0
   99 CONTINUE

C Main loop ... each pass round this loop is similar to Dijkstra's
C algorithm for solving the single source shortest path problem

    DO 100 JORD = 1,N

      IF (JPERM(JORD).NE.0) GO TO 100
C JORD is next unmatched column
C DMIN is the length of shortest path in the tree
      DMIN = RINF
      QLEN = 0
      LOW = N + 1
      UP = N + 1
C CSP is the cost of the shortest augmenting path to unassigned row
C IRN(ISP). The corresponding column index is JSP.
      CSP = RINF
C Build shortest path tree starting from unassigned column (root) JORD
      J = JORD
      PR(J) = -1

C Scan column J
      DO 115 K = IP(J),IP(J+1)-1
      I = IRN(K)
      DNEW = A(K) - U(I)
      IF (DNEW.GE.CSP) GO TO 115
      IF (IPERM(I).EQ.0) THEN
        CSP = DNEW
        ISP = K
        JSP = J
      ELSE
        IF (DNEW.LT.DMIN) DMIN = DNEW
        D(I) = DNEW
        QLEN = QLEN + 1
        Q(QLEN) = K
      ENDIF
  115   CONTINUE
C Initialize heap Q and Q2 with rows held in Q(1:QLEN)
      Q0 = QLEN
      QLEN = 0
      DO 120 KK = 1,Q0
      K = Q(KK)
      I = IRN(K)
      IF (CSP.LE.D(I)) THEN
        D(I) = RINF
        GO TO 120
      ENDIF
      IF (D(I).LE.DMIN) THEN
        LOW = LOW - 1
        Q(LOW) = I
        L(I) = LOW
      ELSE
        QLEN = QLEN + 1
        L(I) = QLEN
        CALL MC64DD(I,N,Q,D,L,2)
      ENDIF
C Update tree
      JJ = IPERM(I)
      OUT(JJ) = K
      PR(JJ) = J
  120   CONTINUE

      DO 150 JDUM = 1,NUM

C If Q2 is empty, extract rows from Q
      IF (LOW.EQ.UP) THEN
        IF (QLEN.EQ.0) GO TO 160
        I = Q(1)
        IF (D(I).GE.CSP) GO TO 160
        DMIN = D(I)
  152     CALL MC64ED(QLEN,N,Q,D,L,2)
        LOW = LOW - 1
        Q(LOW) = I
        L(I) = LOW
        IF (QLEN.EQ.0) GO TO 153
        I = Q(1)
        IF (D(I).GT.DMIN) GO TO 153
        GO TO 152
      ENDIF
C Q0 is row whose distance D(Q0) to the root is smallest
  153   Q0 = Q(UP-1)
      DQ0 = D(Q0)
C Exit loop if path to Q0 is longer than the shortest augmenting path
      IF (DQ0.GE.CSP) GO TO 160
      UP = UP - 1

C Scan column that matches with row Q0
      J = IPERM(Q0)
      VJ = DQ0 - A(JPERM(J)) + U(Q0)
      DO 155 K = IP(J),IP(J+1)-1
        I = IRN(K)
        IF (L(I).GE.UP) GO TO 155
C DNEW is new cost
        DNEW = VJ + A(K)-U(I)
C Do not update D(I) if DNEW ge cost of shortest path
        IF (DNEW.GE.CSP) GO TO 155
        IF (IPERM(I).EQ.0) THEN
C Row I is unmatched; update shortest path info
          CSP = DNEW
          ISP = K
          JSP = J
        ELSE
C Row I is matched; do not update D(I) if DNEW is larger
          DI = D(I)
          IF (DI.LE.DNEW) GO TO 155
          IF (L(I).GE.LOW) GO TO 155
          D(I) = DNEW
          IF (DNEW.LE.DMIN) THEN
          LPOS = L(I)
          IF (LPOS.NE.0)
   *        CALL MC64FD(LPOS,QLEN,N,Q,D,L,2)
          LOW = LOW - 1
          Q(LOW) = I
          L(I) = LOW
          ELSE
          IF (L(I).EQ.0) THEN
            QLEN = QLEN + 1
            L(I) = QLEN
          ENDIF
          CALL MC64DD(I,N,Q,D,L,2)
          ENDIF
C Update tree
          JJ = IPERM(I)
          OUT(JJ) = K
          PR(JJ) = J
        ENDIF
  155   CONTINUE
  150   CONTINUE

C If CSP = RINF, no augmenting path is found
  160   IF (CSP.EQ.RINF) GO TO 190
C Find augmenting path by tracing backward in PR; update IPERM,JPERM
      NUM = NUM + 1
      I = IRN(ISP)
      IPERM(I) = JSP
      JPERM(JSP) = ISP
      J = JSP
      DO 170 JDUM = 1,NUM
      JJ = PR(J)
      IF (JJ.EQ.-1) GO TO 180
      K = OUT(J)
      I = IRN(K)
      IPERM(I) = JJ
      JPERM(JJ) = K
      J = JJ
  170   CONTINUE
C End of dummy loop; this point is never reached

C Update U for rows in Q(UP:N)
  180   DO 185 KK = UP,N
      I = Q(KK)
      U(I) = U(I) + D(I) - CSP
  185   CONTINUE
  190   DO 191 KK = LOW,N
      I = Q(KK)
      D(I) = RINF
      L(I) = 0
  191   CONTINUE
      DO 193 KK = 1,QLEN
      I = Q(KK)
      D(I) = RINF
      L(I) = 0
  193   CONTINUE

  100 CONTINUE
C End of main loop


C Set dual column variable in D(1:N)
 1000 DO 200 J = 1,N
      K = JPERM(J)
      IF (K.NE.0) THEN
      D(J) = A(K) - U(IRN(K))
      ELSE
      D(J) = ZERO
      ENDIF
      IF (IPERM(J).EQ.0) U(J) = ZERO
  200 CONTINUE

    IF (NUM.EQ.N) GO TO 1100

C The matrix is structurally singular, complete IPERM.
C JPERM, OUT are work arrays
    DO 300 J = 1,N
      JPERM(J) = 0
  300 CONTINUE
    K = 0
    DO 310 I = 1,N
      IF (IPERM(I).EQ.0) THEN
      K = K + 1
      OUT(K) = I
      ELSE
      J = IPERM(I)
      JPERM(J) = I
      ENDIF
  310 CONTINUE
    K = 0
    DO 320 J = 1,N
      IF (JPERM(J).NE.0) GO TO 320
      K = K + 1
      JDUM = OUT(K)
      IPERM(JDUM) = - J
  320 CONTINUE
 1100 RETURN
    END


* COPYRIGHT (c) 1993 AEA Technology and
* Council for the Central Laboratory of the Research Councils
C Original date 18 June 2001
C  June 2001: threadsafe version of MF41
C
C 12th July 2004 Version 1.0.0. Version numbering added.
C 3 December 2009. Version 1.1.0. CMPLX replaced by DCMPLX in D version.

    SUBROUTINE MF71AD(N,KASE,X,EST,W,KEEP)
C
C    MF71AD ESTIMATES THE 1-NORM OF A SQUARE, DOUBLE PRECISION
C    COMPLEX MATRIX A.
C    REVERSE COMMUNICATION IS USED FOR EVALUATING
C    MATRIX-VECTOR PRODUCTS.
C
C
C       N     INTEGER
C           THE ORDER OF THE MATRIX.  N .GE. 1.
C
C       KASE    INTEGER
C           SET INITIALLY TO ZERO . IF N .LE. 0 SET TO -1
C           ON INTERMEDIATE RETURN
C           = 1 OR 2.
C          ON FINAL RETURN
C           =  0  ,IF SUCCESS
C           = -1  ,IF N .LE.0
C
C       X     COMPLEX*16 ARRAY OF DIMENSION (N)
C           IF 1-NORM IS REQUIRED
C           MUST BE OVERWRITTEN BY
C
C              A*X,         IF KASE=1,
C              CTRANS(A)*X,   IF KASE=2,
C
C           WHERE  CTRANS(A)  IS THE CONJUGATE TRANSPOSE,
C           AND MF71 MUST BE RE-CALLED, WITH ALL THE OTHER
C           PARAMETERS UNCHANGED.
C           IF INFINITY-NORM IS REQUIRED
C           MUST BE OVERWRITTEN BY
C
C              CTRANS(A)*X,   IF KASE=1,
C              A*X,         IF KASE=2,
C
C           AND MF71 MUST BE RE-CALLED, WITH ALL THE OTHER
C           PARAMETERS UNCHANGED.
C
C       EST   DOUBLE PRECISION
C           CONTAINS AN ESTIMATE (A LOWER BOUND) FOR NORM(A).
C
C       W     COMPLEX*16 ARRAY OF DIMENSION (N)
C           = A*V,   WHERE  EST = NORM(W)/NORM(V)
C                (V  IS NOT RETURNED).
C
C       KEEP    INTEGER ARRAY LENGTH 5 USED TO PRESERVE PRIVATE
C           DATA, JUMP, ITER, J AND JLAST BETWEEN CALLS,
C           KEEP(5) IS SPARE.
C
C    REFERENCE
C    N.J. HIGHAM (1987) FORTRAN CODES FOR ESTIMATING
C    THE ONE-NORM OF A
C    REAL OR COMPLEX MATRIX, WITH APPLICATIONS
C    TO CONDITION  ESTIMATION, ACM TRANS MATH SOFTWARE,VOL 14,
C    1988, PP381-396.
C
C   .. Parameters ..
    INTEGER ITMAX
    PARAMETER (ITMAX=5)
    DOUBLE PRECISION ZERO,ONE,TWO
    PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
    COMPLEX*16 CZERO,CONE
    PARAMETER (CZERO= (0.0D0,0.0D0),CONE= (1.0D0,0.0D0))
C   ..
C   .. Scalar Arguments ..
    DOUBLE PRECISION EST
    INTEGER KASE,N
C   ..
C   .. Array Arguments ..
    COMPLEX*16 W(N),X(N)
    INTEGER KEEP(5)
C   ..
C   .. Local Scalars ..
    DOUBLE PRECISION ALTSGN,TEMP,XMAX
    INTEGER I,ITER,J,JLAST,JUMP
C   ..
C   .. Intrinsic Functions ..
    INTRINSIC ABS,DCMPLX,DBLE
C   ..
C
    IF (N.LE.0) THEN
       KASE = -1
       RETURN

    END IF

    IF (KASE.EQ.0) THEN
       DO 10 I = 1,N
        X(I) = DCMPLX(ONE/DBLE(N))
   10    CONTINUE
       KASE = 1
       JUMP = 1
       KEEP(1) = JUMP
       KEEP(2) = 0
       KEEP(3) = 0
       KEEP(4) = 0
       RETURN

    END IF
C
    JUMP  = KEEP(1)
    ITER  = KEEP(2)
    J   = KEEP(3)
    JLAST = KEEP(4)
C
    GO TO (20,40,80,110,150) JUMP
C
C    ................ ENTRY   (JUMP = 1)
C
   20 CONTINUE
    IF (N.EQ.1) THEN
       W(1) = X(1)
       EST = ABS(W(1))
C       ... QUIT
       GO TO 180

    END IF
C
    DO 30 I = 1,N
       IF (ABS(X(I)).EQ.ZERO) THEN
        X(I) = CONE
       ELSE
        X(I) = X(I)/DCMPLX(ABS(X(I)))
       END IF
   30 CONTINUE
    KASE = 2
    JUMP = 2
    GO TO 1010
C
C    ................ ENTRY   (JUMP = 2)
C   FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
C
   40 CONTINUE
C FIND MAX ENTRY IN X
    J = 1
    XMAX = ABS(X(J))
    DO 50 I = 2,N
       IF (ABS(X(I)).LE.XMAX) GO TO 50
       J = I
       XMAX = ABS(X(J))
   50 CONTINUE
    ITER = 2
C
C    MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
C
   60 CONTINUE
    DO 70 I = 1,N
       X(I) = CZERO
   70 CONTINUE
    X(J) = CONE
    KASE = 1
    JUMP = 3
    GO TO 1010
C
C    ................ ENTRY   (JUMP = 3)
C   X HAS BEEN OVERWRITTEN BY A*X.
C
   80 CONTINUE
C
C   COPY X INTO W
    DO 90 I = 1,N
       W(I) = X(I)
   90 CONTINUE
C
    DO 100 I = 1,N
       IF (ABS(X(I)).EQ.ZERO) THEN
        X(I) = CONE
       ELSE
        X(I) = X(I)/DCMPLX(ABS(X(I)))
       END IF
  100 CONTINUE
    KASE = 2
    JUMP = 4
    GO TO 1010
C
C    ................ ENTRY   (JUMP = 4)
C
  110 CONTINUE
    JLAST = J
C FIND MAX ENTRY IN X
    J = 1
    XMAX = ABS(X(J))
    DO 120 I = 2,N
       IF (ABS(X(I)).LE.XMAX) GO TO 120
       J = I
       XMAX = ABS(X(J))
  120 CONTINUE
    IF (ABS(X(JLAST)).NE.ABS(X(J)) .AND.
   +    (ITER.LT.ITMAX)) THEN
       ITER = ITER + 1
       GO TO 60

    END IF
C
C    ITERATION COMPLETE.  FINAL STAGE.
C
    EST = ZERO
    DO 130 I = 1,N
       EST = EST + ABS(W(I))
  130 CONTINUE
C
    ALTSGN = ONE
    DO 140 I = 1,N
       X(I) = DCMPLX(ALTSGN* (ONE+DBLE(I-1)/DBLE(N-1)))
       ALTSGN = -ALTSGN
  140 CONTINUE
    KASE = 1
    JUMP = 5
    GO TO 1010
C
C    ................ ENTRY   (JUMP = 5)
C
  150 CONTINUE
    TEMP = ZERO
    DO 160 I = 1,N
       TEMP = TEMP + ABS(X(I))
  160 CONTINUE
    TEMP = TWO*TEMP/DBLE(3*N)
    IF (TEMP.GT.EST) THEN
C
C    COPY X INTO W
C
       DO 170 I = 1,N
        W(I) = X(I)
  170    CONTINUE
       EST = TEMP
    END IF
C
  180 KASE = 0
C
 1010 CONTINUE
    KEEP(1) = JUMP
    KEEP(2) = ITER
    KEEP(3) = J
    KEEP(4) = JLAST
    RETURN
C
    END
    SUBROUTINE METIS_NODEND_ME57(N,IPTR,IRN,METFTN,METOPT,INVPRM,PERM)
C Dummy routine that is called if MeTiS is not linked.
    INTEGER N
    INTEGER IPTR(N+1),IRN(*),METFTN,METOPT(8),INVPRM(N),PERM(N)
    PERM(1) = -1
    RETURN
    END
