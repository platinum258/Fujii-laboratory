! COPYRIGHT (c) 1976 AEA Technology and
! Council for the Central Laboratory of the Research Councils
!
! Version 1.2.0
! See ChangeLog for version history.
!
      SUBROUTINE EA22ID(ICNTL,KEEP,RKEEP)
      INTEGER ICNTL(10)
      INTEGER KEEP(30)
      DOUBLE PRECISION RKEEP(20)
      INTEGER I

C  initialize optional user controls
      ICNTL(1) = 6
      ICNTL(2) = 0
      ICNTL(3) = 0
      DO 10 I = 4,10
        ICNTL(I) = 0
   10 CONTINUE

C  initialize persistent data to avoid undefined assignment
      DO 20 I = 1,30
        KEEP(I) = 0
   20 CONTINUE
      DO 30 I = 1,20
        RKEEP(I) = 0.0
   30 CONTINUE

C  initialize random number generator
      CALL FA14ID(KEEP(21))

      RETURN
      END
      SUBROUTINE EA22AD(N,NUMEIG,NUMCOL,EPS,X,NX,D,U,W,LW,IW,KMULT,IPOS,
     +                 ICNTL,KEEP,RKEEP)
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS
      INTEGER IPOS,KMULT,LW,N,NUMCOL,NUMEIG,NX
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION D(NUMCOL),U(*),W(LW),X(NX,NUMCOL)
      INTEGER IW(NUMCOL)
      INTEGER ICNTL(10),KEEP(30)
      DOUBLE PRECISION RKEEP(20)
C     ..
C     .. Local Scalars ..
      INTEGER KOUNT,N1,N2,N3,PSQ,LP
C     ..
C     .. External Subroutines ..
      EXTERNAL EA22BD
C     ..
C     .. Executable Statements ..
C
C  Restore persistent values
      KOUNT = KEEP(12)
C
      LP = ICNTL(1)

      IF (IPOS.EQ.1) THEN
        KOUNT = 0
        IF (N.LT.1) THEN
          IPOS = -1
          IF (LP.GE.0) THEN
            WRITE (LP,FMT=991) IPOS
            WRITE (LP,FMT=992) N
          ENDIF
          GO TO 1000

        END IF

        IF (NUMEIG.LT.1) THEN
          IPOS = -2
          IF (LP.GE.0) THEN
            WRITE (LP,FMT=991) IPOS
            WRITE (LP,FMT=993) NUMEIG
          ENDIF
          GO TO 1000

        END IF

      END IF

      KOUNT = KOUNT + 1
      N1 = N + 1
      PSQ = NUMCOL*NUMCOL
      N2 = N1 + PSQ
      N3 = N2 + PSQ
      CALL EA22BD(N,NUMEIG,NUMCOL,EPS,X,NX,D,U,W,W(N1),W(N2),W(N3),IW,
     +            KMULT,IPOS,ICNTL,KEEP,RKEEP)
      IF (IPOS.GT.1 .AND. KOUNT.EQ.KMULT) THEN
        IPOS = -5
        IF (LP.GE.0) THEN
          WRITE (LP,FMT=991) IPOS
          WRITE (LP,FMT=994) KMULT
        ENDIF
      END IF

  991 FORMAT (/,' ERROR RETURN FROM EA22A/AD. IPOS = ',I2)
  992 FORMAT (' VALUE OF N LESS THAN ONE =',I10)
  993 FORMAT (' VALUE OF NUMEIG LESS THAN ONE =',I10)
  994 FORMAT (' NUMBER MATRIX-VECTOR MULTIPLICATIONS MORE',/,
     +       ' THAN SPECIFIED BY KMULT =',I10)
C
C  Save persistent data and return to caller
 1000 CONTINUE
      KEEP(12) = KOUNT
      RETURN

      END
      SUBROUTINE EA22BD(N,G,P,EPS,X,NX,D,U,W,B,Q,V,R,KMULT,IPOS,ICNTL,
     +                 KEEP,RKEEP)
C **  SUBROUTINES CALLED BY EA22B/BD ARE THE FOLLOWING  **
C FA14A/AD(-1) PRODUCES A PSEUDO-RANDOM NUMBER IN THE RANGE (-1,1).
C EA22E/ED ORTHONORMALIZES COLUMN OF X WRT ALL EARLIER COLUMNS.
C EA22F/FD CALCULATES THE EIGENSOLUTION OF A MATRIX USING CLASSICAL
C     JACOBI.
C EA22G/D AND EA22H/D ARE INNER-PRODUCT ROUTINES.
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS
      INTEGER G,IPOS,KMULT,N,NX,P
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION B(P,P),D(P),Q(P,P),U(N),V(N),W(N),X(NX,P)
      INTEGER R(P)
      INTEGER ICNTL(10),KEEP(30)
      DOUBLE PRECISION RKEEP(20)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION LRES,MACHEP,S,S1,T,ZW
      DOUBLE PRECISION E,E1,EMAX,LEIG,RES
      INTEGER I,II,IMG,IMGP1,IROT,J,JMG,KMG,L,LEFT,LMG,PM1
      INTEGER EM,GP1,IFLAG,IMULT,K,KK,KO,KS,LP,M,PMG,Z
C     ..
C     .. External Functions ..
      DOUBLE PRECISION EA22GD,EA22HD,FA14AD
      EXTERNAL EA22GD,EA22HD,FA14AD
C     ..
C     .. External Subroutines ..
      EXTERNAL EA22ED,EA22FD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,LOG,MAX,MIN,INT,SQRT
C     ..
C     .. Executable Statements ..
C
      LP = ICNTL(1)
C
C  Restore persistent values
      Z      = KEEP(1)
      KO     = KEEP(2)
      IMULT  = KEEP(3)
      K      = KEEP(4)
      KK     = KEEP(5)
      M      = KEEP(6)
      EM     = KEEP(7)
      PMG    = KEEP(8)
      GP1    = KEEP(9)
      KS     = KEEP(10)
      IFLAG  = KEEP(11)
      IMG    = KEEP(12)
      IMGP1  = KEEP(13)
      IROT   = KEEP(14)
      JMG    = KEEP(15)
      KMG    = KEEP(16)
      L      = KEEP(17)
      LEFT   = KEEP(18)
      LMG    = KEEP(19)
      PM1    = KEEP(20)
C
      RES    = RKEEP(1)
      EMAX   = RKEEP(2)
      E      = RKEEP(3)
      E1     = RKEEP(4)
      LEIG   = RKEEP(5)
      LRES   = RKEEP(6)
      MACHEP = RKEEP(7)
      S      = RKEEP(8)
      S1     = RKEEP(9)
      T      = RKEEP(10)
      ZW     = RKEEP(11)
C
      IF (IPOS.NE.1) IMULT = IMULT + 1
      GO TO (10,70,210,400) IPOS
C EMAX IS AN ESTIMATE OF ABS(LARGEST EIGENVALUE)
C EM HOLDS THE NUMBER OF EIGENVALUES REQUIRED.
C G HOLDS THE NUMBER OF EIGENVALUES FOUND SO FAR.
C KO HOLDS THE NUMBER OF INNER-PRODUCTS PERFORMED.
C IMULT HOLDS THE NUMBER OF MATRIX BY VECTOR MULTIPLICATIONS.
C KS IS THE NUMBER OF MATRIX BY MATRIX (EIGENVECTOR ESTIMATES)
C     MULTIPLICATIONS.
C Z COUNTS THE NUMBER OF TIMES THE MAIN LOOP HAS BEEN EXECUTED.
C (-E,E) IS THE RANGE USED FOR CHEBYSHEV ACCELERATION.
C IFLAG IS A SWITCH TO INDICATE WHEN THE ESTIMATE OF THE MAXIMUM
C     UNCONVERGED EIGENVALUE CEASES TO BE MONOTONIC.
   10 EMAX = 0.0D0
      EM = G
      G = 0
      KO = 0
      M = 0
      IMULT = 0
      KS = 0
      E = 0.0D0
      IFLAG = 2
      Z = 1
      MACHEP = EPSILON(MACHEP)
C CHECKS RELATIVE SIZES OF THE PARAMETERS EM,P,NX, AND N.
      IF (NX.LT.N) THEN
        IPOS = -3
        IF (LP.GE.0) THEN
          WRITE (LP,FMT=465) IPOS
          WRITE (LP,FMT=470)
        ENDIF
        GO TO 1000

      END IF

      IF (P.GT.N) THEN
        IPOS = -4
        IF (LP.GE.0) THEN
          WRITE (LP,FMT=465) IPOS
          WRITE (LP,FMT=480)
        ENDIF
        GO TO 1000

      END IF

      IF (EM.GT.P .OR. (EM.EQ.P.AND.ICNTL(2).NE.1)) THEN
        IPOS = -4
        IF (LP.GE.0) THEN
          WRITE (LP,FMT=465) IPOS
          IF (EM.EQ.P .AND. ICNTL(2).NE.1) WRITE (LP,FMT=490)
          IF (P.LT.EM) WRITE (LP,FMT=500)
        ENDIF
        GO TO 1000

      END IF

      IF (ICNTL(2).EQ.1) E = D(1)
      IF (ICNTL(3).EQ.1) GO TO 30
      DO 25 J = 1,P
        DO 20 I = 1,N
          X(I,J) = FA14AD(KEEP(21),-1)
   20   CONTINUE
   25 CONTINUE
C
C HERE IS A LOOP EQUIVALENT TO  DO 430 Z=1,KMULT
C     AN ACTUAL DO LOOP CANNOT BE USED BECAUSE OF RETURNS TO THE
C     CALLING PROGRAM.
   30 GP1 = G + 1
      PMG = P - G
C
C *******  JACOBI STEP  ******
C FIRST ORTHONORMALIZE UNCONVERGED COLUMNS OF X WRT EARLIER COLUMNS.
      DO 40 J = GP1,P
        CALL EA22ED(J,N,X,NX,IPOS,KEEP(21))
   40 CONTINUE
      IF (IPOS.EQ.-7) THEN
        IF (LP.GE.0) THEN
          WRITE (LP,FMT=465) IPOS
          WRITE (LP,FMT=540)
        ENDIF
        GO TO 1000

      END IF

      KO = KO + (P*P+P-G-G*G)/2
      IPOS = 2
      K = GP1
C HERE IS A LOOP EQUIVALENT TO  DO 90 K=GP1,P
C THIS LOOP COMPUTES  XTRANSPOSE * A * X  FOR THE P-G UNCONVERGED
C     VECTORS OF X AND PUTS VECTOR  A * X  IN X.
   50 DO 60 I = 1,N
        W(I) = X(I,K)
   60 CONTINUE
C EXITS TO COMPUTE A * X(.,K) RETURNING WITH RESULT IN U.
      GO TO 1000

   70 KMG = K - G
C B = XTRANSPOSE * A * X.
      DO 80 L = K,P
        LMG = L - G
        B(KMG,LMG) = EA22GD(N,U,X(1,L))
   80 CONTINUE
      KO = KO + P - K + 1
C PUTS  A * X(.,K) INTO X(.,K).
      DO 90 I = 1,N
        X(I,K) = U(I)
   90 CONTINUE
      K = K + 1
      IF (K.LE.P) GO TO 50
C
C COMPUTES EIGENVALUES AND EIGENVECTORS OF B.
      IROT = 20*PMG*PMG
      CALL EA22FD(PMG,B,P,D(GP1),Q,.TRUE.,IROT,R)
C
C REORDERS THE EIGENVALUES AND VECTORS IN DESCENDING ORDER.
C COMPARES D(G+2),D(G+1); THEN D(G+3),D(G+2),D(G+1), THEN ...
C     STOPPING WITHIN EACH SET OF COMPARISONS WHENEVER NO INTERCHANGE
C     IS REQUIRED.
      PM1 = P - 1
      IF (G.EQ.PM1) GO TO 130
      DO 120 J = GP1,PM1
        DO 110 II = GP1,J
          I = GP1 + J - II
          IF (ABS(D(I)).GE.ABS(D(I+1))) GO TO 120
C INTERCHANGE PERFORMED.
          S = D(I)
          D(I) = D(I+1)
          D(I+1) = S
          IMG = I - G
          IMGP1 = IMG + 1
          DO 100 L = 1,PMG
            S = Q(L,IMG)
            Q(L,IMG) = Q(L,IMGP1)
            Q(L,IMGP1) = S
  100     CONTINUE
  110   CONTINUE
  120 CONTINUE
C
C POSTMULTIPLY X BY Q AND PUT RESULT IN X.
C COMPUTATION PERFORMED ON X (UNCONVERGED COLUMNS ONLY) A ROW (I)
C     AT A TIME.
  130 DO 170 I = 1,N
        DO 150 J = GP1,P
          JMG = J - G
          V(J) = EA22HD(PMG,X(I,GP1),NX,Q(1,JMG),1)
  150   CONTINUE
        DO 160 J = GP1,P
          X(I,J) = V(J)
  160   CONTINUE
  170 CONTINUE
      KO = KO + PMG*PMG
C
      KS = KS + 1
C
C ******  CHECK THE RESIDUALS  ******
C FIRST SEE IF ESTIMATE OF MAX. EVALUE HAS INCREASED.  IF NOT SET FLAG.
      IF (IFLAG.GT.0) GO TO 174
      IF (ABS(D(GP1)).GT.LEIG) GO TO 176
      IFLAG = 1
  174 IF (IFLAG.EQ.2) IFLAG = 0
  176 LEIG = ABS(D(GP1))
      DO 180 I = 1,N
        V(I) = X(I,GP1)
  180 CONTINUE
      IPOS = 3
      K = GP1
C HERE IS A LOOP EQUIVALENT TO  DO 230 K=GP1,EM
  190 CALL EA22ED(K,N,X,NX,IPOS,KEEP(21))
      IF (IPOS.EQ.-7) THEN
        IF (LP.GE.0) THEN
          WRITE (LP,FMT=465) IPOS
          WRITE (LP,FMT=540)
        ENDIF
        GO TO 1000

      END IF

      KO = KO + K
      DO 200 I = 1,N
        W(I) = X(I,K)
  200 CONTINUE
C EXITS TO COMPUTE  A * X(.,K)  AND RETURNS WITH RESULT IN  U.
      GO TO 1000

  210 S = 0.0D0
      DO 220 I = 1,N
        S = S + (U(I)-D(K)*X(I,K))**2
  220 CONTINUE
      S = SQRT(S)
      IF (S.GT.EPS) GO TO 240
C 230 CONTINUE
      K = K + 1
      IF (K.LE.EM) GO TO 190
  240 LRES = RES
      RES = S
      IF (K.GT.GP1) IFLAG = 2
      IF (IFLAG.NE.1) GO TO 248
      IF (LRES.GT.RES) GO TO 248
      S = DBLE(N)*ABS(D(GP1))*MACHEP
      IF (RES.GT.S) GO TO 248
      K = K + 1
      IFLAG = 2
      EPS = RES
C
C ******  COMPUTE RANGE OF CHEBYCHEV POLYNOMIAL  ******
C WE NORMALLY USE THE PTH EIGENVALUE OF A, BUT IN THE FIRST THREE
C     ITERATIONS WE USE THE SQUARE ROOT OF THE PTH EIGENVALUE OF
C     A**2.  NOTE THAT X CURRENTLY HOLDS A*X1*Q WHERE X1 IS THE
C     VALUE OF X AT THE START OF THE ITERATION AND IS ORTHONORMAL.
C
  248 IF (ICNTL(2).NE.1 .AND. ABS(D(P)).GT.E) E = ABS(D(P))
      IF (K.GT.GP1 .OR. ICNTL(2).EQ.1) GO TO 300
      DO 250 I = 1,N
        X(I,K) = V(I)
  250 CONTINUE
C B = XTRANSPOSE*X.
      DO 270 I = GP1,P
        IMG = I - G
        DO 260 J = I,P
          JMG = J - G
          B(IMG,JMG) = EA22GD(N,X(1,I),X(1,J))
  260   CONTINUE
  270 CONTINUE
      KO = KO + (PMG* (PMG+1))/2
      IROT = 20*PMG*PMG
      CALL EA22FD(PMG,B,P,W,Q,.FALSE.,IROT,R)
      EMAX = W(1)
      S = EMAX
      IF (PMG.LT.2) GO TO 290
      DO 280 I = 2,PMG
        S = MIN(S,W(I))
        EMAX = MAX(EMAX,W(I))
  280 CONTINUE
C EMAX IS AN ESTIMATE OF THE 2-NORM OF A ... I.E. SPECTRAL NORM.
  290 EMAX = SQRT(EMAX)
      IF (S.GT.E*E) E = SQRT(S)
C
C THIS IS THE ONLY POINT IN THE PROGRAM AT WHICH G IS ALTERED.
  300 G = K - 1
      GP1 = G + 1
      PMG = P - G
C HAVE WE OBTAINED ALL THE EIGENVALUES AND EIGENVECTORS REQUESTED...
      IF (G.EQ.EM) GO TO 450
C
C
C ******  CHEBYSHEV ACCELERATION STEP ******
      IF (Z.GE.3) GO TO 320
C ADD IN RANDOM VECTOR AT A LEVEL THAT DOES NOT SIGNIFICANTLY
C     ALTER RESIDUAL.
      S = RES/ (D(GP1)*SQRT(DBLE(N))*2.0D0)
      DO 310 I = 1,N
        ZW = FA14AD(KEEP(21),-1)
        X(I,GP1) = X(I,GP1) + S*ZW
  310 CONTINUE
      CALL EA22ED(GP1,N,X,NX,IPOS,KEEP(21))
      KO = KO + GP1
      IF (IPOS.EQ.-7) THEN
        IF (LP.GE.0) THEN
          WRITE (LP,FMT=465) IPOS
          WRITE (LP,FMT=540)
        ENDIF
        GO TO 1000

      END IF
C
C FIND ORDER FOR CHEBYSHEV POLYNOMIAL.
  320 IF (E.LE.0.0D0) GO TO 440
      S = ABS(D(GP1)/E)
      IF (G.EQ.0 .AND. S.LT.EMAX/E) S = EMAX/E
      IF (S.LE.1.0D0) S = 1.0D0 + MACHEP
      T = LOG(2.0D0*RES/ (EPS*S))/LOG(SQRT(S*S-1.0D0)+S)
      IF (T.LT.0.0D0) T = T - 1.0D0
      LEFT = (KMULT-IMULT)/PMG - 2
      M = MAX(0.0D0,MIN(T+1,KS+3.0D0,DBLE(LEFT)))
C
C COMPUTE REORTHOGONALIZATION PERIODS, WHICH KEEP 2-NORMS OF
C     COMPONENTS OF CONVERGED VECTORS LESS THAN ONE.  L INDICATES
C     THE LAST VECTOR REQUIRING ANY REORTHOGONALIZATION.
      L = 0
      IF (G.EQ.0) GO TO 340
      DO 330 I = 1,G
        S = ABS(D(I)-D(GP1))
        IF (S.LE.0.D0) T = 0.D0
        IF (S.GT.0.0D0) T = LOG(0.5D0*EPS/S)
        S = ABS(D(I)/E)
        IF (S.LE.1.0D0) THEN
          IPOS = -6
          IF (LP.GE.0) THEN
            WRITE (LP,FMT=465) IPOS
            WRITE (LP,FMT=520)
          ENDIF
          GO TO 1000

        END IF

        R(I) = MAX(1,INT(-T/LOG(S+SQRT(S*S-1.0D0))))
        IF (R(I).LT.M) L = I
  330 CONTINUE
  340 IPOS = 4
      IF (M.EQ.0) GO TO 440
C
C FORM TM(A)X.
      KK = GP1
C HERE IS A LOOP EQUIVALENT TO  DO 430 KK=GP1,P
  350 DO 360 I = 1,N
        V(I) = 0.0D0
        W(I) = X(I,KK)
  360 CONTINUE
      E1 = 1.0D0/E
      K = 1
C HERE IS A LOOP EQUIVALENT TO  DO 420 K=1,M
C
C IF R(J) IS A DIVISOR OF K, THEN V AND W ARE ORTHOGONALIZED WRT J TH
C     COLUMN OF X.
  370 IF (L.EQ.0) GO TO 1000
      DO 390 J = 1,L
        IF ((K/R(J))*R(J).NE.K) GO TO 390
        S = EA22GD(N,V,X(1,J))
        S1 = EA22GD(N,W,X(1,J))
        KO = KO + 2
        DO 380 I = 1,N
          W(I) = W(I) - S1*X(I,J)
          V(I) = V(I) - S*X(I,J)
  380   CONTINUE
  390 CONTINUE
C EXITS TO COMPUTE A * W RETURNING WITH RESULT IN U.
      GO TO 1000

  400 DO 410 I = 1,N
        S = W(I)
        W(I) = E1*U(I) - V(I)
        V(I) = S
  410 CONTINUE
C 420 CONTINUE
      E1 = 2.0D0/E
      K = K + 1
      IF (K.LE.M) GO TO 370
      DO 430 I = 1,N
        X(I,KK) = W(I)
  430 CONTINUE
      KK = KK + 1
      IF (KK.LE.P) GO TO 350
      KS = KS + M
C
  440 Z = Z + 1
C ANOTHER TIME ROUND THE MAJOR LOOP.
      GO TO 30
C
C
C REQUIRED EIGENSOLUTION HAS BEEN FOUND.
  450 IPOS = 1
C
C  Save persistent data and return to caller
C
 1000 CONTINUE
      KEEP(1)  = Z
      KEEP(2)  = KO
      KEEP(3)  = IMULT
      KEEP(4)  = K
      KEEP(5)  = KK
      KEEP(6)  = M
      KEEP(7)  = EM
      KEEP(8)  = PMG
      KEEP(9)  = GP1
      KEEP(10) = KS
      KEEP(11) = IFLAG
      KEEP(12) = IMG
      KEEP(13) = IMGP1
      KEEP(14) = IROT
      KEEP(15) = JMG
      KEEP(16) = KMG
      KEEP(17) = L
      KEEP(18) = LEFT
      KEEP(19) = LMG
      KEEP(20) = PM1
C
      RKEEP(1) = RES
      RKEEP(2) = EMAX
      RKEEP(3) = E
      RKEEP(4) = E1
      RKEEP(5) = LEIG
      RKEEP(6) = LRES
      RKEEP(7) = MACHEP
      RKEEP(8) = S
      RKEEP(9) = S1
      RKEEP(10) = T
      RKEEP(11) = ZW
      RETURN
C
C *******   FORMAT STATEMENTS  ******
  465 FORMAT (/,1X,'ERROR RETURN FROM EA22AD. IPOS = ',I2)
  470 FORMAT (1X,'FIRST DIMENSION OF EIGENVECTOR ARRAY IS LESS THAN ',/,
     +       1X,'DIMENSION OF PROBLEM')
  480 FORMAT (1X,'THE NUMBER OF SIMULTANEOUS ITERATES IS GREATER THAN',
     +       /,1X,'THE DIMENSION')
  490 FORMAT (' VALUE SHOULD BE SPECIFIED FOR (NUMCOL+1) TH EIGENVALUE')
  500 FORMAT (1X,'NUMBER OF EIGENVALUES REQUIRED IS GREATER THAN NUMBER'
     +       ,/,1X,'OF SIMULTANEOUS ITERATES')
  520 FORMAT (' POOR CONVERGENCE BECAUSE OF CLOSE EIGENVALUES')
  540 FORMAT (' UNEXPECTED ERROR IN EA22AD(ORTHON)')

      END
C****************************************************************
C ORTHOGONALIZE K TH COLUMN OF X WITH RESPECT TO EARLIER COLUMNS.
C     IF THE RESULT IS A ZERO VECTOR THEN REPLACE IT BY A RANDOM
C     VECTOR AND REPEAT THE ORTHOGONALIZATION. FINALLY NORMALIZE
C     COLUMN.
      SUBROUTINE EA22ED(K,N,X,NX,IPOS,ISEED)
C     .. Parameters ..
      DOUBLE PRECISION FACT
      PARAMETER (FACT=4.D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER IPOS,K,N,NX,ISEED
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(NX,K)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION S,XKN
      INTEGER I,IRAN,J
C     ..
C     .. External Functions ..
      DOUBLE PRECISION EA22GD,FA14AD
      EXTERNAL EA22GD,FA14AD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
C     .. Executable Statements ..

      DO 60 IRAN = 1,10
        XKN = 0.D0
        DO 20 J = 1,K
          S = EA22GD(N,X(1,K),X(1,J))
          IF (J.EQ.K) GO TO 20
          XKN = XKN + S*S
          DO 10 I = 1,N
            X(I,K) = X(I,K) - S*X(I,J)
   10     CONTINUE
   20   CONTINUE
        IF (S.EQ.0.0D0) GO TO 40
        IF (S*FACT.LE.XKN) GO TO 60
        S = SQRT(S)
        DO 30 I = 1,N
          X(I,K) = X(I,K)/S
   30   CONTINUE
        GO TO 70

   40   DO 50 I = 1,N
          X(I,K) = FA14AD(ISEED,-1)
   50   CONTINUE
   60 CONTINUE
      IPOS = -7
   70 RETURN

      END
      SUBROUTINE EA22FD(N,A,ND,D,V,EIGVEC,IROT,IW)
C THIS SUBROUTINE COMPUTES THE EIGENSYSTEM OF MATRIX A USING THE
C     CLASSICAL JACOBI METHOD.  IT IS IDENTICAL TO SUBROUTINE EA13AD
C     IN THE HARWELL SUBROUTINE LIBRARY EXCEPT THAT THE ERROR
C     RETURN IS SUPPRESSED.
C IW(I) HOLDS THE ROW INDEX OF THE LARGEST ELEMENT ABOVE THE
C     DIAGONAL IN COLUMN I.
C     .. Scalar Arguments ..
      INTEGER IROT,N,ND
      LOGICAL EIGVEC
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(ND,N),D(N),V(ND,N)
      INTEGER IW(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION BMAX,C,DIFF,G,H,HIGH,S,T,THETA
      INTEGER I,IR,IWI,IWJ,J,KROT,M,P,PM1,PP1,Q,QM1,QP1
C     ..
C     .. External Functions ..
      INTEGER EA22JD
      EXTERNAL EA22JD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
C     ..
C     .. Executable Statements ..
C CHECK INPUT PARAMETERS
      IF (.NOT.EIGVEC) GO TO 30
      DO 20 J = 1,N
        DO 10 I = 1,N
          V(I,J) = 0.0
   10   CONTINUE
        V(J,J) = 1.0
   20 CONTINUE
C INITIALIZATION OF IW AND D.
   30 KROT = 1
      D(1) = A(1,1)
      IF (N.EQ.1) GO TO 220
      DO 40 J = 2,N
        IW(J) = EA22JD(J-1,A(1,J))
        D(J) = A(J,J)
   40 CONTINUE
      DO 200 KROT = 1,IROT
C FIND THE LARGEST OFF-DIAGONAL ELEMENT.
        P = 1
        Q = 2
        BMAX = ABS(A(1,2))
        IF (N.EQ.2) GO TO 60
        DO 50 J = 3,N
          M = IW(J)
          HIGH = ABS(A(M,J))
          IF (HIGH.LE.BMAX) GO TO 50
          BMAX = HIGH
          P = M
          Q = J
   50   CONTINUE
   60   IF (BMAX.EQ.0.0D0) GO TO 220
C
C COMPARE THE VALUE OF THE LARGEST OFF-DIAGONAL ELEMENT WITH THE
C     CORRESPONDING DIAGONAL ENTRIES.
        G = BMAX*1.0D02
        T = ABS(D(P))+G
        IF (T.NE.ABS(D(P))) GO TO 70
        T = ABS(D(Q))+G
        IF (T.NE.ABS(D(Q))) GO TO 70
C PUT OFF-DIAGONAL ELEMENT TO ZERO IF IT IS NEGLIGIBLE COMPARED
C     WITH BOTH DIAGONAL ELEMENTS.
        A(P,Q) = 0.0D0
        IW(Q) = EA22JD(Q-1,A(1,Q))
        GO TO 200
C
C CALCULATE THE ROTATIONS.
C COT(2*PHI) = THETA, WHERE PHI IS THE ROTATION ANGLE.
   70   DIFF = D(Q) - D(P)
        IF (G+ABS(DIFF).EQ.ABS(DIFF)) GO TO 80
        THETA = DIFF/ (2.0D0*A(P,Q))
        T = 1.0D0/ (ABS(THETA)+SQRT(1.0D0+THETA**2))
        IF (THETA.LT.0.0D0) T = -T
        GO TO 90

   80   T = A(P,Q)/DIFF
   90   C = 1.0D0/SQRT(1.0D0+T*T)
        S = T*C
        H = T*A(P,Q)
        D(P) = D(P) - H
        D(Q) = D(Q) + H
        A(P,Q) = 0.0D0
        PM1 = P - 1
C
C PERFORM ROTATIONS ON A .. AT THE SAME TIME UPDATING THE IW ARRAY.
        IF (PM1.LT.1) GO TO 110
        DO 100 I = 1,PM1
          G = A(I,P)
          H = A(I,Q)
          A(I,P) = C*G - S*H
          A(I,Q) = S*G + C*H
  100   CONTINUE
  110   PP1 = P + 1
        QM1 = Q - 1
        IF (PP1.GT.QM1) GO TO 140
        DO 130 I = PP1,QM1
          G = A(P,I)
          H = A(I,Q)
          A(P,I) = C*G - S*H
          A(I,Q) = S*G + C*H
          IF (P.EQ.IW(I)) GO TO 120
          IWI = IW(I)
          IF (ABS(A(P,I)).GT.ABS(A(IWI,I))) IW(I) = P
          GO TO 130

  120     IF (ABS(G).GE.ABS(A(P,I))) IW(I) = EA22JD(I-1,A(1,I))
  130   CONTINUE
  140   QP1 = Q + 1
        IF (QP1.GT.N) GO TO 180
        DO 170 J = QP1,N
          G = A(P,J)
          H = A(Q,J)
          A(P,J) = C*G - S*H
          A(Q,J) = S*G + C*H
          IR = P
          IF (ABS(A(P,J)).GE.ABS(A(Q,J))) GO TO 150
          IR = Q
          G = H
  150     IF (IR.EQ.IW(J)) GO TO 160
          IWJ = IW(J)
          IF (ABS(A(IR,J)).GT.ABS(A(IWJ,J))) IW(J) = IR
          GO TO 170

  160     IF (ABS(G).GE.ABS(A(IR,J))) IW(J) = EA22JD(J-1,A(1,J))
  170   CONTINUE
  180   IF (P.NE.1) IW(P) = EA22JD(P-1,A(1,P))
        IW(Q) = EA22JD(Q-1,A(1,Q))
        IF (.NOT.EIGVEC) GO TO 200
C UPDATE EIGENVECTOR ESTIMATES.
        DO 190 I = 1,N
          G = V(I,P)
          H = V(I,Q)
          V(I,P) = C*G - S*H
          V(I,Q) = S*G + C*H
  190   CONTINUE
  200 CONTINUE
C     WRITE(LP,210)
C 210 FORMAT(36H TOO MANY ITERATIONS TAKEN BY EA13AD//)
      RETURN

  220 IROT = KROT - 1
      RETURN

      END
      INTEGER FUNCTION EA22JD(N,V)
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION V(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION BMAX,VI
      INTEGER I,M
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Executable Statements ..
      M = 1
      BMAX = ABS(V(1))
      IF (N.LE.1) GO TO 20
      DO 10 I = 2,N
        VI = ABS(V(I))
        IF (VI.LE.BMAX) GO TO 10
        BMAX = VI
        M = I
   10 CONTINUE
   20 EA22JD = M
      RETURN

      END
      DOUBLE PRECISION FUNCTION EA22GD(N,A,B)
C EA22GD COMPUTES THE INNER PRODUCT OF TWO DOUBLE PRECISION VECTORS.
C     IT IS IDENTICAL TO THE HARWELL SUBROUTINE FM01A/AD.
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),B(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION R1
      INTEGER I
C     ..
C     .. Executable Statements ..
C
C    N   THE LENGTH OF THE VECTORS (IF N<= 0  EA22GD = 0)
C    A   THE FIRST VECTOR
C    B   THE SECOND VECTOR
C    EA22GD  THE RESULT (DOUBLE PRECISION)
C
      R1 = 0.D0
      IF (N.LE.0) GO TO 2
      DO 1 I = 1,N
        R1 = R1 + A(I)*B(I)
    1 CONTINUE
    2 EA22GD = R1
      RETURN

      END
      DOUBLE PRECISION FUNCTION EA22HD(N,A,IA,B,IB)
C EA22HD COMPUTES THE INNER PRODUCT OF TWO DOUBLE PRECISION VECTORS.
C      IT IS IDENTICAL TO THE HARWELL SUBROUTINE FM02A/AD.
C     .. Scalar Arguments ..
      INTEGER IA,IB,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),B(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION R1
      INTEGER I,JA,JB
C     ..
C     .. Executable Statements ..
C
C    N   THE LENGTH OF THE VECTORS (IF N<= 0  EA22HD = 0)
C    A   THE FIRST VECTOR
C    IA  SUBSCRIPT DISPLACEMENT BETWEEN ELEMENTS OF A
C    B   THE SECOND VECTOR
C    IB  SUBSCRIPT DISPLACEMENT BETWEEN ELEMENTS OF B
C    EA22HD  THE RESULT
C
      R1 = 0D0
      IF (N.LE.0) GO TO 2
      JA = 1
      IF (IA.LT.0) JA = 1 - (N-1)*IA
      JB = 1
      IF (IB.LT.0) JB = 1 - (N-1)*IB
      I = 0
    1 I = I + 1
      R1 = R1 + A(JA)*B(JB)
      JA = JA + IA
      JB = JB + IB
      IF (I.LT.N) GO TO 1
    2 EA22HD = R1
      RETURN

      END
! COPYRIGHT (c) 1976 AEA Technology and
! Council for the Central Laboratory of the Research Councils
!
! Version 1.2.0
! See ChangeLog for version history.
!
      SUBROUTINE EA22I(ICNTL,KEEP,RKEEP)
      INTEGER ICNTL(10)
      INTEGER KEEP(30)
      REAL RKEEP(20)
      INTEGER I

C  initialize optional user controls
      ICNTL(1) = 6
      ICNTL(2) = 0
      ICNTL(3) = 0
      DO 10 I = 4,10
        ICNTL(I) = 0
   10 CONTINUE

C  initialize persistent data to avoid undefined assignment
      DO 20 I = 1,30
        KEEP(I) = 0
   20 CONTINUE
      DO 30 I = 1,20
        RKEEP(I) = 0.0
   30 CONTINUE

C  initialize random number generator
      CALL FA14I(KEEP(21))

      RETURN
      END
      SUBROUTINE EA22A(N,NUMEIG,NUMCOL,EPS,X,NX,D,U,W,LW,IW,KMULT,IPOS,
     +                 ICNTL,KEEP,RKEEP)
C     .. Scalar Arguments ..
      REAL EPS
      INTEGER IPOS,KMULT,LW,N,NUMCOL,NUMEIG,NX
C     ..
C     .. Array Arguments ..
      REAL D(NUMCOL),U(*),W(LW),X(NX,NUMCOL)
      INTEGER IW(NUMCOL)
      INTEGER ICNTL(10),KEEP(30)
      REAL RKEEP(20)
C     ..
C     .. Local Scalars ..
      INTEGER KOUNT,N1,N2,N3,PSQ,LP
C     ..
C     .. External Subroutines ..
      EXTERNAL EA22B
C     ..
C     .. Executable Statements ..
C
C  Restore persistent values
      KOUNT = KEEP(12)
C
      LP = ICNTL(1)

      IF (IPOS.EQ.1) THEN
        KOUNT = 0
        IF (N.LT.1) THEN
          IPOS = -1
          IF (LP.GE.0) THEN
            WRITE (LP,FMT=991) IPOS
            WRITE (LP,FMT=992) N
          ENDIF
          GO TO 1000

        END IF

        IF (NUMEIG.LT.1) THEN
          IPOS = -2
          IF (LP.GE.0) THEN
            WRITE (LP,FMT=991) IPOS
            WRITE (LP,FMT=993) NUMEIG
          ENDIF
          GO TO 1000

        END IF

      END IF

      KOUNT = KOUNT + 1
      N1 = N + 1
      PSQ = NUMCOL*NUMCOL
      N2 = N1 + PSQ
      N3 = N2 + PSQ
      CALL EA22B(N,NUMEIG,NUMCOL,EPS,X,NX,D,U,W,W(N1),W(N2),W(N3),IW,
     +           KMULT,IPOS,ICNTL,KEEP,RKEEP)
      IF (IPOS.GT.1 .AND. KOUNT.EQ.KMULT) THEN
        IPOS = -5
        IF (LP.GE.0) THEN
          WRITE (LP,FMT=991) IPOS
          WRITE (LP,FMT=994) KMULT
        ENDIF
      END IF

  991 FORMAT (/,' ERROR RETURN FROM EA22A/AD. IPOS = ',I2)
  992 FORMAT (' VALUE OF N LESS THAN ONE =',I10)
  993 FORMAT (' VALUE OF NUMEIG LESS THAN ONE =',I10)
  994 FORMAT (' NUMBER MATRIX-VECTOR MULTIPLICATIONS MORE',/,
     +       ' THAN SPECIFIED BY KMULT =',I10)
C
C  Save persistent data and return to caller
 1000 CONTINUE
      KEEP(12) = KOUNT
      RETURN

      END
      SUBROUTINE EA22B(N,G,P,EPS,X,NX,D,U,W,B,Q,V,R,KMULT,IPOS,ICNTL,
     +                 KEEP,RKEEP)
C **  SUBROUTINES CALLED BY EA22B/BD ARE THE FOLLOWING  **
C FA14A/AD(-1) PRODUCES A PSEUDO-RANDOM NUMBER IN THE RANGE (-1,1).
C EA22E/ED ORTHONORMALIZES COLUMN OF X WRT ALL EARLIER COLUMNS.
C EA22F/FD CALCULATES THE EIGENSOLUTION OF A MATRIX USING CLASSICAL
C     JACOBI.
C EA22G/D AND EA22H/D ARE INNER-PRODUCT ROUTINES.
C     .. Scalar Arguments ..
      REAL EPS
      INTEGER G,IPOS,KMULT,N,NX,P
C     ..
C     .. Array Arguments ..
      REAL B(P,P),D(P),Q(P,P),U(N),V(N),W(N),X(NX,P)
      INTEGER R(P)
      INTEGER ICNTL(10),KEEP(30)
      REAL RKEEP(20)
C     ..
C     .. Local Scalars ..
      REAL LRES,MACHEP,S,S1,T,ZW
      REAL E,E1,EMAX,LEIG,RES
      INTEGER I,II,IMG,IMGP1,IROT,J,JMG,KMG,L,LEFT,LMG,PM1
      INTEGER EM,GP1,IFLAG,IMULT,K,KK,KO,KS,LP,M,PMG,Z
C     ..
C     .. External Functions ..
      REAL EA22G,EA22H,FA14A
      EXTERNAL EA22G,EA22H,FA14A
C     ..
C     .. External Subroutines ..
      EXTERNAL EA22E,EA22F
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,REAL,LOG,MAX,MIN,INT,SQRT
C     ..
C     .. Executable Statements ..
C
      LP = ICNTL(1)
C
C  Restore persistent values
      Z      = KEEP(1)
      KO     = KEEP(2)
      IMULT  = KEEP(3)
      K      = KEEP(4)
      KK     = KEEP(5)
      M      = KEEP(6)
      EM     = KEEP(7)
      PMG    = KEEP(8)
      GP1    = KEEP(9)
      KS     = KEEP(10)
      IFLAG  = KEEP(11)
      IMG    = KEEP(12)
      IMGP1  = KEEP(13)
      IROT   = KEEP(14)
      JMG    = KEEP(15)
      KMG    = KEEP(16)
      L      = KEEP(17)
      LEFT   = KEEP(18)
      LMG    = KEEP(19)
      PM1    = KEEP(20)
C
      RES    = RKEEP(1)
      EMAX   = RKEEP(2)
      E      = RKEEP(3)
      E1     = RKEEP(4)
      LEIG   = RKEEP(5)
      LRES   = RKEEP(6)
      MACHEP = RKEEP(7)
      S      = RKEEP(8)
      S1     = RKEEP(9)
      T      = RKEEP(10)
      ZW     = RKEEP(11)
C
      IF (IPOS.NE.1) IMULT = IMULT + 1
      GO TO (10,70,210,400) IPOS
C EMAX IS AN ESTIMATE OF ABS(LARGEST EIGENVALUE)
C EM HOLDS THE NUMBER OF EIGENVALUES REQUIRED.
C G HOLDS THE NUMBER OF EIGENVALUES FOUND SO FAR.
C KO HOLDS THE NUMBER OF INNER-PRODUCTS PERFORMED.
C IMULT HOLDS THE NUMBER OF MATRIX BY VECTOR MULTIPLICATIONS.
C KS IS THE NUMBER OF MATRIX BY MATRIX (EIGENVECTOR ESTIMATES)
C     MULTIPLICATIONS.
C Z COUNTS THE NUMBER OF TIMES THE MAIN LOOP HAS BEEN EXECUTED.
C (-E,E) IS THE RANGE USED FOR CHEBYSHEV ACCELERATION.
C IFLAG IS A SWITCH TO INDICATE WHEN THE ESTIMATE OF THE MAXIMUM
C     UNCONVERGED EIGENVALUE CEASES TO BE MONOTONIC.
   10 EMAX = 0.0
      EM = G
      G = 0
      KO = 0
      M = 0
      IMULT = 0
      KS = 0
      E = 0.0
      IFLAG = 2
      Z = 1
      MACHEP = EPSILON(MACHEP)
C CHECKS RELATIVE SIZES OF THE PARAMETERS EM,P,NX, AND N.
      IF (NX.LT.N) THEN
        IPOS = -3
        IF (LP.GE.0) THEN
          WRITE (LP,FMT=465) IPOS
          WRITE (LP,FMT=470)
        ENDIF
        GO TO 1000

      END IF

      IF (P.GT.N) THEN
        IPOS = -4
        IF (LP.GE.0) THEN
          WRITE (LP,FMT=465) IPOS
          WRITE (LP,FMT=480)
        ENDIF
        GO TO 1000

      END IF

      IF (EM.GT.P .OR. (EM.EQ.P.AND.ICNTL(2).NE.1)) THEN
        IPOS = -4
        IF (LP.GE.0) THEN
          WRITE (LP,FMT=465) IPOS
          IF (EM.EQ.P .AND. ICNTL(2).NE.1) WRITE (LP,FMT=490)
          IF (P.LT.EM) WRITE (LP,FMT=500)
        ENDIF
        GO TO 1000

      END IF

      IF (ICNTL(2).EQ.1) E = D(1)
      IF (ICNTL(3).EQ.1) GO TO 30
      DO 25 J = 1,P
        DO 20 I = 1,N
          X(I,J) = FA14A(KEEP(21),-1)
   20   CONTINUE
   25 CONTINUE
C
C HERE IS A LOOP EQUIVALENT TO  DO 430 Z=1,KMULT
C     AN ACTUAL DO LOOP CANNOT BE USED BECAUSE OF RETURNS TO THE
C     CALLING PROGRAM.
   30 GP1 = G + 1
      PMG = P - G
C
C *******  JACOBI STEP  ******
C FIRST ORTHONORMALIZE UNCONVERGED COLUMNS OF X WRT EARLIER COLUMNS.
      DO 40 J = GP1,P
        CALL EA22E(J,N,X,NX,IPOS,KEEP(21))
   40 CONTINUE
      IF (IPOS.EQ.-7) THEN
        IF (LP.GE.0) THEN
          WRITE (LP,FMT=465) IPOS
          WRITE (LP,FMT=540)
        ENDIF
        GO TO 1000

      END IF

      KO = KO + (P*P+P-G-G*G)/2
      IPOS = 2
      K = GP1
C HERE IS A LOOP EQUIVALENT TO  DO 90 K=GP1,P
C THIS LOOP COMPUTES  XTRANSPOSE * A * X  FOR THE P-G UNCONVERGED
C     VECTORS OF X AND PUTS VECTOR  A * X  IN X.
   50 DO 60 I = 1,N
        W(I) = X(I,K)
   60 CONTINUE
C EXITS TO COMPUTE A * X(.,K) RETURNING WITH RESULT IN U.
      GO TO 1000

   70 KMG = K - G
C B = XTRANSPOSE * A * X.
      DO 80 L = K,P
        LMG = L - G
        B(KMG,LMG) = EA22G(N,U,X(1,L))
   80 CONTINUE
      KO = KO + P - K + 1
C PUTS  A * X(.,K) INTO X(.,K).
      DO 90 I = 1,N
        X(I,K) = U(I)
   90 CONTINUE
      K = K + 1
      IF (K.LE.P) GO TO 50
C
C COMPUTES EIGENVALUES AND EIGENVECTORS OF B.
      IROT = 20*PMG*PMG
      CALL EA22F(PMG,B,P,D(GP1),Q,.TRUE.,IROT,R)
C
C REORDERS THE EIGENVALUES AND VECTORS IN DESCENDING ORDER.
C COMPARES D(G+2),D(G+1); THEN D(G+3),D(G+2),D(G+1), THEN ...
C     STOPPING WITHIN EACH SET OF COMPARISONS WHENEVER NO INTERCHANGE
C     IS REQUIRED.
      PM1 = P - 1
      IF (G.EQ.PM1) GO TO 130
      DO 120 J = GP1,PM1
        DO 110 II = GP1,J
          I = GP1 + J - II
          IF (ABS(D(I)).GE.ABS(D(I+1))) GO TO 120
C INTERCHANGE PERFORMED.
          S = D(I)
          D(I) = D(I+1)
          D(I+1) = S
          IMG = I - G
          IMGP1 = IMG + 1
          DO 100 L = 1,PMG
            S = Q(L,IMG)
            Q(L,IMG) = Q(L,IMGP1)
            Q(L,IMGP1) = S
  100     CONTINUE
  110   CONTINUE
  120 CONTINUE
C
C POSTMULTIPLY X BY Q AND PUT RESULT IN X.
C COMPUTATION PERFORMED ON X (UNCONVERGED COLUMNS ONLY) A ROW (I)
C     AT A TIME.
  130 DO 170 I = 1,N
        DO 150 J = GP1,P
          JMG = J - G
          V(J) = EA22H(PMG,X(I,GP1),NX,Q(1,JMG),1)
  150   CONTINUE
        DO 160 J = GP1,P
          X(I,J) = V(J)
  160   CONTINUE
  170 CONTINUE
      KO = KO + PMG*PMG
C
      KS = KS + 1
C
C ******  CHECK THE RESIDUALS  ******
C FIRST SEE IF ESTIMATE OF MAX. EVALUE HAS INCREASED.  IF NOT SET FLAG.
      IF (IFLAG.GT.0) GO TO 174
      IF (ABS(D(GP1)).GT.LEIG) GO TO 176
      IFLAG = 1
  174 IF (IFLAG.EQ.2) IFLAG = 0
  176 LEIG = ABS(D(GP1))
      DO 180 I = 1,N
        V(I) = X(I,GP1)
  180 CONTINUE
      IPOS = 3
      K = GP1
C HERE IS A LOOP EQUIVALENT TO  DO 230 K=GP1,EM
  190 CALL EA22E(K,N,X,NX,IPOS,KEEP(21))
      IF (IPOS.EQ.-7) THEN
        IF (LP.GE.0) THEN
          WRITE (LP,FMT=465) IPOS
          WRITE (LP,FMT=540)
        ENDIF
        GO TO 1000

      END IF

      KO = KO + K
      DO 200 I = 1,N
        W(I) = X(I,K)
  200 CONTINUE
C EXITS TO COMPUTE  A * X(.,K)  AND RETURNS WITH RESULT IN  U.
      GO TO 1000

  210 S = 0.0
      DO 220 I = 1,N
        S = S + (U(I)-D(K)*X(I,K))**2
  220 CONTINUE
      S = SQRT(S)
      IF (S.GT.EPS) GO TO 240
C 230 CONTINUE
      K = K + 1
      IF (K.LE.EM) GO TO 190
  240 LRES = RES
      RES = S
      IF (K.GT.GP1) IFLAG = 2
      IF (IFLAG.NE.1) GO TO 248
      IF (LRES.GT.RES) GO TO 248
      S = REAL(N)*ABS(D(GP1))*MACHEP
      IF (RES.GT.S) GO TO 248
      K = K + 1
      IFLAG = 2
      EPS = RES
C
C ******  COMPUTE RANGE OF CHEBYCHEV POLYNOMIAL  ******
C WE NORMALLY USE THE PTH EIGENVALUE OF A, BUT IN THE FIRST THREE
C     ITERATIONS WE USE THE SQUARE ROOT OF THE PTH EIGENVALUE OF
C     A**2.  NOTE THAT X CURRENTLY HOLDS A*X1*Q WHERE X1 IS THE
C     VALUE OF X AT THE START OF THE ITERATION AND IS ORTHONORMAL.
C
  248 IF (ICNTL(2).NE.1 .AND. ABS(D(P)).GT.E) E = ABS(D(P))
      IF (K.GT.GP1 .OR. ICNTL(2).EQ.1) GO TO 300
      DO 250 I = 1,N
        X(I,K) = V(I)
  250 CONTINUE
C B = XTRANSPOSE*X.
      DO 270 I = GP1,P
        IMG = I - G
        DO 260 J = I,P
          JMG = J - G
          B(IMG,JMG) = EA22G(N,X(1,I),X(1,J))
  260   CONTINUE
  270 CONTINUE
      KO = KO + (PMG* (PMG+1))/2
      IROT = 20*PMG*PMG
      CALL EA22F(PMG,B,P,W,Q,.FALSE.,IROT,R)
      EMAX = W(1)
      S = EMAX
      IF (PMG.LT.2) GO TO 290
      DO 280 I = 2,PMG
        S = MIN(S,W(I))
        EMAX = MAX(EMAX,W(I))
  280 CONTINUE
C EMAX IS AN ESTIMATE OF THE 2-NORM OF A ... I.E. SPECTRAL NORM.
  290 EMAX = SQRT(EMAX)
      IF (S.GT.E*E) E = SQRT(S)
C
C THIS IS THE ONLY POINT IN THE PROGRAM AT WHICH G IS ALTERED.
  300 G = K - 1
      GP1 = G + 1
      PMG = P - G
C HAVE WE OBTAINED ALL THE EIGENVALUES AND EIGENVECTORS REQUESTED...
      IF (G.EQ.EM) GO TO 450
C
C
C ******  CHEBYSHEV ACCELERATION STEP ******
      IF (Z.GE.3) GO TO 320
C ADD IN RANDOM VECTOR AT A LEVEL THAT DOES NOT SIGNIFICANTLY
C     ALTER RESIDUAL.
      S = RES/ (D(GP1)*SQRT(REAL(N))*2.0)
      DO 310 I = 1,N
        ZW = FA14A(KEEP(21),-1)
        X(I,GP1) = X(I,GP1) + S*ZW
  310 CONTINUE
      CALL EA22E(GP1,N,X,NX,IPOS,KEEP(21))
      KO = KO + GP1
      IF (IPOS.EQ.-7) THEN
        IF (LP.GE.0) THEN
          WRITE (LP,FMT=465) IPOS
          WRITE (LP,FMT=540)
        ENDIF
        GO TO 1000

      END IF
C
C FIND ORDER FOR CHEBYSHEV POLYNOMIAL.
  320 IF (E.LE.0.0) GO TO 440
      S = ABS(D(GP1)/E)
      IF (G.EQ.0 .AND. S.LT.EMAX/E) S = EMAX/E
      IF (S.LE.1.0) S = 1.0 + MACHEP
      T = LOG(2.0*RES/ (EPS*S))/LOG(SQRT(S*S-1.0)+S)
      IF (T.LT.0.0) T = T - 1.0
      LEFT = (KMULT-IMULT)/PMG - 2
      M = MAX(0.,MIN(T+1,KS+3.,REAL(LEFT)))
C
C COMPUTE REORTHOGONALIZATION PERIODS, WHICH KEEP 2-NORMS OF
C     COMPONENTS OF CONVERGED VECTORS LESS THAN ONE.  L INDICATES
C     THE LAST VECTOR REQUIRING ANY REORTHOGONALIZATION.
      L = 0
      IF (G.EQ.0) GO TO 340
      DO 330 I = 1,G
        S = ABS(D(I)-D(GP1))
        IF (S.LE.0.0) T = 0.0
        IF (S.GT.0.0) T = LOG(0.5*EPS/S)
        S = ABS(D(I)/E)
        IF (S.LE.1.0) THEN
          IPOS = -6
          IF (LP.GE.0) THEN
            WRITE (LP,FMT=465) IPOS
            WRITE (LP,FMT=520)
          ENDIF
          GO TO 1000

        END IF

        R(I) = MAX(1,INT(-T/LOG(S+SQRT(S*S-1.0))))
        IF (R(I).LT.M) L = I
  330 CONTINUE
  340 IPOS = 4
      IF (M.EQ.0) GO TO 440
C
C FORM TM(A)X.
      KK = GP1
C HERE IS A LOOP EQUIVALENT TO  DO 430 KK=GP1,P
  350 DO 360 I = 1,N
        V(I) = 0.0
        W(I) = X(I,KK)
  360 CONTINUE
      E1 = 1.0/E
      K = 1
C HERE IS A LOOP EQUIVALENT TO  DO 420 K=1,M
C
C IF R(J) IS A DIVISOR OF K, THEN V AND W ARE ORTHOGONALIZED WRT J TH
C     COLUMN OF X.
  370 IF (L.EQ.0) GO TO 1000
      DO 390 J = 1,L
        IF ((K/R(J))*R(J).NE.K) GO TO 390
        S = EA22G(N,V,X(1,J))
        S1 = EA22G(N,W,X(1,J))
        KO = KO + 2
        DO 380 I = 1,N
          W(I) = W(I) - S1*X(I,J)
          V(I) = V(I) - S*X(I,J)
  380   CONTINUE
  390 CONTINUE
C EXITS TO COMPUTE A * W RETURNING WITH RESULT IN U.
      GO TO 1000

  400 DO 410 I = 1,N
        S = W(I)
        W(I) = E1*U(I) - V(I)
        V(I) = S
  410 CONTINUE
C 420 CONTINUE
      E1 = 2.0/E
      K = K + 1
      IF (K.LE.M) GO TO 370
      DO 430 I = 1,N
        X(I,KK) = W(I)
  430 CONTINUE
      KK = KK + 1
      IF (KK.LE.P) GO TO 350
      KS = KS + M
C
  440 Z = Z + 1
C ANOTHER TIME ROUND THE MAJOR LOOP.
      GO TO 30
C
C
C REQUIRED EIGENSOLUTION HAS BEEN FOUND.
  450 IPOS = 1
C
C  Save persistent data and return to caller
C
 1000 CONTINUE
      KEEP(1)  = Z
      KEEP(2)  = KO
      KEEP(3)  = IMULT
      KEEP(4)  = K
      KEEP(5)  = KK
      KEEP(6)  = M
      KEEP(7)  = EM
      KEEP(8)  = PMG
      KEEP(9)  = GP1
      KEEP(10) = KS
      KEEP(11) = IFLAG
      KEEP(12) = IMG
      KEEP(13) = IMGP1
      KEEP(14) = IROT
      KEEP(15) = JMG
      KEEP(16) = KMG
      KEEP(17) = L
      KEEP(18) = LEFT
      KEEP(19) = LMG
      KEEP(20) = PM1
C
      RKEEP(1) = RES
      RKEEP(2) = EMAX
      RKEEP(3) = E
      RKEEP(4) = E1
      RKEEP(5) = LEIG
      RKEEP(6) = LRES
      RKEEP(7) = MACHEP
      RKEEP(8) = S
      RKEEP(9) = S1
      RKEEP(10) = T
      RKEEP(11) = ZW
      RETURN
C
C *******   FORMAT STATEMENTS  ******
  465 FORMAT (/,1X,'ERROR RETURN FROM EA22A. IPOS = ',I2)
  470 FORMAT (1X,'FIRST DIMENSION OF EIGENVECTOR ARRAY IS LESS THAN ',/,
     +       1X,'DIMENSION OF PROBLEM')
  480 FORMAT (1X,'THE NUMBER OF SIMULTANEOUS ITERATES IS GREATER THAN',
     +       /,1X,'THE DIMENSION')
  490 FORMAT (' VALUE SHOULD BE SPECIFIED FOR (NUMCOL+1) TH EIGENVALUE')
  500 FORMAT (1X,'NUMBER OF EIGENVALUES REQUIRED IS GREATER THAN NUMBER'
     +       ,/,1X,'OF SIMULTANEOUS ITERATES')
  520 FORMAT (' POOR CONVERGENCE BECAUSE OF CLOSE EIGENVALUES')
  540 FORMAT (' UNEXPECTED ERROR IN EA22A(ORTHON)')

      END
C****************************************************************
C ORTHOGONALIZE K TH COLUMN OF X WITH RESPECT TO EARLIER COLUMNS.
C     IF THE RESULT IS A ZERO VECTOR THEN REPLACE IT BY A RANDOM
C     VECTOR AND REPEAT THE ORTHOGONALIZATION. FINALLY NORMALIZE
C     COLUMN.
      SUBROUTINE EA22E(K,N,X,NX,IPOS,ISEED)
C     .. Parameters ..
      REAL FACT
      PARAMETER (FACT=4.0)
C     ..
C     .. Scalar Arguments ..
      INTEGER IPOS,K,N,NX,ISEED
C     ..
C     .. Array Arguments ..
      REAL X(NX,K)
C     ..
C     .. Local Scalars ..
      REAL S,XKN
      INTEGER I,IRAN,J
C     ..
C     .. External Functions ..
      REAL EA22G,FA14A
      EXTERNAL EA22G,FA14A
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
C     .. Executable Statements ..

      DO 60 IRAN = 1,10
        XKN = 0.0
        DO 20 J = 1,K
          S = EA22G(N,X(1,K),X(1,J))
          IF (J.EQ.K) GO TO 20
          XKN = XKN + S*S
          DO 10 I = 1,N
            X(I,K) = X(I,K) - S*X(I,J)
   10     CONTINUE
   20   CONTINUE
        IF (S.EQ.0.0) GO TO 40
        IF (S*FACT.LE.XKN) GO TO 60
        S = SQRT(S)
        DO 30 I = 1,N
          X(I,K) = X(I,K)/S
   30   CONTINUE
        GO TO 70

   40   DO 50 I = 1,N
          X(I,K) = FA14A(ISEED,-1)
   50   CONTINUE
   60 CONTINUE
      IPOS = -7
   70 RETURN

      END
      SUBROUTINE EA22F(N,A,ND,D,V,EIGVEC,IROT,IW)
C THIS SUBROUTINE COMPUTES THE EIGENSYSTEM OF MATRIX A USING THE
C     CLASSICAL JACOBI METHOD.  IT IS IDENTICAL TO SUBROUTINE EA13AD
C     IN THE HARWELL SUBROUTINE LIBRARY EXCEPT THAT THE ERROR
C     RETURN IS SUPPRESSED.
C IW(I) HOLDS THE ROW INDEX OF THE LARGEST ELEMENT ABOVE THE
C     DIAGONAL IN COLUMN I.
C     .. Scalar Arguments ..
      INTEGER IROT,N,ND
      LOGICAL EIGVEC
C     ..
C     .. Array Arguments ..
      REAL A(ND,N),D(N),V(ND,N)
      INTEGER IW(N)
C     ..
C     .. Local Scalars ..
      REAL BMAX,C,DIFF,G,H,HIGH,S,T,THETA
      INTEGER I,IR,IWI,IWJ,J,KROT,M,P,PM1,PP1,Q,QM1,QP1
C     ..
C     .. External Functions ..
      INTEGER EA22J
      EXTERNAL EA22J
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
C     ..
C     .. Executable Statements ..
C CHECK INPUT PARAMETERS
      IF (.NOT.EIGVEC) GO TO 30
      DO 20 J = 1,N
        DO 10 I = 1,N
          V(I,J) = 0.0
   10   CONTINUE
        V(J,J) = 1.0
   20 CONTINUE
C INITIALIZATION OF IW AND D.
   30 KROT = 1
      D(1) = A(1,1)
      IF (N.EQ.1) GO TO 220
      DO 40 J = 2,N
        IW(J) = EA22J(J-1,A(1,J))
        D(J) = A(J,J)
   40 CONTINUE
      DO 200 KROT = 1,IROT
C FIND THE LARGEST OFF-DIAGONAL ELEMENT.
        P = 1
        Q = 2
        BMAX = ABS(A(1,2))
        IF (N.EQ.2) GO TO 60
        DO 50 J = 3,N
          M = IW(J)
          HIGH = ABS(A(M,J))
          IF (HIGH.LE.BMAX) GO TO 50
          BMAX = HIGH
          P = M
          Q = J
   50   CONTINUE
   60   IF (BMAX.EQ.0.0) GO TO 220
C
C COMPARE THE VALUE OF THE LARGEST OFF-DIAGONAL ELEMENT WITH THE
C     CORRESPONDING DIAGONAL ENTRIES.
        G = BMAX*1.0E02
        T = ABS(D(P))+G
        IF (T.NE.ABS(D(P))) GO TO 70
        T = ABS(D(Q))+G
        IF (T.NE.ABS(D(Q))) GO TO 70
C PUT OFF-DIAGONAL ELEMENT TO ZERO IF IT IS NEGLIGIBLE COMPARED
C     WITH BOTH DIAGONAL ELEMENTS.
        A(P,Q) = 0.0
        IW(Q) = EA22J(Q-1,A(1,Q))
        GO TO 200
C
C CALCULATE THE ROTATIONS.
C COT(2*PHI) = THETA, WHERE PHI IS THE ROTATION ANGLE.
   70   DIFF = D(Q) - D(P)
        IF (G+ABS(DIFF).EQ.ABS(DIFF)) GO TO 80
        THETA = DIFF/ (2.0*A(P,Q))
        T = 1.0/ (ABS(THETA)+SQRT(1.0+THETA**2))
        IF (THETA.LT.0.0) T = -T
        GO TO 90

   80   T = A(P,Q)/DIFF
   90   C = 1.0/SQRT(1.0+T*T)
        S = T*C
        H = T*A(P,Q)
        D(P) = D(P) - H
        D(Q) = D(Q) + H
        A(P,Q) = 0.0
        PM1 = P - 1
C
C PERFORM ROTATIONS ON A .. AT THE SAME TIME UPDATING THE IW ARRAY.
        IF (PM1.LT.1) GO TO 110
        DO 100 I = 1,PM1
          G = A(I,P)
          H = A(I,Q)
          A(I,P) = C*G - S*H
          A(I,Q) = S*G + C*H
  100   CONTINUE
  110   PP1 = P + 1
        QM1 = Q - 1
        IF (PP1.GT.QM1) GO TO 140
        DO 130 I = PP1,QM1
          G = A(P,I)
          H = A(I,Q)
          A(P,I) = C*G - S*H
          A(I,Q) = S*G + C*H
          IF (P.EQ.IW(I)) GO TO 120
          IWI = IW(I)
          IF (ABS(A(P,I)).GT.ABS(A(IWI,I))) IW(I) = P
          GO TO 130

  120     IF (ABS(G).GE.ABS(A(P,I))) IW(I) = EA22J(I-1,A(1,I))
  130   CONTINUE
  140   QP1 = Q + 1
        IF (QP1.GT.N) GO TO 180
        DO 170 J = QP1,N
          G = A(P,J)
          H = A(Q,J)
          A(P,J) = C*G - S*H
          A(Q,J) = S*G + C*H
          IR = P
          IF (ABS(A(P,J)).GE.ABS(A(Q,J))) GO TO 150
          IR = Q
          G = H
  150     IF (IR.EQ.IW(J)) GO TO 160
          IWJ = IW(J)
          IF (ABS(A(IR,J)).GT.ABS(A(IWJ,J))) IW(J) = IR
          GO TO 170

  160     IF (ABS(G).GE.ABS(A(IR,J))) IW(J) = EA22J(J-1,A(1,J))
  170   CONTINUE
  180   IF (P.NE.1) IW(P) = EA22J(P-1,A(1,P))
        IW(Q) = EA22J(Q-1,A(1,Q))
        IF (.NOT.EIGVEC) GO TO 200
C UPDATE EIGENVECTOR ESTIMATES.
        DO 190 I = 1,N
          G = V(I,P)
          H = V(I,Q)
          V(I,P) = C*G - S*H
          V(I,Q) = S*G + C*H
  190   CONTINUE
  200 CONTINUE
C     WRITE(LP,210)
C 210 FORMAT(36H TOO MANY ITERATIONS TAKEN BY EA13AD//)
      RETURN

  220 IROT = KROT - 1
      RETURN

      END
      INTEGER FUNCTION EA22J(N,V)
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      REAL V(N)
C     ..
C     .. Local Scalars ..
      REAL BMAX,VI
      INTEGER I,M
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Executable Statements ..
      M = 1
      BMAX = ABS(V(1))
      IF (N.LE.1) GO TO 20
      DO 10 I = 2,N
        VI = ABS(V(I))
        IF (VI.LE.BMAX) GO TO 10
        BMAX = VI
        M = I
   10 CONTINUE
   20 EA22J = M
      RETURN

      END
      REAL FUNCTION EA22G(N,A,B)
C EA22G COMPUTES THE INNER PRODUCT OF TWO REAL VECTORS.
C     IT IS IDENTICAL TO THE HARWELL SUBROUTINE FM01A/AD.
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      REAL A(*),B(*)
C     ..
C     .. Local Scalars ..
      REAL R1
      INTEGER I
C     ..
C     .. Executable Statements ..
C
C    N   THE LENGTH OF THE VECTORS (IF N<= 0  EA22G = 0)
C    A   THE FIRST VECTOR
C    B   THE SECOND VECTOR
C    EA22G   THE RESULT (REAL)
C
      R1 = 0.0
      IF (N.LE.0) GO TO 2
      DO 1 I = 1,N
        R1 = R1 + A(I)*B(I)
    1 CONTINUE
    2 EA22G = R1
      RETURN

      END
      REAL FUNCTION EA22H(N,A,IA,B,IB)
C EA22H COMPUTES THE INNER PRODUCT OF TWO REAL VECTORS.
C      IT IS IDENTICAL TO THE HARWELL SUBROUTINE FM02A/AD.
C     .. Scalar Arguments ..
      INTEGER IA,IB,N
C     ..
C     .. Array Arguments ..
      REAL A(*),B(*)
C     ..
C     .. Local Scalars ..
      REAL R1
      INTEGER I,JA,JB
C     ..
C     .. Executable Statements ..
C
C    N   THE LENGTH OF THE VECTORS (IF N<= 0  EA22H = 0)
C    A   THE FIRST VECTOR
C    IA  SUBSCRIPT DISPLACEMENT BETWEEN ELEMENTS OF A
C    B   THE SECOND VECTOR
C    IB  SUBSCRIPT DISPLACEMENT BETWEEN ELEMENTS OF B
C    EA22H   THE RESULT
C
      R1 = 0.0
      IF (N.LE.0) GO TO 2
      JA = 1
      IF (IA.LT.0) JA = 1 - (N-1)*IA
      JB = 1
      IF (IB.LT.0) JB = 1 - (N-1)*IB
      I = 0
    1 I = I + 1
      R1 = R1 + A(JA)*B(JB)
      JA = JA + IA
      JB = JB + IB
      IF (I.LT.N) GO TO 1
    2 EA22H = R1
      RETURN

      END
! COPYRIGHT (c) 1979 AEA Technology and
! Council for the Central Laboratory of the Research Councils
!
! Version 1.0.1
! See ChangeLog for version history
!
      DOUBLE PRECISION FUNCTION FA14AD(IX,I)
C         NEARLY PORTABLE RANDOM NUMBER GENERATOR USING THE RECURSION
C                       IX=IX*A MOD P
C
C    WHERE A=7**5
C    AND P=2**31-1.
C
C         THIS FUNCTION DOES NOT ADHERE TO THE ANSI STANDARD 1966 IN
C    TWO RESPECTS:
C      1) IT ASSUMES AN INTEGER WORD LENGTH OF AT LEAST 32 BITS (I.E.
C    INTEGERS WHICH LIE IN THE RANGE 1-2**31 TO 2**31-1 INCLUSIVE MUST
C    BE REPRESENTABLE);
C      2) IT ASSUMES THAT A POSITIVE INTEGER LESS THAN 2**16 MAY BE
C    FLOATED WITHOUT LOSS OF DIGITS.
C
C         THIS CODE IS BASED ON CODE PUBLISHED BY LINUS SCHRAGE IN
C    T.O.M.S. VOL.5 NO.2 JUNE 1979 (PP 132-138)
C
C
C
C       THE FUNCTION IS USED AS FOLLOWS:
C
C                      R=FA14AD(IX,I)
C
C       WHERE IX IS THE GENERATOR WORD
C             I IS AN INTEGER SET BY THE USER.
C
C
C       THE VALUE RETURNED BY FA14A/AD WILL LIE IN THE RANGE
C                (0.,1.)  IF I IS NON-NEGATIVE
C                (-1.,1.) IF I IS NEGATIVE.
C
C       THE METHOD EMPLOYED IS A MULTIPLICATIVE CONGRUENTIAL
C   ONE USING A MULTIPLIER OF 7**5 AND TAKING THE MODULO TO
C   2**31-1, I.E. THE GENERATOR NUMBER , G = IX, IS UPDATED ON
C   EACH CALL TO THE VALUE
C
C                  5          31
C               G*7  MODULO (2  -1)
C
C       THE RESULT RETURNED IS CALCULATED AS A DOUBLE
C  PRECISION NUMBER HAVING THE VALUE
C
C                      31
C                  G/(2   -1)    IF THE ARGUMENT IS
C                                NON-NEGATIVE
C           OR
C                      31
C                2*G/(2   -1)-1  IF THE ARGUMENT IS NEGATIVE
C
C
C 7**5, 2**15, 2**16, 2**31-1
C     .. Parameters ..
      INTEGER A,B15,B16,P
      PARAMETER (A=16807,B15=32768,B16=65536,P=2147483647)
C     ..
C     .. Scalar Arguments ..
      INTEGER IX,I
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION X
      INTEGER FHI,K,LEFTLO,XALO,XHI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC FLOAT
C     ..
C     .. Executable Statements ..
C
C GET 15 HI ORDER BITS OF IX
      XHI = IX/B16
C GET 16 LO BITS OF IX AND FORM LO PRODUCT
      XALO = (IX-XHI*B16)*A
C GET 15 HI ORDER BITS OF LO PRODUCT
      LEFTLO = XALO/B16
C     FORM THE 31 HIGHEST BITS OF FULL PRODUCT
      FHI = XHI*A + LEFTLO
C GET OVERFLOPAST 31ST BIT OF FULL PRODUCT
      K = FHI/B15
C ASSEMBLE ALL THE PARTS AND PRESUBTRACT P
C THE PARENTHESES ARE ESSENTIAL
      IX = (((XALO-LEFTLO*B16)-P)+ (FHI-K*B15)*B16) + K
C ADD P BACK IN IF NECCESSARY
      IF (IX.LT.0) IX = IX + P
C MULTIPLY BY 1/(2**31-1)
      XHI = IX/B16
      X = (FLOAT(XHI)*65536.0D0) + FLOAT(IX-XHI*B16)
      IF (I.GE.0) FA14AD = X*4.6566128752457969241D-10
      IF (I.LT.0) FA14AD = X*9.3132257504915938482D-10 - 1.0D0
      RETURN

      END
      SUBROUTINE FA14BD(IX,MAX,NRAND)
C         NEARLY PORTABLE RANDOM NUMBER GENERATOR USING THE RECURSION
C                       IX=IX*A MOD P
C
C    WHERE A=7**5
C    AND P=2**31-1.
C
C         THIS SUBROUTINE DOES NOT ADHERE TO THE ANSI STANDARD 1966
C    IN ONE RESPECT:
C         IT ASSUMES AN INTEGER WORD LENGTH OF AT LEAST 32 BITS (I.E.
C    INTEGERS WHICH LIE IN THE RANGE 1-2**31 TO 2**31-1 INCLUSIVE MUST
C    BE REPRESENTABLE).
C
C         THIS CODE IS BASED ON CODE PUBLISHED BY LINUS SCHRAGE IN
C    T.O.M.S. VOL.5 NO.2 JUNE 1979 (PP 132-138)
C
C
C       THE FUNCTION IS USED AS FOLLOWS:
C
C                  CALL FA14BD(IX,MAX,NRAND)
C
C       WHERE IX    IS THE GENERATOR WORD
C             MAX   IS AN INTEGER SET BY THE USER AND
C             NRAND IS AN INTEGER SET BY FA14B/BD.
C
C
C       THE VALUE OF NRAND RETURNED BY FA14B/BD WILL LIE IN THE
C   RANGE
C                        (1,MAX)
C
C       THE METHOD EMPLOYED IS A MULTIPLICATIVE CONGRUENTIAL
C   ONE USING A MULTIPLIER OF 7**5 AND TAKING THE MODULO TO
C   2**31-1, I.E. THE GENERATOR NUMBER , G = IX, IS UPDATED ON
C   EACH CALL TO THE VALUE
C
C                  5          31
C               G*7  MODULO (2  -1)
C
C       THE RESULT RETURNED IS AN INTEGER NUMBER
C   HAVING THE VALUE
C
C                        31
C   INT. PART( (MAX*G)/(2   -1) ) + 1
C
C
C 7**5, 2**15, 2**16, 2**31-1
C 2**30,  2**30-1
C     .. Parameters ..
      INTEGER A,B15,B16,P
      PARAMETER (A=16807,B15=32768,B16=65536,P=2147483647)
      INTEGER B30,Q
      PARAMETER (B30=1073741824,Q=1073741823)
C     ..
C     .. Scalar Arguments ..
      INTEGER IX,MAX,NRAND
C     ..
C     .. Local Scalars ..
      INTEGER BE1,BE2,C,D,F,FHI,G,K,LEFTLO,MHI,MLO,MU,NU,XALO,XHI,XLO
C     ..
C     .. Executable Statements ..
C
C GET 15 HI ORDER BITS OF IX
      XHI = IX/B16
C GET 16 LO BITS OF IX AND FORM LO PRODUCT
      XALO = (IX-XHI*B16)*A
C GET 15 HI ORDER BITS OF LO PRODUCT
      LEFTLO = XALO/B16
C     FORM THE 31 HIGHEST BITS OF FULL PRODUCT
      FHI = XHI*A + LEFTLO
C GET OVERFLOPAST 31ST BIT OF FULL PRODUCT
      K = FHI/B15
C ASSEMBLE ALL THE PARTS AND PRESUBTRACT P
C THE PARENTHESES ARE ESSENTIAL
      IX = (((XALO-LEFTLO*B16)-P)+ (FHI-K*B15)*B16) + K
C ADD P BACK IN IF NECCESSARY
      IF (IX.LT.0) IX = IX + P
C MULTIPLY BY MAX AND DIVIDE BY 2**31-1 IN INTEGER ARITHMETIC
C SPLIT IX AND MAX INTO HI AND LO PARTS
      XHI = IX/B15
      XLO = IX - B15*XHI
      MHI = MAX/B15
      MLO = MAX - B15*MHI
C CALCULATE INTERMEDIATE PRODUCT AND SPLIT INTO HI AND LO PARTS
C PRESUBTRACT P
      F = (XHI*MLO-P) + XLO*MHI
C F IS > 0 IF INTERMEDIATE PRODUCT WOULD HAVE OVERFLOWED
      IF (F.GT.0) GO TO 1
      F = F + P
      BE1 = F/B15
      BE2 = F - BE1*B15
      GO TO 2

    1 F = F - 1
      BE1 = F/B15
      BE2 = F - BE1*B15
      BE1 = BE1 + B16
C FORM PRODUCT OF LO PARTS AND ADD IN LO PART OF INTERMEDIATE PRODUCT
C TO GET LO PART OF COMPLETE PRODUCT
    2 G = B15*BE2 + XLO*MLO
C REPRESENT LO PART OF FULL PRODUCT IN BASE 2**30
      D = G/B30
      C = XHI/2
C CALCULATE FULL PRODUCT DIVIDED BY 2**30
      F = ((2* (C*MHI-Q)-1)+MHI* (XHI-2*C)) + D + BE1
C GET FULL PRODUCT DIVIDED IN BASE 2**31
      IF (F.GT.0) GO TO 3
      F = F + P
      NU = F/2
      MU = F - NU*2
      GO TO 4

    3 F = F - 1
      NU = F/2
      MU = F - 2*NU
      NU = NU + B30
C CALCULATE REMAINDER OF PRODUCT DIVIDED BY 2**31
    4 F = (B30*MU-P) + NU + (G-B30*D)
      NRAND = NU + 1
C  ADD ONE IF REMAINDER IS NOT < 2**31-1
      IF (F.GE.0) NRAND = NRAND + 1
      RETURN

      END
      SUBROUTINE FA14CD(IX,IGEN)
C        FA14CD IS A SUBROUTINE USED IN CONJUNCTION WITH FA14AD OR
C   FA14BD. IT PROVIDES THE USER WITH THE FACILITY OF SAVING THE
C   CURRENT VALUE OF THE GENERATOR NUMBER USED BY FA14AD AND FA14BD.
C
C        USE OF THE ROUTINE IS AS FOLLOWS:
C
C                       CALL FA14CD(IGEN)
C
C     WHERE IX   IS THE GENERATOR WORD
C           IGEN IS AN INTEGER WHICH IS SET BY FA14C/CD TO THE CURRENT
C                VALUE OF THE GENERATOR.
C
C
C     .. Scalar Arguments ..
      INTEGER IX,IGEN
C     ..
C     .. Executable Statements ..
      IGEN = IX
      RETURN

      END
      SUBROUTINE FA14DD(IX,IGEN)
C        FA14DD IS A SUBROUTINE USED IN CONJUNCTION WITH FA14AD OR
C   FA14BD. IT PROVIDES THE USER WITH THE FACILITY OF SETTING THE
C   CURRENT VALUE OF THE GENERATOR NUMBER USED BY FA14AD AND FA14BD.
C
C        USE OF THE ROUTINE IS AS FOLLOWS:
C
C                       CALL FA14DD(IGEN)
C
C    WHERE IX   IS THE GENERATOR WORD
C          IGEN IS AN INTEGER, SET BY THE USER TO THE VALUE TO WHICH
C               THE GENERATOR IS TO BE SET. IT IS RECOMMENDED THAT THIS
C               VALUE BE OBTAINED BY A PREVIOUS CALL TO FA14C/CD.
C
C     .. Scalar Arguments ..
      INTEGER IX,IGEN
C     ..
C     .. Executable Statements ..
      IX = IGEN
      RETURN

      END
      SUBROUTINE FA14ID(IX)
C        FA14ID IS A SUBROUTINE USED TO INITIALIZE THE GENERATOR WORD US
C   BY FA14AD AND FA14BD. IT MUST BE CALLED FIRST BEFORE ANY OF THE OTHE
C   ENTRIES ARE CALLED.
C
C        USE OF THE ROUTINE IS AS FOLLOWS:
C
C                       CALL FA14ID(IX)
C
C    WHERE IX   IS THE GENERATOR WORD
C
C     .. Scalar Arguments ..
      INTEGER IX
C     ..
C     .. Executable Statements ..
      IX = 1
      RETURN

      END
! COPYRIGHT (c) 1979 AEA Technology and
! Council for the Central Laboratory of the Research Councils
!
! Version 1.0.1
! See ChangeLog for version history
!
      REAL FUNCTION FA14A(IX,I)
C         NEARLY PORTABLE RANDOM NUMBER GENERATOR USING THE RECURSION
C                       IX=IX*A MOD P
C
C    WHERE A=7**5
C    AND P=2**31-1.
C
C         THIS FUNCTION DOES NOT ADHERE TO THE ANSI STANDARD 1966 IN
C    TWO RESPECTS:
C      1) IT ASSUMES AN INTEGER WORD LENGTH OF AT LEAST 32 BITS (I.E.
C    INTEGERS WHICH LIE IN THE RANGE 1-2**31 TO 2**31-1 INCLUSIVE MUST
C    BE REPRESENTABLE);
C      2) IT ASSUMES THAT A POSITIVE INTEGER LESS THAN 2**16 MAY BE
C    FLOATED WITHOUT LOSS OF DIGITS.
C
C         THIS CODE IS BASED ON CODE PUBLISHED BY LINUS SCHRAGE IN
C    T.O.M.S. VOL.5 NO.2 JUNE 1979 (PP 132-138)
C
C
C
C       THE FUNCTION IS USED AS FOLLOWS:
C
C                      R=FA14A(IX,I)
C
C       WHERE IX IS THE GENERATOR WORD
C             I IS AN INTEGER SET BY THE USER.
C
C
C       THE VALUE RETURNED BY FA14A/AD WILL LIE IN THE RANGE
C                (0.,1.)  IF I IS NON-NEGATIVE
C                (-1.,1.) IF I IS NEGATIVE.
C
C       THE METHOD EMPLOYED IS A MULTIPLICATIVE CONGRUENTIAL
C   ONE USING A MULTIPLIER OF 7**5 AND TAKING THE MODULO TO
C   2**31-1, I.E. THE GENERATOR NUMBER , G = IX, IS UPDATED ON
C   EACH CALL TO THE VALUE
C
C                  5          31
C               G*7  MODULO (2  -1)
C
C       THE RESULT RETURNED IS CALCULATED AS A DOUBLE
C  PRECISION NUMBER HAVING THE VALUE
C
C                      31
C                  G/(2   -1)    IF THE ARGUMENT IS
C                                NON-NEGATIVE
C           OR
C                      31
C                2*G/(2   -1)-1  IF THE ARGUMENT IS NEGATIVE
C
C
C 7**5, 2**15, 2**16, 2**31-1
C     .. Parameters ..
      INTEGER A,B15,B16,P
      PARAMETER (A=16807,B15=32768,B16=65536,P=2147483647)
C     ..
C     .. Scalar Arguments ..
      INTEGER IX,I
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION X
      INTEGER FHI,K,LEFTLO,XALO,XHI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC FLOAT
C     ..
C     .. Executable Statements ..
C
C GET 15 HI ORDER BITS OF IX
      XHI = IX/B16
C GET 16 LO BITS OF IX AND FORM LO PRODUCT
      XALO = (IX-XHI*B16)*A
C GET 15 HI ORDER BITS OF LO PRODUCT
      LEFTLO = XALO/B16
C     FORM THE 31 HIGHEST BITS OF FULL PRODUCT
      FHI = XHI*A + LEFTLO
C GET OVERFLOPAST 31ST BIT OF FULL PRODUCT
      K = FHI/B15
C ASSEMBLE ALL THE PARTS AND PRESUBTRACT P
C THE PARENTHESES ARE ESSENTIAL
      IX = (((XALO-LEFTLO*B16)-P)+ (FHI-K*B15)*B16) + K
C ADD P BACK IN IF NECCESSARY
      IF (IX.LT.0) IX = IX + P
C MULTIPLY BY 1/(2**31-1)
      XHI = IX/B16
      X = (FLOAT(XHI)*65536.0D0) + FLOAT(IX-XHI*B16)
      IF (I.GE.0) FA14A = X*4.6566128752457969241D-10
      IF (I.LT.0) FA14A = X*9.3132257504915938482D-10 - 1.0D0
      RETURN

      END
      SUBROUTINE FA14B(IX,MAX,NRAND)
C         NEARLY PORTABLE RANDOM NUMBER GENERATOR USING THE RECURSION
C                       IX=IX*A MOD P
C
C    WHERE A=7**5
C    AND P=2**31-1.
C
C         THIS SUBROUTINE DOES NOT ADHERE TO THE ANSI STANDARD 1966
C    IN ONE RESPECT:
C         IT ASSUMES AN INTEGER WORD LENGTH OF AT LEAST 32 BITS (I.E.
C    INTEGERS WHICH LIE IN THE RANGE 1-2**31 TO 2**31-1 INCLUSIVE MUST
C    BE REPRESENTABLE).
C
C         THIS CODE IS BASED ON CODE PUBLISHED BY LINUS SCHRAGE IN
C    T.O.M.S. VOL.5 NO.2 JUNE 1979 (PP 132-138)
C
C
C       THE FUNCTION IS USED AS FOLLOWS:
C
C                  CALL FA14B(IX,MAX,NRAND)
C
C       WHERE IX    IS THE GENERATOR WORD
C             MAX   IS AN INTEGER SET BY THE USER AND
C             NRAND IS AN INTEGER SET BY FA14B/BD.
C
C
C       THE VALUE OF NRAND RETURNED BY FA14B/BD WILL LIE IN THE
C   RANGE
C                        (1,MAX)
C
C       THE METHOD EMPLOYED IS A MULTIPLICATIVE CONGRUENTIAL
C   ONE USING A MULTIPLIER OF 7**5 AND TAKING THE MODULO TO
C   2**31-1, I.E. THE GENERATOR NUMBER , G = IX, IS UPDATED ON
C   EACH CALL TO THE VALUE
C
C                  5          31
C               G*7  MODULO (2  -1)
C
C       THE RESULT RETURNED IS AN INTEGER NUMBER
C   HAVING THE VALUE
C
C                        31
C   INT. PART( (MAX*G)/(2   -1) ) + 1
C
C
C 7**5, 2**15, 2**16, 2**31-1
C 2**30,  2**30-1
C     .. Parameters ..
      INTEGER A,B15,B16,P
      PARAMETER (A=16807,B15=32768,B16=65536,P=2147483647)
      INTEGER B30,Q
      PARAMETER (B30=1073741824,Q=1073741823)
C     ..
C     .. Scalar Arguments ..
      INTEGER IX,MAX,NRAND
C     ..
C     .. Local Scalars ..
      INTEGER BE1,BE2,C,D,F,FHI,G,K,LEFTLO,MHI,MLO,MU,NU,XALO,XHI,XLO
C     ..
C     .. Executable Statements ..
C
C GET 15 HI ORDER BITS OF IX
      XHI = IX/B16
C GET 16 LO BITS OF IX AND FORM LO PRODUCT
      XALO = (IX-XHI*B16)*A
C GET 15 HI ORDER BITS OF LO PRODUCT
      LEFTLO = XALO/B16
C     FORM THE 31 HIGHEST BITS OF FULL PRODUCT
      FHI = XHI*A + LEFTLO
C GET OVERFLOPAST 31ST BIT OF FULL PRODUCT
      K = FHI/B15
C ASSEMBLE ALL THE PARTS AND PRESUBTRACT P
C THE PARENTHESES ARE ESSENTIAL
      IX = (((XALO-LEFTLO*B16)-P)+ (FHI-K*B15)*B16) + K
C ADD P BACK IN IF NECCESSARY
      IF (IX.LT.0) IX = IX + P
C MULTIPLY BY MAX AND DIVIDE BY 2**31-1 IN INTEGER ARITHMETIC
C SPLIT IX AND MAX INTO HI AND LO PARTS
      XHI = IX/B15
      XLO = IX - B15*XHI
      MHI = MAX/B15
      MLO = MAX - B15*MHI
C CALCULATE INTERMEDIATE PRODUCT AND SPLIT INTO HI AND LO PARTS
C PRESUBTRACT P
      F = (XHI*MLO-P) + XLO*MHI
C F IS > 0 IF INTERMEDIATE PRODUCT WOULD HAVE OVERFLOWED
      IF (F.GT.0) GO TO 1
      F = F + P
      BE1 = F/B15
      BE2 = F - BE1*B15
      GO TO 2

    1 F = F - 1
      BE1 = F/B15
      BE2 = F - BE1*B15
      BE1 = BE1 + B16
C FORM PRODUCT OF LO PARTS AND ADD IN LO PART OF INTERMEDIATE PRODUCT
C TO GET LO PART OF COMPLETE PRODUCT
    2 G = B15*BE2 + XLO*MLO
C REPRESENT LO PART OF FULL PRODUCT IN BASE 2**30
      D = G/B30
      C = XHI/2
C CALCULATE FULL PRODUCT DIVIDED BY 2**30
      F = ((2* (C*MHI-Q)-1)+MHI* (XHI-2*C)) + D + BE1
C GET FULL PRODUCT DIVIDED IN BASE 2**31
      IF (F.GT.0) GO TO 3
      F = F + P
      NU = F/2
      MU = F - NU*2
      GO TO 4

    3 F = F - 1
      NU = F/2
      MU = F - 2*NU
      NU = NU + B30
C CALCULATE REMAINDER OF PRODUCT DIVIDED BY 2**31
    4 F = (B30*MU-P) + NU + (G-B30*D)
      NRAND = NU + 1
C  ADD ONE IF REMAINDER IS NOT < 2**31-1
      IF (F.GE.0) NRAND = NRAND + 1
      RETURN

      END
      SUBROUTINE FA14C(IX,IGEN)
C        FA14C IS A SUBROUTINE USED IN CONJUNCTION WITH FA14A OR
C   FA14B. IT PROVIDES THE USER WITH THE FACILITY OF SAVING THE
C   CURRENT VALUE OF THE GENERATOR NUMBER USED BY FA14A AND FA14B.
C
C        USE OF THE ROUTINE IS AS FOLLOWS:
C
C                       CALL FA14C(IX,IGEN)
C
C     WHERE IX   IS THE GENERATOR WORD
C           IGEN IS AN INTEGER WHICH IS SET BY FA14C/CD TO THE CURRENT
C                VALUE OF THE GENERATOR.
C
C
C     .. Scalar Arguments ..
      INTEGER IX,IGEN
C     ..
C     .. Executable Statements ..
      IGEN = IX
      RETURN

      END
      SUBROUTINE FA14D(IX,IGEN)
C        FA14D IS A SUBROUTINE USED IN CONJUNCTION WITH FA14A OR
C   FA14B. IT PROVIDES THE USER WITH THE FACILITY OF SETTING THE
C   CURRENT VALUE OF THE GENERATOR NUMBER USED BY FA14A AND FA14B.
C
C        USE OF THE ROUTINE IS AS FOLLOWS:
C
C                       CALL FA14D(IX,IGEN)
C
C    WHERE IX   IS THE GENERATOR WORD
C          IGEN IS AN INTEGER, SET BY THE USER TO THE VALUE TO WHICH
C               THE GENERATOR IS TO BE SET. IT IS RECOMMENDED THAT THIS
C               VALUE BE OBTAINED BY A PREVIOUS CALL TO FA14C/CD.
C
C     .. Scalar Arguments ..
      INTEGER IX,IGEN
C     ..
C     .. Executable Statements ..
      IX = IGEN
      RETURN

      END
      SUBROUTINE FA14I(IX)
C        FA14I IS A SUBROUTINE USED TO INITIALIZE THE GENERATOR WORD USE
C   BY FA14A AND FA14B. IT MUST BE CALLED FIRST BEFORE ANY OF THE OTHER
C   ENTRIES ARE CALLED.
C
C        USE OF THE ROUTINE IS AS FOLLOWS:
C
C                       CALL FA14I(IX)
C
C    WHERE IX   IS THE GENERATOR WORD
C
C     .. Scalar Arguments ..
      INTEGER IX
C     ..
C     .. Executable Statements ..
      IX = 1
      RETURN

      END
