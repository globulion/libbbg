C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE PMLOCA(NATOMS,NBASIS,NMOS,MAPI,SAO,VECIN,
     &                  MAXIT,CVGLOC,N2,NAE,LPRINT,TRAN)
C
C          Pipek-Mezey molecular orbital localization scheme
C          -------------------------------------------------
C          The algorithm is taken from GAMESS package
C
C                                                 19.08.2013
C
C   Variables:
C     MAPI   - list of atomic indices in bfs order
C     VECIN  - input canonical orbitals
C     CVGLOC - convergence criterion
C     NAE    - no of alpha electrons (assumed closed shell)
C     N2     - no of triangular matrix elements (N+1)N/2
C     
C   Returns:
C     TRAN   - transformation matrix
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION VECIN(NMOS,NBASIS), SAO(NBASIS,NBASIS),
     &          MAPI(NBASIS),RIJ(N2,NATOMS),QPIX(NBASIS),
     &          QPJX(NBASIS),IORD(NBASIS),IIR(NBASIS)
      PARAMETER (ZERO=0.D0,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00,
     &           TENM3=1.0D-03,TENM8=1.0D-08,TENM10=1.0D-10)
      DOUBLE PRECISION TRAN(NMOS,NMOS), LMORND
      LOGICAL LPRINT
Cf2py INTENT(OUT) TRAN
C
      INDX(I,J) = ((MAX(I,J)*(MAX(I,J)-1))/2+MIN(I,J))
C
      NREDO = 0
  110 CONTINUE
      NREDO=NREDO+1
C
C...construct initial atomic populations
      CALL VCLR(RIJ,1,N2*NATOMS)
      IJ = 0
      DO 280 I = 1,NMOS
      DO 280 J = 1,I
         IJ = IJ+1
         DO 260 K = 1,NBASIS
            KK = MAPI(K)
            DO 260 L = 1,NBASIS
               LL = MAPI(L)
C                KL = INDX(K,L)
               SUM = VECIN(I,K)*VECIN(J,L)*SAO(K,L)/TWO
               RIJ(IJ,KK) = RIJ(IJ,KK) + SUM
               RIJ(IJ,LL) = RIJ(IJ,LL) + SUM
  260       CONTINUE
  280 CONTINUE
C
C...seed the random function
      ONEPT0 = ONE
      IF (NREDO.EQ.1) XX = LMORND(ONEPT0,VECIN,NBASIS,NAE,NBASIS)
C
C...initialize transformation matrix
      CALL VCLR(TRAN,1,NMOS*NMOS)
      DO 340 I = 1,NMOS
         TRAN(I,I) = ONE
  340 CONTINUE
C        -------------------------
C        BEGIN LOCALIZATION CYCLES
C        -------------------------
      ITER = 0
      SHIFT = DATAN(ONEPT0)
C
  360 CONTINUE
      CHANGE = ZERO
      ITER = ITER+1
      DO 380 I = 1,NMOS
         IIR(I) = I
  380 CONTINUE
      NNN = NMOS
      DO 400 I = 1,NMOS
         XX = LMORND(CHANGE,VECIN,NBASIS,NAE,NBASIS)
         III = INT(XX*NNN+ONE)
         IORD(I) = IIR(III)
         IIR(III) = IIR(NNN)
         NNN = NNN-1
  400 CONTINUE
C
C        FOR EACH PAIR OF ORBITALS A TWO DIMENSIONAL UNITARY
C        TRANSFORMATION IS PERFORMED. THE TRANSFORMATION IS
C
C           PSI'(I) =  COS(T)*PSI(I) + SIN(T)*PSI(J)  AND
C           PSI'(J) = -SIN(T)*PSI(I) + COS(T)*PSI(J).
C
C        LOCALIZATION REQUIRES THAT T BE SUCH AS TO MAXIMIZE
C        THE SUM OF THE SQUARES OF THE ATOMIC POPULATIONS.
C
      DO 920 III = 1,NMOS
         I  = IORD(III)
         II = INDX(I,I)
         JM = 1
         RM = ZERO
         TM = ZERO
         SM = ZERO
         CM = ONE
         DO 580 J = 1,NMOS
            IF(I.EQ.J) GO TO 580
            IJ = INDX(I,J)
            JJ = INDX(J,J)
            T = ZERO
            TX = ZERO
            DO 480 KK = 1,NATOMS
               T= T + FOUR*RIJ(IJ,KK)**2 - RIJ(II,KK)**2 - RIJ(JJ,KK)**2
     &         + TWO*RIJ(II,KK)*RIJ(JJ,KK)
               TX = TX + RIJ(IJ,KK)*(RIJ(JJ,KK) - RIJ(II,KK))
  480       CONTINUE
            IF ((DABS(T).LE.TENM10).AND.(DABS(TX).LE.TENM10)) GO TO 580
            TX = FOUR*TX
            T = DATAN2(TX,T)/FOUR
            SIGN = ONE
            IF (T.GT.ZERO) SIGN = -ONE
            T = T+SIGN*SHIFT
            ITIM = 0
  500       ITIM = ITIM+1
            S = DSIN(T)
            C = DCOS(T)
            RIN = ZERO
            DO 520 KK = 1,NATOMS
               QPI = C*C*RIJ(II,KK)+S*S*RIJ(JJ,KK)+TWO*C*S*RIJ(IJ,KK)
               QPJ = C*C*RIJ(JJ,KK)+S*S*RIJ(II,KK)-TWO*C*S*RIJ(IJ,KK)
               RIN = RIN+QPI*QPI+QPJ*QPJ-RIJ(II,KK)**2-RIJ(JJ,KK)**2
  520       CONTINUE
            TTEST = DABS(T)-SHIFT
            IF ((DABS(T).LE.TENM8).OR.(DABS(TTEST).LE.TENM8)) GO TO 560
            IF (RIN .GE. -TENM8) GO TO 560
            IF (ITIM .LE. 1) GO TO 540
            RETURN
C
  540       SIGN = ONE
            IF (T .GT. ZERO) SIGN = -ONE
            T = T+SHIFT*SIGN
            GO TO 500
C
  560       IF (RIN .LE. RM) GO TO 580
            RM = RIN
            TM = T
            SM = S
            CM = C
            JM = J
  580    CONTINUE
C
         RIN = RM
         T = TM
         S = SM
         C = CM
         J = JM
         IJ = INDX(I,J)
         JJ = INDX(J,J)
C
C...accumulate the 2x2 rotation
         CHANGE = CHANGE+T*T
         CALL DROT(NMOS,TRAN(1,I),1,TRAN(1,J),1,C,S)
C
C...update the atomic populations
         DO 880 KK = 1,NATOMS
            QPI = C*C*RIJ(II,KK)+S*S*RIJ(JJ,KK)+TWO*C*S*RIJ(IJ,KK)
            QPJ = C*C*RIJ(JJ,KK)+S*S*RIJ(II,KK)-TWO*C*S*RIJ(IJ,KK)
            QPIJ = (C*C-S*S)*RIJ(IJ,KK)+C*S*(RIJ(JJ,KK)-RIJ(II,KK))
            DO 720 K = 1,NMOS
               IF (I.EQ.K.OR.J.EQ.K) GO TO 720
               IK = INDX(I,K)
               JK = INDX(J,K)
               QPIX(K) = C*RIJ(IK,KK)+S*RIJ(JK,KK)
               QPJX(K) = C*RIJ(JK,KK)-S*RIJ(IK,KK)
               RIJ(IK,KK) = QPIX(K)
               RIJ(JK,KK) = QPJX(K)
  720       CONTINUE
            RIN = RIN+QPI+QPJ-RIJ(II,KK)-RIJ(JJ,KK)
            RIJ(II,KK) = QPI
            RIJ(JJ,KK) = QPJ
            RIJ(IJ,KK) = QPIJ
  880    CONTINUE
  920 CONTINUE
C
C        CONVERGED?
C
      CHANGE = DSQRT(TWO*CHANGE/(NMOS*(NMOS-1)))
      IF(ITER.LT.MAXIT  .AND.  CHANGE.GT.TENM3*CVGLOC) GO TO 360
      IF(CHANGE.LE.CVGLOC) GO TO 1000
         IF(NREDO.LE.2) GO TO 110
            RETURN
C          ---------------------------------
C          FINISHED WITH LOCALIZATION CYCLES
C          ---------------------------------
 1000 CONTINUE
      IF (LPRINT) THEN
          WRITE (6,9080) ITER
      ENDIF
 9080 FORMAT(10X,'LOCALIZATION CONVERGED IN',I4,' ITERATIONS')
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE VCLR(C,IC,N)
C
C     Return clear vector of zeros (dp) of size N in increment IC
C
      INTEGER IC,N,KK,M
      DOUBLE PRECISION C(1)
      IF (N.LE.0) GO TO 12
      KK = 1
      DO 10 M=1,N
        C(KK) = 0.D0
        KK = KK + IC
10    CONTINUE
12    RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      subroutine  drot (n,dx,incx,dy,incy,c,s)
C
C     applies a plane rotation.
C     jack dongarra, linpack, 3/11/78.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
      double precision dx(n),dy(n),dtemp,c,s
      integer i,incx,incy,ix,iy,n
C
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
C
C       code for unequal increments or equal increments not equal
C         to 1
C
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = c*dx(ix) + s*dy(iy)
        dy(iy) = c*dy(iy) - s*dx(ix)
        dx(ix) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
C
C       code for both increments equal to 1
C
   20 do 30 i = 1,n
        dtemp = c*dx(i) + s*dy(i)
        dy(i) = c*dy(i) - s*dx(i)
        dx(i) = dtemp
   30 continue
      return
      end
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      DOUBLE PRECISION FUNCTION LMORND(XX,D,L1,NAE,NBASIS)
C
C         Return some random
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATM=2000)
      DIMENSION D(L1,L1),U(1)
      SAVE U
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
C
      PI = DACOS(-ONE)
      IF (XX .EQ. ZERO) GO TO 100
         N = ABS(NAE-NBASIS)+1
         M = N+5
         XY = D(N,M)*DATAN(ONE)
         U(1) = (PI+XY)**5
         IU1 = INT(U(1))
         XY = IU1
         U(1) = U(1)-XY
         LMORND = U(1)
         RETURN
C
  100 CONTINUE
      U(1) = (PI+U(1))**5
      IU1 = INT(U(1))
      XY = IU1
      U(1) = U(1)-XY
      LMORND = U(1)
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|
