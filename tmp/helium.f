      PROGRAM QDMC
C
C    -----------------------------------------------------------
C                MONOATOMIC AND DIATOMIC MOLECULES
C                               ***
C           QUANTUM DIFFIUSION MONTE CARLO SIMULATION 
C                     WITH IMPORTANCE SAMPLING
C    
C      Bartosz Błasiak                   11 Dec 2013
C    -----------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.00D+00,ONE=1.00D+00,HALF=0.5000D+00,
     &           DCONST=0.500D+00,MWALK=900000,NBIN=200)
      DIMENSION X1(MWALK,3),X2(MWALk,3),E(MWALK),IMULT(MWALK),
     &          IBINPOP(NBIN,NBIN,NBIN)
C
C     --- INPUT PARAMETERS
C
      DATA XMIN/-8.0D+00/,XMAX/8.0D+00/
      DATA MAXWLK/MWALK/,MINWLK/10000/,NORMWLK/100000/
      DATA NBLOCKS/2000/,NSTEPS/100/,ET/-2.88700D+00/,NEQUIL/100/
      DATA ISEED/77741/
      DATA TAU/0.01D+00/
C
C     INITIALIZE WALKER POSITIONS
C
      WRITE(*,*) " INITIALIZATION RUNNING ... "
C
      CALL XINIT(X1,X2,XMIN,XMAX,MWALK,MAXWLK,ISEED)
C
      DO K=1,MAXWLK
         E(K) = ZERO
         IMULT(K) = 0
      ENDDO
C
      DO 1234 I=1,NBIN
      DO 1234 J=1,NBIN
      DO 1234 K=1,NBIN
         IBINPOP(I,J,K) = 0
 1234 CONTINUE
      OPEN(901,FILE='energy.dat')
C
C     ITERATE OVER BLOCKS
C
      WRITE(*,*) " RUNNING OVER BLOCKS ... "
C
      NWALKERS = NORMWLK
C
      EAVE = ZERO
      EAV2 = ZERO
      EGAV = ZERO
      EGV2 = ZERO
      NBDONE = 0
      NPOPT = 0
C
      DO 900 IBLOCK=1,NBLOCKS
         IF (IBLOCK.GT.NEQUIL) NBDONE = NBDONE + 1
         IF (IBLOCK-1.EQ.NEQUIL) THEN
             WRITE(*,*) " --- EQUILIBRATION FINISHED.
     & GATHERING STATISTICS STARTED ---"
         ENDIF
C
         NOLD = NWALKERS
         CALL BLOCK(IBLOCK,NWALKERS,NSTEPS,ISEED,E,X1,X2,
     &              EAVE,EAV2,ET,NBDONE,NEQUIL,MWALK,
     &              MINWLK,MAXWLK,NORMWLK,IMULT,TAU,
     &              XMIN,XMAX,NBIN,IBINPOP,EBLOCK,EG,
     &              NPOPT,EGAV,EGV2)
C
c         RATIO = DFLOAT(NOLD)/DFLOAT(NWALKERS)
c         ET = ET + DLOG(RATIO)/(TAU*4.D00)
         WRITE(901,1020) IBLOCK,EBLOCK,EG,NWALKERS
 1020    FORMAT(I8,D22.8,D22.8,I20)
C
 900  CONTINUE
C
      EAVE = EAVE / NBDONE
      EAV2 = EAV2 / NBDONE
      EGAV = EGAV / NBDONE
      EGV2 = EGV2 / NBDONE
      VAR = DSQRT(DABS(EAVE*EAVE-EAV2)/NBDONE)
      VRG = DSQRT(DABS(EGAV*EGAV-EGV2)/NBDONE)
C
      WRITE(*,1030) EAVE,VAR,EGAV,VRG
 1030 FORMAT(/,' ======= GRAND AVERAGES ====== ',/,
     &  ' GROUND STATE ENERGY',T37,D12.6,1X,'+/-',1X,D12.6,/,
     &  ' GROUND STATE ENERGY',T37,D12.6,1X,'+/-',1X,D12.6)
C
      WRITE(*,1040) 
 1040 FORMAT(/,' SAVING WAVEFUNCTION 3D HISTOGRAM',/)
c      CALL SAVBIN(IBINPOP,NBIN,XMIN,XMAX,NPOPT)
C
      STOP
      END
C ==============================================================

      SUBROUTINE BLOCK(IBLOCK,NWALKERS,NSTEPS,ISEED,E,X1,X2,
     &                 EAVE,EAV2,ET,NBDONE,NEQUIL,MWALK,
     &                 MINWLK,MAXWLK,NORMWLK,IMULT,TAU,
     &                 XMIN,XMAX,NBIN,IBINPOP,EBLOCK,EG,
     &                 NPOPT,EGAV,EGV2)
C
C     NSTEPS ITERATIONS FOR ONE BLOCK
C
C---------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.00D+00,ONE=1.00D+00,HALF=0.5000D+00,
     &           TWO=2.00D+00,DCONST=0.500D+00)
      DIMENSION X1(MWALK,3),X2(MWALK,3),Y1(3),Y2(3),FQX1(3),FQY1(3),
     &          E(MWALK),IMULT(MWALK),IBINPOP(NBIN,NBIN,NBIN),
     &          FQX2(3),FQY2(3)
c      DATA Y1/3*0.0D+00/,FQX1/3*0.0D+00/,FQY/3*0.0D+00/
C
      EBLOCK = ZERO
      EG  = ZERO
      DNEL= TWO
      VAR = TAU
      DTAU=DSQRT(TAU)
      DC  = DCONST*TAU*HALF
      DNSTEPS = DFLOAT(NSTEPS)
C
      DO 100 NSTEP=1,NSTEPS
         M = NWALKERS
C
C        --- MOVE WALKERS ---
C
         DO 110 K=1,M
            ACCEPT = ZERO
C
            CALL TRIALF(X1(K,1),X1(K,2),X1(K,3),
     &                  X2(K,1),X2(K,2),X2(K,3),
     &                  PSIX,ELX,FQX1,FQX2)
C
C           --- TRANSITION DUE TO DIFFUSION AND QUANTUM FORCE DRIFT ---
C
            GRX1 = BOXMLR(ISEED,VAR)
            GRY1 = BOXMLR(ISEED,VAR)
            GRZ1 = BOXMLR(ISEED,VAR)
            GRX2 = BOXMLR(ISEED,VAR)
            GRY2 = BOXMLR(ISEED,VAR)
            GRZ2 = BOXMLR(ISEED,VAR)
C
 330        CONTINUE
C
            Y1(1) = X1(K,1) + DCONST*TAU*FQX1(1) + GRX1
            Y1(2) = X1(K,2) + DCONST*TAU*FQX1(2) + GRY1
            Y1(3) = X1(K,3) + DCONST*TAU*FQX1(3) + GRZ1
            Y2(1) = X2(K,1) + DCONST*TAU*FQX2(1) + GRX2
            Y2(2) = X2(K,2) + DCONST*TAU*FQX2(2) + GRY2
            Y2(3) = X2(K,3) + DCONST*TAU*FQX2(3) + GRZ2
C
            CALL TRIALF(Y1(1),Y1(2),Y1(3),
     &                  Y2(1),Y2(2),Y2(2),
     &                  PSIY,ELY,FQY1,FQY2)
C
C           --- PERFORM GENERALIZED METROPOLIS SAMPLING ---
C
            Q = ZERO
C
            Q = Q + HALF*(FQX1(1)+FQY1(1))*(DC
     &                                *(FQX1(1)-FQY1(1))
     &                                -(Y1(1)-X1(K,1))) +
C
     &              HALF*(FQX1(2)+FQY1(2))*(DC
     &                                *(FQX1(2)-FQY1(2))
     &                                -(Y1(2)-X1(K,2))) +
C
     &              HALF*(FQX1(3)+FQY1(3))*(DC
     &                                *(FQX1(3)-FQY1(3))
     &                                -(Y1(3)-X1(K,3))) +
C
     &              HALF*(FQX2(1)+FQY2(1))*(DC
     &                                *(FQX2(1)-FQY2(1))
     &                                -(Y2(1)-X2(K,1))) +
C
     &              HALF*(FQX2(2)+FQY2(2))*(DC
     &                                *(FQX2(2)-FQY2(2))
     &                                -(Y2(2)-X2(K,2))) +
C
     &              HALF*(FQX2(3)+FQY2(3))*(DC
     &                                *(FQX2(3)-FQY2(3))
     &                                -(Y2(3)-X2(K,3)))
            PP = PSIY/PSIX
            A = MIN(ONE,PP*PP*DEXP(Q))
C
c            IF(A.GT.ROULET(ISEED)) THEN
            IF(A.GT.URAN(0)) THEN
               X1(K,1) = Y1(1)
               X1(K,2) = Y1(2)
               X1(K,3) = Y1(3)
               X2(K,1) = Y2(1)
               X2(K,2) = Y2(2)
               X2(K,3) = Y2(3)
               ACCEPT = ACCEPT + ONE
               E(K) = ELY
C   
C              ESTIMATE THE BRANCHING WEIGHT
C
               GB = DEXP(-(HALF*(ELX+ELY)-ET)*TAU*ACCEPT)
c              IMULT(K) = INT( GB + ROULET(ISEED) )
               IMULT(K) = INT( GB + URAN(0) )
            ELSE 
c                GO TO 330
               IMULT(K) = 1
            ENDIF
C
C           END K-LOOP OVER WALKERS
C
 110     CONTINUE
C
C        --- DO THE BRANCHING ---
C
         K = 0
 120     K = K + 1
         IF(K.GT.NWALKERS) GOTO 150
 130     IF(IMULT(K).GT.0) GOTO 120
C
C        DELETE DEAD WALKERS
C
         IF(K.EQ.NWALKERS) GOTO 140
         IMULT(K) = IMULT(NWALKERS)
         E(K) = E(NWALKERS)
         X1(K,1) = X1(NWALKERS,1)
         X1(K,2) = X1(NWALKERS,2)
         X1(K,3) = X1(NWALKERS,3)
         X2(K,1) = X2(NWALKERS,1)
         X2(K,2) = X2(NWALKERS,2)
         X2(K,3) = X2(NWALKERS,3)
         NWALKERS = NWALKERS - 1
         GOTO 130
 140     NWALKERS = NWALKERS - 1
 150     CONTINUE
         NEWPOP=0 
         DO 160 IW=1,NWALKERS
            NEWPOP = NEWPOP + IMULT(IW)
 160     CONTINUE
C
C        REPLICATE AND ASSIGN X() AND E() TO NEWLY CREATED WALKERS
C
         K = NWALKERS
         IW = 0
 170     IW = IW + 1
         IF(IW.GT.NWALKERS) GOTO 190
         JJ = 0
 180     JJ = JJ + 1
         IF(JJ.GT.IMULT(IW)-1) GOTO 170
         K = K + 1
C
C        TAKE CARE OF OVERFLOWING THE ENSEMBLE ARRAY
C
         IF(K.GT.MAXWLK) THEN
            NWALKERS = NORMWLK 
            GO TO 200
         ENDIF
         E(K) = E(IW)
         X1(K,1) = X1(IW,1)
         X1(K,2) = X1(IW,2)
         X1(K,3) = X1(IW,3)
         X2(K,1) = X2(IW,1)
         X2(K,2) = X2(IW,2)
         X2(K,3) = X2(IW,3)
         GO TO 180
 190     CONTINUE
         NWALKERS = NEWPOP
C
C        TAKE CARE OF UNDERFLOWS
C 
         IF(NWALKERS.LT.1) THEN
            NWALKERS = MINWLK
         ENDIF
200      CONTINUE
C
         RATIO = DFLOAT(M)/DFLOAT(NWALKERS)
         EG = EG + (DLOG(RATIO)/TAU + ET)
C         RATIO = (10000.0D00/DFLOAT(NWALKERS))
C         ET = ET + DLOG(RATIO)/10.D00
C         ET = ET + DLOG(ONE/RATIO)/(TAU*4.000D+00)
C
         DO K=1,NWALKERS
            EBLOCK = EBLOCK + E(K)
         ENDDO 
C
C        CONTINUE THE NEXT ITERATION
C
 100  CONTINUE
C
C     SUMMARIZE BLOCK RUN AND SAVE PROPERTIES
C
      EBLOCK = EBLOCK / NSTEPS / NWALKERS
      EG = EG / NSTEPS
      IF (IBLOCK.LE.NEQUIL) THEN
          ET = HALF*(ET+EG)
      ELSE
          ET = ET + HALF*(EG-ET)/NBDONE 
      ENDIF
      IF (IBLOCK.GT.NEQUIL) THEN
          EAVE = EAVE + EBLOCK
          EAV2 = EAV2 + EBLOCK*EBLOCK
          EGAV = EGAV + EG
          EGV2 = EGV2 + EG*EG
C
C         ACCUMULATE WAVEFUNCTION 3D-HISTOGRAM
C
          BINWIDTH = (XMAX-XMIN)/NBIN
C
          DO 1111 IW=1,NWALKERS
C
          XX1 = X1(IW,1)
          YY1 = X1(IW,2)
          ZZ1 = X1(IW,3)
          XX2 = X2(IW,1)
          YY2 = X2(IW,2)
          ZZ2 = X2(IW,3)
C
          AMINX = XMIN
          AMINY = XMIN
          AMINZ = XMIN
C
          DO IBX = 1,NBIN
             AMAXX = AMINX + BINWIDTH
             IF(XX1.GE.AMINX.AND.XX1.LT.AMAXX) THEN
             DO IBY = 1,NBIN
                AMAXY = AMINY + BINWIDTH
                IF(YY1.GE.AMINY.AND.YY1.LT.AMAXY) THEN
                DO IBZ = 1,NBIN
                   AMAXZ = AMINZ + BINWIDTH
                   IF(ZZ1.GE.AMINZ.AND.ZZ1.LT.AMAXZ) THEN
                      IBINPOP(IBX,IBY,IBZ)=IBINPOP(IBX,IBY,IBZ)+1
                      NPOPT = NPOPT + 1
                      GOTO 1111
                   ENDIF
                   AMINZ = AMAXZ
                ENDDO
                ENDIF
                AMINY = AMAXY
             ENDDO
             ENDIF
             AMINX = AMAXX
          ENDDO
 1111     CONTINUE
C
      ENDIF
      WRITE(*,*)IBLOCK,"TRIAL ENERGY= ",ET," EBLOCK= ",EBLOCK," EG= ",
     &          EG, NWALKERS
C
      RETURN
      END
C ==============================================================

      SUBROUTINE TRIALF(X1,Y1,Z1,X2,Y2,Z2,PSI,EL,FQ1,FQ2)
C
C     EVALUATE TRIAL PSI VALUE, QUANTUM FORCE AND LOCAL ENERGY
C     FOR ELECTRON I AND WALKER K 
C
C---------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.00D+00,ONE=1.00D+00,TWO=2.00D+00)
      DIMENSION FQ1(3),FQ2(3)
C
      R1SQ = X1*X1+Y1*Y1+Z1*Z1
      R2SQ = X2*X2+Y2*Y2+Z2*Z2
      R1  = DSQRT(R1SQ)
      R2  = DSQRT(R2SQ)
      X12 = X1-X2
      Y12 = Y1-Y2
      Z12 = Z1-Z2 
      R12SQ = X12*x12+y12*y12+z12*z12
      R12 = DSQRT(R12SQ)
C
C     WAVE FUNCTION
C
      V = -TWO*(R1+R2) + 0.500D+00*R12
      PSI = DEXP(V)
C
C     QUANTUM FORCE
C
      F1 = 4.0D+00/R1
      F2 = 4.0D+00/R2
      FQ1(1) = X12/R12 - F1*X1
      FQ1(2) = Y12/R12 - F1*Y1
      FQ1(3) = Z12/R12 - F1*Z1
      FQ2(1) =-X12/R12 - F2*X2
      FQ2(2) =-Y12/R12 - F2*Y2
      FQ2(3) =-Z12/R12 - F2*Z2
C
C     LOCAL ENERGY
C
      T1 = X12*X1+Y12*Y1+Z12*Z1
      T2 = X12*X2+Y12*Y2+Z12*Z2
      EL  = -17.0D+00/4.0D+00 + (T1/R1-T2/R2)/R12
C
      RETURN
      END
C ==============================================================

      SUBROUTINE XINIT(X1,X2,XMIN,XMAX,MWALK,MAXWLK,ISEED)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X1(MWALK,3),X2(MWALK,3)
C
C     INITIALIZE RANDOM NUMBER GENERATOR AND WALKER POSITIONS
C
C---------------------------------------------------------------
      DUM=URAN(ISEED)
      XRANGE = XMAX-XMIN
      DO I=1,3
         DO K=1,MAXWLK
c            X(K,I) = ROULET(ISEED) * XRANGE + XMIN
            X1(K,I) = URAN(0) * XRANGE + XMIN
            X2(K,I) = URAN(0) * XRANGE + XMIN
         ENDDO
      ENDDO
      RETURN
      END
C ==============================================================

      FUNCTION BOXMLR(ISEED,VAR)
C
C     1D BOX-MULLER TRANSFORM
C
C---------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TWO=2.00D+00,TWOPI=6.283185307D+00)
c      U1 = ROULET(ISEED)+0.1D-10
c      U2 = ROULET(ISEED)+0.1D-10
      U1 = URAN(0)
      U2 = URAN(0)
      BOXMLR = DSQRT(-TWO*VAR*DLOG(U1))
      BOXMLR = BOXMLR * DSIN(TWOPI*U2)
      RETURN
      END
C ==============================================================

      SUBROUTINE SAVBIN(IBINPOP,NBIN,XMIN,XMAX,NPOPT)
C
C          Save 3D wavefunction histogram
C
C---------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IBINPOP(NBIN,NBIN,NBIN)
C
      OPEN(111,FILE='wfn.dat')
      TOT = DFLOAT(NPOPT)
      BINWIDTH = (XMAX-XMIN)/NBIN
C
      DO 120 I=1,NBIN
         RI = XMIN + (I-1)*BINWIDTH
         DO 130 J=1,NBIN
            RJ = XMIN + (J-1)*BINWIDTH
            DO 140 K=1,NBIN
               RK = XMIN + (K-1)*BINWIDTH
               WRITE(111,1050) RI,RJ,RK,DFLOAT(IBINPOP(I,J,K))/TOT
 140        CONTINUE
 130     CONTINUE
 120  CONTINUE
C
 1050 FORMAT(3(D20.8),D20.8)
C
      RETURN
      END
C ==============================================================

      DOUBLE PRECISION FUNCTION ROULET(ISEED)
C
C         Create a pseudo-random number in the range [0:1]
C
C---------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RR(97)
      INTEGER ISEED
      SAVE IX1,IX2,IX3,RR
C
C... parameter values for overflows at 2**24
      DATA M1, IA1, IC1/ 31104,   625,   6571/
      DATA M2, IA2, IC2/ 12960,  1741,   2731/
      DATA M3, IA3, IC3/ 14000,  1541,   2957/
      RM1 = 1./M1
      RM2 = 1./M2
C
C... initialize the shuffling vector RR to hold the random numbers 
      IF (ISEED.LT.0) THEN
          IX1 = MOD(IC1-ISEED  ,M1)
          IX1 = MOD(IA1*IX1+IC1,M1)
          IX2 = MOD(IX1, M2)
          IX1 = MOD(IA1*IX1+IC1,M1)
          IX3 = MOD(IX1, M3)
C
C...load vector RR
          DO 11 J=1, 97
             IX1 = MOD(IA1*IX1+IC1,M1)
             IX2 = MOD(IA2*IX2+IC2,M2)
             RR(J)= (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
  11      CONTINUE
          ISEED=1
      ENDIF
C
C...randomly sample vector RR
      IX3 = MOD(IA3*IX3+IC3,M3)
      J = 1 + (97*IX3)/M3
C
      IF (J.GT.97 .OR. J.LT.1)  WRITE(6,99)
  99  FORMAT(//5X,'Array size for RR violated in ROULET'/)
C
C...change interval from [0,1] to [-1,1]
c      ROULET = 2.D0*RR(J)-1.D0
      ROULET = RR(J)
C
C...replace this RR(J) with the next value in the sequence
      IX1 = MOD(IA1*IX1+IC1,M1)
      IX2 = MOD(IA2*IX2+IC2,M2)
      RR(J)= (FLOAT(IX1) + FLOAT(IX2)*RM2)*RM1
C
      RETURN
      END
C ==============================================================

      FUNCTION URAN(I)
C
C     PSEUDO RANDOM NUMBERIUS
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IA=16807,IM=2147483647,IQ=127773,IR=2836)
      SAVE IS
      AM=1.D00/DFLOAT(IM)
      IF (I.EQ.0) THEN
          IK1=IS/IQ
          IS=IA*(IS-IK1*IQ)-IR*IK1
          IF (IS.LT.0) IS=IS+IM 
          URAN=AM*IS
      ELSEIF(I.GT.0) THEN
          IS=I
      ELSE
          I=IS
      ENDIF
      RETURN
      END
C ==============================================================

