C-----|--|---------|---------|---------|---------|---------|---------|--|------|
C
C      filename: gentcf.f
C
C                   Compute a generic autocorrelation function
C                   From 'Molecular Dynamics Simulation' J.M.Haile (app.J)
C                   Modified version
C
C                      version 3.2a    3 Aug 2013    Bartosz Błasiak
C -----------------------------------------------------------------------------
C      The autocorrelation function here is normalized:
C    
C                   < A(T) A(T + t) >_T
C         C(t) =   ---------------------
C                     < A(T) A(T) >_T
C
C -----------------------------------------------------------------------------
C
      SUBROUTINE GENTCF(DFILE,NMAX,NSKIP,NORGNS,LPRINT,
     &                                                TCFDTA,NDELS)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL LPRINT
      REAL*8  TCFDTA(1000000)
Cf2py INTENT(OUT) TCFDTA
Cf2py INTENT(INOUT) NDELS
      CHARACTER*32 Dfile
C
C...enter filename that has data
      OPEN (10,FILE=DFILE,ACCESS='sequential',FORM='formatted')
C
C...set parameters
      CALL TCPARM(NMAX,NSKIP,NDELS)
C
C...read down disk thru number of initial time-steps to skip
      IF (NSKIP.GT.0) THEN 
          DO 20 I = 1, NSKIP
             READ(10,91) NDUMMY
 91          FORMAT(I6)
 20       CONTINUE
      ENDIF
C
C...loop ovel time-origins
      CALL TCLOOP(NORGNS,NDELS)
      CLOSE(10)
C
C...normalize and print time correlation function
      CALL TCNORM(NDELS,LPRINT,TCFDTA)
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE TCPARM(NMAX,NSKIP,
     &                             NDELS)
C
C          calculational parameters
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /DATER / X(1000000),NN(1000000)
      COMMON /SCALRS/ RL(3), IG(4)
      EQUIVALENCE (RL(1), TSTEP)
      EQUIVALENCE (IG(2), NLEFT), (IG(4), KINTVL)
C
C...check for sufficient time steps on disk to get to requested time
C   delay
      NLEFT = NMAX - NSKIP
      IF (NDELS.GT.NLEFT) NDELS = NLEFT-1
C
C...move (NDELS+1) sets of data from disc to fast memory
      NDELP1 = NDELS+1
      DO 100 I=1, NDELP1
         READ(10,93) NN(I), X(I)
  93     FORMAT(I6,D13.5)
 100  CONTINUE
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|
 
      SUBROUTINE TCLOOP(NORGNS,NDELS)
C
C          Loop over time_origins for autocorrelation
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /DATER / X(1000000),NN(1000000)
      COMMON /ANS   / TCF(1000000), NORM(1000000)
      COMMON /SCALRS/ RL(3), IG(4)
      EQUIVALENCE (IG(2), NLEFT)
      EQUIVALENCE (RL(2), XAVE),  (RL(3), XQSAV)
C
      FSUM = 0.
      FSUMSQ=0.
C
C...loop over time_origins
      DO 300 K=1, NORGNS
         XX     = X(1)
         FSUM  = FSUM + XX
         FSUMSQ= FSUMSQ + XX*XX
C
C...at intervals print reassurascence to user
         IF (MOD(NN(1),1000000) .EQ. 0) WRITE(6,95) NN(1)
  95     FORMAT(5X,'gentcf at time_origin ',I6)
C
C...check number of  time_steps remaining on disk
         NCHECK = NLEFT-K
         IF (NDELS.LE.NCHECK) NDELAY = NDELS
         IF (NDELS.GT.NCHECK) NDELAY = NCHECK
         IF (NDELAY.GT.0) THEN
C
C...loop over delay time from origin
            DO 200 J=1,NDELAY
               JP1 = J+1
               NORM(J) = NORM(J) + 1
               TCF(J) = TCF(J) + XX*X(JP1)
 200        CONTINUE
C
C...realling data in fast memory
            CALL TCSHFT(NDELS,NDELAY,NCHECK)
         ENDIF
 300  CONTINUE
C
      XAVE = FSUM/DFLOAT(NORGNS)
      XQSAV=FSUMSQ/DFLOAT(NORGNS)
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE TCSHFT(NDELS,NDELAY,NCHECK)
C
C          realling data in fast memory
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /DATER / X(1000000),NN(1000000)
C
C...sjift contents of fast memory up one slot
      DO 400 J=1,NDELAY
         JP1 = J+1
         NN(J)=NN(JP1)
         X(J) =X(JP1)
 400  CONTINUE
C
C...if there are still time_steps on disc, read in the next one
      IF (NDELS.LT.NCHECK) THEN
          NP1 = NDELAY+1
          READ(10,93) NN(NP1), X(NP1)
  93      FORMAT(I6,D13.5)
      ENDIF
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE TCNORM(NDELS,LPRINT,
     &                               TCFDTA)
C
C          Normalize and print results
C
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL LPRINT
      REAL*8  TCFDTA(1000000)
      COMMON /ANS   / TCF(1000000), NORM(1000000)
      COMMON /SCALRS/ RL(3), IG(4)
      EQUIVALENCE (IG(4), KINTVL)
      EQUIVALENCE (RL(1), TSTEP), (RL(2), XAVE),  (RL(3), XQSAV)
C
      WRITE(6,97) XAVE
  97  FORMAT(/5X,"Average = ", F8.4/)
C
C...loop over delay times
      IF (LPRINT) THEN
        DO 600 J=1,NDELS
           TCFNM = TCF(J)/ (XQSAV*DFLOAT(NORM(J)))
           TCFDTA(J) = TCFNM
           WRITE(6,99) J, TCFNM
 600    CONTINUE
      ELSE
         DO 700 J=1,NDELS
               TCFNM = TCF(J)/ (XQSAV*DFLOAT(NORM(J)))
               TCFDTA(J) = TCFNM
 700     CONTINUE
      ENDIF
C
  99  FORMAT(10X,I3,3X,1PE13.3)
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      BLOCK DATA
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /ANS   / TCF(1000000), NORM(1000000)
      DATA TCF/1000000*0.D0/
      DATA NORM/1000000*0/
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|
