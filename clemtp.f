C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE CLEMTP(NORGA,RA,CA,DA,QA,OA,
     &                  NORGB,RB,CB,DB,QB,OB,
     &                  ELMTP,A,B,C,D,E,
     &                  CCTOT,CDTOT,CQTOT,
     &                  CTTOT,DDTOT,DQTOT)
C
C -----------------------------------------------------------------------------
C
C         ELECTROSTATIC INTERACTION ENERGY FROM MULTIPOLE EXPANSION
C
C                                                     23.08.2013
C
C -----------------------------------------------------------------------------
C   Variables:
C     NORG?   - number of distributed sites
C     X?      - positions, charges, dipoles, quadrupoles and octupoles
C     ELMTP   - interaction energy
C     A-E     - R-1 to R-5 terms
C     ??TOT   - total contributions to ELMTP
C     
C   Returns:
C     ELMTP, A-E, {??TOT}
C -----------------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RA(NORGA,3),CA(NORGA),DA(NORGA,3),QA(NORGA,3,3),
     &          OA(NORGA,3,3,3),RB(NORGB,3),CB(NORGB),DB(NORGB,3),
     &          QB(NORGB,3,3),OB(NORGB,3,3,3)
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,THREE=3.0D+00,
     &           FOUR=4.0D+00,FIVE=5.0D+00)
Cf2py INTENT(OUT) ELMTP,A,B,C,D,E,CCTOT,CDTOT,CQTOT,CTTOT,DDTOT,DQTOT
C
      CC = ZERO
      CD = ZERO
      DC = ZERO
      DD = ZERO
      CQ = ZERO
      QC = ZERO 
      CT = ZERO
      TC = ZERO
      DQ = ZERO
      QD = ZERO
C
C.....loop over all pairs of origins
      DO 9000 I=1,NORGA
      DO 9000 J=1,NORGB
C
C........accumulate important sums
         SUM1 = ZERO
         SUM2 = ZERO
         SUM3 = ZERO
         SUM4 = ZERO
         SUM5 = ZERO
         SUM6 = ZERO
         SUM7 = ZERO
         SUM8 = ZERO
         SUM9 = ZERO
         RAB = ZERO
         DO 100 K=1,3
            RK = RB(J,K) - RA(I,K)
            RAB = RAB + RK * RK
            DBJK = DB(J,K)
            DAIK = DA(I,K)
            SUM1 = SUM1 + DBJK * RK
            SUM2 = SUM2 + DAIK * RK
            SUM3 = SUM3 + DAIK * DBJK
            DO 200 L=1,3
               RR = RK * (RB(J,L) - RA(I,L))
               QBJKL = QB(J,K,L)
               QAIKL = QA(I,K,L)
               SUM4 = SUM4 + QBJKL * RR
               SUM5 = SUM5 + QAIKL * RR
               SUM8 = SUM8 + QBJKL * DA(I,L) * RK
               SUM9 = SUM9 + QAIKL * DB(J,L) * RK
               DO 300 M=1,3
                  RRR = RR * (RB(J,M) - RA(I,M))
                  SUM6 = SUM6 + OB(J,K,L,M) * RRR
                  SUM7 = SUM7 + OA(I,K,L,M) * RRR
  300          CONTINUE
  200       CONTINUE
  100    CONTINUE
         RAB  = DSQRT(RAB)
         RAB3 = RAB**3
         RAB5 = RAB**5
         RAB7 = RAB**7
C
C........evaluate respective terms
         CC = CC + ( CA(I) * CB(J) ) / RAB
C
         CDIJ = - CA(I) * SUM1 / RAB3
         DCIJ =   CB(J) * SUM2 / RAB3
         CD = CD + CDIJ
         DC = DC + DCIJ
C
         DD = DD + ( SUM3 / RAB3 )
         DD = DD - (( THREE * SUM1 * SUM2 ) / RAB5)
C
         CQIJ = CA(I) * SUM4 / RAB5
         QCIJ = CB(J) * SUM5 / RAB5
         CQ = CQ + CQIJ
         QC = QC + QCIJ
C
         CTIJ = - CA(I) * SUM6 / RAB7
         TCIJ =   CB(J) * SUM7 / RAB7
         CT = CT + CTIJ
         TC = TC + TCIJ
C
         DQIJ = (- TWO * SUM8 / RAB5) + (FIVE * SUM2 * SUM4 / RAB7)
         QDIJ = (  TWO * SUM9 / RAB5) - (FIVE * SUM1 * SUM5 / RAB7)
         DQ = DQ + DQIJ
         QD = QD + QDIJ
C
 9000 CONTINUE
C
C.....collect the results
      CCTOT = CC
      CDTOT = CD + DC
      DDTOT = DD
      CQTOT = CQ + QC
      CTTOT = CT + TC
      DQTOT = DQ + QD
C
      A = CCTOT
      B = A + CDTOT
      C = B + CQTOT + DDTOT
      D = C + CTTOT + DQTOT
      E = ZERO
C
      ELMTP = D
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|
