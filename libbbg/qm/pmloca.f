C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE PMLOCA(NATOMS,NBASIS,NMOS,MAPI,SAO,VECIN,
     &                  MAXIT,CVGLOC,LPRINT,TRAN)
C
C -----------------------------------------------------------------------------
C              PIPEK-MEZEY MOLECULAR ORBITAL LOCALIZATION SCHEME
C          J. PIPEK AND P. G. MEZEY  J. CHEM. PHYS. 90, 4916 (1989),
C 
C                The algorithm is taken from Psi4 package
C                  Code adapted from C++ to Fortran 77
C
C                                    Created:      19.08.2013 Seoul
C                                    Revised:      17.07.2019 Gundelfingen
C -----------------------------------------------------------------------------
C   Variables:
C     MAPI   - list of atomic indices in bfs order
C     VECIN  - input canonical orbitals
C     CVGLOC - convergence criterion (population)
C     
C   Returns:
C     TRAN   - transformation matrix
C -----------------------------------------------------------------------------
C   Original implementation license and credit (algorithm):
C
C     Psi4: an open-source quantum chemistry software package
C
C     Copyright (c) 2007-2018 The Psi4 Developers.
C
C     Psi4 is free software; you can redistribute it and/or modify                   
C     it under the terms of the GNU Lesser General Public License as published by
C     the Free Software Foundation, version 3.
C                                                                                    
C     Psi4 is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU Lesser General Public License for more details.
C                                                                                    
C     You should have received a copy of the GNU Lesser General Public License along
C     with Psi4; if not, write to the Free Software Foundation, Inc.,
C     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
C
C -----------------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION VECIN(NMOS,NBASIS), SAO(NBASIS,NBASIS),
     &          MAPI(NBASIS)
      PARAMETER (T45DEG=3.141592653589793D+00/4.0D+00)
      PARAMETER (ZERO=0.D0,ONE=1.0D+00,TWO=2.0D+00,HALF=0.5D+00,
     &           TENM8=1.0D-08)
      LOGICAL LPRINT
c
c.....internal
      DIMENSION CS(NBASIS,NMOS)
      DIMENSION VECOUT(NBASIS,NMOS)
c     DIMENSION IORD(NMOS)
      DIMENSION NNATL(NATOMS+1)
      DIMENSION TRAN(NMOS,NMOS)
Cf2py INTENT(OUT) TRAN
c
c...initialize transformation matrix
      DO I = 1, NMOS
         TRAN(I,I) = ONE
         DO J = 1, I-1            
              TRAN(I,J) = ZERO
              TRAN(J,I) = ZERO
         END DO
      END DO
c
c.....copy input vectors C to initialize C_loc
      DO K = 1, NBASIS
      DO I = 1, NMOS
         VECOUT(K, I) = VECIN(I, K)
      END DO
      END DO
c
c.....Compute CS matrix
      DO K = 1, NBASIS
      DO I = 1, NMOS
         CS_KI = ZERO
         DO L = 1, NBASIS
            CS_KI = CS_KI + SAO(K,L)*VECOUT(L,I)
         END DO
         CS(K,I) = CS_KI
      END DO
      END DO
c
c.....Populate basis function number per atom
      CALL FILLNA(NNATL, MAPI, NBASIS, NATOMS)
c
c.....Population metric (initial)
      DMET = ZERO
      DO I = 1, NMOS
         DO K = 1, NATOMS
            KOF= NNATL(K)
            KM = NNATL(K+1) - KOF
            PA = CMYDOT(CS,VECOUT,KOF,KM,I,I,NBASIS,NMOS)
            DMET = DMET + PA*PA
         END DO
      END DO
c
      DMET_OLD = DMET
      write(*,*) "Init Metric: ", DMET
c
c.....Iteration cycles
      DO ITER = 1, MAXIT
         DO I2 = 1, NMOS - 1
         DO J2 = I2 + 1, NMOS 
c           I = IORD(I2)
c           J = IORD(J2)
            I = I2
            J = J2
c...........Determine rotation angle
            A = ZERO 
            B = ZERO 
            C = ZERO 
            DO K = 1, NATOMS
               KOF= NNATL(K)
               KM = NNATL(K+1) - KOF
               AII= CMYDOT(CS,VECOUT,KOF,KM,I,I,NBASIS,NMOS)
               AJJ= CMYDOT(CS,VECOUT,KOF,KM,J,J,NBASIS,NMOS)
               AIJ= CMYDOT(CS,VECOUT,KOF,KM,I,J,NBASIS,NMOS)
     &            + CMYDOT(CS,VECOUT,KOF,KM,J,I,NBASIS,NMOS)
               AIJ= HALF * AIJ
               AD = AII - AJJ
               AO = TWO * AIJ
               A  = A + AD * AD
               B  = B + AO * AO
               C  = C + AD * AO
            END DO
c
            HD = A - B
            HO = TWO * C
            T  = HALF * DATAN2(HO, HD + DSQRT(HD*HD + HO*HO))
c
c...........check for trivial rotation angle
            IF (DABS(T).LT.TENM8) THEN
                O0 = ZERO
                O1 = ZERO
                DO K = 1, NATOMS
                   KOF= NNATL(K)
                   KM = NNATL(K+1) - KOF
                   AII= CMYDOT(CS,VECOUT,KOF,KM,I,I,NBASIS,NMOS) 
                   AJJ= CMYDOT(CS,VECOUT,KOF,KM,J,J,NBASIS,NMOS)
                   AIJ= CMYDOT(CS,VECOUT,KOF,KM,I,J,NBASIS,NMOS)
     &                + CMYDOT(CS,VECOUT,KOF,KM,J,I,NBASIS,NMOS)
                   AIJ= HALF * AIJ
                   O0 = O0 + AIJ * AIJ
                   DA = AJJ - AII
                   O1 = O1 + DA*DA
                END DO
                O1 = O1 * 0.25D+00
                IF (O1.LT.O0) T = T45DEG
            END IF
c
c...........Plane rotation
            CC = DCOS(T)
            SS = DSIN(T)
            CALL CMYROT(CS,CC,SS,NBASIS,NMOS,I,J,.TRUE.)
            CALL CMYROT(VECOUT,CC,SS,NBASIS,NMOS,I,J,.TRUE.)
            CALL CMYROT(TRAN,CC,SS,NMOS,NMOS,I,J,.TRUE.)
         END DO 
         END DO
c
c........Population metric
         DMET = ZERO
         DO I = 1, NMOS                                    
            DO K = 1, NATOMS
               KOF= NNATL(K)
               KM = NNATL(K+1) - KOF
               PA = CMYDOT(CS,VECOUT,KOF,KM,I,I,NBASIS,NMOS)
               DMET = DMET + PA*PA
            END DO
         END DO
         write(*,*) "Metric: ", DMET, ITER
c
c........Convergence criterion
         CONV = DABS(DMET - DMET_OLD) / DABS(DMET_OLD)
         DMET_OLD = DMET
c
c........Converged?
         IF (CONV.LT.CVGLOC) GO TO 500
      END DO
C
 500  CONTINUE
c
c.....Print out
      IF (LPRINT) THEN
          WRITE (6,1000) ITER
      ENDIF
C
 1000 FORMAT(10X,'LOCALIZATION CONVERGED IN',I6,' ITERATIONS')
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE CMYROT(U,C,S,NA,NB,I,J,COL)
C
C     Internal adaptation of CROT
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL COL
      DIMENSION U(NA,NB)
c.....Rotation of column pair
      IF (COL) THEN
          DO K = 1, NA                 
             TEMP   = C*U(K,I) + S*U(K,J) 
             U(K,J) =-S*U(K,I) + C*U(K,J)
             U(K,I) = TEMP
          END DO
c.....Rotation of row pair
      ELSE
          DO K = 1, NB
             TEMP   = C*U(I,K) + S*U(J,K) 
             U(J,K) =-S*U(I,K) + C*U(J,K)
             U(I,K) = TEMP
          END DO
      END IF
c
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      DOUBLE PRECISION FUNCTION CMYDOT(C1,C2,NSTART,NGO,I,J,NBASIS,NMOS)
C
C     Internal adaptation of CDOT
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION C1(NBASIS, NMOS)
      DIMENSION C2(NBASIS, NMOS)
      CMYDOT = 0.0D+00
C
      DO N = NSTART, NSTART+NGO-1
         CMYDOT = CMYDOT + C1(N, I) * C2(N, J)
      END DO
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE FILLNA(NNATL, MAPI, NBASIS, NATOMS)
C
C     Fills in the starting values of basis function numbers
C     Per atom. Last element is NBASIS+1.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION NNATL(NATOMS+1), MAPI(NBASIS)
C
      K = 1
      IPREV = -1
C
      DO N = 1, NBASIS
         INEXT = MAPI(N)
         IF (INEXT.NE.IPREV) THEN
             NNATL(K) = N
             K = K + 1
         END IF
         IPREV = INEXT
      END DO
C
      NNATL(NATOMS+1) = NBASIS+1
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|
