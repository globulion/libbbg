C-----|--|---------|---------|---------|---------|---------|---------|--|------|
C
C      filename: efprot.f
C
C                              THE ROTATION 
C                OF EFFECTIVE FRAGMENT POTENTIAL PARAMETERS
C
C                               version 0.0a    28 Oct 2013    Bartosz Błasiak
C NOTES:
C
C            R - Rotation matrix
C      
C                XX XY XZ     1  2  3
C                YX YY YZ     4  5  6 
C                ZX ZY ZZ     7  8  9
C
C            Turn of indices in AOs (PyQuante convention):
C
C            P - shell:       X  Y  Z
C            D - shell:       XX YY ZZ XY XZ YZ
C -----------------------------------------------------------------------------
      SUBROUTINE VECROT(NMOS,NBASIS,VEC,ROT,ITYP)
C
C        ROTATE THE WAVE FUNCTION PARAMETERS
C
C -----------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION VEC(NMOS,NBASIS),ROT(3,3),ITYP(NBASIS)
Cf2py INTENT(IN,OUT) VEC
      RXX = ROT(1,1)
      RYY = ROT(2,2)
      RZZ = ROT(3,3)
      RXY = ROT(1,2)
      RYX = ROT(2,1)
      RXZ = ROT(1,3)
      RZX = ROT(3,1)
      RYZ = ROT(2,3)
      RZY = ROT(3,2)
C
      NP = 0
      ND = 0
      DO 1010 K=1,NBASIS
         ITYPK = ITYP(K)
C
C        P-SHELL
C
         IF (ITYPK.EQ.1) THEN
             NP = NP + 1
             IF ((MOD((NP+2),3)).EQ.0) THEN
                DO 2010 I=1,NMOS
                   XO = VEC(I,K  )
                   YO = VEC(I,K+1)
                   ZO = VEC(I,K+2)
C
                   XN = RXX * XO + RYX * YO + RZX * ZO
                   YN = RXY * XO + RYY * YO + RZY * ZO
                   ZN = RXZ * XO + RYZ * YO + RZZ * ZO
C
                   VEC(I,K  ) = XN
                   VEC(I,K+1) = YN
                   VEC(I,K+2) = ZN
 2010           CONTINUE
             ENDIF
C
C        D-SHELL
C
         ELSE IF (ITYPK.EQ.2) THEN
             ND = ND + 1
             IF ((MOD((NP+5),6)).EQ.0) THEN
                DO 3010 I=1,NMOS
                   XXO = VEC(I,K  )
                   YYO = VEC(I,K+1)
                   ZZO = VEC(I,K+2)
                   XYO = VEC(I,K+3)
                   XZO = VEC(I,K+4)
                   YZO = VEC(I,K+5)
C
                   XXN = RXX * RXX * XXO + 
     &                   RXX * RYX * XYO + 
     &                   RXX * RZX * XZO +
     &                   RYX * RXX * XYO +
     &                   RYX * RYX * YYO +
     &                   RYX * RZX * YZO +
     &                   RZX * RXX * XZO +
     &                   RZX * RYX * YZO +
     &                   RZX * RZX * ZZO
C
                   YYN = RXY * RXY * XXO + 
     &                   RXY * RYY * XYO + 
     &                   RXY * RZY * XZO +
     &                   RYY * RXY * XYO +
     &                   RYY * RYY * YYO +
     &                   RYY * RZY * YZO +
     &                   RZY * RXY * XZO +
     &                   RZY * RYY * YZO +
     &                   RZY * RZY * ZZO
C
                   ZZN = RXZ * RXZ * XXO + 
     &                   RXZ * RYZ * XYO + 
     &                   RXZ * RZZ * XZO +
     &                   RYZ * RXZ * XYO +
     &                   RYZ * RYZ * YYO +
     &                   RYZ * RZZ * YZO +
     &                   RZZ * RXZ * XZO +
     &                   RZZ * RYZ * YZO +
     &                   RZZ * RZZ * ZZO
C
                   XYN = RXX * RXY * XXO + 
     &                   RXX * RYY * XYO + 
     &                   RXX * RZY * XZO +
     &                   RYX * RXY * XYO +
     &                   RYX * RYY * YYO +
     &                   RYX * RZY * YZO +
     &                   RZX * RXY * XZO +
     &                   RZX * RYY * YZO +
     &                   RZX * RZY * ZZO
C
                   XZN = RXX * RXZ * XXO + 
     &                   RXX * RYZ * XYO + 
     &                   RXX * RZZ * XZO +
     &                   RYX * RXZ * XYO +
     &                   RYX * RYZ * YYO +
     &                   RYX * RZZ * YZO +
     &                   RZX * RXZ * XZO +
     &                   RZX * RYZ * YZO +
     &                   RZX * RZZ * ZZO
C
                   YZN = RXY * RXZ * XXO + 
     &                   RXY * RYZ * XYO + 
     &                   RXY * RZZ * XZO +
     &                   RYY * RXZ * XYO +
     &                   RYY * RYZ * YYO +
     &                   RYY * RZZ * YZO +
     &                   RZY * RXZ * XZO +
     &                   RZY * RYZ * YZO +
     &                   RZY * RZZ * ZZO
C
                   VEC(I,K  ) = XXN
                   VEC(I,K+1) = YYN
                   VEC(I,K+2) = ZZN
                   VEC(I,K+3) = XYN
                   VEC(I,K+4) = XZN
                   VEC(I,K+5) = YZN
 3010           CONTINUE
             ENDIF
         ENDIF
C
 1010 CONTINUE
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|
      SUBROUTINE VC1ROT(NMODES,NMOS,NBASIS,VEC1,ROT,ITYP)
C
C          ROTATE VECL1 TENSORS
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION VEC1(NMODES,NMOS,NBASIS),ROT(3,3),ITYP(NBASIS)
Cf2py INTENT(IN,OUT) VEC1
      RXX = ROT(1,1)
      RYY = ROT(2,2)
      RZZ = ROT(3,3)
      RXY = ROT(1,2)
      RYX = ROT(2,1)
      RXZ = ROT(1,3)
      RZX = ROT(3,1)
      RYZ = ROT(2,3)
      RZY = ROT(3,2)
C
      NP = 0
      ND = 0
      DO 1010 K=1,NBASIS
         ITYPK = ITYP(K)
C
C        P-SHELL
C
         IF (ITYPK.EQ.1) THEN
             NP = NP + 1
             IF ((MOD((NP+2),3)).EQ.0) THEN
                DO 2010 M=1,NMODES
                DO 2010 I=1,NMOS
                   XO = VEC1(M,I,K  )
                   YO = VEC1(M,I,K+1)
                   ZO = VEC1(M,I,K+2)
C
                   XN = RXX * XO + RYX * YO + RZX * ZO
                   YN = RXY * XO + RYY * YO + RZY * ZO
                   ZN = RXZ * XO + RYZ * YO + RZZ * ZO
C
                   VEC1(M,I,K  ) = XN
                   VEC1(M,I,K+1) = YN
                   VEC1(M,I,K+2) = ZN
 2010           CONTINUE
             ENDIF
C
C        D-SHELL
C
         ELSE IF (ITYPK.EQ.2) THEN
             ND = ND + 1
             IF ((MOD((NP+5),6)).EQ.0) THEN
                DO 3010 M=1,NMODES
                DO 3010 I=1,NMOS
                   XXO = VEC1(M,I,K  )
                   YYO = VEC1(M,I,K+1)
                   ZZO = VEC1(M,I,K+2)
                   XYO = VEC1(M,I,K+3)
                   XZO = VEC1(M,I,K+4)
                   YZO = VEC1(M,I,K+5)
C
                   XXN = RXX * RXX * XXO + 
     &                   RXX * RYX * XYO + 
     &                   RXX * RZX * XZO +
     &                   RYX * RXX * XYO +
     &                   RYX * RYX * YYO +
     &                   RYX * RZX * YZO +
     &                   RZX * RXX * XZO +
     &                   RZX * RYX * YZO +
     &                   RZX * RZX * ZZO
C
                   YYN = RXY * RXY * XXO + 
     &                   RXY * RYY * XYO + 
     &                   RXY * RZY * XZO +
     &                   RYY * RXY * XYO +
     &                   RYY * RYY * YYO +
     &                   RYY * RZY * YZO +
     &                   RZY * RXY * XZO +
     &                   RZY * RYY * YZO +
     &                   RZY * RZY * ZZO
C
                   ZZN = RXZ * RXZ * XXO + 
     &                   RXZ * RYZ * XYO + 
     &                   RXZ * RZZ * XZO +
     &                   RYZ * RXZ * XYO +
     &                   RYZ * RYZ * YYO +
     &                   RYZ * RZZ * YZO +
     &                   RZZ * RXZ * XZO +
     &                   RZZ * RYZ * YZO +
     &                   RZZ * RZZ * ZZO
C
                   XYN = RXX * RXY * XXO + 
     &                   RXX * RYY * XYO + 
     &                   RXX * RZY * XZO +
     &                   RYX * RXY * XYO +
     &                   RYX * RYY * YYO +
     &                   RYX * RZY * YZO +
     &                   RZX * RXY * XZO +
     &                   RZX * RYY * YZO +
     &                   RZX * RZY * ZZO
C
                   XZN = RXX * RXZ * XXO + 
     &                   RXX * RYZ * XYO + 
     &                   RXX * RZZ * XZO +
     &                   RYX * RXZ * XYO +
     &                   RYX * RYZ * YYO +
     &                   RYX * RZZ * YZO +
     &                   RZX * RXZ * XZO +
     &                   RZX * RYZ * YZO +
     &                   RZX * RZZ * ZZO
C
                   YZN = RXY * RXZ * XXO + 
     &                   RXY * RYZ * XYO + 
     &                   RXY * RZZ * XZO +
     &                   RYY * RXZ * XYO +
     &                   RYY * RYZ * YYO +
     &                   RYY * RZZ * YZO +
     &                   RZY * RXZ * XZO +
     &                   RZY * RYZ * YZO +
     &                   RZY * RZZ * ZZO
C
                   VEC1(M,I,K  ) = XXN
                   VEC1(M,I,K+1) = YYN
                   VEC1(M,I,K+2) = ZZN
                   VEC1(M,I,K+3) = XYN
                   VEC1(M,I,K+4) = XZN
                   VEC1(M,I,K+5) = YZN
 3010           CONTINUE
             ENDIF
         ENDIF
C
 1010 CONTINUE
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE ROTDMA(DIP,QAD,OCT,ROT,NDMA)
C
C -----------------------------------------------------------------------------
C
C               ROTATE THE DISTRIBUTED MULTIPOLE MOMENT TENSORS
C 
C              Bartosz Blasiak                        08.11.2013
C
C -----------------------------------------------------------------------------
C
C   Description:
C     Performs the unitary transformations:
C
C       D_a   = D_a'     * R_a'a
C       Q_ab  = Q_a'b'   * R_a'a * R_b'b
C       O_abc = O_a'b'c' * R_a'a * R_b'b * R_c'c
C
C   Input variables:
C     NDMA    - array of numbers of distributed electrostatic sites
C     DIP,
C     QAD,OCT - non-scalar distributed multipoles
C     ROT     - unitary rotation matrix
C     
C   Returns:
C     DIP,QAD,OCT - rotated moments
C
C   Notes:
C     The reduced format of tensor storage is used:
C
C     DIP(i) X   Y   Z
C            1   2   3
C     QAD(i) XX  YY  ZZ  XY  XZ  YZ
C            1   2   3   4   5   6
C     OCT(i) XXX YYY ZZZ XXY XXZ XYY YYZ XZZ YZZ XYZ
C            1   2   3   4   5   6   7   8   9   10
C -----------------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DIP(NDMA,3),QAD(NDMA,6),OCT(NDMA,10),ROT(3,3)
      PARAMETER (ZERO=0.0D+00)
Cf2py INTENT(IN,OUT) DIP, QAD, OCT
C
      RXX = ROT(1,1)
      RYY = ROT(2,2)
      RZZ = ROT(3,3)
      RXY = ROT(1,2)
      RYX = ROT(2,1)
      RXZ = ROT(1,3)
      RZX = ROT(3,1)
      RYZ = ROT(2,3)
      RZY = ROT(3,2)
C
      DO 11 N=1,NDMA
C
C       UNPACK
C
         DX   = DIP(N,1)
         DY   = DIP(N,2)
         DZ   = DIP(N,3)
C
         QXX  = QAD(N,1)
         QYY  = QAD(N,2)
         QZZ  = QAD(N,3)
         QXY  = QAD(N,4)
         QXZ  = QAD(N,5)
         QYZ  = QAD(N,6)
C
         OXXX = OCT(N,1)
         OYYY = OCT(N,2)
         OZZZ = OCT(N,3)
         OXXY = OCT(N,4)
         OXXZ = OCT(N,5)
         OXYY = OCT(N,6)
         OYYZ = OCT(N,7)
         OXZZ = OCT(N,8)
         OYZZ = OCT(N,9)
         OXYZ = OCT(N,10)
C
C        ROTATE DIPOLES
C
         rDX = RXX * DX + RYX * DY + RZX * DZ
         rDY = RXY * DX + RYY * DY + RZY * DZ
         rDZ = RXZ * DX + RYZ * DY + RZZ * DZ
C
C        ROTATE QUADRUPOLES
C
         rQXX = RXX * RXX * QXX +
     &          RXX * RYX * QXY +
     &          RXX * RZX * QXZ +
     &          RYX * RXX * QXY +
     &          RYX * RYX * QYY +
     &          RYX * RZX * QYZ +
     &          RZX * RXX * QXZ +
     &          RZX * RYX * QYZ +
     &          RZX * RZX * QZZ
C
         rQYY = RXY * RXY * QXX +
     &          RXY * RYY * QXY +
     &          RXY * RZY * QXZ +
     &          RYY * RXY * QXY +
     &          RYY * RYY * QYY +
     &          RYY * RZY * QYZ +
     &          RZY * RXY * QXZ +
     &          RZY * RYY * QYZ +
     &          RZY * RZY * QZZ
C
         rQZZ = RXZ * RXZ * QXX +
     &          RXZ * RYZ * QXY +
     &          RXZ * RZZ * QXZ +
     &          RYZ * RXZ * QXY +
     &          RYZ * RYZ * QYY +
     &          RYZ * RZZ * QYZ +
     &          RZZ * RXZ * QXZ +
     &          RZZ * RYZ * QYZ +
     &          RZZ * RZZ * QZZ
C
         rQXY = RXX * RXY * QXX +
     &          RXX * RYY * QXY +
     &          RXX * RZY * QXZ +
     &          RYX * RXY * QXY +
     &          RYX * RYY * QYY +
     &          RYX * RZY * QYZ +
     &          RZX * RXY * QXZ +
     &          RZX * RYY * QYZ +
     &          RZX * RZY * QZZ
C
         rQXZ = RXX * RXZ * QXX +
     &          RXX * RYZ * QXY +
     &          RXX * RZZ * QXZ +
     &          RYX * RXZ * QXY +
     &          RYX * RYZ * QYY +
     &          RYX * RZZ * QYZ +
     &          RZX * RXZ * QXZ +
     &          RZX * RYZ * QYZ +
     &          RZX * RZZ * QZZ
C
         rQYZ = RXY * RXZ * QXX + 
     &          RXY * RYZ * QXY + 
     &          RXY * RZZ * QXZ +
     &          RYY * RXZ * QXY +
     &          RYY * RYZ * QYY +
     &          RYY * RZZ * QYZ +
     &          RZY * RXZ * QXZ +
     &          RZY * RYZ * QYZ +
     &          RZY * RZZ * QZZ
C
C        ROTATE OCTUPOLES
C
         rOXXX = RXX * RXX * RXX * OXXX +
     &           RXX * RXX * RYX * OXXY +
     &           RXX * RXX * RZX * OXXZ +
     &           RXX * RYX * RXX * OXXY +
     &           RXX * RYX * RYX * OXYY +
     &           RXX * RYX * RZX * OXYZ +
     &           RXX * RZX * RXX * OXXZ +
     &           RXX * RZX * RYX * OXYZ +
     &           RXX * RZX * RZX * OXZZ +
     &           RYX * RXX * RXX * OXXY +
     &           RYX * RXX * RYX * OXYY +
     &           RYX * RXX * RZX * OXYZ +
     &           RYX * RYX * RXX * OXYY +
     &           RYX * RYX * RYX * OYYY +
     &           RYX * RYX * RZX * OYYZ +
     &           RYX * RZX * RXX * OXYZ +
     &           RYX * RZX * RYX * OYYZ +
     &           RYX * RZX * RZX * OYZZ +
     &           RZX * RXX * RXX * OXXZ +
     &           RZX * RXX * RYX * OXYZ +
     &           RZX * RXX * RZX * OXZZ +
     &           RZX * RYX * RXX * OXYZ +
     &           RZX * RYX * RYX * OYYZ +
     &           RZX * RYX * RZX * OYZZ +
     &           RZX * RZX * RXX * OXZZ +
     &           RZX * RZX * RYX * OYZZ +
     &           RZX * RZX * RZX * OZZZ
C
         rOYYY = RXY * RXY * RXY * OXXX +
     &           RXY * RXY * RYY * OXXY +
     &           RXY * RXY * RZY * OXXZ +
     &           RXY * RYY * RXY * OXXY +
     &           RXY * RYY * RYY * OXYY +
     &           RXY * RYY * RZY * OXYZ +
     &           RXY * RZY * RXY * OXXZ +
     &           RXY * RZY * RYY * OXYZ +
     &           RXY * RZY * RZY * OXZZ +
     &           RYY * RXY * RXY * OXXY +
     &           RYY * RXY * RYY * OXYY +
     &           RYY * RXY * RZY * OXYZ +
     &           RYY * RYY * RXY * OXYY +
     &           RYY * RYY * RYY * OYYY +
     &           RYY * RYY * RZY * OYYZ +
     &           RYY * RZY * RXY * OXYZ +
     &           RYY * RZY * RYY * OYYZ +
     &           RYY * RZY * RZY * OYZZ +
     &           RZY * RXY * RXY * OXXZ +
     &           RZY * RXY * RYY * OXYZ +
     &           RZY * RXY * RZY * OXZZ +
     &           RZY * RYY * RXY * OXYZ +
     &           RZY * RYY * RYY * OYYZ +
     &           RZY * RYY * RZY * OYZZ +
     &           RZY * RZY * RXY * OXZZ +
     &           RZY * RZY * RYY * OYZZ +
     &           RZY * RZY * RZY * OZZZ
C
         rOZZZ = RXZ * RXZ * RXZ * OXXX +
     &           RXZ * RXZ * RYZ * OXXY +
     &           RXZ * RXZ * RZZ * OXXZ +
     &           RXZ * RYZ * RXZ * OXXY +
     &           RXZ * RYZ * RYZ * OXYY +
     &           RXZ * RYZ * RZZ * OXYZ +
     &           RXZ * RZZ * RXZ * OXXZ +
     &           RXZ * RZZ * RYZ * OXYZ +
     &           RXZ * RZZ * RZZ * OXZZ +
     &           RYZ * RXZ * RXZ * OXXY +
     &           RYZ * RXZ * RYZ * OXYY +
     &           RYZ * RXZ * RZZ * OXYZ +
     &           RYZ * RYZ * RXZ * OXYY +
     &           RYZ * RYZ * RYZ * OYYY +
     &           RYZ * RYZ * RZZ * OYYZ +
     &           RYZ * RZZ * RXZ * OXYZ +
     &           RYZ * RZZ * RYZ * OYYZ +
     &           RYZ * RZZ * RZZ * OYZZ +
     &           RZZ * RXZ * RXZ * OXXZ +
     &           RZZ * RXZ * RYZ * OXYZ +
     &           RZZ * RXZ * RZZ * OXZZ +
     &           RZZ * RYZ * RXZ * OXYZ +
     &           RZZ * RYZ * RYZ * OYYZ +
     &           RZZ * RYZ * RZZ * OYZZ +
     &           RZZ * RZZ * RXZ * OXZZ +
     &           RZZ * RZZ * RYZ * OYZZ +
     &           RZZ * RZZ * RZZ * OZZZ
C
         rOXXY = RXX * RXX * RXY * OXXX +
     &           RXX * RXX * RYY * OXXY +
     &           RXX * RXX * RZY * OXXZ +
     &           RXX * RYX * RXY * OXXY +
     &           RXX * RYX * RYY * OXYY +
     &           RXX * RYX * RZY * OXYZ +
     &           RXX * RZX * RXY * OXXZ +
     &           RXX * RZX * RYY * OXYZ +
     &           RXX * RZX * RZY * OXZZ +
     &           RYX * RXX * RXY * OXXY +
     &           RYX * RXX * RYY * OXYY +
     &           RYX * RXX * RZY * OXYZ +
     &           RYX * RYX * RXY * OXYY +
     &           RYX * RYX * RYY * OYYY +
     &           RYX * RYX * RZY * OYYZ +
     &           RYX * RZX * RXY * OXYZ +
     &           RYX * RZX * RYY * OYYZ +
     &           RYX * RZX * RZY * OYZZ +
     &           RZX * RXX * RXY * OXXZ +
     &           RZX * RXX * RYY * OXYZ +
     &           RZX * RXX * RZY * OXZZ +
     &           RZX * RYX * RXY * OXYZ +
     &           RZX * RYX * RYY * OYYZ +
     &           RZX * RYX * RZY * OYZZ +
     &           RZX * RZX * RXY * OXZZ +
     &           RZX * RZX * RYY * OYZZ +
     &           RZX * RZX * RZY * OZZZ
C
         rOXXZ = RXX * RXX * RXZ * OXXX +
     &           RXX * RXX * RYZ * OXXY +
     &           RXX * RXX * RZZ * OXXZ +
     &           RXX * RYX * RXZ * OXXY +
     &           RXX * RYX * RYZ * OXYY +
     &           RXX * RYX * RZZ * OXYZ +
     &           RXX * RZX * RXZ * OXXZ +
     &           RXX * RZX * RYZ * OXYZ +
     &           RXX * RZX * RZZ * OXZZ +
     &           RYX * RXX * RXZ * OXXY +
     &           RYX * RXX * RYZ * OXYY +
     &           RYX * RXX * RZZ * OXYZ +
     &           RYX * RYX * RXZ * OXYY +
     &           RYX * RYX * RYZ * OYYY +
     &           RYX * RYX * RZZ * OYYZ +
     &           RYX * RZX * RXZ * OXYZ +
     &           RYX * RZX * RYZ * OYYZ +
     &           RYX * RZX * RZZ * OYZZ +
     &           RZX * RXX * RXZ * OXXZ +
     &           RZX * RXX * RYZ * OXYZ +
     &           RZX * RXX * RZZ * OXZZ +
     &           RZX * RYX * RXZ * OXYZ +
     &           RZX * RYX * RYZ * OYYZ +
     &           RZX * RYX * RZZ * OYZZ +
     &           RZX * RZX * RXZ * OXZZ +
     &           RZX * RZX * RYZ * OYZZ +
     &           RZX * RZX * RZZ * OZZZ
C
         rOXYY = RXX * RXY * RXY * OXXX +
     &           RXX * RXY * RYY * OXXY +
     &           RXX * RXY * RZY * OXXZ +
     &           RXX * RYY * RXY * OXXY +
     &           RXX * RYY * RYY * OXYY +
     &           RXX * RYY * RZY * OXYZ +
     &           RXX * RZY * RXY * OXXZ +
     &           RXX * RZY * RYY * OXYZ +
     &           RXX * RZY * RZY * OXZZ +
     &           RYX * RXY * RXY * OXXY +
     &           RYX * RXY * RYY * OXYY +
     &           RYX * RXY * RZY * OXYZ +
     &           RYX * RYY * RXY * OXYY +
     &           RYX * RYY * RYY * OYYY +
     &           RYX * RYY * RZY * OYYZ +
     &           RYX * RZY * RXY * OXYZ +
     &           RYX * RZY * RYY * OYYZ +
     &           RYX * RZY * RZY * OYZZ +
     &           RZX * RXY * RXY * OXXZ +
     &           RZX * RXY * RYY * OXYZ +
     &           RZX * RXY * RZY * OXZZ +
     &           RZX * RYY * RXY * OXYZ +
     &           RZX * RYY * RYY * OYYZ +
     &           RZX * RYY * RZY * OYZZ +
     &           RZX * RZY * RXY * OXZZ +
     &           RZX * RZY * RYY * OYZZ +
     &           RZX * RZY * RZY * OZZZ
C
         rOYYZ = RXY * RXY * RXZ * OXXX +
     &           RXY * RXY * RYZ * OXXY +
     &           RXY * RXY * RZZ * OXXZ +
     &           RXY * RYY * RXZ * OXXY +
     &           RXY * RYY * RYZ * OXYY +
     &           RXY * RYY * RZZ * OXYZ +
     &           RXY * RZY * RXZ * OXXZ +
     &           RXY * RZY * RYZ * OXYZ +
     &           RXY * RZY * RZZ * OXZZ +
     &           RYY * RXY * RXZ * OXXY +
     &           RYY * RXY * RYZ * OXYY +
     &           RYY * RXY * RZZ * OXYZ +
     &           RYY * RYY * RXZ * OXYY +
     &           RYY * RYY * RYZ * OYYY +
     &           RYY * RYY * RZZ * OYYZ +
     &           RYY * RZY * RXZ * OXYZ +
     &           RYY * RZY * RYZ * OYYZ +
     &           RYY * RZY * RZZ * OYZZ +
     &           RZY * RXY * RXZ * OXXZ +
     &           RZY * RXY * RYZ * OXYZ +
     &           RZY * RXY * RZZ * OXZZ +
     &           RZY * RYY * RXZ * OXYZ +
     &           RZY * RYY * RYZ * OYYZ +
     &           RZY * RYY * RZZ * OYZZ +
     &           RZY * RZY * RXZ * OXZZ +
     &           RZY * RZY * RYZ * OYZZ +
     &           RZY * RZY * RZZ * OZZZ
C
         rOXZZ = RXX * RXZ * RXZ * OXXX +
     &           RXX * RXZ * RYZ * OXXY +
     &           RXX * RXZ * RZZ * OXXZ +
     &           RXX * RYZ * RXZ * OXXY +
     &           RXX * RYZ * RYZ * OXYY +
     &           RXX * RYZ * RZZ * OXYZ +
     &           RXX * RZZ * RXZ * OXXZ +
     &           RXX * RZZ * RYZ * OXYZ +
     &           RXX * RZZ * RZZ * OXZZ +
     &           RYX * RXZ * RXZ * OXXY +
     &           RYX * RXZ * RYZ * OXYY +
     &           RYX * RXZ * RZZ * OXYZ +
     &           RYX * RYZ * RXZ * OXYY +
     &           RYX * RYZ * RYZ * OYYY +
     &           RYX * RYZ * RZZ * OYYZ +
     &           RYX * RZZ * RXZ * OXYZ +
     &           RYX * RZZ * RYZ * OYYZ +
     &           RYX * RZZ * RZZ * OYZZ +
     &           RZX * RXZ * RXZ * OXXZ +
     &           RZX * RXZ * RYZ * OXYZ +
     &           RZX * RXZ * RZZ * OXZZ +
     &           RZX * RYZ * RXZ * OXYZ +
     &           RZX * RYZ * RYZ * OYYZ +
     &           RZX * RYZ * RZZ * OYZZ +
     &           RZX * RZZ * RXZ * OXZZ +
     &           RZX * RZZ * RYZ * OYZZ +
     &           RZX * RZZ * RZZ * OZZZ
C
         rOYZZ = RXY * RXZ * RXZ * OXXX +
     &           RXY * RXZ * RYZ * OXXY +
     &           RXY * RXZ * RZZ * OXXZ +
     &           RXY * RYZ * RXZ * OXXY +
     &           RXY * RYZ * RYZ * OXYY +
     &           RXY * RYZ * RZZ * OXYZ +
     &           RXY * RZZ * RXZ * OXXZ +
     &           RXY * RZZ * RYZ * OXYZ +
     &           RXY * RZZ * RZZ * OXZZ +
     &           RYY * RXZ * RXZ * OXXY +
     &           RYY * RXZ * RYZ * OXYY +
     &           RYY * RXZ * RZZ * OXYZ +
     &           RYY * RYZ * RXZ * OXYY +
     &           RYY * RYZ * RYZ * OYYY +
     &           RYY * RYZ * RZZ * OYYZ +
     &           RYY * RZZ * RXZ * OXYZ +
     &           RYY * RZZ * RYZ * OYYZ +
     &           RYY * RZZ * RZZ * OYZZ +
     &           RZY * RXZ * RXZ * OXXZ +
     &           RZY * RXZ * RYZ * OXYZ +
     &           RZY * RXZ * RZZ * OXZZ +
     &           RZY * RYZ * RXZ * OXYZ +
     &           RZY * RYZ * RYZ * OYYZ +
     &           RZY * RYZ * RZZ * OYZZ +
     &           RZY * RZZ * RXZ * OXZZ +
     &           RZY * RZZ * RYZ * OYZZ +
     &           RZY * RZZ * RZZ * OZZZ
C
         rOXYZ = RXX * RXY * RXZ * OXXX +
     &           RXX * RXY * RYZ * OXXY +
     &           RXX * RXY * RZZ * OXXZ +
     &           RXX * RYY * RXZ * OXXY +
     &           RXX * RYY * RYZ * OXYY +
     &           RXX * RYY * RZZ * OXYZ +
     &           RXX * RZY * RXZ * OXXZ +
     &           RXX * RZY * RYZ * OXYZ +
     &           RXX * RZY * RZZ * OXZZ +
     &           RYX * RXY * RXZ * OXXY +
     &           RYX * RXY * RYZ * OXYY +
     &           RYX * RXY * RZZ * OXYZ +
     &           RYX * RYY * RXZ * OXYY +
     &           RYX * RYY * RYZ * OYYY +
     &           RYX * RYY * RZZ * OYYZ +
     &           RYX * RZY * RXZ * OXYZ +
     &           RYX * RZY * RYZ * OYYZ +
     &           RYX * RZY * RZZ * OYZZ +
     &           RZX * RXY * RXZ * OXXZ +
     &           RZX * RXY * RYZ * OXYZ +
     &           RZX * RXY * RZZ * OXZZ +
     &           RZX * RYY * RXZ * OXYZ +
     &           RZX * RYY * RYZ * OYYZ +
     &           RZX * RYY * RZZ * OYZZ +
     &           RZX * RZY * RXZ * OXZZ +
     &           RZX * RZY * RYZ * OYZZ +
     &           RZX * RZY * RZZ * OZZZ
C
C        SAVE
C
         DIP(N,1) = rDX
         DIP(N,2) = rDY
         DIP(N,3) = rDZ
C
         QAD(N,1) = rQXX
         QAD(N,2) = rQYY
         QAD(N,3) = rQZZ
         QAD(N,4) = rQXY
         QAD(N,5) = rQXZ
         QAD(N,6) = rQYZ
C
         OCT(N,1) = rOXXX
         OCT(N,2) = rOYYY
         OCT(N,3) = rOZZZ
         OCT(N,4) = rOXXY
         OCT(N,5) = rOXXZ
         OCT(N,6) = rOXYY
         OCT(N,7) = rOYYZ
         OCT(N,8) = rOXZZ
         OCT(N,9) = rOYZZ
         OCT(N,10)= rOXYZ
C
 11   CONTINUE
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE ROTDM1(DIP,QAD,OCT,ROT,NMODES,NDMA)
C
C -----------------------------------------------------------------------------
C
C               ROTATE THE DISTRIBUTED MULTIPOLE MOMENT TENSORS
C            WHICH ARE THE FIRST DERIVATIVES WRT NORMAL COORDINATE
C 
C              Bartosz Blasiak                        28.11.2013
C
C -----------------------------------------------------------------------------
C
C   Description:
C     Performs the unitary transformations:
C
C       D_a   = D_a'     * R_a'a
C       Q_ab  = Q_a'b'   * R_a'a * R_b'b
C       O_abc = O_a'b'c' * R_a'a * R_b'b * R_c'c
C
C   Input variables:
C     NDMA    - array of numbers of distributed electrostatic sites
C     DIP,
C     QAD,OCT - non-scalar distributed multipoles
C     ROT     - unitary rotation matrix
C     
C   Returns:
C     DIP,QAD,OCT - rotated moments
C
C   Notes:
C     The reduced format of tensor storage is used:
C
C     DIP(i) X   Y   Z
C            1   2   3
C     QAD(i) XX  YY  ZZ  XY  XZ  YZ
C            1   2   3   4   5   6
C     OCT(i) XXX YYY ZZZ XXY XXZ XYY YYZ XZZ YZZ XYZ
C            1   2   3   4   5   6   7   8   9   10
C -----------------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DIP(NMODES*NDMA*3),QAD(NMODES*NDMA*6),
     &          OCT(NMODES*NDMA*10),ROT(3,3)
      PARAMETER (ZERO=0.0D+00)
Cf2py INTENT(IN,OUT) DIP, QAD, OCT
C
      RXX = ROT(1,1)
      RYY = ROT(2,2)
      RZZ = ROT(3,3)
      RXY = ROT(1,2)
      RYX = ROT(2,1)
      RXZ = ROT(1,3)
      RZX = ROT(3,1)
      RYZ = ROT(2,3)
      RZY = ROT(3,2)
C
      NDMA3 = NDMA*3
      NDMA6 = NDMA*6
      NDMA10= NDMA*10
C
      DO N=1,NDMA
      DO M=1,NMODES
         N3 = NDMA3 *(M-1) + 3*(N-1) + 1
         N6 = NDMA6 *(M-1) + 6*(N-1) + 1
         N10= NDMA10*(M-1) +10*(N-1) + 1
C
C       UNPACK
C
         DX   = DIP(N3  )
         DY   = DIP(N3+1)
         DZ   = DIP(N3+2)
C
         QXX  = QAD(N6  )
         QYY  = QAD(N6+1)
         QZZ  = QAD(N6+2)
         QXY  = QAD(N6+3)
         QXZ  = QAD(N6+4)
         QYZ  = QAD(N6+5)
C
         OXXX = OCT(N10  )
         OYYY = OCT(N10+1)
         OZZZ = OCT(N10+2)
         OXXY = OCT(N10+3)
         OXXZ = OCT(N10+4)
         OXYY = OCT(N10+5)
         OYYZ = OCT(N10+6)
         OXZZ = OCT(N10+7)
         OYZZ = OCT(N10+8)
         OXYZ = OCT(N10+9)
C
C        ROTATE DIPOLES
C
         rDX = RXX * DX + RYX * DY + RZX * DZ
         rDY = RXY * DX + RYY * DY + RZY * DZ
         rDZ = RXZ * DX + RYZ * DY + RZZ * DZ
C
C        ROTATE QUADRUPOLES
C
         rQXX = RXX * RXX * QXX +
     &          RXX * RYX * QXY +
     &          RXX * RZX * QXZ +
     &          RYX * RXX * QXY +
     &          RYX * RYX * QYY +
     &          RYX * RZX * QYZ +
     &          RZX * RXX * QXZ +
     &          RZX * RYX * QYZ +
     &          RZX * RZX * QZZ
C
         rQYY = RXY * RXY * QXX +
     &          RXY * RYY * QXY +
     &          RXY * RZY * QXZ +
     &          RYY * RXY * QXY +
     &          RYY * RYY * QYY +
     &          RYY * RZY * QYZ +
     &          RZY * RXY * QXZ +
     &          RZY * RYY * QYZ +
     &          RZY * RZY * QZZ
C
         rQZZ = RXZ * RXZ * QXX +
     &          RXZ * RYZ * QXY +
     &          RXZ * RZZ * QXZ +
     &          RYZ * RXZ * QXY +
     &          RYZ * RYZ * QYY +
     &          RYZ * RZZ * QYZ +
     &          RZZ * RXZ * QXZ +
     &          RZZ * RYZ * QYZ +
     &          RZZ * RZZ * QZZ
C
         rQXY = RXX * RXY * QXX +
     &          RXX * RYY * QXY +
     &          RXX * RZY * QXZ +
     &          RYX * RXY * QXY +
     &          RYX * RYY * QYY +
     &          RYX * RZY * QYZ +
     &          RZX * RXY * QXZ +
     &          RZX * RYY * QYZ +
     &          RZX * RZY * QZZ
C
         rQXZ = RXX * RXZ * QXX +
     &          RXX * RYZ * QXY +
     &          RXX * RZZ * QXZ +
     &          RYX * RXZ * QXY +
     &          RYX * RYZ * QYY +
     &          RYX * RZZ * QYZ +
     &          RZX * RXZ * QXZ +
     &          RZX * RYZ * QYZ +
     &          RZX * RZZ * QZZ
C
         rQYZ = RXY * RXZ * QXX + 
     &          RXY * RYZ * QXY + 
     &          RXY * RZZ * QXZ +
     &          RYY * RXZ * QXY +
     &          RYY * RYZ * QYY +
     &          RYY * RZZ * QYZ +
     &          RZY * RXZ * QXZ +
     &          RZY * RYZ * QYZ +
     &          RZY * RZZ * QZZ
C
C        ROTATE OCTUPOLES
C
         rOXXX = RXX * RXX * RXX * OXXX +
     &           RXX * RXX * RYX * OXXY +
     &           RXX * RXX * RZX * OXXZ +
     &           RXX * RYX * RXX * OXXY +
     &           RXX * RYX * RYX * OXYY +
     &           RXX * RYX * RZX * OXYZ +
     &           RXX * RZX * RXX * OXXZ +
     &           RXX * RZX * RYX * OXYZ +
     &           RXX * RZX * RZX * OXZZ +
     &           RYX * RXX * RXX * OXXY +
     &           RYX * RXX * RYX * OXYY +
     &           RYX * RXX * RZX * OXYZ +
     &           RYX * RYX * RXX * OXYY +
     &           RYX * RYX * RYX * OYYY +
     &           RYX * RYX * RZX * OYYZ +
     &           RYX * RZX * RXX * OXYZ +
     &           RYX * RZX * RYX * OYYZ +
     &           RYX * RZX * RZX * OYZZ +
     &           RZX * RXX * RXX * OXXZ +
     &           RZX * RXX * RYX * OXYZ +
     &           RZX * RXX * RZX * OXZZ +
     &           RZX * RYX * RXX * OXYZ +
     &           RZX * RYX * RYX * OYYZ +
     &           RZX * RYX * RZX * OYZZ +
     &           RZX * RZX * RXX * OXZZ +
     &           RZX * RZX * RYX * OYZZ +
     &           RZX * RZX * RZX * OZZZ
C
         rOYYY = RXY * RXY * RXY * OXXX +
     &           RXY * RXY * RYY * OXXY +
     &           RXY * RXY * RZY * OXXZ +
     &           RXY * RYY * RXY * OXXY +
     &           RXY * RYY * RYY * OXYY +
     &           RXY * RYY * RZY * OXYZ +
     &           RXY * RZY * RXY * OXXZ +
     &           RXY * RZY * RYY * OXYZ +
     &           RXY * RZY * RZY * OXZZ +
     &           RYY * RXY * RXY * OXXY +
     &           RYY * RXY * RYY * OXYY +
     &           RYY * RXY * RZY * OXYZ +
     &           RYY * RYY * RXY * OXYY +
     &           RYY * RYY * RYY * OYYY +
     &           RYY * RYY * RZY * OYYZ +
     &           RYY * RZY * RXY * OXYZ +
     &           RYY * RZY * RYY * OYYZ +
     &           RYY * RZY * RZY * OYZZ +
     &           RZY * RXY * RXY * OXXZ +
     &           RZY * RXY * RYY * OXYZ +
     &           RZY * RXY * RZY * OXZZ +
     &           RZY * RYY * RXY * OXYZ +
     &           RZY * RYY * RYY * OYYZ +
     &           RZY * RYY * RZY * OYZZ +
     &           RZY * RZY * RXY * OXZZ +
     &           RZY * RZY * RYY * OYZZ +
     &           RZY * RZY * RZY * OZZZ
C
         rOZZZ = RXZ * RXZ * RXZ * OXXX +
     &           RXZ * RXZ * RYZ * OXXY +
     &           RXZ * RXZ * RZZ * OXXZ +
     &           RXZ * RYZ * RXZ * OXXY +
     &           RXZ * RYZ * RYZ * OXYY +
     &           RXZ * RYZ * RZZ * OXYZ +
     &           RXZ * RZZ * RXZ * OXXZ +
     &           RXZ * RZZ * RYZ * OXYZ +
     &           RXZ * RZZ * RZZ * OXZZ +
     &           RYZ * RXZ * RXZ * OXXY +
     &           RYZ * RXZ * RYZ * OXYY +
     &           RYZ * RXZ * RZZ * OXYZ +
     &           RYZ * RYZ * RXZ * OXYY +
     &           RYZ * RYZ * RYZ * OYYY +
     &           RYZ * RYZ * RZZ * OYYZ +
     &           RYZ * RZZ * RXZ * OXYZ +
     &           RYZ * RZZ * RYZ * OYYZ +
     &           RYZ * RZZ * RZZ * OYZZ +
     &           RZZ * RXZ * RXZ * OXXZ +
     &           RZZ * RXZ * RYZ * OXYZ +
     &           RZZ * RXZ * RZZ * OXZZ +
     &           RZZ * RYZ * RXZ * OXYZ +
     &           RZZ * RYZ * RYZ * OYYZ +
     &           RZZ * RYZ * RZZ * OYZZ +
     &           RZZ * RZZ * RXZ * OXZZ +
     &           RZZ * RZZ * RYZ * OYZZ +
     &           RZZ * RZZ * RZZ * OZZZ
C
         rOXXY = RXX * RXX * RXY * OXXX +
     &           RXX * RXX * RYY * OXXY +
     &           RXX * RXX * RZY * OXXZ +
     &           RXX * RYX * RXY * OXXY +
     &           RXX * RYX * RYY * OXYY +
     &           RXX * RYX * RZY * OXYZ +
     &           RXX * RZX * RXY * OXXZ +
     &           RXX * RZX * RYY * OXYZ +
     &           RXX * RZX * RZY * OXZZ +
     &           RYX * RXX * RXY * OXXY +
     &           RYX * RXX * RYY * OXYY +
     &           RYX * RXX * RZY * OXYZ +
     &           RYX * RYX * RXY * OXYY +
     &           RYX * RYX * RYY * OYYY +
     &           RYX * RYX * RZY * OYYZ +
     &           RYX * RZX * RXY * OXYZ +
     &           RYX * RZX * RYY * OYYZ +
     &           RYX * RZX * RZY * OYZZ +
     &           RZX * RXX * RXY * OXXZ +
     &           RZX * RXX * RYY * OXYZ +
     &           RZX * RXX * RZY * OXZZ +
     &           RZX * RYX * RXY * OXYZ +
     &           RZX * RYX * RYY * OYYZ +
     &           RZX * RYX * RZY * OYZZ +
     &           RZX * RZX * RXY * OXZZ +
     &           RZX * RZX * RYY * OYZZ +
     &           RZX * RZX * RZY * OZZZ
C
         rOXXZ = RXX * RXX * RXZ * OXXX +
     &           RXX * RXX * RYZ * OXXY +
     &           RXX * RXX * RZZ * OXXZ +
     &           RXX * RYX * RXZ * OXXY +
     &           RXX * RYX * RYZ * OXYY +
     &           RXX * RYX * RZZ * OXYZ +
     &           RXX * RZX * RXZ * OXXZ +
     &           RXX * RZX * RYZ * OXYZ +
     &           RXX * RZX * RZZ * OXZZ +
     &           RYX * RXX * RXZ * OXXY +
     &           RYX * RXX * RYZ * OXYY +
     &           RYX * RXX * RZZ * OXYZ +
     &           RYX * RYX * RXZ * OXYY +
     &           RYX * RYX * RYZ * OYYY +
     &           RYX * RYX * RZZ * OYYZ +
     &           RYX * RZX * RXZ * OXYZ +
     &           RYX * RZX * RYZ * OYYZ +
     &           RYX * RZX * RZZ * OYZZ +
     &           RZX * RXX * RXZ * OXXZ +
     &           RZX * RXX * RYZ * OXYZ +
     &           RZX * RXX * RZZ * OXZZ +
     &           RZX * RYX * RXZ * OXYZ +
     &           RZX * RYX * RYZ * OYYZ +
     &           RZX * RYX * RZZ * OYZZ +
     &           RZX * RZX * RXZ * OXZZ +
     &           RZX * RZX * RYZ * OYZZ +
     &           RZX * RZX * RZZ * OZZZ
C
         rOXYY = RXX * RXY * RXY * OXXX +
     &           RXX * RXY * RYY * OXXY +
     &           RXX * RXY * RZY * OXXZ +
     &           RXX * RYY * RXY * OXXY +
     &           RXX * RYY * RYY * OXYY +
     &           RXX * RYY * RZY * OXYZ +
     &           RXX * RZY * RXY * OXXZ +
     &           RXX * RZY * RYY * OXYZ +
     &           RXX * RZY * RZY * OXZZ +
     &           RYX * RXY * RXY * OXXY +
     &           RYX * RXY * RYY * OXYY +
     &           RYX * RXY * RZY * OXYZ +
     &           RYX * RYY * RXY * OXYY +
     &           RYX * RYY * RYY * OYYY +
     &           RYX * RYY * RZY * OYYZ +
     &           RYX * RZY * RXY * OXYZ +
     &           RYX * RZY * RYY * OYYZ +
     &           RYX * RZY * RZY * OYZZ +
     &           RZX * RXY * RXY * OXXZ +
     &           RZX * RXY * RYY * OXYZ +
     &           RZX * RXY * RZY * OXZZ +
     &           RZX * RYY * RXY * OXYZ +
     &           RZX * RYY * RYY * OYYZ +
     &           RZX * RYY * RZY * OYZZ +
     &           RZX * RZY * RXY * OXZZ +
     &           RZX * RZY * RYY * OYZZ +
     &           RZX * RZY * RZY * OZZZ
C
         rOYYZ = RXY * RXY * RXZ * OXXX +
     &           RXY * RXY * RYZ * OXXY +
     &           RXY * RXY * RZZ * OXXZ +
     &           RXY * RYY * RXZ * OXXY +
     &           RXY * RYY * RYZ * OXYY +
     &           RXY * RYY * RZZ * OXYZ +
     &           RXY * RZY * RXZ * OXXZ +
     &           RXY * RZY * RYZ * OXYZ +
     &           RXY * RZY * RZZ * OXZZ +
     &           RYY * RXY * RXZ * OXXY +
     &           RYY * RXY * RYZ * OXYY +
     &           RYY * RXY * RZZ * OXYZ +
     &           RYY * RYY * RXZ * OXYY +
     &           RYY * RYY * RYZ * OYYY +
     &           RYY * RYY * RZZ * OYYZ +
     &           RYY * RZY * RXZ * OXYZ +
     &           RYY * RZY * RYZ * OYYZ +
     &           RYY * RZY * RZZ * OYZZ +
     &           RZY * RXY * RXZ * OXXZ +
     &           RZY * RXY * RYZ * OXYZ +
     &           RZY * RXY * RZZ * OXZZ +
     &           RZY * RYY * RXZ * OXYZ +
     &           RZY * RYY * RYZ * OYYZ +
     &           RZY * RYY * RZZ * OYZZ +
     &           RZY * RZY * RXZ * OXZZ +
     &           RZY * RZY * RYZ * OYZZ +
     &           RZY * RZY * RZZ * OZZZ
C
         rOXZZ = RXX * RXZ * RXZ * OXXX +
     &           RXX * RXZ * RYZ * OXXY +
     &           RXX * RXZ * RZZ * OXXZ +
     &           RXX * RYZ * RXZ * OXXY +
     &           RXX * RYZ * RYZ * OXYY +
     &           RXX * RYZ * RZZ * OXYZ +
     &           RXX * RZZ * RXZ * OXXZ +
     &           RXX * RZZ * RYZ * OXYZ +
     &           RXX * RZZ * RZZ * OXZZ +
     &           RYX * RXZ * RXZ * OXXY +
     &           RYX * RXZ * RYZ * OXYY +
     &           RYX * RXZ * RZZ * OXYZ +
     &           RYX * RYZ * RXZ * OXYY +
     &           RYX * RYZ * RYZ * OYYY +
     &           RYX * RYZ * RZZ * OYYZ +
     &           RYX * RZZ * RXZ * OXYZ +
     &           RYX * RZZ * RYZ * OYYZ +
     &           RYX * RZZ * RZZ * OYZZ +
     &           RZX * RXZ * RXZ * OXXZ +
     &           RZX * RXZ * RYZ * OXYZ +
     &           RZX * RXZ * RZZ * OXZZ +
     &           RZX * RYZ * RXZ * OXYZ +
     &           RZX * RYZ * RYZ * OYYZ +
     &           RZX * RYZ * RZZ * OYZZ +
     &           RZX * RZZ * RXZ * OXZZ +
     &           RZX * RZZ * RYZ * OYZZ +
     &           RZX * RZZ * RZZ * OZZZ
C
         rOYZZ = RXY * RXZ * RXZ * OXXX +
     &           RXY * RXZ * RYZ * OXXY +
     &           RXY * RXZ * RZZ * OXXZ +
     &           RXY * RYZ * RXZ * OXXY +
     &           RXY * RYZ * RYZ * OXYY +
     &           RXY * RYZ * RZZ * OXYZ +
     &           RXY * RZZ * RXZ * OXXZ +
     &           RXY * RZZ * RYZ * OXYZ +
     &           RXY * RZZ * RZZ * OXZZ +
     &           RYY * RXZ * RXZ * OXXY +
     &           RYY * RXZ * RYZ * OXYY +
     &           RYY * RXZ * RZZ * OXYZ +
     &           RYY * RYZ * RXZ * OXYY +
     &           RYY * RYZ * RYZ * OYYY +
     &           RYY * RYZ * RZZ * OYYZ +
     &           RYY * RZZ * RXZ * OXYZ +
     &           RYY * RZZ * RYZ * OYYZ +
     &           RYY * RZZ * RZZ * OYZZ +
     &           RZY * RXZ * RXZ * OXXZ +
     &           RZY * RXZ * RYZ * OXYZ +
     &           RZY * RXZ * RZZ * OXZZ +
     &           RZY * RYZ * RXZ * OXYZ +
     &           RZY * RYZ * RYZ * OYYZ +
     &           RZY * RYZ * RZZ * OYZZ +
     &           RZY * RZZ * RXZ * OXZZ +
     &           RZY * RZZ * RYZ * OYZZ +
     &           RZY * RZZ * RZZ * OZZZ
C
         rOXYZ = RXX * RXY * RXZ * OXXX +
     &           RXX * RXY * RYZ * OXXY +
     &           RXX * RXY * RZZ * OXXZ +
     &           RXX * RYY * RXZ * OXXY +
     &           RXX * RYY * RYZ * OXYY +
     &           RXX * RYY * RZZ * OXYZ +
     &           RXX * RZY * RXZ * OXXZ +
     &           RXX * RZY * RYZ * OXYZ +
     &           RXX * RZY * RZZ * OXZZ +
     &           RYX * RXY * RXZ * OXXY +
     &           RYX * RXY * RYZ * OXYY +
     &           RYX * RXY * RZZ * OXYZ +
     &           RYX * RYY * RXZ * OXYY +
     &           RYX * RYY * RYZ * OYYY +
     &           RYX * RYY * RZZ * OYYZ +
     &           RYX * RZY * RXZ * OXYZ +
     &           RYX * RZY * RYZ * OYYZ +
     &           RYX * RZY * RZZ * OYZZ +
     &           RZX * RXY * RXZ * OXXZ +
     &           RZX * RXY * RYZ * OXYZ +
     &           RZX * RXY * RZZ * OXZZ +
     &           RZX * RYY * RXZ * OXYZ +
     &           RZX * RYY * RYZ * OYYZ +
     &           RZX * RYY * RZZ * OYZZ +
     &           RZX * RZY * RXZ * OXZZ +
     &           RZX * RZY * RYZ * OYZZ +
     &           RZX * RZY * RZZ * OZZZ
C
C        SAVE
C
         DIP(N3  ) = rDX
         DIP(N3+1) = rDY
         DIP(N3+2) = rDZ
C
         QAD(N6  ) = rQXX
         QAD(N6+1) = rQYY
         QAD(N6+2) = rQZZ
         QAD(N6+3) = rQXY
         QAD(N6+4) = rQXZ
         QAD(N6+5) = rQYZ
C
         OCT(N10  ) = rOXXX
         OCT(N10+1) = rOYYY
         OCT(N10+2) = rOZZZ
         OCT(N10+3) = rOXXY
         OCT(N10+4) = rOXXZ
         OCT(N10+5) = rOXYY
         OCT(N10+6) = rOYYZ
         OCT(N10+7) = rOXZZ
         OCT(N10+8) = rOYZZ
         OCT(N10+9) = rOXYZ
C
      ENDDO
      ENDDO
C
      RETURN
      END

C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE TRACLS(QAD,OCT,NDMA)
C
C -----------------------------------------------------------------------------
C
C               FORM TRACELESS QUADRUPOLE AND OCTUPOLE MOMENTS
C                     ACCORDING TO BUCKINGHAM CONVENTION
C 
C              Bartosz Blasiak                        08.11.2013
C
C -----------------------------------------------------------------------------
C
C   Description:
C     Performs the transformations:
C
C       Q_ab  = Q_a'b'   * R_a'a * R_b'b
C       O_abc = O_a'b'c' * R_a'a * R_b'b * R_c'c
C
C   Input variables:
C     NDMA    - array of numbers of distributed electrostatic sites
C     QAD,OCT - distributed quadrupoles and octupoles
C     
C   Returns:
C     QAD,OCT - tensors in traceless forms
C
C   Notes:
C     The reduced format of tensor storage is used:
C
C     QAD(i) XX  YY  ZZ  XY  XZ  YZ
C            1   2   3   4   5   6
C     OCT(i) XXX YYY ZZZ XXY XXZ XYY YYZ XZZ YZZ XYZ
C            1   2   3   4   5   6   7   8   9   10
C -----------------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION QAD(NDMA,6),OCT(NDMA,10)
      PARAMETER (POLTRA=1.500000000D+00,HALF=0.50000000000D+00,
     &           ONEHLF=2.50000D+00,THREE=3.00000D+00)
Cf2py INTENT(IN,OUT) QAD,OCT
C
      DO 11 N=1,NDMA
C
C       UNPACK
C
         QXX  = QAD(N,1)
         QYY  = QAD(N,2)
         QZZ  = QAD(N,3)
         QXY  = QAD(N,4)
         QXZ  = QAD(N,5)
         QYZ  = QAD(N,6)
C
         OXXX = OCT(N,1)
         OYYY = OCT(N,2)
         OZZZ = OCT(N,3)
         OXXY = OCT(N,4)
         OXXZ = OCT(N,5)
         OXYY = OCT(N,6)
         OYYZ = OCT(N,7)
         OXZZ = OCT(N,8)
         OYZZ = OCT(N,9)
         OXYZ = OCT(N,10)
C
C        TRACELESS QUADRUPOLES
C
         TRACE =(QXX+QYY+QZZ)*HALF
C
         tQXX  = QXX * POLTRA - TRACE
         tQYY  = QYY * POLTRA - TRACE
         tQZZ  = QZZ * POLTRA - TRACE
         tQXY  = QXY * POLTRA
         tQXZ  = QXZ * POLTRA
         tQYZ  = QYZ * POLTRA
C
C        TRACELESS OCTUPOLES
C
         TX =(OXXX+OXYY+OXZZ)*HALF
         TY =(OXXY+OYYY+OYZZ)*HALF
         TZ =(OXXZ+OYYZ+OZZZ)*HALF
C
         tOXXX = OXXX * ONEHLF - TX * THREE
         tOYYY = OYYY * ONEHLF - TY * THREE
         tOZZZ = OZZZ * ONEHLF - TZ * THREE
         tOXXY = OXXY * ONEHLF - TY
         tOXXZ = OXXZ * ONEHLF - TZ
         tOXYY = OXYY * ONEHLF - TX
         tOYYZ = OYYZ * ONEHLF - TZ
         tOXZZ = OXZZ * ONEHLF - TX
         tOYZZ = OYZZ * ONEHLF - TY
         tOXYZ = OXYZ * ONEHLF
C
C        SAVE
C
         QAD(N,1) = tQXX
         QAD(N,2) = tQYY
         QAD(N,3) = tQZZ
         QAD(N,4) = tQXY
         QAD(N,5) = tQXZ
         QAD(N,6) = tQYZ
C
         OCT(N,1) = tOXXX
         OCT(N,2) = tOYYY
         OCT(N,3) = tOZZZ
         OCT(N,4) = tOXXY
         OCT(N,5) = tOXXZ
         OCT(N,6) = tOXYY
         OCT(N,7) = tOYYZ
         OCT(N,8) = tOXZZ
         OCT(N,9) = tOYZZ
         OCT(N,10)= tOXYZ
C
 11   CONTINUE
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE TRACL1(QAD,OCT,NMODES,NDMA)
C
C -----------------------------------------------------------------------------
C
C               FORM TRACELESS QUADRUPOLE AND OCTUPOLE MOMENTS
C                     ACCORDING TO BUCKINGHAM CONVENTION
C                       TENSORS ARE FIRST DERIVATIVES
C                           WRT NORMAL COORDINATES
C 
C              Bartosz Blasiak                        28.11.2013
C
C -----------------------------------------------------------------------------
C
C   Description:
C     Performs the transformations:
C
C       Q_ab  = Q_a'b'   * R_a'a * R_b'b
C       O_abc = O_a'b'c' * R_a'a * R_b'b * R_c'c
C
C   Input variables:
C     NDMA    - array of numbers of distributed electrostatic sites
C     QAD,OCT - distributed quadrupoles and octupoles
C     
C   Returns:
C     QAD,OCT - tensors in traceless forms
C
C   Notes:
C     The reduced format of tensor storage is used:
C
C     QAD(i) XX  YY  ZZ  XY  XZ  YZ
C            1   2   3   4   5   6
C     OCT(i) XXX YYY ZZZ XXY XXZ XYY YYZ XZZ YZZ XYZ
C            1   2   3   4   5   6   7   8   9   10
C -----------------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION QAD(NMODES*NDMA*6),OCT(NMODES*NDMA*10)
      PARAMETER (POLTRA=1.500000000D+00,HALF=0.50000000000D+00,
     &           ONEHLF=2.50000D+00,THREE=3.00000D+00)
Cf2py INTENT(IN,OUT) QAD,OCT
C
      NDMA6 = NDMA*6
      NDMA10= NDMA*10
C
      DO N=1,NDMA
      DO M=1,NMODES
         N6 = NDMA6 *(M-1) + 6*(N-1) + 1
         N10= NDMA10*(M-1) +10*(N-1) + 1
C
C       UNPACK
C
         QXX  = QAD(N6  )
         QYY  = QAD(N6+1)
         QZZ  = QAD(N6+2)
         QXY  = QAD(N6+3)
         QXZ  = QAD(N6+4)
         QYZ  = QAD(N6+5)
C
         OXXX = OCT(N10  )
         OYYY = OCT(N10+1)
         OZZZ = OCT(N10+2)
         OXXY = OCT(N10+3)
         OXXZ = OCT(N10+4)
         OXYY = OCT(N10+5)
         OYYZ = OCT(N10+6)
         OXZZ = OCT(N10+7)
         OYZZ = OCT(N10+8)
         OXYZ = OCT(N10+9)
C
C        TRACELESS QUADRUPOLES
C
         TRACE =(QXX+QYY+QZZ)*HALF
C
         tQXX  = QXX * POLTRA - TRACE
         tQYY  = QYY * POLTRA - TRACE
         tQZZ  = QZZ * POLTRA - TRACE
         tQXY  = QXY * POLTRA
         tQXZ  = QXZ * POLTRA
         tQYZ  = QYZ * POLTRA
C
C        TRACELESS OCTUPOLES
C
         TX =(OXXX+OXYY+OXZZ)*HALF
         TY =(OXXY+OYYY+OYZZ)*HALF
         TZ =(OXXZ+OYYZ+OZZZ)*HALF
C
         tOXXX = OXXX * ONEHLF - TX * THREE
         tOYYY = OYYY * ONEHLF - TY * THREE
         tOZZZ = OZZZ * ONEHLF - TZ * THREE
         tOXXY = OXXY * ONEHLF - TY
         tOXXZ = OXXZ * ONEHLF - TZ
         tOXYY = OXYY * ONEHLF - TX
         tOYYZ = OYYZ * ONEHLF - TZ
         tOXZZ = OXZZ * ONEHLF - TX
         tOYZZ = OYZZ * ONEHLF - TY
         tOXYZ = OXYZ * ONEHLF
C
C        SAVE
C
         QAD(N6  ) = tQXX
         QAD(N6+1) = tQYY
         QAD(N6+2) = tQZZ
         QAD(N6+3) = tQXY
         QAD(N6+4) = tQXZ
         QAD(N6+5) = tQYZ
C
         OCT(N10  )= tOXXX
         OCT(N10+1)= tOYYY
         OCT(N10+2)= tOZZZ
         OCT(N10+3)= tOXXY
         OCT(N10+4)= tOXXZ
         OCT(N10+5)= tOXYY
         OCT(N10+6)= tOYYZ
         OCT(N10+7)= tOXZZ
         OCT(N10+8)= tOYZZ
         OCT(N10+9)= tOXYZ
C
      ENDDO
      ENDDO
C
      RETURN
      END

C-----|--|---------|---------|---------|---------|---------|---------|--|------|
