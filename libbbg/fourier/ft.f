C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE FFT(FR,FI,NF,M)
C
C     FAST FOURIER TRANSFORM
C     IMPLEMENTS THE SIMPLEST VERSION OF COOLEY-TUKEY ALGORITHM
C     THE EFFICIENCY SCALES AS NF*LOG_2(NF)
C
C                                         Oct 9 2014
C                                         Nov11 2014 - added title
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION FR(0:NF-1), FI(0:NF-1)
      PARAMETER (ZERO=0.00D+00, ONE=1.00D+00, TWO=2.00D+00)
      PARAMETER (PI=3.141592653589793D+00)
Cf2py INTENT(IN,OUT) FR,FI
C
      NH = NF/2
      NP = 2**M
      K  = 1
C
      DO I=0,NF-2
         IF(I.LT.K-1) THEN
            F1 = FR(K-1)
            F2 = FI(K-1)
            FR(K-1) = FR(I)
            FI(K-1) = FI(I)
            FR(I) = F1
            FI(I) = F2
         ENDIF
         J = NH
         DO WHILE (J.LT.K) 
            K = K - J
            J = J / 2
         ENDDO
         K = K + J
      ENDDO
C
      K = 1
      DO I=0,M-1
         W = ZERO
         J = K
         K = J*2
         DO IP=0,J-1
            U = DCOS(W)
            V =-DSIN(W)
            W = W + PI/DFLOAT(J)
            DO IQ=IP,NF-1,K
               II = IQ + J
               F1 = FR(II)*U - FI(II)*V
               F2 = FR(II)*V + FI(II)*U
               FR(II) = FR(IQ) - F1
               FR(IQ) = FR(IQ) + F1
               FI(II) = FI(IQ) - F2
               FI(IQ) = FI(IQ) + F2
            ENDDO      
         ENDDO
      ENDDO
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE DFT(FR,FI,NF,GR,GI)
C
C     IMPLEMENTS DIRECT FOURIER TRANSFORM ALGORITHM
C     THE EFFICIENCY SCALES AS NF^2
C
C                                Oct 8 2014
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      DIMENSION FR(NF), FI(NF), GR(NF), GI(NF)
      PARAMETER (ZERO=0.00D+00, ONE=1.00D+00, TWO=2.00D+00)
      PARAMETER (PI=3.141592653589793D+00)
Cf2py INTENT(IN,OUT) GR,GI
C
      RNF = DFLOAT(NF)
      X   = TWO*PI/RNF
      DO I=1,NF
         IM = I-1
         GRV = ZERO
         GIV = ZERO
         DO J=1,NF
             JM = J-1
             Q  = X*IM*JM 
             FRJ= FR(J)
             FIJ= FI(J)
             DCOSQ = DCOS(Q)
             DSINQ = DSIN(Q)
             GRV = GRV + FRJ*DCOSQ + FIJ*DSINQ
             GIV = GIV + FIJ*DCOSQ - FRJ*DSINQ
         ENDDO
         GR(I) = GRV
         GI(I) = GIV
      ENDDO
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE FFTINV(GR,GI,NG,M)
C
C     FAST INVERSE FOURIER TRANSFORM
C     IMPLEMENTS THE SIMPLEST VERSION OF COOLEY-TUKEY ALGORITHM
C     THE EFFICIENCY SCALES AS NF*LOG_2(NF)
C
C                                         Nov 18 2014
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GR(0:NG-1), GI(0:NG-1)
      PARAMETER (ZERO=0.00D+00, ONE=1.00D+00, TWO=2.00D+00)
      PARAMETER (PI=3.141592653589793D+00)
Cf2py INTENT(IN,OUT) GR,GI
C
      NH = NG/2
      NP = 2**M
      K  = 1
C
      DO I=0,NG-2
         IF(I.LT.K-1) THEN
            G1 = GR(K-1)
            G2 = GI(K-1)
            GR(K-1) = GR(I)
            GI(K-1) = GI(I)
            GR(I) = G1
            GI(I) = G2
         ENDIF
         J = NH
         DO WHILE (J.LT.K) 
            K = K - J
            J = J / 2
         ENDDO
         K = K + J
      ENDDO
C
      K = 1
      DO I=0,M-1
         W = ZERO
         J = K
         K = J*2
         DO IP=0,J-1
            U = DCOS(W)
            V = DSIN(W)
            W = W + PI/DFLOAT(J)
            DO IQ=IP,NG-1,K
               II = IQ + J
               F1 = GR(II)*U - GI(II)*V
               F2 = GR(II)*V + GI(II)*U
               GR(II) = GR(IQ) - F1
               GR(IQ) = GR(IQ) + F1
               GI(II) = GI(IQ) - F2
               GI(IQ) = GI(IQ) + F2
            ENDDO      
         ENDDO
      ENDDO
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE DFTINV(GR,GI,NG,FR,FI)
C
C     IMPLEMENTS DIRECT INVERSE FOURIER TRANSFORM ALGORITHM
C     THE EFFICIENCY SCALES AS NF^2
C
C                                Nov 18 2014
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      DIMENSION FR(NG), FI(NG), GR(NG), GI(NG)
      PARAMETER (ZERO=0.00D+00, ONE=1.00D+00, TWO=2.00D+00)
      PARAMETER (PI=3.141592653589793D+00)
Cf2py INTENT(IN,OUT) FR,FI
C
      RNG = DFLOAT(NF)
      X   = TWO*PI/RNG
      DO I=1,NG
         IM = I-1
         FRV = ZERO
         FIV = ZERO
         DO J=1,NG
             JM = J-1
             Q  = X*IM*JM 
             GRJ= GR(J)
             GIJ= GI(J)
             DCOSQ = DCOS(Q)
             DSINQ = DSIN(Q)
             FRV = FRV + GRJ*DCOSQ - GIJ*DSINQ
             FIV = FIV + GIJ*DCOSQ + GRJ*DSINQ
         ENDDO
         FR(I) = FRV
         FI(I) = FIV
      ENDDO
C
      RETURN
      END
