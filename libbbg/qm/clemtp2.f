C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE TDISP6(RPOL,RPOL1,DPOL,NPOL,DPOL1,GIJJ,REDMSS,FREQ,
     &                  NMOLS,NPOLC,NMODES,MPOL,MRPOL,MODE,DISP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DPOL(MPOL),DPOL1(NMODES*NPOLC*12*9),NPOL(NMOLS),
     &          RPOL(MRPOL),RPOL1(NMODES*NPOLC*3),
     &          GIJJ(NMODES),REDMSS(NMODES),FREQ(NMODES),
     &          GIVEC(30),WEIGHT(12),ABSICA(12)
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,THREE=3.0D+00,
     &           FOUR=4.0D+00,FIVE=5.0D+00,SIX=6.0D+00,
     &           NINE=9.0D+00,V24=24.0D+00,V90=90.0D+00,
     &           TOCMRC=219474.63067873946D+00,PI=3.141592653589793D+00)
C
Cf2py INTENT(OUT) DISP
      DATA GIVEC/30*0.0D+00/,XN0/0.3000D+00/
C
C     12-POINT GAUSS-LEGENDRE QUADRATURE FOR [-1: 1]
C
      DATA WEIGHT/0.0471753363865118D+00,0.1069393259953184D+00,
     &            0.1600783285433462D+00,0.2031674267230659D+00,
     &            0.2334925365383548D+00,0.2491470458134028D+00,
     &            0.2491470458134028D+00,0.2334925365383548D+00,
     &            0.2031674267230659D+00,0.1600783285433462D+00,
     &            0.1069393259953184D+00,0.0471753363865118D+00/
C
      DATA ABSICA/-0.9815606342467192D+00,-0.9041172563704749D+00,
     &            -0.7699026741943047D+00,-0.5873179542866175D+00,
     &            -0.3678314989981802D+00,-0.1252334085114689D+00,
     &             0.1252334085114689D+00, 0.3678314989981802D+00,
     &             0.5873179542866175D+00, 0.7699026741943047D+00,
     &             0.9041172563704749D+00, 0.9815606342467192D+00/
C
      DISP = ZERO
C
C     COMPUTE VIBRATIONAL WEIGHTS
C
      DO 999 M=1, NMODES
         FREQM = FREQ(M)
         GIVEC(M) = GIJJ(M)/(REDMSS(M)*FREQM*FREQM)
 999  CONTINUE
C                                                                          
C     ITERATE OVER IR MOLECULE AND ITS POLARIZABLE CENTERS
C                                                                       
      NPOLI = 0
      DO IMOL=1,1
      NIM  = NPOL(IMOL)
      NPOLI = NPOLI + NIM
      DO I=1,NIM
         NIX3 = 3*(NPOLI-NIM) + 3*(I-1) + 1
C
C        EXTRACT CENTERS OF I MOLECULE
C
         RPOLIX = RPOL(NIX3  )
         RPOLIY = RPOL(NIX3+1)
         RPOLIZ = RPOL(NIX3+2)
C
C        ITERATE OVER IMAGINARY FREQUENCIES
C
         DO N=1,12
C
            NNIX9 = 12*9*(NPOLI-NIM) + 12*9*(I-1) + 9*(N-1) + 1
C
C           COMPUTE GAUSS-LEGENDRE WEIGHT
C
            TX = ONE - ABSICA(N)
            G12W = WEIGHT(N) * TWO * XN0 / (TX*TX)
C
C           EXTRACT ALL POLARIZABILITIES FOR A GIVEN N (9 numbers)
C
            VAIXX = DPOL(NNIX9  )
            VAIXY = DPOL(NNIX9+1)
            VAIXZ = DPOL(NNIX9+2)
            VAIYX = DPOL(NNIX9+3)
            VAIYY = DPOL(NNIX9+4)
            VAIYZ = DPOL(NNIX9+5)
            VAIZX = DPOL(NNIX9+6)
            VAIZY = DPOL(NNIX9+7)
            VAIZZ = DPOL(NNIX9+8)
C
C           ITERATE OVER NORMAL MODES OF IR MOLECULE
C 
            DO M=1, NMODES                      
C                                                                        
               GRF = GIVEC(M)
C
               MNIX3  = NPOLC*   3*(M-1) +    3*(I-1)           + 1
               MNNIX9 = NPOLC*12*9*(M-1) + 12*9*(I-1) + 9*(N-1) + 1
C                                                                         
C              EXTRACT DERIVATIVES OF LMO CENTROIDS
C
               RX1 = RPOL1(MNIX3  )
               RY1 = RPOL1(MNIX3+1)
               RZ1 = RPOL1(MNIX3+2)
C
C              EXTRACT ALL MODE POLARIZABILITY DERIVATIVES FOR A GIVEN N (9 numbers)
C
               VAIMXX = DPOL1(MNNIX9  )
               VAIMXY = DPOL1(MNNIX9+1)
               VAIMXZ = DPOL1(MNNIX9+2)
               VAIMYX = DPOL1(MNNIX9+3)
               VAIMYY = DPOL1(MNNIX9+4)
               VAIMYZ = DPOL1(MNNIX9+5)
               VAIMZX = DPOL1(MNNIX9+6)
               VAIMZY = DPOL1(MNNIX9+7)
               VAIMZZ = DPOL1(MNNIX9+8)
C                                                                            
C              ITERATE OVER SOLVENT MOLECULE POLARIZABLE CENTERS         
C                                                                     
               NPOLJ = NPOL(1)
               DO JMOL=2, NMOLS
                  NJM = NPOL(JMOL)
                  NPOLJ = NPOLJ + NJM
                  DO J=1,NJM
C                                                                        
                     NJX3 = 3*(NPOLJ-NJM) + 3*(J-1) + 1
C                                                                        
                     NNJX9 = 12*9*(NPOLJ-NJM) + 12*9*(J-1) + 9*(N-1) + 1
C                                                                        
C                    COMPUTE DISTANCES BETWEEN I AND J CENTER
C
                     RIJX = RPOLIX - RPOL(NJX3  )                    
                     RIJY = RPOLIY - RPOL(NJX3+1)
                     RIJZ = RPOLIZ - RPOL(NJX3+2)
C                                                                    
                     RIJ  = ONE/DSQRT(RIJX*RIJX+RIJY*RIJY+RIJZ*RIJZ)
                     RIJ2 = RIJ  * RIJ
                     RIJ4 = RIJ2 * RIJ2
                     RIJ6 = RIJ4 * RIJ2
                     RIJ8 = RIJ6 * RIJ2
                     RIJ10= RIJ8 * RIJ2
                     RIJ12= RIJ10* RIJ2
C                                                                        
                     RMR = RIJX * RX1 + RIJY * RY1 + RIJZ * RZ1
C                                                                           
C                    EXTRACT ALL POLARIZABILITIES FOR CENTER J (9 numbers)
C
                     VAJXX = DPOL(NNJX9  )
                     VAJXY = DPOL(NNJX9+1)
                     VAJXZ = DPOL(NNJX9+2)
                     VAJYX = DPOL(NNJX9+3)
                     VAJYY = DPOL(NNJX9+4)
                     VAJYZ = DPOL(NNJX9+5)
                     VAJZX = DPOL(NNJX9+6)
                     VAJZY = DPOL(NNJX9+7)
                     VAJZZ = DPOL(NNJX9+8)
C
C ---------------------------------------------------------------------------------
C                    QUADRUPLE SUM OVER X, Y and Z TO COMPUTE TERMS I AND II
C
C ---------------------------------------------------------------------------------
C
C                    EVALUATE ALL I_IJKL ELEMENTS
C
                     VIXXXX = VAIXX * VAJXX
                     VIXXXY = VAIXX * VAJXY
                     VIXXXZ = VAIXX * VAJXZ
                     VIXXYX = VAIXY * VAJXX
                     VIXXYY = VAIXY * VAJXY
                     VIXXYZ = VAIXY * VAJXZ
                     VIXXZX = VAIXZ * VAJXX
                     VIXXZY = VAIXZ * VAJXY
                     VIXXZZ = VAIXZ * VAJXZ
                     VIXYXX = VAIXX * VAJYX
                     VIXYXY = VAIXX * VAJYY
                     VIXYXZ = VAIXX * VAJYZ
                     VIXYYX = VAIXY * VAJYX
                     VIXYYY = VAIXY * VAJYY
                     VIXYYZ = VAIXY * VAJYZ
                     VIXYZX = VAIXZ * VAJYX
                     VIXYZY = VAIXZ * VAJYY
                     VIXYZZ = VAIXZ * VAJYZ
                     VIXZXX = VAIXX * VAJZX
                     VIXZXY = VAIXX * VAJZY
                     VIXZXZ = VAIXX * VAJZZ
                     VIXZYX = VAIXY * VAJZX
                     VIXZYY = VAIXY * VAJZY
                     VIXZYZ = VAIXY * VAJZZ
                     VIXZZX = VAIXZ * VAJZX
                     VIXZZY = VAIXZ * VAJZY
                     VIXZZZ = VAIXZ * VAJZZ
                     VIYXXX = VAIYX * VAJXX
                     VIYXXY = VAIYX * VAJXY
                     VIYXXZ = VAIYX * VAJXZ
                     VIYXYX = VAIYY * VAJXX
                     VIYXYY = VAIYY * VAJXY
                     VIYXYZ = VAIYY * VAJXZ
                     VIYXZX = VAIYZ * VAJXX
                     VIYXZY = VAIYZ * VAJXY
                     VIYXZZ = VAIYZ * VAJXZ
                     VIYYXX = VAIYX * VAJYX
                     VIYYXY = VAIYX * VAJYY
                     VIYYXZ = VAIYX * VAJYZ
                     VIYYYX = VAIYY * VAJYX
                     VIYYYY = VAIYY * VAJYY
                     VIYYYZ = VAIYY * VAJYZ
                     VIYYZX = VAIYZ * VAJYX
                     VIYYZY = VAIYZ * VAJYY
                     VIYYZZ = VAIYZ * VAJYZ
                     VIYZXX = VAIYX * VAJZX
                     VIYZXY = VAIYX * VAJZY
                     VIYZXZ = VAIYX * VAJZZ
                     VIYZYX = VAIYY * VAJZX
                     VIYZYY = VAIYY * VAJZY
                     VIYZYZ = VAIYY * VAJZZ
                     VIYZZX = VAIYZ * VAJZX
                     VIYZZY = VAIYZ * VAJZY
                     VIYZZZ = VAIYZ * VAJZZ
                     VIZXXX = VAIZX * VAJXX
                     VIZXXY = VAIZX * VAJXY
                     VIZXXZ = VAIZX * VAJXZ
                     VIZXYX = VAIZY * VAJXX
                     VIZXYY = VAIZY * VAJXY
                     VIZXYZ = VAIZY * VAJXZ
                     VIZXZX = VAIZZ * VAJXX
                     VIZXZY = VAIZZ * VAJXY
                     VIZXZZ = VAIZZ * VAJXZ
                     VIZYXX = VAIZX * VAJYX
                     VIZYXY = VAIZX * VAJYY
                     VIZYXZ = VAIZX * VAJYZ
                     VIZYYX = VAIZY * VAJYX
                     VIZYYY = VAIZY * VAJYY
                     VIZYYZ = VAIZY * VAJYZ
                     VIZYZX = VAIZZ * VAJYX
                     VIZYZY = VAIZZ * VAJYY
                     VIZYZZ = VAIZZ * VAJYZ
                     VIZZXX = VAIZX * VAJZX
                     VIZZXY = VAIZX * VAJZY
                     VIZZXZ = VAIZX * VAJZZ
                     VIZZYX = VAIZY * VAJZX
                     VIZZYY = VAIZY * VAJZY
                     VIZZYZ = VAIZY * VAJZZ
                     VIZZZX = VAIZZ * VAJZX
                     VIZZZY = VAIZZ * VAJZY
                     VIZZZZ = VAIZZ * VAJZZ
C
C                    EVALUATE ALL Q_IJKL ELEMENTS
C
                     VQXXXX = VAIMXX * VAJXX
                     VQXXXY = VAIMXX * VAJXY
                     VQXXXZ = VAIMXX * VAJXZ
                     VQXXYX = VAIMXY * VAJXX
                     VQXXYY = VAIMXY * VAJXY
                     VQXXYZ = VAIMXY * VAJXZ
                     VQXXZX = VAIMXZ * VAJXX
                     VQXXZY = VAIMXZ * VAJXY
                     VQXXZZ = VAIMXZ * VAJXZ
                     VQXYXX = VAIMXX * VAJYX
                     VQXYXY = VAIMXX * VAJYY
                     VQXYXZ = VAIMXX * VAJYZ
                     VQXYYX = VAIMXY * VAJYX
                     VQXYYY = VAIMXY * VAJYY
                     VQXYYZ = VAIMXY * VAJYZ
                     VQXYZX = VAIMXZ * VAJYX
                     VQXYZY = VAIMXZ * VAJYY
                     VQXYZZ = VAIMXZ * VAJYZ
                     VQXZXX = VAIMXX * VAJZX
                     VQXZXY = VAIMXX * VAJZY
                     VQXZXZ = VAIMXX * VAJZZ
                     VQXZYX = VAIMXY * VAJZX
                     VQXZYY = VAIMXY * VAJZY
                     VQXZYZ = VAIMXY * VAJZZ
                     VQXZZX = VAIMXZ * VAJZX
                     VQXZZY = VAIMXZ * VAJZY
                     VQXZZZ = VAIMXZ * VAJZZ
                     VQYXXX = VAIMYX * VAJXX
                     VQYXXY = VAIMYX * VAJXY
                     VQYXXZ = VAIMYX * VAJXZ
                     VQYXYX = VAIMYY * VAJXX
                     VQYXYY = VAIMYY * VAJXY
                     VQYXYZ = VAIMYY * VAJXZ
                     VQYXZX = VAIMYZ * VAJXX
                     VQYXZY = VAIMYZ * VAJXY
                     VQYXZZ = VAIMYZ * VAJXZ
                     VQYYXX = VAIMYX * VAJYX
                     VQYYXY = VAIMYX * VAJYY
                     VQYYXZ = VAIMYX * VAJYZ
                     VQYYYX = VAIMYY * VAJYX
                     VQYYYY = VAIMYY * VAJYY
                     VQYYYZ = VAIMYY * VAJYZ
                     VQYYZX = VAIMYZ * VAJYX
                     VQYYZY = VAIMYZ * VAJYY
                     VQYYZZ = VAIMYZ * VAJYZ
                     VQYZXX = VAIMYX * VAJZX
                     VQYZXY = VAIMYX * VAJZY
                     VQYZXZ = VAIMYX * VAJZZ
                     VQYZYX = VAIMYY * VAJZX
                     VQYZYY = VAIMYY * VAJZY
                     VQYZYZ = VAIMYY * VAJZZ
                     VQYZZX = VAIMYZ * VAJZX
                     VQYZZY = VAIMYZ * VAJZY
                     VQYZZZ = VAIMYZ * VAJZZ
                     VQZXXX = VAIMZX * VAJXX
                     VQZXXY = VAIMZX * VAJXY
                     VQZXXZ = VAIMZX * VAJXZ
                     VQZXYX = VAIMZY * VAJXX
                     VQZXYY = VAIMZY * VAJXY
                     VQZXYZ = VAIMZY * VAJXZ
                     VQZXZX = VAIMZZ * VAJXX
                     VQZXZY = VAIMZZ * VAJXY
                     VQZXZZ = VAIMZZ * VAJXZ
                     VQZYXX = VAIMZX * VAJYX
                     VQZYXY = VAIMZX * VAJYY
                     VQZYXZ = VAIMZX * VAJYZ
                     VQZYYX = VAIMZY * VAJYX
                     VQZYYY = VAIMZY * VAJYY
                     VQZYYZ = VAIMZY * VAJYZ
                     VQZYZX = VAIMZZ * VAJYX
                     VQZYZY = VAIMZZ * VAJYY
                     VQZYZZ = VAIMZZ * VAJYZ
                     VQZZXX = VAIMZX * VAJZX
                     VQZZXY = VAIMZX * VAJZY
                     VQZZXZ = VAIMZX * VAJZZ
                     VQZZYX = VAIMZY * VAJZX
                     VQZZYY = VAIMZY * VAJZY
                     VQZZYZ = VAIMZY * VAJZZ
                     VQZZZX = VAIMZZ * VAJZX
                     VQZZZY = VAIMZZ * VAJZY
                     VQZZZZ = VAIMZZ * VAJZZ
C
C                    TERM I
C
C  TERM1 R10
           TERM1 =   (RIJX * RIJX * RIJX * RIJX * VQXXXX +
     &                 RIJX * RIJX * RIJX * RIJY * VQXXXY +
     &                 RIJX * RIJX * RIJX * RIJZ * VQXXXZ +
     &                 RIJX * RIJX * RIJY * RIJX * VQXXYX +
     &                 RIJX * RIJX * RIJY * RIJY * VQXXYY +
     &                 RIJX * RIJX * RIJY * RIJZ * VQXXYZ +
     &                 RIJX * RIJX * RIJZ * RIJX * VQXXZX +
     &                 RIJX * RIJX * RIJZ * RIJY * VQXXZY +
     &                 RIJX * RIJX * RIJZ * RIJZ * VQXXZZ +
     &                 RIJX * RIJY * RIJX * RIJX * VQXYXX +
     &                 RIJX * RIJY * RIJX * RIJY * VQXYXY +
     &                 RIJX * RIJY * RIJX * RIJZ * VQXYXZ +
     &                 RIJX * RIJY * RIJY * RIJX * VQXYYX +
     &                 RIJX * RIJY * RIJY * RIJY * VQXYYY +
     &                 RIJX * RIJY * RIJY * RIJZ * VQXYYZ +
     &                 RIJX * RIJY * RIJZ * RIJX * VQXYZX +
     &                 RIJX * RIJY * RIJZ * RIJY * VQXYZY +
     &                 RIJX * RIJY * RIJZ * RIJZ * VQXYZZ +
     &                 RIJX * RIJZ * RIJX * RIJX * VQXZXX +
     &                 RIJX * RIJZ * RIJX * RIJY * VQXZXY +
     &                 RIJX * RIJZ * RIJX * RIJZ * VQXZXZ +
     &                 RIJX * RIJZ * RIJY * RIJX * VQXZYX +
     &                 RIJX * RIJZ * RIJY * RIJY * VQXZYY +
     &                 RIJX * RIJZ * RIJY * RIJZ * VQXZYZ +
     &                 RIJX * RIJZ * RIJZ * RIJX * VQXZZX +
     &                 RIJX * RIJZ * RIJZ * RIJY * VQXZZY +
     &                 RIJX * RIJZ * RIJZ * RIJZ * VQXZZZ +
     &                 RIJY * RIJX * RIJX * RIJX * VQYXXX +
     &                 RIJY * RIJX * RIJX * RIJY * VQYXXY +
     &                 RIJY * RIJX * RIJX * RIJZ * VQYXXZ +
     &                 RIJY * RIJX * RIJY * RIJX * VQYXYX +
     &                 RIJY * RIJX * RIJY * RIJY * VQYXYY +
     &                 RIJY * RIJX * RIJY * RIJZ * VQYXYZ +
     &                 RIJY * RIJX * RIJZ * RIJX * VQYXZX +
     &                 RIJY * RIJX * RIJZ * RIJY * VQYXZY +
     &                 RIJY * RIJX * RIJZ * RIJZ * VQYXZZ +
     &                 RIJY * RIJY * RIJX * RIJX * VQYYXX +
     &                 RIJY * RIJY * RIJX * RIJY * VQYYXY +
     &                 RIJY * RIJY * RIJX * RIJZ * VQYYXZ +
     &                 RIJY * RIJY * RIJY * RIJX * VQYYYX +
     &                 RIJY * RIJY * RIJY * RIJY * VQYYYY +
     &                 RIJY * RIJY * RIJY * RIJZ * VQYYYZ +
     &                 RIJY * RIJY * RIJZ * RIJX * VQYYZX +
     &                 RIJY * RIJY * RIJZ * RIJY * VQYYZY +
     &                 RIJY * RIJY * RIJZ * RIJZ * VQYYZZ +
     &                 RIJY * RIJZ * RIJX * RIJX * VQYZXX +
     &                 RIJY * RIJZ * RIJX * RIJY * VQYZXY +
     &                 RIJY * RIJZ * RIJX * RIJZ * VQYZXZ +
     &                 RIJY * RIJZ * RIJY * RIJX * VQYZYX +
     &                 RIJY * RIJZ * RIJY * RIJY * VQYZYY +
     &                 RIJY * RIJZ * RIJY * RIJZ * VQYZYZ +
     &                 RIJY * RIJZ * RIJZ * RIJX * VQYZZX +
     &                 RIJY * RIJZ * RIJZ * RIJY * VQYZZY +
     &                 RIJY * RIJZ * RIJZ * RIJZ * VQYZZZ +
     &                 RIJZ * RIJX * RIJX * RIJX * VQZXXX +
     &                 RIJZ * RIJX * RIJX * RIJY * VQZXXY +
     &                 RIJZ * RIJX * RIJX * RIJZ * VQZXXZ +
     &                 RIJZ * RIJX * RIJY * RIJX * VQZXYX +
     &                 RIJZ * RIJX * RIJY * RIJY * VQZXYY +
     &                 RIJZ * RIJX * RIJY * RIJZ * VQZXYZ +
     &                 RIJZ * RIJX * RIJZ * RIJX * VQZXZX +
     &                 RIJZ * RIJX * RIJZ * RIJY * VQZXZY +
     &                 RIJZ * RIJX * RIJZ * RIJZ * VQZXZZ +
     &                 RIJZ * RIJY * RIJX * RIJX * VQZYXX +
     &                 RIJZ * RIJY * RIJX * RIJY * VQZYXY +
     &                 RIJZ * RIJY * RIJX * RIJZ * VQZYXZ +
     &                 RIJZ * RIJY * RIJY * RIJX * VQZYYX +
     &                 RIJZ * RIJY * RIJY * RIJY * VQZYYY +
     &                 RIJZ * RIJY * RIJY * RIJZ * VQZYYZ +
     &                 RIJZ * RIJY * RIJZ * RIJX * VQZYZX +
     &                 RIJZ * RIJY * RIJZ * RIJY * VQZYZY +
     &                 RIJZ * RIJY * RIJZ * RIJZ * VQZYZZ +
     &                 RIJZ * RIJZ * RIJX * RIJX * VQZZXX +
     &                 RIJZ * RIJZ * RIJX * RIJY * VQZZXY +
     &                 RIJZ * RIJZ * RIJX * RIJZ * VQZZXZ +
     &                 RIJZ * RIJZ * RIJY * RIJX * VQZZYX +
     &                 RIJZ * RIJZ * RIJY * RIJY * VQZZYY +
     &                 RIJZ * RIJZ * RIJY * RIJZ * VQZZYZ +
     &                 RIJZ * RIJZ * RIJZ * RIJX * VQZZZX +
     &                 RIJZ * RIJZ * RIJZ * RIJY * VQZZZY +
     &                 RIJZ * RIJZ * RIJZ * RIJZ * VQZZZZ )
     &                 * NINE / RIJ10
     &             - ((RIJX * RIJX + RIJX * RIJX ) * VQXXXX +
     &                 RIJX * RIJY * VQXXXY +
     &                 RIJX * RIJZ * VQXXXZ +
     &                 RIJY * RIJX * VQXXYX +
     &                (RIJX * RIJX + RIJY * RIJY ) * VQXXYY +
     &                 RIJY * RIJZ * VQXXYZ +
     &                 RIJZ * RIJX * VQXXZX +
     &                 RIJZ * RIJY * VQXXZY +
     &                (RIJX * RIJX + RIJZ * RIJZ ) * VQXXZZ +
     &                 RIJX * RIJY * VQXYXX +
     &                 RIJX * RIJY * VQXYYY +
     &                 RIJX * RIJY * VQXYZZ +
     &                 RIJX * RIJZ * VQXZXX +
     &                 RIJX * RIJZ * VQXZYY +
     &                 RIJX * RIJZ * VQXZZZ +
     &                 RIJY * RIJX * VQYXXX +
     &                 RIJY * RIJX * VQYXYY +
     &                 RIJY * RIJX * VQYXZZ +
     &                (RIJY * RIJY + RIJX * RIJX ) * VQYYXX +
     &                 RIJX * RIJY * VQYYXY +
     &                 RIJX * RIJZ * VQYYXZ +
     &                 RIJY * RIJX * VQYYYX +
     &                (RIJY * RIJY + RIJY * RIJY ) * VQYYYY +
     &                 RIJY * RIJZ * VQYYYZ +
     &                 RIJZ * RIJX * VQYYZX +
     &                 RIJZ * RIJY * VQYYZY +
     &                (RIJY * RIJY + RIJZ * RIJZ ) * VQYYZZ +
     &                 RIJY * RIJZ * VQYZXX +
     &                 RIJY * RIJZ * VQYZYY +
     &                 RIJY * RIJZ * VQYZZZ +
     &                 RIJZ * RIJX * VQZXXX +
     &                 RIJZ * RIJX * VQZXYY +
     &                 RIJZ * RIJX * VQZXZZ +
     &                 RIJZ * RIJY * VQZYXX +
     &                 RIJZ * RIJY * VQZYYY +
     &                 RIJZ * RIJY * VQZYZZ +
     &                (RIJZ * RIJZ + RIJX * RIJX ) * VQZZXX +
     &                 RIJX * RIJY * VQZZXY +
     &                 RIJX * RIJZ * VQZZXZ +
     &                 RIJY * RIJX * VQZZYX +
     &                (RIJZ * RIJZ + RIJY * RIJY ) * VQZZYY +
     &                 RIJY * RIJZ * VQZZYZ +
     &                 RIJZ * RIJX * VQZZZX +
     &                 RIJZ * RIJY * VQZZZY +
     &                (RIJZ * RIJZ + RIJZ * RIJZ ) * VQZZZZ )
     &                * THREE / RIJ8
     &             + (VQXXXX + VQXXYY + VQXXZZ + VQYYXX +
     &                VQYYYY + VQYYZZ + VQZZXX + VQZZYY +
     &                VQZZZZ ) / RIJ6
C
C                    TERM II
C
                     TERM2 = ZERO
                     ...
C ---------------------------------------------------------------------------------
C
C                    ACCUMULATE DISPERSION FREQUENCY SHIFT
C
                     DISP = DISP + G12W * GRF * (TERM1 + TERM2)
C                                                                           
                  ENDDO
C                                                                           
               ENDDO
C
            ENDDO
         ENDDO
C
      ENDDO
      ENDDO
C
      DISP = DISP / (TWO * PI)
      DISP = DISP / (TWO * REDMSS(MODE) * FREQ(MODE))
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

