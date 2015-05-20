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
           TERM1 =   (RIJX * RIJX * RIJX * RIJX * VQXXXX +
     &                RIJX * RIJX * RIJX * RIJY * VQXXXY +
     &                RIJX * RIJX * RIJX * RIJZ * VQXXXZ +
     &                RIJX * RIJX * RIJY * RIJX * VQXXYX +
     &                RIJX * RIJX * RIJY * RIJY * VQXXYY +
     &                RIJX * RIJX * RIJY * RIJZ * VQXXYZ +
     &                RIJX * RIJX * RIJZ * RIJX * VQXXZX +
     &                RIJX * RIJX * RIJZ * RIJY * VQXXZY +
     &                RIJX * RIJX * RIJZ * RIJZ * VQXXZZ +
     &                RIJX * RIJY * RIJX * RIJX * VQXYXX +
     &                RIJX * RIJY * RIJX * RIJY * VQXYXY +
     &                RIJX * RIJY * RIJX * RIJZ * VQXYXZ +
     &                RIJX * RIJY * RIJY * RIJX * VQXYYX +
     &                RIJX * RIJY * RIJY * RIJY * VQXYYY +
     &                RIJX * RIJY * RIJY * RIJZ * VQXYYZ +
     &                RIJX * RIJY * RIJZ * RIJX * VQXYZX +
     &                RIJX * RIJY * RIJZ * RIJY * VQXYZY +
     &                RIJX * RIJY * RIJZ * RIJZ * VQXYZZ +
     &                RIJX * RIJZ * RIJX * RIJX * VQXZXX +
     &                RIJX * RIJZ * RIJX * RIJY * VQXZXY +
     &                RIJX * RIJZ * RIJX * RIJZ * VQXZXZ +
     &                RIJX * RIJZ * RIJY * RIJX * VQXZYX +
     &                RIJX * RIJZ * RIJY * RIJY * VQXZYY +
     &                RIJX * RIJZ * RIJY * RIJZ * VQXZYZ +
     &                RIJX * RIJZ * RIJZ * RIJX * VQXZZX +
     &                RIJX * RIJZ * RIJZ * RIJY * VQXZZY +
     &                RIJX * RIJZ * RIJZ * RIJZ * VQXZZZ +
     &                RIJY * RIJX * RIJX * RIJX * VQYXXX +
     &                RIJY * RIJX * RIJX * RIJY * VQYXXY +
     &                RIJY * RIJX * RIJX * RIJZ * VQYXXZ +
     &                RIJY * RIJX * RIJY * RIJX * VQYXYX +
     &                RIJY * RIJX * RIJY * RIJY * VQYXYY +
     &                RIJY * RIJX * RIJY * RIJZ * VQYXYZ +
     &                RIJY * RIJX * RIJZ * RIJX * VQYXZX +
     &                RIJY * RIJX * RIJZ * RIJY * VQYXZY +
     &                RIJY * RIJX * RIJZ * RIJZ * VQYXZZ +
     &                RIJY * RIJY * RIJX * RIJX * VQYYXX +
     &                RIJY * RIJY * RIJX * RIJY * VQYYXY +
     &                RIJY * RIJY * RIJX * RIJZ * VQYYXZ +
     &                RIJY * RIJY * RIJY * RIJX * VQYYYX +
     &                RIJY * RIJY * RIJY * RIJY * VQYYYY +
     &                RIJY * RIJY * RIJY * RIJZ * VQYYYZ +
     &                RIJY * RIJY * RIJZ * RIJX * VQYYZX +
     &                RIJY * RIJY * RIJZ * RIJY * VQYYZY +
     &                RIJY * RIJY * RIJZ * RIJZ * VQYYZZ +
     &                RIJY * RIJZ * RIJX * RIJX * VQYZXX +
     &                RIJY * RIJZ * RIJX * RIJY * VQYZXY +
     &                RIJY * RIJZ * RIJX * RIJZ * VQYZXZ +
     &                RIJY * RIJZ * RIJY * RIJX * VQYZYX +
     &                RIJY * RIJZ * RIJY * RIJY * VQYZYY +
     &                RIJY * RIJZ * RIJY * RIJZ * VQYZYZ +
     &                RIJY * RIJZ * RIJZ * RIJX * VQYZZX +
     &                RIJY * RIJZ * RIJZ * RIJY * VQYZZY +
     &                RIJY * RIJZ * RIJZ * RIJZ * VQYZZZ +
     &                RIJZ * RIJX * RIJX * RIJX * VQZXXX +
     &                RIJZ * RIJX * RIJX * RIJY * VQZXXY +
     &                RIJZ * RIJX * RIJX * RIJZ * VQZXXZ +
     &                RIJZ * RIJX * RIJY * RIJX * VQZXYX +
     &                RIJZ * RIJX * RIJY * RIJY * VQZXYY +
     &                RIJZ * RIJX * RIJY * RIJZ * VQZXYZ +
     &                RIJZ * RIJX * RIJZ * RIJX * VQZXZX +
     &                RIJZ * RIJX * RIJZ * RIJY * VQZXZY +
     &                RIJZ * RIJX * RIJZ * RIJZ * VQZXZZ +
     &                RIJZ * RIJY * RIJX * RIJX * VQZYXX +
     &                RIJZ * RIJY * RIJX * RIJY * VQZYXY +
     &                RIJZ * RIJY * RIJX * RIJZ * VQZYXZ +
     &                RIJZ * RIJY * RIJY * RIJX * VQZYYX +
     &                RIJZ * RIJY * RIJY * RIJY * VQZYYY +
     &                RIJZ * RIJY * RIJY * RIJZ * VQZYYZ +
     &                RIJZ * RIJY * RIJZ * RIJX * VQZYZX +
     &                RIJZ * RIJY * RIJZ * RIJY * VQZYZY +
     &                RIJZ * RIJY * RIJZ * RIJZ * VQZYZZ +
     &                RIJZ * RIJZ * RIJX * RIJX * VQZZXX +
     &                RIJZ * RIJZ * RIJX * RIJY * VQZZXY +
     &                RIJZ * RIJZ * RIJX * RIJZ * VQZZXZ +
     &                RIJZ * RIJZ * RIJY * RIJX * VQZZYX +
     &                RIJZ * RIJZ * RIJY * RIJY * VQZZYY +
     &                RIJZ * RIJZ * RIJY * RIJZ * VQZZYZ +
     &                RIJZ * RIJZ * RIJZ * RIJX * VQZZZX +
     &                RIJZ * RIJZ * RIJZ * RIJY * VQZZZY +
     &                RIJZ * RIJZ * RIJZ * RIJZ * VQZZZZ )
     &                * NINE / RIJ10
     &            - ((RIJX * RIJX + RIJX * RIJX ) * VQXXXX +
     &                RIJX * RIJY * VQXXXY +
     &                RIJX * RIJZ * VQXXXZ +
     &                RIJY * RIJX * VQXXYX +
     &               (RIJX * RIJX + RIJY * RIJY ) * VQXXYY +
     &                RIJY * RIJZ * VQXXYZ +
     &                RIJZ * RIJX * VQXXZX +
     &                RIJZ * RIJY * VQXXZY +
     &               (RIJX * RIJX + RIJZ * RIJZ ) * VQXXZZ +
     &                RIJX * RIJY * VQXYXX +
     &                RIJX * RIJY * VQXYYY +
     &                RIJX * RIJY * VQXYZZ +
     &                RIJX * RIJZ * VQXZXX +
     &                RIJX * RIJZ * VQXZYY +
     &                RIJX * RIJZ * VQXZZZ +
     &                RIJY * RIJX * VQYXXX +
     &                RIJY * RIJX * VQYXYY +
     &                RIJY * RIJX * VQYXZZ +
     &               (RIJY * RIJY + RIJX * RIJX ) * VQYYXX +
     &                RIJX * RIJY * VQYYXY +
     &                RIJX * RIJZ * VQYYXZ +
     &                RIJY * RIJX * VQYYYX +
     &               (RIJY * RIJY + RIJY * RIJY ) * VQYYYY +
     &                RIJY * RIJZ * VQYYYZ +
     &                RIJZ * RIJX * VQYYZX +
     &                RIJZ * RIJY * VQYYZY +
     &               (RIJY * RIJY + RIJZ * RIJZ ) * VQYYZZ +
     &                RIJY * RIJZ * VQYZXX +
     &                RIJY * RIJZ * VQYZYY +
     &                RIJY * RIJZ * VQYZZZ +
     &                RIJZ * RIJX * VQZXXX +
     &                RIJZ * RIJX * VQZXYY +
     &                RIJZ * RIJX * VQZXZZ +
     &                RIJZ * RIJY * VQZYXX +
     &                RIJZ * RIJY * VQZYYY +
     &                RIJZ * RIJY * VQZYZZ +
     &               (RIJZ * RIJZ + RIJX * RIJX ) * VQZZXX +
     &                RIJX * RIJY * VQZZXY +
     &                RIJX * RIJZ * VQZZXZ +
     &                RIJY * RIJX * VQZZYX +
     &               (RIJZ * RIJZ + RIJY * RIJY ) * VQZZYY +
     &                RIJY * RIJZ * VQZZYZ +
     &                RIJZ * RIJX * VQZZZX +
     &                RIJZ * RIJY * VQZZZY +
     &               (RIJZ * RIJZ + RIJZ * RIJZ ) * VQZZZZ )
     &               * THREE / RIJ8
     &            + (VQXXXX + VQXXYY + VQXXZZ + VQYYXX +
     &               VQYYYY + VQYYZZ + VQZZXX + VQZZYY +
     &               VQZZZZ ) / RIJ6
C
C                    TERM II
C
           TERM2 =   ((RX1 * RIJX * RIJX * RIJX +
     &                 RIJX * RX1 * RIJX * RIJX +
     &                 RIJX * RIJX * RX1 * RIJX +
     &                 RIJX * RIJX * RIJX * RX1 ) * VIXXXX +
     &                (RX1 * RIJX * RIJX * RIJY +
     &                 RIJX * RX1 * RIJX * RIJY +
     &                 RIJX * RIJX * RX1 * RIJY +
     &                 RIJX * RIJX * RIJX * RY1 ) * VIXXXY +
     &                (RX1 * RIJX * RIJX * RIJZ +
     &                 RIJX * RX1 * RIJX * RIJZ +
     &                 RIJX * RIJX * RX1 * RIJZ +
     &                 RIJX * RIJX * RIJX * RZ1 ) * VIXXXZ +
     &                (RX1 * RIJX * RIJY * RIJX +
     &                 RIJX * RX1 * RIJY * RIJX +
     &                 RIJX * RIJX * RY1 * RIJX +
     &                 RIJX * RIJX * RIJY * RX1 ) * VIXXYX +
     &                (RX1 * RIJX * RIJY * RIJY +
     &                 RIJX * RX1 * RIJY * RIJY +
     &                 RIJX * RIJX * RY1 * RIJY +
     &                 RIJX * RIJX * RIJY * RY1 ) * VIXXYY +
     &                (RX1 * RIJX * RIJY * RIJZ +
     &                 RIJX * RX1 * RIJY * RIJZ +
     &                 RIJX * RIJX * RY1 * RIJZ +
     &                 RIJX * RIJX * RIJY * RZ1 ) * VIXXYZ +
     &                (RX1 * RIJX * RIJZ * RIJX +
     &                 RIJX * RX1 * RIJZ * RIJX +
     &                 RIJX * RIJX * RZ1 * RIJX +
     &                 RIJX * RIJX * RIJZ * RX1 ) * VIXXZX +
     &                (RX1 * RIJX * RIJZ * RIJY +
     &                 RIJX * RX1 * RIJZ * RIJY +
     &                 RIJX * RIJX * RZ1 * RIJY +
     &                 RIJX * RIJX * RIJZ * RY1 ) * VIXXZY +
     &                (RX1 * RIJX * RIJZ * RIJZ +
     &                 RIJX * RX1 * RIJZ * RIJZ +
     &                 RIJX * RIJX * RZ1 * RIJZ +
     &                 RIJX * RIJX * RIJZ * RZ1 ) * VIXXZZ +
     &                (RX1 * RIJY * RIJX * RIJX +
     &                 RIJX * RY1 * RIJX * RIJX +
     &                 RIJX * RIJY * RX1 * RIJX +
     &                 RIJX * RIJY * RIJX * RX1 ) * VIXYXX +
     &                (RX1 * RIJY * RIJX * RIJY +
     &                 RIJX * RY1 * RIJX * RIJY +
     &                 RIJX * RIJY * RX1 * RIJY +
     &                 RIJX * RIJY * RIJX * RY1 ) * VIXYXY +
     &                (RX1 * RIJY * RIJX * RIJZ +
     &                 RIJX * RY1 * RIJX * RIJZ +
     &                 RIJX * RIJY * RX1 * RIJZ +
     &                 RIJX * RIJY * RIJX * RZ1 ) * VIXYXZ +
     &                (RX1 * RIJY * RIJY * RIJX +
     &                 RIJX * RY1 * RIJY * RIJX +
     &                 RIJX * RIJY * RY1 * RIJX +
     &                 RIJX * RIJY * RIJY * RX1 ) * VIXYYX +
     &                (RX1 * RIJY * RIJY * RIJY +
     &                 RIJX * RY1 * RIJY * RIJY +
     &                 RIJX * RIJY * RY1 * RIJY +
     &                 RIJX * RIJY * RIJY * RY1 ) * VIXYYY +
     &                (RX1 * RIJY * RIJY * RIJZ +
     &                 RIJX * RY1 * RIJY * RIJZ +
     &                 RIJX * RIJY * RY1 * RIJZ +
     &                 RIJX * RIJY * RIJY * RZ1 ) * VIXYYZ +
     &                (RX1 * RIJY * RIJZ * RIJX +
     &                 RIJX * RY1 * RIJZ * RIJX +
     &                 RIJX * RIJY * RZ1 * RIJX +
     &                 RIJX * RIJY * RIJZ * RX1 ) * VIXYZX +
     &                (RX1 * RIJY * RIJZ * RIJY +
     &                 RIJX * RY1 * RIJZ * RIJY +
     &                 RIJX * RIJY * RZ1 * RIJY +
     &                 RIJX * RIJY * RIJZ * RY1 ) * VIXYZY +
     &                (RX1 * RIJY * RIJZ * RIJZ +
     &                 RIJX * RY1 * RIJZ * RIJZ +
     &                 RIJX * RIJY * RZ1 * RIJZ +
     &                 RIJX * RIJY * RIJZ * RZ1 ) * VIXYZZ +
     &                (RX1 * RIJZ * RIJX * RIJX +
     &                 RIJX * RZ1 * RIJX * RIJX +
     &                 RIJX * RIJZ * RX1 * RIJX +
     &                 RIJX * RIJZ * RIJX * RX1 ) * VIXZXX +
     &                (RX1 * RIJZ * RIJX * RIJY +
     &                 RIJX * RZ1 * RIJX * RIJY +
     &                 RIJX * RIJZ * RX1 * RIJY +
     &                 RIJX * RIJZ * RIJX * RY1 ) * VIXZXY +
     &                (RX1 * RIJZ * RIJX * RIJZ +
     &                 RIJX * RZ1 * RIJX * RIJZ +
     &                 RIJX * RIJZ * RX1 * RIJZ +
     &                 RIJX * RIJZ * RIJX * RZ1 ) * VIXZXZ +
     &                (RX1 * RIJZ * RIJY * RIJX +
     &                 RIJX * RZ1 * RIJY * RIJX +
     &                 RIJX * RIJZ * RY1 * RIJX +
     &                 RIJX * RIJZ * RIJY * RX1 ) * VIXZYX +
     &                (RX1 * RIJZ * RIJY * RIJY +
     &                 RIJX * RZ1 * RIJY * RIJY +
     &                 RIJX * RIJZ * RY1 * RIJY +
     &                 RIJX * RIJZ * RIJY * RY1 ) * VIXZYY +
     &                (RX1 * RIJZ * RIJY * RIJZ +
     &                 RIJX * RZ1 * RIJY * RIJZ +
     &                 RIJX * RIJZ * RY1 * RIJZ +
     &                 RIJX * RIJZ * RIJY * RZ1 ) * VIXZYZ +
     &                (RX1 * RIJZ * RIJZ * RIJX +
     &                 RIJX * RZ1 * RIJZ * RIJX +
     &                 RIJX * RIJZ * RZ1 * RIJX +
     &                 RIJX * RIJZ * RIJZ * RX1 ) * VIXZZX +
     &                (RX1 * RIJZ * RIJZ * RIJY +
     &                 RIJX * RZ1 * RIJZ * RIJY +
     &                 RIJX * RIJZ * RZ1 * RIJY +
     &                 RIJX * RIJZ * RIJZ * RY1 ) * VIXZZY +
     &                (RX1 * RIJZ * RIJZ * RIJZ +
     &                 RIJX * RZ1 * RIJZ * RIJZ +
     &                 RIJX * RIJZ * RZ1 * RIJZ +
     &                 RIJX * RIJZ * RIJZ * RZ1 ) * VIXZZZ +
     &                (RY1 * RIJX * RIJX * RIJX +
     &                 RIJY * RX1 * RIJX * RIJX +
     &                 RIJY * RIJX * RX1 * RIJX +
     &                 RIJY * RIJX * RIJX * RX1 ) * VIYXXX +
     &                (RY1 * RIJX * RIJX * RIJY +
     &                 RIJY * RX1 * RIJX * RIJY +
     &                 RIJY * RIJX * RX1 * RIJY +
     &                 RIJY * RIJX * RIJX * RY1 ) * VIYXXY +
     &                (RY1 * RIJX * RIJX * RIJZ +
     &                 RIJY * RX1 * RIJX * RIJZ +
     &                 RIJY * RIJX * RX1 * RIJZ +
     &                 RIJY * RIJX * RIJX * RZ1 ) * VIYXXZ +
     &                (RY1 * RIJX * RIJY * RIJX +
     &                 RIJY * RX1 * RIJY * RIJX +
     &                 RIJY * RIJX * RY1 * RIJX +
     &                 RIJY * RIJX * RIJY * RX1 ) * VIYXYX +
     &                (RY1 * RIJX * RIJY * RIJY +
     &                 RIJY * RX1 * RIJY * RIJY +
     &                 RIJY * RIJX * RY1 * RIJY +
     &                 RIJY * RIJX * RIJY * RY1 ) * VIYXYY +
     &                (RY1 * RIJX * RIJY * RIJZ +
     &                 RIJY * RX1 * RIJY * RIJZ +
     &                 RIJY * RIJX * RY1 * RIJZ +
     &                 RIJY * RIJX * RIJY * RZ1 ) * VIYXYZ +
     &                (RY1 * RIJX * RIJZ * RIJX +
     &                 RIJY * RX1 * RIJZ * RIJX +
     &                 RIJY * RIJX * RZ1 * RIJX +
     &                 RIJY * RIJX * RIJZ * RX1 ) * VIYXZX +
     &                (RY1 * RIJX * RIJZ * RIJY +
     &                 RIJY * RX1 * RIJZ * RIJY +
     &                 RIJY * RIJX * RZ1 * RIJY +
     &                 RIJY * RIJX * RIJZ * RY1 ) * VIYXZY +
     &                (RY1 * RIJX * RIJZ * RIJZ +
     &                 RIJY * RX1 * RIJZ * RIJZ +
     &                 RIJY * RIJX * RZ1 * RIJZ +
     &                 RIJY * RIJX * RIJZ * RZ1 ) * VIYXZZ +
     &                (RY1 * RIJY * RIJX * RIJX +
     &                 RIJY * RY1 * RIJX * RIJX +
     &                 RIJY * RIJY * RX1 * RIJX +
     &                 RIJY * RIJY * RIJX * RX1 ) * VIYYXX +
     &                (RY1 * RIJY * RIJX * RIJY +
     &                 RIJY * RY1 * RIJX * RIJY +
     &                 RIJY * RIJY * RX1 * RIJY +
     &                 RIJY * RIJY * RIJX * RY1 ) * VIYYXY +
     &                (RY1 * RIJY * RIJX * RIJZ +
     &                 RIJY * RY1 * RIJX * RIJZ +
     &                 RIJY * RIJY * RX1 * RIJZ +
     &                 RIJY * RIJY * RIJX * RZ1 ) * VIYYXZ +
     &                (RY1 * RIJY * RIJY * RIJX +
     &                 RIJY * RY1 * RIJY * RIJX +
     &                 RIJY * RIJY * RY1 * RIJX +
     &                 RIJY * RIJY * RIJY * RX1 ) * VIYYYX +
     &                (RY1 * RIJY * RIJY * RIJY +
     &                 RIJY * RY1 * RIJY * RIJY +
     &                 RIJY * RIJY * RY1 * RIJY +
     &                 RIJY * RIJY * RIJY * RY1 ) * VIYYYY +
     &                (RY1 * RIJY * RIJY * RIJZ +
     &                 RIJY * RY1 * RIJY * RIJZ +
     &                 RIJY * RIJY * RY1 * RIJZ +
     &                 RIJY * RIJY * RIJY * RZ1 ) * VIYYYZ +
     &                (RY1 * RIJY * RIJZ * RIJX +
     &                 RIJY * RY1 * RIJZ * RIJX +
     &                 RIJY * RIJY * RZ1 * RIJX +
     &                 RIJY * RIJY * RIJZ * RX1 ) * VIYYZX +
     &                (RY1 * RIJY * RIJZ * RIJY +
     &                 RIJY * RY1 * RIJZ * RIJY +
     &                 RIJY * RIJY * RZ1 * RIJY +
     &                 RIJY * RIJY * RIJZ * RY1 ) * VIYYZY +
     &                (RY1 * RIJY * RIJZ * RIJZ +
     &                 RIJY * RY1 * RIJZ * RIJZ +
     &                 RIJY * RIJY * RZ1 * RIJZ +
     &                 RIJY * RIJY * RIJZ * RZ1 ) * VIYYZZ +
     &                (RY1 * RIJZ * RIJX * RIJX +
     &                 RIJY * RZ1 * RIJX * RIJX +
     &                 RIJY * RIJZ * RX1 * RIJX +
     &                 RIJY * RIJZ * RIJX * RX1 ) * VIYZXX +
     &                (RY1 * RIJZ * RIJX * RIJY +
     &                 RIJY * RZ1 * RIJX * RIJY +
     &                 RIJY * RIJZ * RX1 * RIJY +
     &                 RIJY * RIJZ * RIJX * RY1 ) * VIYZXY +
     &                (RY1 * RIJZ * RIJX * RIJZ +
     &                 RIJY * RZ1 * RIJX * RIJZ +
     &                 RIJY * RIJZ * RX1 * RIJZ +
     &                 RIJY * RIJZ * RIJX * RZ1 ) * VIYZXZ +
     &                (RY1 * RIJZ * RIJY * RIJX +
     &                 RIJY * RZ1 * RIJY * RIJX +
     &                 RIJY * RIJZ * RY1 * RIJX +
     &                 RIJY * RIJZ * RIJY * RX1 ) * VIYZYX +
     &                (RY1 * RIJZ * RIJY * RIJY +
     &                 RIJY * RZ1 * RIJY * RIJY +
     &                 RIJY * RIJZ * RY1 * RIJY +
     &                 RIJY * RIJZ * RIJY * RY1 ) * VIYZYY +
     &                (RY1 * RIJZ * RIJY * RIJZ +
     &                 RIJY * RZ1 * RIJY * RIJZ +
     &                 RIJY * RIJZ * RY1 * RIJZ +
     &                 RIJY * RIJZ * RIJY * RZ1 ) * VIYZYZ +
     &                (RY1 * RIJZ * RIJZ * RIJX +
     &                 RIJY * RZ1 * RIJZ * RIJX +
     &                 RIJY * RIJZ * RZ1 * RIJX +
     &                 RIJY * RIJZ * RIJZ * RX1 ) * VIYZZX +
     &                (RY1 * RIJZ * RIJZ * RIJY +
     &                 RIJY * RZ1 * RIJZ * RIJY +
     &                 RIJY * RIJZ * RZ1 * RIJY +
     &                 RIJY * RIJZ * RIJZ * RY1 ) * VIYZZY +
     &                (RY1 * RIJZ * RIJZ * RIJZ +
     &                 RIJY * RZ1 * RIJZ * RIJZ +
     &                 RIJY * RIJZ * RZ1 * RIJZ +
     &                 RIJY * RIJZ * RIJZ * RZ1 ) * VIYZZZ +
     &                (RZ1 * RIJX * RIJX * RIJX +
     &                 RIJZ * RX1 * RIJX * RIJX +
     &                 RIJZ * RIJX * RX1 * RIJX +
     &                 RIJZ * RIJX * RIJX * RX1 ) * VIZXXX +
     &                (RZ1 * RIJX * RIJX * RIJY +
     &                 RIJZ * RX1 * RIJX * RIJY +
     &                 RIJZ * RIJX * RX1 * RIJY +
     &                 RIJZ * RIJX * RIJX * RY1 ) * VIZXXY +
     &                (RZ1 * RIJX * RIJX * RIJZ +
     &                 RIJZ * RX1 * RIJX * RIJZ +
     &                 RIJZ * RIJX * RX1 * RIJZ +
     &                 RIJZ * RIJX * RIJX * RZ1 ) * VIZXXZ +
     &                (RZ1 * RIJX * RIJY * RIJX +
     &                 RIJZ * RX1 * RIJY * RIJX +
     &                 RIJZ * RIJX * RY1 * RIJX +
     &                 RIJZ * RIJX * RIJY * RX1 ) * VIZXYX +
     &                (RZ1 * RIJX * RIJY * RIJY +
     &                 RIJZ * RX1 * RIJY * RIJY +
     &                 RIJZ * RIJX * RY1 * RIJY +
     &                 RIJZ * RIJX * RIJY * RY1 ) * VIZXYY +
     &                (RZ1 * RIJX * RIJY * RIJZ +
     &                 RIJZ * RX1 * RIJY * RIJZ +
     &                 RIJZ * RIJX * RY1 * RIJZ +
     &                 RIJZ * RIJX * RIJY * RZ1 ) * VIZXYZ +
     &                (RZ1 * RIJX * RIJZ * RIJX +
     &                 RIJZ * RX1 * RIJZ * RIJX +
     &                 RIJZ * RIJX * RZ1 * RIJX +
     &                 RIJZ * RIJX * RIJZ * RX1 ) * VIZXZX +
     &                (RZ1 * RIJX * RIJZ * RIJY +
     &                 RIJZ * RX1 * RIJZ * RIJY +
     &                 RIJZ * RIJX * RZ1 * RIJY +
     &                 RIJZ * RIJX * RIJZ * RY1 ) * VIZXZY +
     &                (RZ1 * RIJX * RIJZ * RIJZ +
     &                 RIJZ * RX1 * RIJZ * RIJZ +
     &                 RIJZ * RIJX * RZ1 * RIJZ +
     &                 RIJZ * RIJX * RIJZ * RZ1 ) * VIZXZZ +
     &                (RZ1 * RIJY * RIJX * RIJX +
     &                 RIJZ * RY1 * RIJX * RIJX +
     &                 RIJZ * RIJY * RX1 * RIJX +
     &                 RIJZ * RIJY * RIJX * RX1 ) * VIZYXX +
     &                (RZ1 * RIJY * RIJX * RIJY +
     &                 RIJZ * RY1 * RIJX * RIJY +
     &                 RIJZ * RIJY * RX1 * RIJY +
     &                 RIJZ * RIJY * RIJX * RY1 ) * VIZYXY +
     &                (RZ1 * RIJY * RIJX * RIJZ +
     &                 RIJZ * RY1 * RIJX * RIJZ +
     &                 RIJZ * RIJY * RX1 * RIJZ +
     &                 RIJZ * RIJY * RIJX * RZ1 ) * VIZYXZ +
     &                (RZ1 * RIJY * RIJY * RIJX +
     &                 RIJZ * RY1 * RIJY * RIJX +
     &                 RIJZ * RIJY * RY1 * RIJX +
     &                 RIJZ * RIJY * RIJY * RX1 ) * VIZYYX +
     &                (RZ1 * RIJY * RIJY * RIJY +
     &                 RIJZ * RY1 * RIJY * RIJY +
     &                 RIJZ * RIJY * RY1 * RIJY +
     &                 RIJZ * RIJY * RIJY * RY1 ) * VIZYYY +
     &                (RZ1 * RIJY * RIJY * RIJZ +
     &                 RIJZ * RY1 * RIJY * RIJZ +
     &                 RIJZ * RIJY * RY1 * RIJZ +
     &                 RIJZ * RIJY * RIJY * RZ1 ) * VIZYYZ +
     &                (RZ1 * RIJY * RIJZ * RIJX +
     &                 RIJZ * RY1 * RIJZ * RIJX +
     &                 RIJZ * RIJY * RZ1 * RIJX +
     &                 RIJZ * RIJY * RIJZ * RX1 ) * VIZYZX +
     &                (RZ1 * RIJY * RIJZ * RIJY +
     &                 RIJZ * RY1 * RIJZ * RIJY +
     &                 RIJZ * RIJY * RZ1 * RIJY +
     &                 RIJZ * RIJY * RIJZ * RY1 ) * VIZYZY +
     &                (RZ1 * RIJY * RIJZ * RIJZ +
     &                 RIJZ * RY1 * RIJZ * RIJZ +
     &                 RIJZ * RIJY * RZ1 * RIJZ +
     &                 RIJZ * RIJY * RIJZ * RZ1 ) * VIZYZZ +
     &                (RZ1 * RIJZ * RIJX * RIJX +
     &                 RIJZ * RZ1 * RIJX * RIJX +
     &                 RIJZ * RIJZ * RX1 * RIJX +
     &                 RIJZ * RIJZ * RIJX * RX1 ) * VIZZXX +
     &                (RZ1 * RIJZ * RIJX * RIJY +
     &                 RIJZ * RZ1 * RIJX * RIJY +
     &                 RIJZ * RIJZ * RX1 * RIJY +
     &                 RIJZ * RIJZ * RIJX * RY1 ) * VIZZXY +
     &                (RZ1 * RIJZ * RIJX * RIJZ +
     &                 RIJZ * RZ1 * RIJX * RIJZ +
     &                 RIJZ * RIJZ * RX1 * RIJZ +
     &                 RIJZ * RIJZ * RIJX * RZ1 ) * VIZZXZ +
     &                (RZ1 * RIJZ * RIJY * RIJX +
     &                 RIJZ * RZ1 * RIJY * RIJX +
     &                 RIJZ * RIJZ * RY1 * RIJX +
     &                 RIJZ * RIJZ * RIJY * RX1 ) * VIZZYX +
     &                (RZ1 * RIJZ * RIJY * RIJY +
     &                 RIJZ * RZ1 * RIJY * RIJY +
     &                 RIJZ * RIJZ * RY1 * RIJY +
     &                 RIJZ * RIJZ * RIJY * RY1 ) * VIZZYY +
     &                (RZ1 * RIJZ * RIJY * RIJZ +
     &                 RIJZ * RZ1 * RIJY * RIJZ +
     &                 RIJZ * RIJZ * RY1 * RIJZ +
     &                 RIJZ * RIJZ * RIJY * RZ1 ) * VIZZYZ +
     &                (RZ1 * RIJZ * RIJZ * RIJX +
     &                 RIJZ * RZ1 * RIJZ * RIJX +
     &                 RIJZ * RIJZ * RZ1 * RIJX +
     &                 RIJZ * RIJZ * RIJZ * RX1 ) * VIZZZX +
     &                (RZ1 * RIJZ * RIJZ * RIJY +
     &                 RIJZ * RZ1 * RIJZ * RIJY +
     &                 RIJZ * RIJZ * RZ1 * RIJY +
     &                 RIJZ * RIJZ * RIJZ * RY1 ) * VIZZZY +
     &                (RZ1 * RIJZ * RIJZ * RIJZ +
     &                 RIJZ * RZ1 * RIJZ * RIJZ +
     &                 RIJZ * RIJZ * RZ1 * RIJZ +
     &                 RIJZ * RIJZ * RIJZ * RZ1 ) * VIZZZZ )
     &               * NINE / RIJ10
     &             - (RIJX * RIJX * RIJX * RIJX * VIXXXX +
     &                RIJX * RIJX * RIJX * RIJY * VIXXXY +
     &                RIJX * RIJX * RIJX * RIJZ * VIXXXZ +
     &                RIJX * RIJX * RIJY * RIJX * VIXXYX +
     &                RIJX * RIJX * RIJY * RIJY * VIXXYY +
     &                RIJX * RIJX * RIJY * RIJZ * VIXXYZ +
     &                RIJX * RIJX * RIJZ * RIJX * VIXXZX +
     &                RIJX * RIJX * RIJZ * RIJY * VIXXZY +
     &                RIJX * RIJX * RIJZ * RIJZ * VIXXZZ +
     &                RIJX * RIJY * RIJX * RIJX * VIXYXX +
     &                RIJX * RIJY * RIJX * RIJY * VIXYXY +
     &                RIJX * RIJY * RIJX * RIJZ * VIXYXZ +
     &                RIJX * RIJY * RIJY * RIJX * VIXYYX +
     &                RIJX * RIJY * RIJY * RIJY * VIXYYY +
     &                RIJX * RIJY * RIJY * RIJZ * VIXYYZ +
     &                RIJX * RIJY * RIJZ * RIJX * VIXYZX +
     &                RIJX * RIJY * RIJZ * RIJY * VIXYZY +
     &                RIJX * RIJY * RIJZ * RIJZ * VIXYZZ +
     &                RIJX * RIJZ * RIJX * RIJX * VIXZXX +
     &                RIJX * RIJZ * RIJX * RIJY * VIXZXY +
     &                RIJX * RIJZ * RIJX * RIJZ * VIXZXZ +
     &                RIJX * RIJZ * RIJY * RIJX * VIXZYX +
     &                RIJX * RIJZ * RIJY * RIJY * VIXZYY +
     &                RIJX * RIJZ * RIJY * RIJZ * VIXZYZ +
     &                RIJX * RIJZ * RIJZ * RIJX * VIXZZX +
     &                RIJX * RIJZ * RIJZ * RIJY * VIXZZY +
     &                RIJX * RIJZ * RIJZ * RIJZ * VIXZZZ +
     &                RIJY * RIJX * RIJX * RIJX * VIYXXX +
     &                RIJY * RIJX * RIJX * RIJY * VIYXXY +
     &                RIJY * RIJX * RIJX * RIJZ * VIYXXZ +
     &                RIJY * RIJX * RIJY * RIJX * VIYXYX +
     &                RIJY * RIJX * RIJY * RIJY * VIYXYY +
     &                RIJY * RIJX * RIJY * RIJZ * VIYXYZ +
     &                RIJY * RIJX * RIJZ * RIJX * VIYXZX +
     &                RIJY * RIJX * RIJZ * RIJY * VIYXZY +
     &                RIJY * RIJX * RIJZ * RIJZ * VIYXZZ +
     &                RIJY * RIJY * RIJX * RIJX * VIYYXX +
     &                RIJY * RIJY * RIJX * RIJY * VIYYXY +
     &                RIJY * RIJY * RIJX * RIJZ * VIYYXZ +
     &                RIJY * RIJY * RIJY * RIJX * VIYYYX +
     &                RIJY * RIJY * RIJY * RIJY * VIYYYY +
     &                RIJY * RIJY * RIJY * RIJZ * VIYYYZ +
     &                RIJY * RIJY * RIJZ * RIJX * VIYYZX +
     &                RIJY * RIJY * RIJZ * RIJY * VIYYZY +
     &                RIJY * RIJY * RIJZ * RIJZ * VIYYZZ +
     &                RIJY * RIJZ * RIJX * RIJX * VIYZXX +
     &                RIJY * RIJZ * RIJX * RIJY * VIYZXY +
     &                RIJY * RIJZ * RIJX * RIJZ * VIYZXZ +
     &                RIJY * RIJZ * RIJY * RIJX * VIYZYX +
     &                RIJY * RIJZ * RIJY * RIJY * VIYZYY +
     &                RIJY * RIJZ * RIJY * RIJZ * VIYZYZ +
     &                RIJY * RIJZ * RIJZ * RIJX * VIYZZX +
     &                RIJY * RIJZ * RIJZ * RIJY * VIYZZY +
     &                RIJY * RIJZ * RIJZ * RIJZ * VIYZZZ +
     &                RIJZ * RIJX * RIJX * RIJX * VIZXXX +
     &                RIJZ * RIJX * RIJX * RIJY * VIZXXY +
     &                RIJZ * RIJX * RIJX * RIJZ * VIZXXZ +
     &                RIJZ * RIJX * RIJY * RIJX * VIZXYX +
     &                RIJZ * RIJX * RIJY * RIJY * VIZXYY +
     &                RIJZ * RIJX * RIJY * RIJZ * VIZXYZ +
     &                RIJZ * RIJX * RIJZ * RIJX * VIZXZX +
     &                RIJZ * RIJX * RIJZ * RIJY * VIZXZY +
     &                RIJZ * RIJX * RIJZ * RIJZ * VIZXZZ +
     &                RIJZ * RIJY * RIJX * RIJX * VIZYXX +
     &                RIJZ * RIJY * RIJX * RIJY * VIZYXY +
     &                RIJZ * RIJY * RIJX * RIJZ * VIZYXZ +
     &                RIJZ * RIJY * RIJY * RIJX * VIZYYX +
     &                RIJZ * RIJY * RIJY * RIJY * VIZYYY +
     &                RIJZ * RIJY * RIJY * RIJZ * VIZYYZ +
     &                RIJZ * RIJY * RIJZ * RIJX * VIZYZX +
     &                RIJZ * RIJY * RIJZ * RIJY * VIZYZY +
     &                RIJZ * RIJY * RIJZ * RIJZ * VIZYZZ +
     &                RIJZ * RIJZ * RIJX * RIJX * VIZZXX +
     &                RIJZ * RIJZ * RIJX * RIJY * VIZZXY +
     &                RIJZ * RIJZ * RIJX * RIJZ * VIZZXZ +
     &                RIJZ * RIJZ * RIJY * RIJX * VIZZYX +
     &                RIJZ * RIJZ * RIJY * RIJY * VIZZYY +
     &                RIJZ * RIJZ * RIJY * RIJZ * VIZZYZ +
     &                RIJZ * RIJZ * RIJZ * RIJX * VIZZZX +
     &                RIJZ * RIJZ * RIJZ * RIJY * VIZZZY +
     &                RIJZ * RIJZ * RIJZ * RIJZ * VIZZZZ )
     &              * RMR * V90 / RIJ12
     &            + ((RIJX * RIJX + RIJX * RIJX ) * VIXXXX +
     &                RIJX * RIJY * VIXXXY +
     &                RIJX * RIJZ * VIXXXZ +
     &                RIJY * RIJX * VIXXYX +
     &               (RIJX * RIJX + RIJY * RIJY ) * VIXXYY +
     &                RIJY * RIJZ * VIXXYZ +
     &                RIJZ * RIJX * VIXXZX +
     &                RIJZ * RIJY * VIXXZY +
     &               (RIJX * RIJX + RIJZ * RIJZ ) * VIXXZZ +
     &                RIJX * RIJY * VIXYXX +
     &                RIJX * RIJY * VIXYYY +
     &                RIJX * RIJY * VIXYZZ +
     &                RIJX * RIJZ * VIXZXX +
     &                RIJX * RIJZ * VIXZYY +
     &                RIJX * RIJZ * VIXZZZ +
     &                RIJY * RIJX * VIYXXX +
     &                RIJY * RIJX * VIYXYY +
     &                RIJY * RIJX * VIYXZZ +
     &               (RIJY * RIJY + RIJX * RIJX ) * VIYYXX +
     &                RIJX * RIJY * VIYYXY +
     &                RIJX * RIJZ * VIYYXZ +
     &                RIJY * RIJX * VIYYYX +
     &               (RIJY * RIJY + RIJY * RIJY ) * VIYYYY +
     &                RIJY * RIJZ * VIYYYZ +
     &                RIJZ * RIJX * VIYYZX +
     &                RIJZ * RIJY * VIYYZY +
     &               (RIJY * RIJY + RIJZ * RIJZ ) * VIYYZZ +
     &                RIJY * RIJZ * VIYZXX +
     &                RIJY * RIJZ * VIYZYY +
     &                RIJY * RIJZ * VIYZZZ +
     &                RIJZ * RIJX * VIZXXX +
     &                RIJZ * RIJX * VIZXYY +
     &                RIJZ * RIJX * VIZXZZ +
     &                RIJZ * RIJY * VIZYXX +
     &                RIJZ * RIJY * VIZYYY +
     &                RIJZ * RIJY * VIZYZZ +
     &               (RIJZ * RIJZ + RIJX * RIJX ) * VIZZXX +
     &                RIJX * RIJY * VIZZXY +
     &                RIJX * RIJZ * VIZZXZ +
     &                RIJY * RIJX * VIZZYX +
     &               (RIJZ * RIJZ + RIJY * RIJY ) * VIZZYY +
     &                RIJY * RIJZ * VIZZYZ +
     &                RIJZ * RIJX * VIZZZX +
     &                RIJZ * RIJY * VIZZZY +
     &               (RIJZ * RIJZ + RIJZ * RIJZ ) * VIZZZZ )
     &              * RMR * V24 / RIJ10
     &             - ((RX1 * RIJX + RIJX * RX1 +
     &                 RX1 * RIJX + RIJX * RX1 ) * VIXXXX +
     &                (RX1 * RIJY + RIJX * RY1 ) * VIXXXY +
     &                (RX1 * RIJZ + RIJX * RZ1 ) * VIXXXZ +
     &                (RY1 * RIJX + RIJY * RX1 ) * VIXXYX +
     &                (RX1 * RIJX + RIJX * RX1 +
     &                 RY1 * RIJY + RIJY * RY1 ) * VIXXYY +
     &                (RY1 * RIJZ + RIJY * RZ1 ) * VIXXYZ +
     &                (RZ1 * RIJX + RIJZ * RX1 ) * VIXXZX +
     &                (RZ1 * RIJY + RIJZ * RY1 ) * VIXXZY +
     &                (RX1 * RIJX + RIJX * RX1 +
     &                 RZ1 * RIJZ + RIJZ * RZ1 ) * VIXXZZ +
     &                (RX1 * RIJY + RIJX * RY1 ) * VIXYXX +
     &                (RX1 * RIJY + RIJX * RY1 ) * VIXYYY +
     &                (RX1 * RIJY + RIJX * RY1 ) * VIXYZZ +
     &                (RX1 * RIJZ + RIJX * RZ1 ) * VIXZXX +
     &                (RX1 * RIJZ + RIJX * RZ1 ) * VIXZYY +
     &                (RX1 * RIJZ + RIJX * RZ1 ) * VIXZZZ +
     &                (RY1 * RIJX + RIJY * RX1 ) * VIYXXX +
     &                (RY1 * RIJX + RIJY * RX1 ) * VIYXYY +
     &                (RY1 * RIJX + RIJY * RX1 ) * VIYXZZ +
     &                (RY1 * RIJY + RIJY * RY1 +
     &                 RX1 * RIJX + RIJX * RX1 ) * VIYYXX +
     &                (RX1 * RIJY + RIJX * RY1 ) * VIYYXY +
     &                (RX1 * RIJZ + RIJX * RZ1 ) * VIYYXZ +
     &                (RY1 * RIJX + RIJY * RX1 ) * VIYYYX +
     &                (RY1 * RIJY + RIJY * RY1 +
     &                 RY1 * RIJY + RIJY * RY1 ) * VIYYYY +
     &                (RY1 * RIJZ + RIJY * RZ1 ) * VIYYYZ +
     &                (RZ1 * RIJX + RIJZ * RX1 ) * VIYYZX +
     &                (RZ1 * RIJY + RIJZ * RY1 ) * VIYYZY +
     &                (RY1 * RIJY + RIJY * RY1 +
     &                 RZ1 * RIJZ + RIJZ * RZ1 ) * VIYYZZ +
     &                (RY1 * RIJZ + RIJY * RZ1 ) * VIYZXX +
     &                (RY1 * RIJZ + RIJY * RZ1 ) * VIYZYY +
     &                (RY1 * RIJZ + RIJY * RZ1 ) * VIYZZZ +
     &                (RZ1 * RIJX + RIJZ * RX1 ) * VIZXXX +
     &                (RZ1 * RIJX + RIJZ * RX1 ) * VIZXYY +
     &                (RZ1 * RIJX + RIJZ * RX1 ) * VIZXZZ +
     &                (RZ1 * RIJY + RIJZ * RY1 ) * VIZYXX +
     &                (RZ1 * RIJY + RIJZ * RY1 ) * VIZYYY +
     &                (RZ1 * RIJY + RIJZ * RY1 ) * VIZYZZ +
     &                (RZ1 * RIJZ + RIJZ * RZ1 +
     &                 RX1 * RIJX + RIJX * RX1 ) * VIZZXX +
     &                (RX1 * RIJY + RIJX * RY1 ) * VIZZXY +
     &                (RX1 * RIJZ + RIJX * RZ1 ) * VIZZXZ +
     &                (RY1 * RIJX + RIJY * RX1 ) * VIZZYX +
     &                (RZ1 * RIJZ + RIJZ * RZ1 +
     &                 RY1 * RIJY + RIJY * RY1 ) * VIZZYY +
     &                (RY1 * RIJZ + RIJY * RZ1 ) * VIZZYZ +
     &                (RZ1 * RIJX + RIJZ * RX1 ) * VIZZZX +
     &                (RZ1 * RIJY + RIJZ * RY1 ) * VIZZZY +
     &                (RZ1 * RIJZ + RIJZ * RZ1 +
     &                 RZ1 * RIJZ + RIJZ * RZ1 ) * VIZZZZ )
     &                * THREE / RIJ8
     &              - (VIXXXX + VIXXYY + VIXXZZ + VIYYXX +
     &                 VIYYYY + VIYYZZ + VIZZXX + VIZZYY + 
     &                 VIZZZZ ) * RMR * SIX / RIJ8
C
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

