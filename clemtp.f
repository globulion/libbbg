C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE DMTCOR(RDMA,NDMA,CHG,DIP,QAD,OCT,
     *                  CHGM,DIPM,QADM,OCTM,
     *                  REDMSS,FREQ,GIJJ,LVEC,
     *                  SHIFT,RF2,RF3,RF4,RK2,RK3,RK4,B,C,D,
     *                  NMOLS,NDMAS,NDMAC,NMODES,MODE,LWRITE)
C
C -----------------------------------------------------------------------------
C
C           ELECTROSTATIC FREQUENCY SHIFT FROM DISTRIBUTED MULTIPOLE
C                   EXPANSION MODEL - CENTRAL MOLECULE MODEL
C                            CORRECTION TERMS
C
C              Bartosz Błasiak                         31 Dec 2013
C
C -----------------------------------------------------------------------------
C
C   Notes:
C     The reduced format of tensor storage is used:
C
C     CHG(i) .
C            1
C     DIP(i) X   Y   Z
C            1   2   3
C     QAD(i) XX  YY  ZZ  XY  XZ  YZ
C            1   2   3   4   5   6
C     OCT(i) XXX YYY ZZZ XXY XXZ XYY YYZ XZZ YZZ XYZ
C            1   2   3   4   5   6   7   8   9
C -----------------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RDMA(NDMAS*3),CHG(NDMAS),DIP(NDMAS*3),QAD(NDMAS*6),
     &          OCT(NDMAS*10),NDMA(NMOLS),
     &          CHGM(NMODES*NDMAC),DIPM(NMODES*NDMAC*3),
     &          QADM(NMODES*NDMAC*6),OCTM(NMODES*NDMAC*10),
     &          REDMSS(NMODES),FREQ(NMODES),GIJJ(NMODES),GIVEC(NMODES)
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,THREE=3.0D+00,
     &           FOUR=4.0D+00,FIVE=5.0D+00,SIX=6.0D+00,
     &           TOCMRC=219474.63067873946D+00)
      DOUBLE PRECISION LX,LY,LZ,LVEC(NMODES*NDMAC*3)
      LOGICAL LWRITE
Cf2py INTENT(OUT) SHIFT,RF2,RF3,RF4,RK2,RK3,RK4,B,C,D
C
      RF2 = ZERO
      RF3 = ZERO
      RF4 = ZERO
C
      RK2 = ZERO
      RK3 = ZERO
      RK4 = ZERO
C
      MODEM1 = MODE - 1
C     
C     EXTRACT THE FIRST MOLECULE
C
      NI  = NDMA(1)
      NI3 = NI*3
      NI6 = NI*6
      NI10= NI*10
C
      DO MM=1,NMODES
         FMM = FREQ(MM)
         GIVEC(MM) = GIJJ(MM) / (REDMSS(MM) * FMM * FMM )
      ENDDO
C
      TMW = - TWO * REDMSS(MODE) * FREQ(MODE)
C
C     --- LOOP OVER SURROUNDING MOLECULES ---
C
      NMOLJ = NI
      DO 190 J=2,NMOLS
         NJ = NDMA(J)
         NMOLJ = NMOLJ + NJ
         DO N=1,NJ
            NJX0 =   (NMOLJ-NJ) +    N       
            NJX3 = 3*(NMOLJ-NJ) + 3*(N-1) + 1
            NJX6 = 6*(NMOLJ-NJ) + 6*(N-1) + 1
            NJX10=10*(NMOLJ-NJ) +10*(N-1) + 1
C                                                     
            NJY3 = NJX3 + 1
            NJZ3 = NJY3 + 1
C                                                     
C           UNPACK FOR MOLECULE J
C                                                     
            CJ = CHG(NJX0)
C
            DJX= DIP(NJX3)
            DJY= DIP(NJY3)
            DJZ= DIP(NJZ3)
C
            QJXX = QAD(NJX6  )
            QJYY = QAD(NJX6+1)
            QJZZ = QAD(NJX6+2)
            QJXY = QAD(NJX6+3)
            QJXZ = QAD(NJX6+4)
            QJYZ = QAD(NJX6+5)
C
            OJXXX = OCT(NJX10  )
            OJYYY = OCT(NJX10+1)
            OJZZZ = OCT(NJX10+2)
            OJXXY = OCT(NJX10+3)
            OJXXZ = OCT(NJX10+4)
            OJXYY = OCT(NJX10+5)
            OJYYZ = OCT(NJX10+6)
            OJXZZ = OCT(NJX10+7)
            OJYZZ = OCT(NJX10+8)
            OJXYZ = OCT(NJX10+9)
C
C           --- ITERATE OVER CENTRAL MOLECULE I
C
            DO M=1,NI
               NIX0 =    M
               NIX3 = 3*(M-1) + 1
               NIX6 = 6*(M-1) + 1
               NIX10=10*(M-1) + 1
C
               NIY3 = NIX3 + 1
               NIZ3 = NIY3 + 1
C
C              UNPACK DISTANCES FOR MOLECULE I
C
               RIX = RDMA(NIX3)
               RIY = RDMA(NIY3)
               RIZ = RDMA(NIZ3)
C
               RX = RIX - RDMA(NJX3)
               RY = RIY - RDMA(NJY3)
               RZ = RIZ - RDMA(NJZ3)
C
               RMN  = DSQRT(RX*RX+RY*RY+RZ*RZ)
               RMN3 = ONE/(RMN*RMN*RMN)
               RMN5 = RMN3/(RMN*RMN)
               RMN7 = RMN5/(RMN*RMN)
C
C              EXTRACT DMA OF CENTRAL MOLECULE
C
               CI  = CHG(NIX0)
C
               DIX = DIP(NIX3)
               DIY = DIP(NIY3)
               DIZ = DIP(NIZ3)
C
               QIXX = QAD(NIX6  )
               QIYY = QAD(NIX6+1)
               QIZZ = QAD(NIX6+2)
               QIXY = QAD(NIX6+3)
               QIXZ = QAD(NIX6+4)
               QIYZ = QAD(NIX6+5)
C
               OIXXX = OCT(NIX10  )
               OIYYY = OCT(NIX10+1)
               OIZZZ = OCT(NIX10+2)
               OIXXY = OCT(NIX10+3)
               OIXXZ = OCT(NIX10+4)
               OIXYY = OCT(NIX10+5)
               OIYYZ = OCT(NIX10+6)
               OIXZZ = OCT(NIX10+7)
               OIYZZ = OCT(NIX10+8)
               OIXYZ = OCT(NIX10+9)
C
C              EXTRACT FIRST DERIVATIVES OF DMA WRT MODE
C
               NIXM0 = NI  *MODEM1 +    M
               NIXM3 = NI3 *MODEM1 + 3*(M-1) + 1
               NIXM6 = NI6 *MODEM1 + 6*(M-1) + 1
               NIXM10= NI10*MODEM1 +10*(M-1) + 1
C
               NIYM3 = NIXM3 + 1
               NIZM3 = NIYM3 + 1
C
               CI1  = CHGM(NIXM0)
C
               DIX1 = DIPM(NIXM3)
               DIY1 = DIPM(NIYM3)
               DIZ1 = DIPM(NIZM3)
C
               QIXX1 = QADM(NIXM6  )
               QIYY1 = QADM(NIXM6+1)
               QIZZ1 = QADM(NIXM6+2)
               QIXY1 = QADM(NIXM6+3)
               QIXZ1 = QADM(NIXM6+4)
               QIYZ1 = QADM(NIXM6+5)
C
               OIXXX1 = OCTM(NIXM10  )
               OIYYY1 = OCTM(NIXM10+1)
               OIZZZ1 = OCTM(NIXM10+2)
               OIXXY1 = OCTM(NIXM10+3)
               OIXXZ1 = OCTM(NIXM10+4)
               OIXYY1 = OCTM(NIXM10+5)
               OIYYZ1 = OCTM(NIXM10+6)
               OIXZZ1 = OCTM(NIXM10+7)
               OIYZZ1 = OCTM(NIXM10+8)
               OIXYZ1 = OCTM(NIXM10+9)
C
C              ELECTRONIC ANAHARMONICITY CORRECTION TERMS
C
               LX = LVEC(NIXM3)
               LY = LVEC(NIYM3)
               LZ = LVEC(NIZM3)
C
               S1 =   RX*LX +   RY*LY +   RZ*LZ
               S2 =  DJX*LX +  DJY*LY +  DJZ*LZ
               S3 =  DJX*RX +  DJY*RY +  DJZ*RZ
               S4 =   LX*LX +   LY*LY +   LZ*LZ
               S5 =  DIX*LX +  DIY*LY +  DIZ*LZ
               S6 =  DIX*RX +  DIY*RY +  DIZ*RZ
               S7 =  DIX*DJX+  DIY*DJY+  DIZ*DJZ
               S8 = DIX1*LX + DIY1*LY + DIZ1*LZ
               S9 = DIX1*RX + DIY1*RY + DIZ1*RZ
               S10= DIX1*DJX+ DIY1*DJY+ DIZ1*DJZ
C
               S11= (QJXX * RX * RX       +  
     &               QJXY * RX * RY * TWO +  
     &               QJXZ * RX * RZ * TWO +  
     &               QJYY * RY * RY       +  
     &               QJYZ * RY * RZ * TWO +  
     &               QJZZ * RZ * RZ)
C
               S14= (QIXX * RX * RX       +  
     &               QIXY * RX * RY * TWO +  
     &               QIXZ * RX * RZ * TWO +  
     &               QIYY * RY * RY       +  
     &               QIYZ * RY * RZ * TWO +  
     &               QIZZ * RZ * RZ)
C
               S12=  QJXX * LX * RX +  
     &               QJXY * LX * RY +
     &               QJXZ * LX * RZ +
     &               QJXY * LY * RX +
     &               QJYY * LY * RY +
     &               QJYZ * LY * RZ +
     &               QJXZ * LZ * RX +
     &               QJYZ * LZ * RY +
     &               QJZZ * LZ * RZ 
C
               S13=  QIXX * LX * RX +  
     &               QIXY * LX * RY +
     &               QIXZ * LX * RZ +
     &               QIXY * LY * RX +
     &               QIYY * LY * RY +
     &               QIYZ * LY * RZ +
     &               QIXZ * LZ * RX +
     &               QIYZ * LZ * RY +
     &               QIZZ * LZ * RZ 
C
               S15= (QIXX1 * RX * RX       +  
     &               QIXY1 * RX * RY * TWO +  
     &               QIXZ1 * RX * RZ * TWO +  
     &               QIYY1 * RY * RY       +  
     &               QIYZ1 * RY * RZ * TWO +  
     &               QIZZ1 * RZ * RZ)
C
               S16=  QIXX1 * LX * RX +  
     &               QIXY1 * LX * RY +
     &               QIXZ1 * LX * RZ +
     &               QIXY1 * LY * RX +
     &               QIYY1 * LY * RY +
     &               QIYZ1 * LY * RZ +
     &               QIXZ1 * LZ * RX +
     &               QIYZ1 * LZ * RY +
     &               QIZZ1 * LZ * RZ 
C
C              ACCUMULATE THE TERMS
C
               RK2 = RK2   - 2.00D+00 * CJ * CI1 * S1 * RMN3
C
               RK3 = RK3 + ( 2.00D+00 * CI1 * S2 - CI * CJ * S4 
     &                     - 2.00D+00 * CJ * S8 ) * RMN3
C
               RK3 = RK3 + ( 3.00D+00 * CI * CJ * S1 
     &                     - 6.00D+00 * CI1 * S1 * S3
     &                     + 6.00D+00 * CJ * S1 * S9 ) * RMN5
C
               RK4 = RK4 + ( 4.00D+00 * CI1 * S12
     &                     + 4.00D+00 * CJ * S16 
     &                     - 6.00D+00 * ( S9*S2 + S1*S10 + S8*S3 )
     &                     - 3.00D+00 * CI * ( 2.00D+00 * S1*S2 + S4*S3)
     &                     + 3.00D+00 * CJ * ( 2.00D+00 * S1*S5 + S4*S6)
     &                     ) * RMN5
C
               RK4 = RK4 + (30.00D+00 * S9 * S1 * S3 
     &                     -10.00D+00 * CI1 * S11 * S1
     &                     -10.00D+00 * CJ * S15 * S1
     &                     +15.00D+00 * CI * S1 * S1 * S3
     &                     -15.00D+00 * CJ * S1 * S1 * S6 ) * RMN7
C
C              --- ITERATE OVER NORMAL COORDINATES OF MOLECULE I ---
C
               DO MM=1,NMODES
                  NIXM3 = NI3 *(MM-1) + 3*(M-1) + 1
                  NIYM3 = NIXM3 + 1
                  NIZM3 = NIYM3 + 1
C
                  GIVECM = GIVEC(MM)
C
C                 UNPACK EIGENVECTOR ELEMENTS
C
                  LX = LVEC(NIXM3)
                  LY = LVEC(NIYM3)
                  LZ = LVEC(NIZM3)
C                                                     
C                 MECHANICAL ANHARMONICITY CORRECTION TERMS
C                                                           
                  S1 =   RX*LX +   RY*LY +   RZ*LZ
                  S2 =  DJX*LX +  DJY*LY +  DJZ*LZ
                  S3 =  DJX*RX +  DJY*RY +  DJZ*RZ
                  S4 =   LX*LX +   LY*LY +   LZ*LZ
                  S5 =  DIX*LX +  DIY*LY +  DIZ*LZ
                  S6 =  DIX*RX +  DIY*RY +  DIZ*RZ
                  S7 =  DIX*DJX+  DIY*DJY+  DIZ*DJZ
C
                  S11= (QJXX * RX * RX       +  
     &                  QJXY * RX * RY * TWO +  
     &                  QJXZ * RX * RZ * TWO +  
     &                  QJYY * RY * RY       +  
     &                  QJYZ * RY * RZ * TWO +  
     &                  QJZZ * RZ * RZ)
C                                                           
                  S14= (QIXX * RX * RX       +  
     &                  QIXY * RX * RY * TWO +  
     &                  QIXZ * RX * RZ * TWO +  
     &                  QIYY * RY * RY       +  
     &                  QIYZ * RY * RZ * TWO +  
     &                  QIZZ * RZ * RZ)
C                                                           
                  S12=  QJXX * LX * RX +  
     &                  QJXY * LX * RY +
     &                  QJXZ * LX * RZ +
     &                  QJXY * LY * RX +
     &                  QJYY * LY * RY +
     &                  QJYZ * LY * RZ +
     &                  QJXZ * LZ * RX +
     &                  QJYZ * LZ * RY +
     &                  QJZZ * LZ * RZ 
C                                                           
                  S13=  QIXX * LX * RX +  
     &                  QIXY * LX * RY +
     &                  QIXZ * LX * RZ +
     &                  QIXY * LY * RX +
     &                  QIYY * LY * RY +
     &                  QIYZ * LY * RZ +
     &                  QIXZ * LZ * RX +
     &                  QIYZ * LZ * RY +
     &                  QIZZ * LZ * RZ 
C                                                     
C                 ACCUMULATE THE TERMS
C                                                           
                  RF2 = RF2 - CI * CJ * S1 * RMN3 * GIVECM
C
                  RF3 = RF3 + ( CI * S2 - CJ * S5 ) * RMN3 * GIVECM
C
                  RF3 = RF3 + ( 3.00D+00 * ( CJ*S1*S6 - CI*S1*S3 )
     &                        ) * RMN5*GIVECM
C
                  RF4 = RF4 + ( 2.00D+00 * ( CJ * S13 + CI * S12 )
     &                        - 3.00D+00 * ( S6*S2 + S1*S7 + S3*S5 )
     &                        ) * RMN5 * GIVECM
C
                  RF4 = RF4 + (15.00D+00 * S6 * S3 * S1
     &                        - 5.00D+00 * ( CI*S11*S1 + CJ*S14*S1 )
     &                        ) * RMN7 * GIVECM
            ENDDO
         ENDDO
      ENDDO
C
 190  CONTINUE
C
      RF2 =  RF2/TMW
      RF3 =  RF3/TMW
      RF4 =  RF4/TMW
C
      RK2 = -RK2/TMW
      RK3 = -RK3/TMW
      RK4 = -RK4/TMW
C
      A = ZERO
      B = A + RF2 + RK2
      C = B + RF3 + RK3
      D = C + RF4 + RK4
      E = ZERO
C
      SHIFT = D
C
      IF (LWRITE) THEN
C
          CRF2 = RF2 * TOCMRC
          CRF3 = RF3 * TOCMRC
          CRF4 = RF4 * TOCMRC
          CRK2 = RK2 * TOCMRC
          CRK3 = RK3 * TOCMRC
          CRK4 = RK4 * TOCMRC
C
          CB = B * TOCMRC
          CC = C * TOCMRC
          CD = D * TOCMRC
C
          WRITE(*,9919) CRF2, CRK2, CRF3, CRK3, CRF4, CRK4
          WRITE(*,9917) CB, CC, CD
      ENDIF
C
 9919 FORMAT(/,14X,'MA',12X,'EA',/,'R2',D14.2,D14.2,/
     &                             'R3',D14.2,D14.2,/
     &                             'R4',D14.2,D14.2,/)
C
 9917 FORMAT(/,' CORRECTIONS TOTAL',/,'B',8X,D10.2,/
     &                                'C',8X,D10.2,/
     &                                'D',8X,D10.2,/)
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE SDMTPE(RDMA,NDMA,CHG,DIP,QAD,OCT,
     *                  CHGM,DIPM,QADM,OCTM,
     *                  REDMSS,FREQ,
     *                  SHIFT,A,B,C,D,E,
     *                  NMOLS,NDMAS,NDMAC,NMODES,MODE,LWRITE)
C
C -----------------------------------------------------------------------------
C
C           ELECTROSTATIC FREQUENCY SHIFT FROM DISTRIBUTED MULTIPOLE
C                   EXPANSION MODEL - CENTRAL MOLECULE MODEL
C                           ELECTRONIC ANHARMONICITY
C                           WITHOUT CORRECTION TERMS
C
C              Bartosz Błasiak                         23 Nov 2013
C
C -----------------------------------------------------------------------------
C
C   Notes:
C     The reduced format of tensor storage is used:
C
C     CHG(i) .
C            1
C     DIP(i) X   Y   Z
C            1   2   3
C     QAD(i) XX  YY  ZZ  XY  XZ  YZ
C            1   2   3   4   5   6
C     OCT(i) XXX YYY ZZZ XXY XXZ XYY YYZ XZZ YZZ XYZ
C            1   2   3   4   5   6   7   8   9
C -----------------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RDMA(NDMAS*3),CHG(NDMAS),DIP(NDMAS*3),QAD(NDMAS*6),
     &          OCT(NDMAS*10),NDMA(NMOLS),
     &          CHGM(NDMAC),DIPM(NDMAC*3),
     &          QADM(NDMAC*6),OCTM(NDMAC*10),
     &          REDMSS(NMODES),FREQ(NMODES)
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,THREE=3.0D+00,
     &           FOUR=4.0D+00,FIVE=5.0D+00,SIX=6.0D+00)
      LOGICAL LWRITE
Cf2py INTENT(OUT) SHIFT,A,B,C,D,E
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
C     EXTRACT THE FIRST MOLECULE
C
      NI = NDMA(1)
      DO 90 M=1,NI
         NIX0 =    M
         NIX3 = 3*(M-1) + 1
         NIX6 = 6*(M-1) + 1
         NIX10=10*(M-1) + 1
C
         NIY3 = NIX3 + 1
         NIZ3 = NIY3 + 1
C
C        UNPACK FOR MOLECULE I
C
         RIX = RDMA(NIX3)
         RIY = RDMA(NIY3)
         RIZ = RDMA(NIZ3)
         CI = CHGM(NIX0)
         DIX= DIPM(NIX3)
         DIY= DIPM(NIY3)
         DIZ= DIPM(NIZ3)
         QIXX = QADM(NIX6  )
         QIYY = QADM(NIX6+1)
         QIZZ = QADM(NIX6+2)
         QIXY = QADM(NIX6+3)
         QIXZ = QADM(NIX6+4)
         QIYZ = QADM(NIX6+5)
         OIXXX = OCTM(NIX10  )
         OIYYY = OCTM(NIX10+1)
         OIZZZ = OCTM(NIX10+2)
         OIXXY = OCTM(NIX10+3)
         OIXXZ = OCTM(NIX10+4)
         OIXYY = OCTM(NIX10+5)
         OIYYZ = OCTM(NIX10+6)
         OIXZZ = OCTM(NIX10+7)
         OIYZZ = OCTM(NIX10+8)
         OIXYZ = OCTM(NIX10+9)
C
C        --- LOOP OVER SURROUNDING MOLECULE ---
C
         NMOLJ = NI
         DO J=2,NMOLS
            NJ = NDMA(J)
            NMOLJ = NMOLJ + NJ
            DO N=1,NJ
               NJX0 =   (NMOLJ-NJ) +    N       
               NJX3 = 3*(NMOLJ-NJ) + 3*(N-1) + 1
               NJX6 = 6*(NMOLJ-NJ) + 6*(N-1) + 1
               NJX10=10*(NMOLJ-NJ) +10*(N-1) + 1
C                                                        
               NJY3 = NJX3 + 1
               NJZ3 = NJY3 + 1
C                                                        
C              UNPACK FOR MOLECULE J
C                                                        
               RX = RDMA(NJX3) - RIX
               RY = RDMA(NJY3) - RIY
               RZ = RDMA(NJZ3) - RIZ
C                                                        
               CJ = CHG(NJX0)
               DJX= DIP(NJX3)
               DJY= DIP(NJY3)
               DJZ= DIP(NJZ3)
               QJXX = QAD(NJX6  )
               QJYY = QAD(NJX6+1)
               QJZZ = QAD(NJX6+2)
               QJXY = QAD(NJX6+3)
               QJXZ = QAD(NJX6+4)
               QJYZ = QAD(NJX6+5)
               OJXXX = OCT(NJX10  )
               OJYYY = OCT(NJX10+1)
               OJZZZ = OCT(NJX10+2)
               OJXXY = OCT(NJX10+3)
               OJXXZ = OCT(NJX10+4)
               OJXYY = OCT(NJX10+5)
               OJYYZ = OCT(NJX10+6)
               OJXZZ = OCT(NJX10+7)
               OJYZZ = OCT(NJX10+8)
               OJXYZ = OCT(NJX10+9)
C                                                        
               RMN= DSQRT(RX*RX+RY*RY+RZ*RZ)
               RMN3 = ONE/(RMN*RMN*RMN)
               RMN5 = RMN3/(RMN*RMN)
               RMN7 = RMN5/(RMN*RMN)
C                                                        
C              TENSORDOTS
C                                                        
               S1 = -(DJX*RX+DJY*RY+DJZ*RZ)
               S2 =  (DIX*RX+DIY*RY+DIZ*RZ)
C                                                        
               S3 =  (DIX*DJX+DIY*DJY+DIZ*DJZ)
               S4 =   S1 * S2
C                                                        
               S5 = (QJXX * RX * RX       +  
     &               QJXY * RX * RY * TWO +  
     &               QJXZ * RX * RZ * TWO +  
     &               QJYY * RY * RY       +  
     &               QJYZ * RY * RZ * TWO +  
     &               QJZZ * RZ * RZ)
C                                                        
               S6 = (QIXX * RX * RX       +  
     &               QIXY * RX * RY * TWO +  
     &               QIXZ * RX * RZ * TWO +  
     &               QIYY * RY * RY       +  
     &               QIYZ * RY * RZ * TWO +  
     &               QIZZ * RZ * RZ)
C                                                        
               S7 =-(QJXX * DIX * RX +  
     &               QJXY * DIX * RY +
     &               QJXZ * DIX * RZ +
     &               QJXY * DIY * RX +
     &               QJYY * DIY * RY +
     &               QJYZ * DIY * RZ +
     &               QJXZ * DIZ * RX +
     &               QJYZ * DIZ * RY +
     &               QJZZ * DIZ * RZ)
C                                                        
               S8 = (QIXX * DJX * RX +  
     &               QIXY * DJX * RY +
     &               QIXZ * DJX * RZ +
     &               QIXY * DJY * RX +
     &               QIYY * DJY * RY +
     &               QIYZ * DJY * RZ +
     &               QIXZ * DJZ * RX +
     &               QIYZ * DJZ * RY +
     &               QIZZ * DJZ * RZ)
C        
               S9 =  S5 * S2
               S10=  S6 * S1
C                                                        
               RXRXRX = RX * RX * RX
               RXRXRY = RX * RX * RY * THREE
               RXRYRY = RX * RY * RY * THREE
               RYRYRY = RY * RY * RY
               RYRYRZ = RY * RY * RZ * THREE
               RYRZRZ = RY * RZ * RZ * THREE
               RZRZRZ = RZ * RZ * RZ
               RXRYRZ = RX * RY * RZ * SIX
               RXRXRZ = RX * RX * RZ * THREE
               RXRZRZ = RX * RZ * RZ * THREE
C                                                        
               S11=-(OJXXX * RXRXRX +
     &               OJXXY * RXRXRY +
     &               OJXYY * RXRYRY +
     &               OJYYY * RYRYRY +
     &               OJYYZ * RYRYRZ +
     &               OJYZZ * RYRZRZ +
     &               OJZZZ * RZRZRZ +
     &               OJXYZ * RXRYRZ +
     &               OJXXZ * RXRXRZ +
     &               OJXZZ * RXRZRZ)
C                                                        
               S12=  OIXXX * RXRXRX +
     &               OIXXY * RXRXRY +
     &               OIXYY * RXRYRY +
     &               OIYYY * RYRYRY +
     &               OIYYZ * RYRYRZ +
     &               OIYZZ * RYRZRZ +
     &               OIZZZ * RZRZRZ +
     &               OIXYZ * RXRYRZ +
     &               OIXXZ * RXRXRZ +
     &               OIXZZ * RXRZRZ 
C                                                        
C              ACCUMULATE THE TERMS
C                                                        
               CC = CC + CI * CJ / RMN 
               CD = CD + S1 * CI * RMN3
               DC = DC + S2 * CJ * RMN3 
               CQ = CQ + S5 * CI * RMN5
               QC = QC + S6 * CJ * RMN5 
               CT = CT + S11 * CI * RMN7
               TC = TC + S12 * CJ * RMN7
               DD = DD + (S3 * RMN3 + THREE * S4 * RMN5)
               DQ = DQ + (TWO * S7 * RMN5 + FIVE * S9  * RMN7 )
               QD = QD + (TWO * S8 * RMN5 + FIVE * S10 * RMN7 )
            ENDDO
         ENDDO
C
 90   CONTINUE
C
      TMF = TWO * REDMSS(MODE) * FREQ(MODE)
C
      CC = CC / TMF
      CD = CD / TMF
      DC = DC / TMF
      CQ = CQ / TMF
      QC = QC / TMF
      CT = CT / TMF
      TC = TC / TMF
      DD = DD / TMF
      DQ = DQ / TMF
      QD = QD / TMF
C
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
      SHIFT = D
C
      IF (LWRITE) THEN
          WRITE(*,91) "CC",CC,"CD",CD,"DC",DC,"CQ",CQ,"QC",QC,
     &                "CO",CT,"OC",TC,"DD",DD,"DQ",DQ,"QD",QD
      ENDIF
 91   FORMAT(10(3X,A," = ",F20.16,/))

      RETURN
      END

C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE SDMTPM(RDMA,NDMA,CHG,DIP,QAD,OCT,
     *                  CHGM,DIPM,QADM,OCTM,
     *                  REDMSS,FREQ,GIJJ,
     *                  SHIFT,A,B,C,D,E,
     *                  NMOLS,NDMAS,NDMAC,NMODES,MODE,LWRITE)
C
C -----------------------------------------------------------------------------
C
C           ELECTROSTATIC FREQUENCY SHIFT FROM DISTRIBUTED MULTIPOLE
C                   EXPANSION MODEL - CENTRAL MOLECULE MODEL
C                           MECHANICAL ANHARMONICITY
C                           WITHOUT CORRECTION TERMS
C
C              Bartosz Błasiak                         22 Dec 2013
C
C -----------------------------------------------------------------------------
C
C   Notes:
C     The reduced format of tensor storage is used:
C
C     CHG(i) .
C            1
C     DIP(i) X   Y   Z
C            1   2   3
C     QAD(i) XX  YY  ZZ  XY  XZ  YZ
C            1   2   3   4   5   6
C     OCT(i) XXX YYY ZZZ XXY XXZ XYY YYZ XZZ YZZ XYZ
C            1   2   3   4   5   6   7   8   9
C -----------------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RDMA(NDMAS*3),CHG(NDMAS),DIP(NDMAS*3),QAD(NDMAS*6),
     &          OCT(NDMAS*10),NDMA(NMOLS),
     &          CHGM(NMODES*NDMAC),DIPM(NMODES*NDMAC*3),
     &          QADM(NMODES*NDMAC*6),OCTM(NMODES*NDMAC*10),
     &          REDMSS(NMODES),FREQ(NMODES),GIJJ(NMODES),GIVEC(NMODES)
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,THREE=3.0D+00,
     &           FOUR=4.0D+00,FIVE=5.0D+00,SIX=6.0D+00)
      LOGICAL LWRITE
Cf2py INTENT(OUT) SHIFT,A,B,C,D,E
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
C     EXTRACT THE FIRST MOLECULE
C
      NI = NDMA(1)
      NI3= NI*3
      NI6= NI*6
      NI10=NI*10
C
      DO MM=1,NMODES
         FMM = FREQ(MM)
         GIVEC(MM) = GIJJ(MM) / (REDMSS(MM) * FMM * FMM )
      ENDDO
C
      TMW = - TWO * REDMSS(MODE) * FREQ(MODE)
C
C     --- LOOP OVER SURROUNDING MOLECULE ---
C
      NMOLJ = NI
      DO 190 J=2,NMOLS
         NJ = NDMA(J)
         NMOLJ = NMOLJ + NJ
         DO N=1,NJ
            NJX0 =   (NMOLJ-NJ) +    N       
            NJX3 = 3*(NMOLJ-NJ) + 3*(N-1) + 1
            NJX6 = 6*(NMOLJ-NJ) + 6*(N-1) + 1
            NJX10=10*(NMOLJ-NJ) +10*(N-1) + 1
C                                                     
            NJY3 = NJX3 + 1
            NJZ3 = NJY3 + 1
C                                                     
C           UNPACK FOR MOLECULE J
C                                                     
            CJ = CHG(NJX0)
            DJX= DIP(NJX3)
            DJY= DIP(NJY3)
            DJZ= DIP(NJZ3)
            QJXX = QAD(NJX6  )
            QJYY = QAD(NJX6+1)
            QJZZ = QAD(NJX6+2)
            QJXY = QAD(NJX6+3)
            QJXZ = QAD(NJX6+4)
            QJYZ = QAD(NJX6+5)
            OJXXX = OCT(NJX10  )
            OJYYY = OCT(NJX10+1)
            OJZZZ = OCT(NJX10+2)
            OJXXY = OCT(NJX10+3)
            OJXXZ = OCT(NJX10+4)
            OJXYY = OCT(NJX10+5)
            OJYYZ = OCT(NJX10+6)
            OJXZZ = OCT(NJX10+7)
            OJYZZ = OCT(NJX10+8)
            OJXYZ = OCT(NJX10+9)
C
C           --- ITERATE OVER CENTRAL MOLECULE I
C
            DO M=1,NI             
               NIX3 = 3*(M-1) + 1
               NIY3 = NIX3 + 1
               NIZ3 = NIY3 + 1
C                                    
C              UNPACK DISTANCES FOR MOLECULE I
C                                    
               RIX = RDMA(NIX3)
               RIY = RDMA(NIY3)
               RIZ = RDMA(NIZ3)
C
               RX = RDMA(NJX3) - RIX
               RY = RDMA(NJY3) - RIY
               RZ = RDMA(NJZ3) - RIZ
C                                                     
               RMN= DSQRT(RX*RX+RY*RY+RZ*RZ)
               RMN3 = ONE/(RMN*RMN*RMN)
               RMN5 = RMN3/(RMN*RMN)
               RMN7 = RMN5/(RMN*RMN)
C
C              --- ITERATE OVER NORMAL COORDINATES OF MOLECULE I ---
C
               DO MM=1,NMODES
                  NIXM0 = NI  *(MM-1) +    M
                  NIXM3 = NI3 *(MM-1) + 3*(M-1) + 1
                  NIXM6 = NI6 *(MM-1) + 6*(M-1) + 1
                  NIXM10= NI10*(MM-1) +10*(M-1) + 1
C
                  NIYM3 = NIXM3 + 1
                  NIZM3 = NIYM3 + 1
C
                  GIVECM = GIVEC(MM)
C
                  CI  = CHGM(NIXM0)       
                  DIX = DIPM(NIXM3)
                  DIY = DIPM(NIYM3)
                  DIZ = DIPM(NIZM3)
                  QIXX = QADM(NIXM6  )
                  QIYY = QADM(NIXM6+1)
                  QIZZ = QADM(NIXM6+2)
                  QIXY = QADM(NIXM6+3)
                  QIXZ = QADM(NIXM6+4)
                  QIYZ = QADM(NIXM6+5)
                  OIXXX = OCTM(NIXM10  )
                  OIYYY = OCTM(NIXM10+1)
                  OIZZZ = OCTM(NIXM10+2)
                  OIXXY = OCTM(NIXM10+3)
                  OIXXZ = OCTM(NIXM10+4)
                  OIXYY = OCTM(NIXM10+5)
                  OIYYZ = OCTM(NIXM10+6)
                  OIXZZ = OCTM(NIXM10+7)
                  OIYZZ = OCTM(NIXM10+8)
                  OIXYZ = OCTM(NIXM10+9)
C                                                     
C                 TENSORDOTS                                 
C                                                           
                  S1 = -(DJX*RX+DJY*RY+DJZ*RZ)
                  S2 =  (DIX*RX+DIY*RY+DIZ*RZ)
C                                                           
                  S3 =  (DIX*DJX+DIY*DJY+DIZ*DJZ)
                  S4 =   S1 * S2
C                                                           
                  S5 = (QJXX * RX * RX       +  
     &                  QJXY * RX * RY * TWO +  
     &                  QJXZ * RX * RZ * TWO +  
     &                  QJYY * RY * RY       +  
     &                  QJYZ * RY * RZ * TWO +  
     &                  QJZZ * RZ * RZ)
C                                                           
                  S6 = (QIXX * RX * RX       +  
     &                  QIXY * RX * RY * TWO +  
     &                  QIXZ * RX * RZ * TWO +  
     &                  QIYY * RY * RY       +  
     &                  QIYZ * RY * RZ * TWO +  
     &                  QIZZ * RZ * RZ)
C                                                           
                  S7 =-(QJXX * DIX * RX +  
     &                  QJXY * DIX * RY +
     &                  QJXZ * DIX * RZ +
     &                  QJXY * DIY * RX +
     &                  QJYY * DIY * RY +
     &                  QJYZ * DIY * RZ +
     &                  QJXZ * DIZ * RX +
     &                  QJYZ * DIZ * RY +
     &                  QJZZ * DIZ * RZ)
C                                                           
                  S8 = (QIXX * DJX * RX +  
     &                  QIXY * DJX * RY +
     &                  QIXZ * DJX * RZ +
     &                  QIXY * DJY * RX +
     &                  QIYY * DJY * RY +
     &                  QIYZ * DJY * RZ +
     &                  QIXZ * DJZ * RX +
     &                  QIYZ * DJZ * RY +
     &                  QIZZ * DJZ * RZ)
C                                                            
                  S9 =  S5 * S2
                  S10=  S6 * S1
C                                                     
                  RXRXRX = RX * RX * RX                            
                  RXRXRY = RX * RX * RY * THREE
                  RXRYRY = RX * RY * RY * THREE
                  RYRYRY = RY * RY * RY
                  RYRYRZ = RY * RY * RZ * THREE
                  RYRZRZ = RY * RZ * RZ * THREE
                  RZRZRZ = RZ * RZ * RZ
                  RXRYRZ = RX * RY * RZ * SIX
                  RXRXRZ = RX * RX * RZ * THREE
                  RXRZRZ = RX * RZ * RZ * THREE
C                                                           
                  S11=-(OJXXX * RXRXRX +
     &                  OJXXY * RXRXRY +
     &                  OJXYY * RXRYRY +
     &                  OJYYY * RYRYRY +
     &                  OJYYZ * RYRYRZ +
     &                  OJYZZ * RYRZRZ +
     &                  OJZZZ * RZRZRZ +
     &                  OJXYZ * RXRYRZ +
     &                  OJXXZ * RXRXRZ +
     &                  OJXZZ * RXRZRZ)
C                                                           
                  S12=  OIXXX * RXRXRX +
     &                  OIXXY * RXRXRY +
     &                  OIXYY * RXRYRY +
     &                  OIYYY * RYRYRY +
     &                  OIYYZ * RYRYRZ +
     &                  OIYZZ * RYRZRZ +
     &                  OIZZZ * RZRZRZ +
     &                  OIXYZ * RXRYRZ +
     &                  OIXXZ * RXRXRZ +
     &                  OIXZZ * RXRZRZ 
C                                                           
C                 ACCUMULATE THE TERMS
C                                                           
                  CC = CC + CI * CJ / RMN * GIVECM
                  CD = CD + S1 * CI * RMN3 * GIVECM
                  DC = DC + S2 * CJ * RMN3 * GIVECM
                  CQ = CQ + S5 * CI * RMN5 * GIVECM
                  QC = QC + S6 * CJ * RMN5 * GIVECM
                  CT = CT + S11 * CI * RMN7 * GIVECM
                  TC = TC + S12 * CJ * RMN7 * GIVECM
                  DD = DD + (S3 * RMN3 + THREE * S4 * RMN5) * GIVECM
                  DQ = DQ + (TWO * S7 * RMN5 + FIVE * S9  * RMN7 ) 
     &                                  * GIVECM
                  QD = QD + (TWO * S8 * RMN5 + FIVE * S10 * RMN7 ) 
     &                                  * GIVECM
            ENDDO
         ENDDO
      ENDDO
C
 190  CONTINUE
C
      CC = CC / TMW
      CD = CD / TMW
      DC = DC / TMW
      CQ = CQ / TMW
      QC = QC / TMW
      CT = CT / TMW
      TC = TC / TMW
      DD = DD / TMW
      DQ = DQ / TMW
      QD = QD / TMW
C
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
      SHIFT = D
C
      IF (LWRITE) THEN
          WRITE(*,91) "CC",CC,"CD",CD,"DC",DC,"CQ",CQ,"QC",QC,
     &                "CO",CT,"OC",TC,"DD",DD,"DQ",DQ,"QD",QD
      ENDIF
 91   FORMAT(10(3X,A," = ",F20.16,/))
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE EDMTPC(RDMA,CHG,DIP,QAD,OCT,NDMA,
     *                  EINT,A,B,C,D,E,NMOLS,NDMAS,LWRITE)
C
C -----------------------------------------------------------------------------
C
C           ELECTROSTATIC INTERACTION ENERGY FROM DISTRIBUTED MULTIPOLE
C                   EXPANSION MODEL - CENTRAL MOLECULE MODEL
C
C              Bartosz Błasiak                         14 Nov 2013
C
C -----------------------------------------------------------------------------
C
C   Notes:
C     The reduced format of tensor storage is used:
C
C     CHG(i) .
C            1
C     DIP(i) X   Y   Z
C            1   2   3
C     QAD(i) XX  YY  ZZ  XY  XZ  YZ
C            1   2   3   4   5   6
C     OCT(i) XXX YYY ZZZ XXY XXZ XYY YYZ XZZ YZZ XYZ
C            1   2   3   4   5   6   7   8   9
C -----------------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RDMA(NDMAS*3),CHG(NDMAS),DIP(NDMAS*3),QAD(NDMAS*6),
     &          OCT(NDMAS*10),NDMA(NMOLS)
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,THREE=3.0D+00,
     &           FOUR=4.0D+00,FIVE=5.0D+00,SIX=6.0D+00)
      LOGICAL LWRITE
Cf2py INTENT(OUT) EINT,A,B,C,D,E
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
C     EXTRACT THE FIRST MOLECULE
C
      NI = NDMA(1)
      DO 90 M=1,NI
         NIX0 =    M
         NIX3 = 3*(M-1) + 1
         NIX6 = 6*(M-1) + 1
         NIX10=10*(M-1) + 1
C
         NIY3 = NIX3 + 1
         NIZ3 = NIY3 + 1
C
C        UNPACK FOR MOLECULE I
C
         RIX = RDMA(NIX3)
         RIY = RDMA(NIY3)
         RIZ = RDMA(NIZ3)
         CI = CHG(NIX0)
         DIX= DIP(NIX3)
         DIY= DIP(NIY3)
         DIZ= DIP(NIZ3)
         QIXX = QAD(NIX6  )
         QIYY = QAD(NIX6+1)
         QIZZ = QAD(NIX6+2)
         QIXY = QAD(NIX6+3)
         QIXZ = QAD(NIX6+4)
         QIYZ = QAD(NIX6+5)
         OIXXX = OCT(NIX10  )
         OIYYY = OCT(NIX10+1)
         OIZZZ = OCT(NIX10+2)
         OIXXY = OCT(NIX10+3)
         OIXXZ = OCT(NIX10+4)
         OIXYY = OCT(NIX10+5)
         OIYYZ = OCT(NIX10+6)
         OIXZZ = OCT(NIX10+7)
         OIYZZ = OCT(NIX10+8)
         OIXYZ = OCT(NIX10+9)
C
C        --- LOOP OVER SURROUNDING MOLECULE ---
C
         NMOLJ = NI
         DO J=2,NMOLS
            NJ = NDMA(J)
            NMOLJ = NMOLJ + NJ
            DO N=1,NJ
               NJX0 =   (NMOLJ-NJ) +    N       
               NJX3 = 3*(NMOLJ-NJ) + 3*(N-1) + 1
               NJX6 = 6*(NMOLJ-NJ) + 6*(N-1) + 1
               NJX10=10*(NMOLJ-NJ) +10*(N-1) + 1
C                                                        
               NJY3 = NJX3 + 1
               NJZ3 = NJY3 + 1
C                                                        
C              UNPACK FOR MOLECULE J
C                                                        
               RX = RDMA(NJX3) - RIX
               RY = RDMA(NJY3) - RIY
               RZ = RDMA(NJZ3) - RIZ
C                                                        
               CJ = CHG(NJX0)
               DJX= DIP(NJX3)
               DJY= DIP(NJY3)
               DJZ= DIP(NJZ3)
               QJXX = QAD(NJX6  )
               QJYY = QAD(NJX6+1)
               QJZZ = QAD(NJX6+2)
               QJXY = QAD(NJX6+3)
               QJXZ = QAD(NJX6+4)
               QJYZ = QAD(NJX6+5)
               OJXXX = OCT(NJX10  )
               OJYYY = OCT(NJX10+1)
               OJZZZ = OCT(NJX10+2)
               OJXXY = OCT(NJX10+3)
               OJXXZ = OCT(NJX10+4)
               OJXYY = OCT(NJX10+5)
               OJYYZ = OCT(NJX10+6)
               OJXZZ = OCT(NJX10+7)
               OJYZZ = OCT(NJX10+8)
               OJXYZ = OCT(NJX10+9)
C                                                        
               RMN= DSQRT(RX*RX+RY*RY+RZ*RZ)
               RMN3 = ONE/(RMN*RMN*RMN)
               RMN5 = RMN3/(RMN*RMN)
               RMN7 = RMN5/(RMN*RMN)
C                                                        
C              TENSORDOTS
C                                                        
               S1 = -(DJX*RX+DJY*RY+DJZ*RZ)
               S2 =  (DIX*RX+DIY*RY+DIZ*RZ)
C                                                        
               S3 =  (DIX*DJX+DIY*DJY+DIZ*DJZ)
               S4 =   S1 * S2
C                                                        
               S5 = (QJXX * RX * RX       +  
     &               QJXY * RX * RY * TWO +  
     &               QJXZ * RX * RZ * TWO +  
     &               QJYY * RY * RY       +  
     &               QJYZ * RY * RZ * TWO +  
     &               QJZZ * RZ * RZ)
C                                                        
               S6 = (QIXX * RX * RX       +  
     &               QIXY * RX * RY * TWO +  
     &               QIXZ * RX * RZ * TWO +  
     &               QIYY * RY * RY       +  
     &               QIYZ * RY * RZ * TWO +  
     &               QIZZ * RZ * RZ)
C                                                        
               S7 =-(QJXX * DIX * RX +  
     &               QJXY * DIX * RY +
     &               QJXZ * DIX * RZ +
     &               QJXY * DIY * RX +
     &               QJYY * DIY * RY +
     &               QJYZ * DIY * RZ +
     &               QJXZ * DIZ * RX +
     &               QJYZ * DIZ * RY +
     &               QJZZ * DIZ * RZ)
C                                                        
               S8 = (QIXX * DJX * RX +  
     &               QIXY * DJX * RY +
     &               QIXZ * DJX * RZ +
     &               QIXY * DJY * RX +
     &               QIYY * DJY * RY +
     &               QIYZ * DJY * RZ +
     &               QIXZ * DJZ * RX +
     &               QIYZ * DJZ * RY +
     &               QIZZ * DJZ * RZ)
C        
               S9 =  S5 * S2
               S10=  S6 * S1
C                                                        
               RXRXRX = RX * RX * RX
               RXRXRY = RX * RX * RY * THREE
               RXRYRY = RX * RY * RY * THREE
               RYRYRY = RY * RY * RY
               RYRYRZ = RY * RY * RZ * THREE
               RYRZRZ = RY * RZ * RZ * THREE
               RZRZRZ = RZ * RZ * RZ
               RXRYRZ = RX * RY * RZ * SIX
               RXRXRZ = RX * RX * RZ * THREE
               RXRZRZ = RX * RZ * RZ * THREE
C                                                        
               S11=-(OJXXX * RXRXRX +
     &               OJXXY * RXRXRY +
     &               OJXYY * RXRYRY +
     &               OJYYY * RYRYRY +
     &               OJYYZ * RYRYRZ +
     &               OJYZZ * RYRZRZ +
     &               OJZZZ * RZRZRZ +
     &               OJXYZ * RXRYRZ +
     &               OJXXZ * RXRXRZ +
     &               OJXZZ * RXRZRZ)
C                                                        
               S12=  OIXXX * RXRXRX +
     &               OIXXY * RXRXRY +
     &               OIXYY * RXRYRY +
     &               OIYYY * RYRYRY +
     &               OIYYZ * RYRYRZ +
     &               OIYZZ * RYRZRZ +
     &               OIZZZ * RZRZRZ +
     &               OIXYZ * RXRYRZ +
     &               OIXXZ * RXRXRZ +
     &               OIXZZ * RXRZRZ 
C                                                        
C              ACCUMULATE THE TERMS
C                                                        
               CC = CC + CI * CJ / RMN 
               CD = CD + S1 * CI * RMN3
               DC = DC + S2 * CJ * RMN3 
               CQ = CQ + S5 * CI * RMN5
               QC = QC + S6 * CJ * RMN5 
               CT = CT + S11 * CI * RMN7
               TC = TC + S12 * CJ * RMN7
               DD = DD + (S3 * RMN3 + THREE * S4 * RMN5)
               DQ = DQ + (TWO * S7 * RMN5 + FIVE * S9  * RMN7 )
               QD = QD + (TWO * S8 * RMN5 + FIVE * S10 * RMN7 )
            ENDDO
         ENDDO
C
 90   CONTINUE
C
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
      EINT = D
C
      IF (LWRITE) THEN
          WRITE(*,91) "CC",CC,"CD",CD,"DC",DC,"CQ",CQ,"QC",QC,
     &                "CO",CT,"OC",TC,"DD",DD,"DQ",DQ,"QD",QD
      ENDIF
 91   FORMAT(10(3X,A," = ",F10.4,/))

      RETURN
      END

C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE EDMTPA(RDMA,CHG,DIP,QAD,OCT,NDMA,
     &                  EINT,NMOLS,NDMAS,LWRITE)
C
C          PAIRWISE ALL ELECTROSTATIC ENERGY FROM DISTRIBUTED MULTIPOLES
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RDMA(NDMAS*3),CHG(NDMAS),DIP(NDMAS*3),QAD(NDMAS*6),
     &          OCT(NDMAS*10),NDMA(NMOLS)
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,THREE=3.0D+00,
     &           FOUR=4.0D+00,FIVE=5.0D+00,SIX=6.0D+00)
      LOGICAL LWRITE
Cf2py INTENT(OUT) EINT
      CC = ZERO
      CD = ZERO
c      DC = ZERO
      DD = ZERO
      CQ = ZERO
c      QC = ZERO
      CT = ZERO
c      TC = ZERO
      DQ = ZERO
c      QD = ZERO
C
C     --- LOOP OVER MOLECULE I ---
C
      NMOLI = 0
      DO 90 I=1,NMOLS
         NI = NDMA(I)
         NMOLI = NMOLI + NI
         DO M=1,NI
            NIX0 =   (NMOLI-NI) +    M
            NIX3 = 3*(NMOLI-NI) + 3*(M-1) + 1
            NIX6 = 6*(NMOLI-NI) + 6*(M-1) + 1
            NIX10=10*(NMOLI-NI) +10*(M-1) + 1
C
            NIY3 = NIX3 + 1
            NIZ3 = NIY3 + 1
C
C           UNPACK FOR MOLECULE I
C
            RIX = RDMA(NIX3)
            RIY = RDMA(NIY3)
            RIZ = RDMA(NIZ3)
            CI = CHG(NIX0)
            DIX= DIP(NIX3)
            DIY= DIP(NIY3)
            DIZ= DIP(NIZ3)
            QIXX = QAD(NIX6  )
            QIYY = QAD(NIX6+1)
            QIZZ = QAD(NIX6+2)
            QIXY = QAD(NIX6+3)
            QIXZ = QAD(NIX6+4)
            QIYZ = QAD(NIX6+5)
            OIXXX = OCT(NIX10  )
            OIYYY = OCT(NIX10+1)
            OIZZZ = OCT(NIX10+2)
            OIXXY = OCT(NIX10+3)
            OIXXZ = OCT(NIX10+4)
            OIXYY = OCT(NIX10+5)
            OIYYZ = OCT(NIX10+6)
            OIXZZ = OCT(NIX10+7)
            OIYZZ = OCT(NIX10+8)
            OIXYZ = OCT(NIX10+9)
C
C           --- LOOP OVER MOLECULE J ---
C
            NMOLJ = 0
            DO J=1,I-1
               NJ = NDMA(J)
               NMOLJ = NMOLJ + NJ
               DO N=1,NJ
                  NJX0 =   (NMOLJ-NJ) +    N       
                  NJX3 = 3*(NMOLJ-NJ) + 3*(N-1) + 1
                  NJX6 = 6*(NMOLJ-NJ) + 6*(N-1) + 1
                  NJX10=10*(NMOLJ-NJ) +10*(N-1) + 1
C
                  NJY3 = NJX3 + 1
                  NJZ3 = NJY3 + 1
C
C                 UNPACK FOR MOLECULE J
C
                  RX = RDMA(NJX3) - RIX
                  RY = RDMA(NJY3) - RIY
                  RZ = RDMA(NJZ3) - RIZ
C
                  CJ = CHG(NJX0)
                  DJX= DIP(NJX3)
                  DJY= DIP(NJY3)
                  DJZ= DIP(NJZ3)
                  QJXX = QAD(NJX6  )
                  QJYY = QAD(NJX6+1)
                  QJZZ = QAD(NJX6+2)
                  QJXY = QAD(NJX6+3)
                  QJXZ = QAD(NJX6+4)
                  QJYZ = QAD(NJX6+5)
                  OJXXX = OCT(NJX10  )
                  OJYYY = OCT(NJX10+1)
                  OJZZZ = OCT(NJX10+2)
                  OJXXY = OCT(NJX10+3)
                  OJXXZ = OCT(NJX10+4)
                  OJXYY = OCT(NJX10+5)
                  OJYYZ = OCT(NJX10+6)
                  OJXZZ = OCT(NJX10+7)
                  OJYZZ = OCT(NJX10+8)
                  OJXYZ = OCT(NJX10+9)
C
                  RMN= DSQRT(RX*RX+RY*RY+RZ*RZ)
                  RMN3 = ONE/(RMN*RMN*RMN)
                  RMN5 = RMN3/(RMN*RMN)
                  RMN7 = RMN5/(RMN*RMN)
C
C                 TENSORDOTS
C
                  S1 = -(DJX*RX+DJY*RY+DJZ*RZ)
                  S2 =  (DIX*RX+DIY*RY+DIZ*RZ)
C
                  S3 =  (DIX*DJX+DIY*DJY+DIZ*DJZ)
                  S4 =   S1 * S2
C
                  S5 = (QJXX * RX * RX       +  
     &                  QJXY * RX * RY * TWO +  
     &                  QJXZ * RX * RZ * TWO +  
     &                  QJYY * RY * RY       +  
     &                  QJYZ * RY * RZ * TWO +  
     &                  QJZZ * RZ * RZ)
C
                  S6 = (QIXX * RX * RX       +  
     &                  QIXY * RX * RY * TWO +  
     &                  QIXZ * RX * RZ * TWO +  
     &                  QIYY * RY * RY       +  
     &                  QIYZ * RY * RZ * TWO +  
     &                  QIZZ * RZ * RZ)
C
                  S7 =-(QJXX * DIX * RX +  
     &                  QJXY * DIX * RY +
     &                  QJXZ * DIX * RZ +
     &                  QJXY * DIY * RX +
     &                  QJYY * DIY * RY +
     &                  QJYZ * DIY * RZ +
     &                  QJXZ * DIZ * RX +
     &                  QJYZ * DIZ * RY +
     &                  QJZZ * DIZ * RZ)
C
                  S8 = (QIXX * DJX * RX +  
     &                  QIXY * DJX * RY +
     &                  QIXZ * DJX * RZ +
     &                  QIXY * DJY * RX +
     &                  QIYY * DJY * RY +
     &                  QIYZ * DJY * RZ +
     &                  QIXZ * DJZ * RX +
     &                  QIYZ * DJZ * RY +
     &                  QIZZ * DJZ * RZ)
C        
                  S9 =  S5 * S2
                  S10=  S6 * S1
C
                  RXRXRX = RX * RX * RX
                  RXRXRY = RX * RX * RY * THREE
                  RXRYRY = RX * RY * RY * THREE
                  RYRYRY = RY * RY * RY
                  RYRYRZ = RY * RY * RZ * THREE
                  RYRZRZ = RY * RZ * RZ * THREE
                  RZRZRZ = RZ * RZ * RZ
                  RXRYRZ = RX * RY * RZ * SIX
                  RXRXRZ = RX * RX * RZ * THREE
                  RXRZRZ = RX * RZ * RZ * THREE
C
                  S11=-(OJXXX * RXRXRX +
     &                  OJXXY * RXRXRY +
     &                  OJXYY * RXRYRY +
     &                  OJYYY * RYRYRY +
     &                  OJYYZ * RYRYRZ +
     &                  OJYZZ * RYRZRZ +
     &                  OJZZZ * RZRZRZ +
     &                  OJXYZ * RXRYRZ +
     &                  OJXXZ * RXRXRZ +
     &                  OJXZZ * RXRZRZ)
C
                  S12=  OIXXX * RXRXRX +
     &                  OIXXY * RXRXRY +
     &                  OIXYY * RXRYRY +
     &                  OIYYY * RYRYRY +
     &                  OIYYZ * RYRYRZ +
     &                  OIYZZ * RYRZRZ +
     &                  OIZZZ * RZRZRZ +
     &                  OIXYZ * RXRYRZ +
     &                  OIXXZ * RXRXRZ +
     &                  OIXZZ * RXRZRZ 
C
C                 ACCUMULATE THE TERMS
C
                  CC = CC + (CI * CJ / RMN )
                  CD = CD + (S1 * CI + S2 * CJ) * RMN3
                  CQ = CQ + (S5 * CI + S6 * CJ) * RMN5
                  CT = CT + (S11 * CI + S12 * CJ) * RMN7
                  DD = DD + (S3 * RMN3 + THREE * S4 * RMN5)
                  DQ = DQ + (TWO  * (S7 + S8 ) * RMN5 + 
     &                       FIVE * (S9 + S10) * RMN7 )
               ENDDO
            ENDDO
         ENDDO
C
 90   CONTINUE
      EINT = CC + CD + CQ + CT + DD + DQ
C
      IF (LWRITE) THEN
          WRITE(*,91) "CC",CC,"CD",CD,"CQ",CQ,"CO",CT,"DD",DD,"DQ",DQ
      ENDIF
 91   FORMAT(6(3X,A," = ",F10.4,/))
C
      RETURN
      END
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
