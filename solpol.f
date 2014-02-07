C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE MOLLST(RC,RM,IC,IP,IE,ICM,IPM,IEM, 
     &                  NAT,CCUT,PCUT,ECUT,
     &                  NMOLS,NCOORD,NC)
C
C -----------------------------------------------------------------------------
C
C         DETERMINE THE MOLECULES LYING WITHIN COULOMB, POLARIZATION 
C            AND EXCHANGE-REPULSION RADII FROM CENTRAL MOLECULE
C
C              Bartosz Błasiak                       12 Nov 2013
C
C -----------------------------------------------------------------------------
C
C   Input variables:
C
C     ** Double precision
C     RC         - coordinates of central molecule (dimension 3,NC)
C     RM         - array of coordinates of other molecules (dimension 3,NCOORD)
C     CCUT       - Coulomb cutoff distance
C     PCUT       - Polarization cutoff distance
C     ECUT       - Exchange-repulsion cutoff
C
C     ** Integer
C     IC         - array of condition numbers for Coulomb sphere
C     IP         - array of condition numbers for Polarization sphere
C     IE         - array of condition numbers for Short-Range sphere
C     NMOLS      - number of other molecules (apart from central one)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RC(NC,3),RM(3*NCOORD),NAT(NMOLS)
      LOGICAL IC(NCOORD),IP(NCOORD),IE(NCOORD),
     &        ICM(NMOLS),IPM(NMOLS),IEM(NMOLS)
      PARAMETER (ZERO=0.0D+00)
Cf2py INTENT(IN,OUT) IC,IP,IE,ICM,IPM,IEM
C
C     CENTER OF GEOMETRY OF CENTRAL MOLECULE
C
      RCX = ZERO
      RCY = ZERO
      RCZ = ZERO
      DNC = DFLOAT(NC)
C
      DO I=1,NC
         RCX = RCX + RC(I,1)
         RCY = RCY + RC(I,2)
         RCZ = RCZ + RC(I,3)
      ENDDO
C
      RCX = RCX / DNC
      RCY = RCY / DNC
      RCZ = RCZ / DNC
C
C     LOOP OVER ALL OTHER MOLECULES
C
      NATSUM = 0
      DO 99 I=1,NMOLS
         NATI = NAT(I)
         DNA = DFLOAT(NATI)
         NATSUM = NATSUM + NATI
C
C        CENTER OF MASS OF A MOLECULE
C
         RMX = ZERO
         RMY = ZERO
         RMZ = ZERO
C
         DO J=1,NATI
            IX = 3*(NATSUM-NATI) + 3*(J-1) + 1
            IY = IX + 1
            IZ = IY + 1
C
            RMX = RMX + RM(IX)
            RMY = RMY + RM(IY) 
            RMZ = RMZ + RM(IZ)
         ENDDO
C
         RMX = RMX / DNA
         RMY = RMY / DNA
         RMZ = RMZ / DNA
C
         DX = RCX - RMX
         DY = RCY - RMY
         DZ = RCZ - RMZ
         RR = DSQRT(DX*DX+DY*DY+DZ*DZ)
C
         IF (RR.LT.CCUT) THEN
             ICM(I) = .TRUE.
             DO J=1,NATI
                IX = NATSUM-NATI + J
                IC(IX) = .TRUE.
             ENDDO
         ENDIF
C
         IF (RR.LT.PCUT) THEN
             IPM(I) = .TRUE.
             DO J=1,NATI
                IX = NATSUM-NATI + J
                IP(IX) = .TRUE.
             ENDDO
         ENDIF
C
         IF (RR.LT.ECUT) THEN
             IEM(I) = .TRUE.
             DO J=1,NATI
                IX = NATSUM-NATI + J
                IE(IX) = .TRUE.
             ENDDO
         ENDIF
C
99    CONTINUE
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE SFTPOL(RDMA,CHG,DIP,QAD,OCT,CHG1,DIP1,QAD1,OCT1,RPOL,
     *                  POL,EPOL,SHIFT,
     *                  DMAT,FLDS,DIPIND,DIMAT,FIVEC,SDIPND,AVEC,VEC1,
     *                  MAT1,REDMSS,FREQ,GIJJ,RPOL1,POL1,LVEC,
     *                  NMOLS,NDMA,NPOL,NDIM,NDMAS,MODE,NMODES,NDMAC,
     *                  NPOLC,MDIP,MQAD,MOCT,MRPOL,MPOL,LWRITE)
C
C -----------------------------------------------------------------------------
C
C         ELECTROSTATIC POLARIZATION FREQUENCY SHIFT FROM MULTIPOLE EXPANSION
C                  DISTRIBUTED DIPOLE POLARIZABILITY MODEL
C 
C              Bartosz Błasiak                        20 Nov 2013
C
C -----------------------------------------------------------------------------
C
C   Input variables:
C
C   Returns:
C     EPOL    - polarization energy
C     SHIFT   - frequency shift
C
C   External:
C     ILAENV  - evaluate optimum block size (lapack)
C     DGETRF  - LU decomposition (lapack) 
C     DGETRI  - rectanqular real matrix inversion (lapack)
C     DDOT    - dot product of two real vectors (lapack)
C     DGEMM   - matrix-matrix multiplication
C
C   Notes:
C     The reduced format of tensor storage is used
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RDMA(MDIP),CHG(NDMAS),DIP(MDIP),QAD(MQAD),OCT(MOCT),
     &          RPOL(MRPOL),POL(MPOL),DMAT(NDIM,NDIM),FLDS(NDIM),
     &          DIP1(NMODES*NDMAC*3),QAD1(NMODES*NDMAC*6),
     &          OCT1(NMODES*NDMAC*10),CHG1(NMODES*NDMAC),
     &          DIMAT(NDIM,NDIM),FIVEC(NDIM),SDIPND(NDIM),AVEC(NDIM),
     &          REDMSS(NMODES),FREQ(NMODES),GIJJ(NMODES),
     &          POL1(NMODES*NPOLC*9),
     &          DIPIND(NDIM),NDMA(NMOLS),NPOL(NMOLS),
     &          RPOL1(NMODES*NPOLC*3),
     &          IPIV(10000),WORK(5000000),FI(40),VEC1(NDIM)
      PARAMETER (ZERO=0.0D+00,HALF=0.50D+00,ONE=1.00D+00)
      DOUBLE PRECISION DDOT,MAT1(NDIM,NDIM)
      DOUBLE PRECISION LVEC((NMODES+6)*NMODES)
      LOGICAL LWRITE
      EXTERNAL ILAENV,DGETRF,DGETRI,DDOT,DGEMM
Cf2py INTENT(OUT) EPOL, SHIFT
      EPOL = ZERO
      SHIFT= ZERO
      DATA IPIV/10000*0/
      DATA WORK/5000000*0.0D+00/
      DATA FI/40*0.0D+00/
C
C     CALCULATE FIELDS AT POLARIZABLE CENTERS AND D-MATRIX
C
      CALL FIELDS(RDMA,CHG,DIP,QAD,OCT,RPOL,POL,DMAT,FLDS,
     ^            NMOLS,NPOL,NDMA,NDIM,NDMAS,
     ^            MDIP,MQAD,MOCT,MRPOL,MPOL)
C
C     COPY OFFDIAGONAL BLOCKS OF DMAT TO MAT1 TO FORM T-TENSOR. 
C     LEAVE THE DIAGONAL QUADRATIC 3x3 BLOCKS TO BE ZERO
C
      DO 10 I=1,NDIM
      DO 10 J=1,(I-3)
         DMIJ = DMAT(I,J)
         MAT1(I,J) = DMIJ
         MAT1(J,I) = DMIJ
 10   CONTINUE
C
C     THEN INVERT D-MATRIX
C
      NB = ILAENV(1,"DGETRI","NALIWKU",NDIM,-1,-1,-1)
      LWORK = NDIM*NB
      CALL DGETRF(NDIM,NDIM,DMAT,NDIM,IPIV,INFO)
      IF (LWRITE) WRITE(*,*) " LWORK= ", LWORK, "INFO= ",INFO
      CALL DGETRI(NDIM,DMAT,NDIM,IPIV,WORK,LWORK,INFO)
      CALL CHECK(INFO)
      IF (INFO.NE.0) GOTO 1123
C
C     CALCULATE INDUCED DIPOLES AND POLARIZATION ENERGY
C
      CALL DGMV('N',DMAT,FLDS,DIPIND,NDIM)
      EPOL = - DDOT(NDIM,FLDS,1,DIPIND,1) * HALF
C
      IF (LWRITE) THEN
          CALL MATWRT(DMAT,NDIM,NDIM,-1,"dmat.dat")
          CALL VECWRT(SDIPND,NDIM,-1,"sdipnd.dat")
      ENDIF
C
C     CALCULATE DIMAT AND FIVEC AND ACCUMULATE THEM TO -AVEC-
C
      CALL AVECEV(RDMA,CHG,DIP,QAD,OCT,CHG1,DIP1,QAD1,OCT1,
     *            RPOL,RPOL1,POL,DMAT,FLDS,
     *            DIMAT,FIVEC,AVEC,VEC1,SDIPND,
     *            GIJJ,REDMSS,FREQ,POL1,LVEC,
     *            NMOLS,NPOL,NDMA,NDIM,NDMAS,
     *            MDIP,MQAD,MOCT,MRPOL,MPOL,MODE,NMODES,
     *            NDMAC,NPOLC)
C
      IF (LWRITE) THEN
          CALL VECWRT(AVEC,NDIM,-1,"avec.dat")
          CALL MATWRT(DIMAT,NDIM,NDIM,-1,"dimat.dat")
      ENDIF
C
C     CALCULATE FREQUENCY SHIFTS
C
      SHIFT = - DDOT(NDIM,FLDS,1,AVEC,1) * HALF
      CALL VECWRT(AVEC,NDIM,-1,"avec.dat")
C
C     EVALUATE [ 1 +  T D-1 ]-1 MATRIX
C
      CALL DGEMM('N','N',NDIM,NDIM,NDIM,ONE,
     &                   MAT1,NDIM,DMAT,NDIM,
     &                   ZERO,DIMAT,NDIM)
C
      DO I=1,NDIM
         DIMAT(I,I) = DIMAT(I,I) + ONE
      ENDDO
C
      CALL DGETRF(NDIM,NDIM,DIMAT,NDIM,IPIV,INFO)
      CALL DGETRI(NDIM,DIMAT,NDIM,IPIV,WORK,LWORK,INFO)
      CALL CHECK(INFO)
      IF (INFO.NE.0) GOTO 1123
C
C     CALCULATE SOLVATOCHROMIC INDUCED DIPOLE MOMENTS
C
      CALL DGMV('N',DIMAT,AVEC,SDIPND,NDIM)
C
      IF (LWRITE) THEN
          CALL VECWRT(SDIPND,NDIM,-1,"sdipnd.dat")
          CALL VECWRT(FLDS,NDIM,-1,"fields.dat")
      ENDIF
      CALL VECWRT(SDIPND,NDIM,-1,"sdipnd.dat")
C
 1123 CONTINUE
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE AVECEV(RDMA,CHG,DIP,QAD,OCT,CHG1,DIP1,QAD1,OCT1,
     *                  RPOL,RPOL1,POL,DMAT,FLDS,
     *                  DIMAT,FIVEC,AVEC,VEC1,SDIPND,
     *                  GIJJ,REDMSS,FREQ,POL1,LVEC,
     *                  NMOLS,NPOL,NDMA,NDIM,NDMAS,
     *                  MDIP,MQAD,MOCT,MRPOL,MPOL,MODE,NMODES,
     *                  NDMAC,NPOLC)
C     -----------------------------------------------------------------------
C     EVALUATE DIMAT AND FIVEC DUE TO EFP FRAGMENTS AND THEN CONSTRUCT -AVEC-
C     
C     -----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RDMA(MDIP),CHG(NDMAS),DIP(MDIP),
     &          QAD(MQAD),OCT(MOCT),CHG1(NMODES*NDMAC),
     &          DIP1(NMODES*NDMAC*3),QAD1(NMODES*NDMAC*6),
     &          OCT1(NMODES*NDMAC*10),
     &          RPOL1(NMODES*NPOLC*3),
     &          RPOL(MRPOL),POL(MPOL),DMAT(NDIM,NDIM),FLDS(NDIM),
     &          DIMAT(NDIM,NDIM),FIVEC(NDIM),AVEC(NDIM),
     &          GIJJ(NMODES),REDMSS(NMODES),FREQ(NMODES),
     &          POL1(NMODES*NPOLC*9),
     &          NPOL(NMOLS),NDMA(NMOLS),VEC1(NDIM),SDIPND(NDIM),
     &          WORKI(30),APOL(3,3),IPIVP(3),PM(3,3),PMT(3,3),GIVEC(30)
      DOUBLE PRECISION MAT2(NDIM,NDIM),MAT3(NDIM,NDIM)
      DOUBLE PRECISION LVEC((NMODES+6)*NMODES)
      EXTERNAL DGETRI,DGETRF,DGEMM
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,THREE=3.0D+00,
     &           FOUR=4.0D+00,FIVE=5.0D+00,SIX=6.0D+00,SEVEN=7.0D+00,
     &           QNINE=9.0D+00,HALF=0.50D+00)
      DATA WORKI/30*0.0D+00/
      DATA APOL /9*0.0D+00/
      DATA PM   /9*0.0D+00/
      DATA PMT  /9*0.0D+00/
      DATA IPIVP/3*0/
      DATA GIVEC/30*0.0D+00/
C
C     AUXILIARY MODE VECTOR
C
      GJ = ONE / (REDMSS(MODE) * FREQ(MODE))
      DO MM=1,NMODES
         FREQM = FREQ(MM)
         GIVEC(MM) = GJ * GIJJ(MM)/(REDMSS(MM)*FREQM*FREQM)
      ENDDO
C
C     ITERATE OVER IR-OTHER MOLECULE PAIRS AND THEIR POLARIZABLE CENTERS
C
      NPOLI = 0 
      DO 7878 IMOL=1,1
      NIM  = NPOL(IMOL)
      NPOLI = NPOLI + NIM
      DO 6767 I=1,NIM
         NIX3 = 3*(NPOLI-NIM) + 3*(I-1) + 1
         NIX9 = 9*(NPOLI-NIM) + 9*(I-1) + 1
C 
         NIY3 = NIX3 + 1
         NIZ3 = NIY3 + 1
C
         RPOLIX = RPOL(NIX3)
         RPOLIY = RPOL(NIY3)
         RPOLIZ = RPOL(NIZ3)
C
C        --- CALCULATE DIMAT DIAGONALS ---
C
         APOL(1,1) = POL(NIX9  )
         APOL(1,2) = POL(NIX9+1)
         APOL(1,3) = POL(NIX9+2)
         APOL(2,1) = POL(NIX9+3)
         APOL(2,2) = POL(NIX9+4)
         APOL(2,3) = POL(NIX9+5)
         APOL(3,1) = POL(NIX9+6)
         APOL(3,2) = POL(NIX9+7)
         APOL(3,3) = POL(NIX9+8)
C
C        INVERT DISTRIBUTED POLARIZABILITY
C
         CALL DGETRF(3,3,APOL,3,IPIVP,INFO)
         CALL DGETRI(3,APOL,3,IPIVP,WORKI,20,INFO)
C
         CALL CHECK(INFO)
         IF (INFO.NE.0) GOTO 19191
C
C        ITERATE OVER NORMAL COORDINATES
C
         DIXIX = ZERO 
         DIYIY = ZERO
         DIZIZ = ZERO
C
         DIXIY = ZERO
         DIXIZ = ZERO
         DIYIZ = ZERO
C
         DIYIX = ZERO
         DIZIX = ZERO
         DIZIY = ZERO
C
         FLDX1 = ZERO
         FLDY1 = ZERO
         FLDZ1 = ZERO
C
         DO MM=1,NMODES
            MNIX9 = 9*NPOLC*(MM-1) + 9*(I-1) + 1
            MNIX3 = 3*NPOLC*(MM-1) + 3*(I-1) + 1
            GRF = GIVEC(MM) / TWO
C
            PM(1,1) = POL1(MNIX9  )
            PM(1,2) = POL1(MNIX9+1)
            PM(1,3) = POL1(MNIX9+2)
            PM(2,1) = POL1(MNIX9+3)
            PM(2,2) = POL1(MNIX9+4)
            PM(2,3) = POL1(MNIX9+5)
            PM(3,1) = POL1(MNIX9+6)
            PM(3,2) = POL1(MNIX9+7)
            PM(3,3) = POL1(MNIX9+8)
C
C           EVALUATE (a-1) da/dQ (a-1)
C
            CALL DGEMM('N','N',3,3,3,ONE,APOL,3,PM,3,ZERO,PMT,3)
            CALL DGEMM('N','N',3,3,3,ONE,PMT,3,APOL,3,ZERO,PM,3)
C
            DIXIX = DIXIX - PM(1,1) * GRF
            DIYIY = DIYIY - PM(2,2) * GRF
            DIZIZ = DIZIZ - PM(3,3) * GRF
C
            DIXIY = DIXIY - PM(1,2) * GRF
            DIXIZ = DIXIZ - PM(1,3) * GRF
            DIYIZ = DIYIZ - PM(2,3) * GRF
C 
            DIYIX = DIYIX - PM(2,1) * GRF
            DIZIX = DIZIX - PM(3,1) * GRF
            DIZIY = DIZIY - PM(3,2) * GRF
C
C           --- CALCULATE FIELD DERIVATIVES ON SOLUTE CENTERS ---
C
            NDMAJ = NDMA(1)
            DO JMOL=2,NMOLS
               NJM  = NDMA(JMOL)
               NDMAJ = NDMAJ + NJM
               DO J=1,NJM
                  NJX0 =    (NDMAJ-NJM) +    J
                  NJX3 =  3*(NDMAJ-NJM) + 3*(J-1) + 1
                  NJX6 =  6*(NDMAJ-NJM) + 6*(J-1) + 1
                  NJX10= 10*(NDMAJ-NJM) +10*(J-1) + 1
C
                  NJY3 = NJX3 + 1
                  NJZ3 = NJY3 + 1
C
                  RIJX = RPOLIX - RDMA(NJX3) 
                  RIJY = RPOLIY - RDMA(NJY3)
                  RIJZ = RPOLIZ - RDMA(NJZ3)
C
                  RIJ  = DSQRT(RIJX*RIJX+
     &                         RIJY*RIJY+
     &                         RIJZ*RIJZ)
C
                  RIJ2 = ONE/(RIJ*RIJ)
                  RIJ3 = RIJ2 / RIJ
                  RIJ4 = RIJ2 * RIJ2
                  RIJ5 = RIJ3 * RIJ2
                  RIJ7 = RIJ5 * RIJ2
                  RIJ9 = RIJ7 * RIJ2
C
C                 UNPACK THE TENSORS
C
                  RX1 = RPOL1(MNIX3  )
                  RY1 = RPOL1(MNIX3+1)
                  RZ1 = RPOL1(MNIX3+2)
C
                  DIPX = DIP(NJX3)
                  DIPY = DIP(NJY3)
                  DIPZ = DIP(NJZ3)
C
                  QXX= QAD(NJX6  )
                  QYY= QAD(NJX6+1)
                  QZZ= QAD(NJX6+2)
                  QXY= QAD(NJX6+3)
                  QXZ= QAD(NJX6+4)
                  QYZ= QAD(NJX6+5)
C
                  OXXX=OCT(NJX10  )
                  OYYY=OCT(NJX10+1)
                  OZZZ=OCT(NJX10+2)
C
                  OXXY=OCT(NJX10+3)
                  OXXZ=OCT(NJX10+4)
                  OXYY=OCT(NJX10+5)
C
                  OYYZ=OCT(NJX10+6)
                  OXZZ=OCT(NJX10+7)
                  OYZZ=OCT(NJX10+8)
                  OXYZ=OCT(NJX10+9)
C
C                 AUXILIARY DOT PRODUCTS
C
                  RIJRI = RIJX * RX1 + 
     &                    RIJY * RY1 +
     &                    RIJZ * RZ1
C
                  DJRIJ = DIPX * RIJX +
     &                    DIPY * RIJY +
     &                    DIPZ * RIJZ
C
                  DJRI  = DIPX * RX1 +
     &                    DIPY * RY1 +
     &                    DIPZ * RZ1
C
                  QJRIJ2 = QXX * RIJX * RIJX       +
     &                     QXY * RIJX * RIJY * TWO +
     &                     QXZ * RIJX * RIJZ * TWO +
     &                     QYY * RIJY * RIJY       +
     &                     QYZ * RIJY * RIJZ * TWO +
     &                     QZZ * RIJZ * RIJZ
C
                  QJRIJI = QXX * RX1 * RIJX +
     &                     QXY * RX1 * RIJY +
     &                     QXY * RY1 * RIJX + 
     &                     QXZ * RX1 * RIJZ +
     &                     QXZ * RZ1 * RIJX +
     &                     QYY * RY1 * RIJY +
     &                     QYZ * RY1 * RIJZ +
     &                     QYZ * RZ1 * RIJY +
     &                     QZZ * RZ1 * RIJZ
C
               SUM1 = OXXX * RIJX * RIJX * RIJX         +
     &                OXXY * RIJX * RIJX * RIJY * THREE +
     &                OXYY * RIJX * RIJY * RIJY * THREE +
     &                OYYY * RIJY * RIJY * RIJY         +
     &                OYYZ * RIJY * RIJY * RIJZ * THREE +
     &                OYZZ * RIJY * RIJZ * RIJZ * THREE +
     &                OZZZ * RIJZ * RIJZ * RIJZ         +
     &                OXYZ * RIJX * RIJY * RIJZ * SIX   +
     &                OXXZ * RIJX * RIJX * RIJZ * THREE +
     &                OXZZ * RIJX * RIJZ * RIJZ * THREE
C
               SUM2 = OXXX * RX1 * RIJX * RIJX         +
     &                OXXY * RX1 * RIJX * RIJY * TWO   +
     &                OXXY * RY1 * RIJX * RIJX         +
     &                OXYY * RX1 * RIJY * RIJY         +
     &                OXYY * RY1 * RIJX * RIJY * TWO   +
     &                OYYY * RY1 * RIJY * RIJY         +
     &                OYYZ * RY1 * RIJY * RIJZ * TWO   +
     &                OYYZ * RZ1 * RIJY * RIJY         +
     &                OYZZ * RY1 * RIJZ * RIJZ         +
     &                OYZZ * RZ1 * RIJY * RIJZ * TWO   +
     &                OZZZ * RZ1 * RIJZ * RIJZ         +
     &                OXYZ * RX1 * RIJY * RIJZ * TWO   +
     &                OXYZ * RY1 * RIJX * RIJZ * TWO   +
     &                OXYZ * RZ1 * RIJX * RIJY * TWO   +
     &                OXXZ * RX1 * RIJX * RIJZ * TWO   +
     &                OXXZ * RZ1 * RIJX * RIJX         +
     &                OXZZ * RX1 * RIJZ * RIJZ         +
     &                OXZZ * RZ1 * RIJX * RIJZ
C
                  VOX  = OXXX * RIJX * RIJX       + 
     &                   OXXY * RIJX * RIJY * TWO +
     &                   OXXZ * RIJX * RIJZ * TWO +
     &                   OXYY * RIJY * RIJY       +
     &                   OXYZ * RIJY * RIJZ * TWO +
     &                   OXZZ * RIJZ * RIJZ
C
                  VOY  = OXXY * RIJX * RIJX       + 
     &                   OXYY * RIJX * RIJY * TWO +
     &                   OXYZ * RIJX * RIJZ * TWO +
     &                   OYYY * RIJY * RIJY       +
     &                   OYYZ * RIJY * RIJZ * TWO +
     &                   OYZZ * RIJZ * RIJZ
C
                  VOZ  = OXXZ * RIJX * RIJX       + 
     &                   OXYZ * RIJX * RIJY * TWO +
     &                   OXZZ * RIJX * RIJZ * TWO +
     &                   OYYZ * RIJY * RIJY       +
     &                   OYZZ * RIJY * RIJZ * TWO +
     &                   OZZZ * RIJZ * RIJZ
C
                  VOX1 = OXXX * RIJX * RX1       + 
     &                   OXXY * RIJX * RY1       +
     &                   OXXY * RIJY * RX1       +
     &                   OXXZ * RIJX * RZ1       +
     &                   OXXZ * RIJZ * RX1       +
     &                   OXYY * RIJY * RY1       +
     &                   OXYZ * RIJY * RZ1       +
     &                   OXYZ * RIJZ * RY1       +
     &                   OXZZ * RIJZ * RZ1
C
                  VOY1 = OXXY * RIJX * RX1       + 
     &                   OXYY * RIJX * RY1       +
     &                   OXYY * RIJY * RX1       +
     &                   OXYZ * RIJX * RZ1       +
     &                   OXYZ * RIJZ * RX1       +
     &                   OYYY * RIJY * RY1       +
     &                   OYYZ * RIJY * RZ1       +
     &                   OYYZ * RIJZ * RY1       +
     &                   OYZZ * RIJZ * RZ1
C
                  VOZ1 = OXXZ * RIJX * RX1       + 
     &                   OXYZ * RIJX * RY1       +
     &                   OXYZ * RIJY * RX1       +
     &                   OXZZ * RIJX * RZ1       +
     &                   OXZZ * RIJZ * RX1       +
     &                   OYYZ * RIJY * RY1       +
     &                   OYZZ * RIJY * RZ1       +
     &                   OYZZ * RIJZ * RY1       +
     &                   OZZZ * RIJZ * RZ1
C
                  VQX  = QXX * RIJX + QXY * RIJY + QXZ * RIJZ
                  VQY  = QXY * RIJX + QYY * RIJY + QYZ * RIJZ
                  VQZ  = QXZ * RIJX + QYZ * RIJY + QZZ * RIJZ
C
                  VQX1 = QXX * RX1 + QXY * RY1 + QXZ * RZ1
                  VQY1 = QXY * RX1 + QYY * RY1 + QYZ * RZ1
                  VQZ1 = QXZ * RX1 + QYZ * RY1 + QZZ * RZ1
C
                  RRR  = THREE * RIJRI * RIJ2
C
C                 CHARGE CONTRIBUTION
C
                  CHGJ = CHG(NJX0) * RIJ3
                  FX = CHGJ * ( RX1 - RRR * RIJX )
                  FY = CHGJ * ( RY1 - RRR * RIJY )
                  FZ = CHGJ * ( RZ1 - RRR * RIJZ )
C
C                 DIPOLE CONTRIBUTION
C
                  FRDI = FIVE * RIJ2 * DJRIJ * RIJRI
C
                  FX = FX + THREE * RIJ5 * ( DJRIJ * RX1 +
     &                 DJRI * RIJX - FRDI * RIJX + DIPX * RIJRI )
C
                  FY = FY + THREE * RIJ5 * ( DJRIJ * RY1 +
     &                 DJRI * RIJY - FRDI * RIJY + DIPY * RIJRI )
C
                  FZ = FZ + THREE * RIJ5 * ( DJRIJ * RZ1 +
     &                 DJRI * RIJZ - FRDI * RIJZ + DIPZ * RIJRI )
C
C                 QUADRUPOLE CONTRIBUTION
C
                  FVV = FIVE * RIJ2
C
                  FX = FX + RIJ5 * ( FVV * ( TWO   * RIJRI  * VQX + 
     &                                       TWO   * QJRIJI * RIJX +
     &                                      QJRIJ2 * RX1 ) - TWO * VQX1
     &                      - 35.00D+00 * RIJRI * QJRIJ2 * RIJ4 * RIJX )
C
                  FY = FY + RIJ5 * ( FVV * ( TWO   * RIJRI  * VQY + 
     &                                       TWO   * QJRIJI * RIJY +
     &                                      QJRIJ2 * RY1 ) - TWO * VQY1
     &                      - 35.00D+00 * RIJRI * QJRIJ2 * RIJ4 * RIJY )
C
                  FZ = FZ + RIJ5 * ( FVV * ( TWO   * RIJRI  * VQZ + 
     &                                       TWO   * QJRIJI * RIJZ +
     &                                      QJRIJ2 * RZ1 ) - TWO * VQZ1
     &                      - 35.00D+00 * RIJRI * QJRIJ2 * RIJ4 * RIJZ )
C
C                 OCTUPOLE CONTRIBUTION
C
                  FX = FX + SEVEN*RIJ9 * (SUM1*RX1 - 
     &                      QNINE*RIJ2*SUM1*RIJRI*RIJX + 
     &                      THREE*SUM2*RIJX +
     &                      THREE*RIJRI*VOX) -
     &                 SIX*RIJ7*VOX1
C
                  FY = FY + SEVEN*RIJ9 * (SUM1*RY1 - 
     &                      QNINE*RIJ2*SUM1*RIJRI*RIJY + 
     &                      THREE*SUM2*RIJY +
     &                      THREE*RIJRI*VOY) -
     &                 SIX*RIJ7*VOY1
C
                  FZ = FZ + SEVEN*RIJ9 * (SUM1*RZ1 - 
     &                      QNINE*RIJ2*SUM1*RIJRI*RIJZ + 
     &                      THREE*SUM2*RIJZ +
     &                      THREE*RIJRI*VOZ) -
     &                 SIX*RIJ7*VOZ1
C
C                 WEIGHT EACH CONTRIBUTION BY MODE COEFFICIENT GRF
C
                  FLDX1 = FLDX1 + FX * GRF
                  FLDY1 = FLDY1 + FY * GRF
                  FLDZ1 = FLDZ1 + FZ * GRF
C
               ENDDO
            ENDDO
         ENDDO
C
C        SAVE THE DIMAT DIAGONALS
C
         DIMAT(NIX3,NIX3) =  DIXIX 
         DIMAT(NIY3,NIY3) =  DIYIY 
         DIMAT(NIZ3,NIZ3) =  DIZIZ 
C
         DIMAT(NIX3,NIY3) =  DIXIY 
         DIMAT(NIX3,NIZ3) =  DIXIZ 
         DIMAT(NIY3,NIZ3) =  DIYIZ 
C
         DIMAT(NIY3,NIX3) =  DIYIX
         DIMAT(NIZ3,NIX3) =  DIZIX
         DIMAT(NIZ3,NIY3) =  DIZIY
C
C        SAVE FIELD DERIVATIVES
C
         FIVEC(NIX3) = FLDX1
         FIVEC(NIY3) = FLDY1
         FIVEC(NIZ3) = FLDZ1
C
C        --- CALCULATE DIMAT OFFDIAGONALS ---
C
         NPOLJ = NPOL(1)
         DO JMOL=2,NMOLS
         NJM   = NPOL(JMOL)
         NPOLJ = NPOLJ + NJM
c         IF (IMOL.EQ.JMOL) GOTO 921
         DO J=1,NJM
            NJX3 = 3*(NPOLJ-NJM) + 3*(J-1) + 1
            NJY3 = NJX3 + 1
            NJZ3 = NJY3 + 1
C
            RIJX = RPOLIX - RPOL(NJX3)
            RIJY = RPOLIY - RPOL(NJY3)
            RIJZ = RPOLIZ - RPOL(NJZ3)
C
            RIJ  = DSQRT(RIJX*RIJX+
     &                   RIJY*RIJY+
     &                   RIJZ*RIJZ)
C
            RIJ3 = ONE/(RIJ*RIJ*RIJ)
            RIJ5 = RIJ3/(RIJ*RIJ)
C
            THR5 = THREE * RIJ5
            FVR2 = FIVE / (RIJ * RIJ)
C
C           ITERATE OVER NORMAL COORDINATES
C
            DMXX = ZERO 
            DMXY = ZERO
            DMXZ = ZERO
            DMYY = ZERO
            DMYZ = ZERO
            DMZZ = ZERO
C
            DO MM=1,NMODES
               MNIX3 = 3*NPOLC*(MM-1) + 3*(I-1) + 1
               GRF = GIVEC(MM)
C
C              UNPACK THE DERIVATIVES OF LMO CENTROIDS
C
               RIX1 = RPOL1(MNIX3  )
               RIY1 = RPOL1(MNIX3+1)
               RIZ1 = RPOL1(MNIX3+2)
C
               R1R  = RIJX * RIX1 + 
     &                RIJY * RIY1 +
     &                RIJZ * RIZ1
C
               DMXX = DMXX + ( R1R * (ONE-FVR2*RIJX*RIJX) 
     &                            + TWO*RIJX*RIX1)         * GRF
               DMXY = DMXY + (-R1R *      FVR2*RIJX*RIJY  
     &                            + RIJX*RIY1 + RIX1*RIJY) * GRF
               DMXZ = DMXZ + (-R1R *      FVR2*RIJX*RIJZ  
     &                            + RIJX*RIZ1 + RIX1*RIJZ) * GRF
               DMYY = DMYY + ( R1R * (ONE-FVR2*RIJY*RIJY) 
     &                            + TWO*RIJY*RIY1)         * GRF
               DMYZ = DMYZ + (-R1R *      FVR2*RIJY*RIJZ
     &                            + RIJY*RIZ1 + RIY1*RIJZ) * GRF
               DMZZ = DMZZ + ( R1R * (ONE-FVR2*RIJZ*RIJZ) 
     &                            + TWO*RIJZ*RIZ1)         * GRF
            ENDDO
C
            THR5 = THR5 / TWO
C
            DMXX = DMXX * THR5
            DMXY = DMXY * THR5
            DMYY = DMYY * THR5
            DMXZ = DMXZ * THR5
            DMYZ = DMYZ * THR5
            DMZZ = DMZZ * THR5
C                                                                  
C           SAVE THE OFFDIAGONALS
C
            DIMAT(NIX3,NJX3) = -DMXX
            DIMAT(NIX3,NJY3) = -DMXY
            DIMAT(NIY3,NJX3) = -DMXY
            DIMAT(NIY3,NJY3) = -DMYY
            DIMAT(NIX3,NJZ3) = -DMXZ
            DIMAT(NIZ3,NJX3) = -DMXZ
            DIMAT(NIY3,NJZ3) = -DMYZ
            DIMAT(NIZ3,NJY3) = -DMYZ
            DIMAT(NIZ3,NJZ3) = -DMZZ
C
            DIMAT(NJX3,NIX3) = -DMXX
            DIMAT(NJY3,NIX3) = -DMXY
            DIMAT(NJX3,NIY3) = -DMXY
            DIMAT(NJY3,NIY3) = -DMYY
            DIMAT(NJZ3,NIX3) = -DMXZ
            DIMAT(NJX3,NIZ3) = -DMXZ
            DIMAT(NJZ3,NIY3) = -DMYZ
            DIMAT(NJY3,NIZ3) = -DMYZ
            DIMAT(NJZ3,NIZ3) = -DMZZ
C
         ENDDO
c 921     CONTINUE
         ENDDO
C
 6767 CONTINUE
 7878 CONTINUE
C
C     --- EVALUATE FIELD DERIVATIVES ON SOLVENT CENTERS ---
C
      NNN   = NMODES+6
      NNNH3 = NNN/3
      NNN6  = NNN*2
      NNN10 = NNNH3*10
c      NDMAC = NDMA(1)
      NPOLJ = NPOL(1)
C
      DO JMOL=2,NMOLS
         NJM   = NPOL(JMOL)
         NPOLJ = NPOLJ + NJM
C
         DO J=1,NJM
            NJX0 =   (NPOLJ-NJM) +    J
            NJX3 = 3*(NPOLJ-NJM) + 3*(J-1) + 1
            NJY3 = NJX3 + 1
            NJZ3 = NJY3 + 1
C
            RPOLJX = RPOL(NJX3)
            RPOLJY = RPOL(NJY3)
            RPOLJZ = RPOL(NJZ3)
C
            FLDX1 = ZERO
            FLDY1 = ZERO
            FLDZ1 = ZERO 
C
C           ITERATE OVER DMA POINTS OF CENTER MOLECULE
C
            DO I=1,NDMAC
               NIX3 = 3*(I-1) + 1
               NIX6 = 6*(I-1) + 1
               NIX10=10*(I-1) + 1
C
               NIY3 = NIX3 + 1
               NIZ3 = NIY3 + 1
C
               RDMAIX = RDMA(NIX3)
               RDMAIY = RDMA(NIY3)
               RDMAIZ = RDMA(NIZ3)
C
               RJIX = RPOLJX - RDMAIX
               RJIY = RPOLJY - RDMAIY
               RJIZ = RPOLJZ - RDMAIZ
C
               RJI  = DSQRT(RJIX*RJIX+
     &                      RJIY*RJIY+
     &                      RJIZ*RJIZ)
C
               RJI2 = ONE/(RJI*RJI)
               RJI3 = RJI2 / RJI
               RJI4 = RJI2 * RJI2
               RJI5 = RJI3 * RJI2
               RJI7 = RJI5 * RJI2
               RJI9 = RJI7 * RJI2
C
C              UNPACK THE DMA TENSORS
C
               CHGI = CHG(I)
C
               DIPX = DIP(NIX3)
               DIPY = DIP(NIY3)
               DIPZ = DIP(NIZ3)
C
               QXX= QAD(NIX6  )
               QYY= QAD(NIX6+1)
               QZZ= QAD(NIX6+2)
               QXY= QAD(NIX6+3)
               QXZ= QAD(NIX6+4)
               QYZ= QAD(NIX6+5)
C
               OXXX=OCT(NIX10  )
               OYYY=OCT(NIX10+1)
               OZZZ=OCT(NIX10+2)
C
               OXXY=OCT(NIX10+3)
               OXXZ=OCT(NIX10+4)
               OXYY=OCT(NIX10+5)
C
               OYYZ=OCT(NIX10+6)
               OXZZ=OCT(NIX10+7)
               OYZZ=OCT(NIX10+8)
               OXYZ=OCT(NIX10+9)
C
C              AUXILIARY DOT PRODUCTS
C
               DJRJI = DIPX * RJIX +
     &                 DIPY * RJIY +
     &                 DIPZ * RJIZ
C
               QJRJI2 = QXX * RJIX * RJIX       +
     &                  QXY * RJIX * RJIY * TWO +
     &                  QXZ * RJIX * RJIZ * TWO +
     &                  QYY * RJIY * RJIY       +
     &                  QYZ * RJIY * RJIZ * TWO +
     &                  QZZ * RJIZ * RJIZ
C
               SUM1 = OXXX * RJIX * RJIX * RJIX         +
     &                OXXY * RJIX * RJIX * RJIY * THREE +
     &                OXYY * RJIX * RJIY * RJIY * THREE +
     &                OYYY * RJIY * RJIY * RJIY         +
     &                OYYZ * RJIY * RJIY * RJIZ * THREE +
     &                OYZZ * RJIY * RJIZ * RJIZ * THREE +
     &                OZZZ * RJIZ * RJIZ * RJIZ         +
     &                OXYZ * RJIX * RJIY * RJIZ * SIX   +
     &                OXXZ * RJIX * RJIX * RJIZ * THREE +
     &                OXZZ * RJIX * RJIZ * RJIZ * THREE
C
                  VOX  = OXXX * RJIX * RJIX       + 
     &                   OXXY * RJIX * RJIY * TWO +
     &                   OXXZ * RJIX * RJIZ * TWO +
     &                   OXYY * RJIY * RJIY       +
     &                   OXYZ * RJIY * RJIZ * TWO +
     &                   OXZZ * RJIZ * RJIZ
C                                         
                  VOY  = OXXY * RJIX * RJIX       + 
     &                   OXYY * RJIX * RJIY * TWO +
     &                   OXYZ * RJIX * RJIZ * TWO +
     &                   OYYY * RJIY * RJIY       +
     &                   OYYZ * RJIY * RJIZ * TWO +
     &                   OYZZ * RJIZ * RJIZ
C                                         
                  VOZ  = OXXZ * RJIX * RJIX       + 
     &                   OXYZ * RJIX * RJIY * TWO +
     &                   OXZZ * RJIX * RJIZ * TWO +
     &                   OYYZ * RJIY * RJIY       +
     &                   OYZZ * RJIY * RJIZ * TWO +
     &                   OZZZ * RJIZ * RJIZ
C                                  
               VQX  = QXX * RJIX + QXY * RJIY + QXZ * RJIZ
               VQY  = QXY * RJIX + QYY * RJIY + QYZ * RJIZ
               VQZ  = QXZ * RJIX + QYZ * RJIY + QZZ * RJIZ
C
C              ITERATE OVER NORMAL MODES
C
               DO M=1,NMODES
                  NM0  = NNNH3*(M-1) +    I
                  NMX3 = NNN  *(M-1) + 3*(I-1) + 1
                  NMX6 = NNN6 *(M-1) + 6*(I-1) + 1
                  NMX10= NNN10*(M-1) +10*(I-1) + 1
                  GRF  = GIVEC(M) / TWO
C
C                 UNPACK THE MASS-WEIGHTED EIGENVECTORS
C
                  RLMX = LVEC(NMX3  )
                  RLMY = LVEC(NMX3+1)
                  RLMZ = LVEC(NMX3+2)
C
C                 UNPACK THE DERIVATIVES OF DMTP
C
                  CHG1I = CHG1(NM0)
C
                  DIPX1 = DIP1(NMX3  )
                  DIPY1 = DIP1(NMX3+1)
                  DIPZ1 = DIP1(NMX3+2)
C              
                  Q1XX= QAD1(NMX6  )
                  Q1YY= QAD1(NMX6+1)
                  Q1ZZ= QAD1(NMX6+2)
                  Q1XY= QAD1(NMX6+3)
                  Q1XZ= QAD1(NMX6+4)
                  Q1YZ= QAD1(NMX6+5)
C              
                  O1XXX=OCT1(NMX10  )
                  O1YYY=OCT1(NMX10+1)
                  O1ZZZ=OCT1(NMX10+2)
C              
                  O1XXY=OCT1(NMX10+3)
                  O1XXZ=OCT1(NMX10+4)
                  O1XYY=OCT1(NMX10+5)
C              
                  O1YYZ=OCT1(NMX10+6)
                  O1XZZ=OCT1(NMX10+7)
                  O1YZZ=OCT1(NMX10+8)
                  O1XYZ=OCT1(NMX10+9)
C
C                 AUXILIARY MODE-DEPENDENT DOT PRODUCTS
C
                  RJIRI = RJIX * RLMX + 
     &                    RJIY * RLMY +
     &                    RJIZ * RLMZ
C
                  RJID1 = RJIX * DIPX1 +
     &                    RJIY * DIPY1 +
     &                    RJIZ * DIPZ1
C
                  DJRI  = DIPX * RLMX +
     &                    DIPY * RLMY +
     &                    DIPZ * RLMZ
C
                  QJRJII = QXX * RLMX * RJIX +
     &                     QXY * RLMX * RJIY +
     &                     QXY * RLMY * RJIX + 
     &                     QXZ * RLMX * RJIZ +
     &                     QXZ * RLMZ * RJIX +
     &                     QYY * RLMY * RJIY +
     &                     QYZ * RLMY * RJIZ +
     &                     QYZ * RLMZ * RJIY +
     &                     QZZ * RLMZ * RJIZ
C
                  Q1RJI2 = Q1XX * RJIX * RJIX       +
     &                     Q1XY * RJIX * RJIY * TWO +
     &                     Q1XZ * RJIX * RJIZ * TWO +
     &                     Q1YY * RJIY * RJIY       +
     &                     Q1YZ * RJIY * RJIZ * TWO +
     &                     Q1ZZ * RJIZ * RJIZ
C
                  VQX1 = QXX * RLMX + QXY * RLMY + QXZ * RLMZ
                  VQY1 = QXY * RLMX + QYY * RLMY + QYZ * RLMZ
                  VQZ1 = QXZ * RLMX + QYZ * RLMY + QZZ * RLMZ
C
                  Q1RJIX = Q1XX * RJIX + Q1XY * RJIY + Q1XZ * RJIZ
                  Q1RJIY = Q1XY * RJIX + Q1YY * RJIY + Q1YZ * RJIZ
                  Q1RJIZ = Q1XZ * RJIX + Q1YZ * RJIY + Q1ZZ * RJIZ
C
               SUM2 = OXXX * RLMX * RJIX * RJIX         +
     &                OXXY * RLMX * RJIX * RJIY * TWO   +
     &                OXXY * RLMY * RJIX * RJIX         +
     &                OXYY * RLMX * RJIY * RJIY         +
     &                OXYY * RLMY * RJIX * RJIY * TWO   +
     &                OYYY * RLMY * RJIY * RJIY         +
     &                OYYZ * RLMY * RJIY * RJIZ * TWO   +
     &                OYYZ * RLMZ * RJIY * RJIY         +
     &                OYZZ * RLMY * RJIZ * RJIZ         +
     &                OYZZ * RLMZ * RJIY * RJIZ * TWO   +
     &                OZZZ * RLMZ * RJIZ * RJIZ         +
     &                OXYZ * RLMX * RJIY * RJIZ * TWO   +
     &                OXYZ * RLMY * RJIX * RJIZ * TWO   +
     &                OXYZ * RLMZ * RJIX * RJIY * TWO   +
     &                OXXZ * RLMX * RJIX * RJIZ * TWO   +
     &                OXXZ * RLMZ * RJIX * RJIX         +
     &                OXZZ * RLMX * RJIZ * RJIZ         +
     &                OXZZ * RLMZ * RJIX * RJIZ
C
               SUM3 = O1XXX * RJIX * RJIX * RJIX         +
     &                O1XXY * RJIX * RJIX * RJIY * THREE +
     &                O1XYY * RJIX * RJIY * RJIY * THREE +
     &                O1YYY * RJIY * RJIY * RJIY         +
     &                O1YYZ * RJIY * RJIY * RJIZ * THREE +
     &                O1YZZ * RJIY * RJIZ * RJIZ * THREE +
     &                O1ZZZ * RJIZ * RJIZ * RJIZ         +
     &                O1XYZ * RJIX * RJIY * RJIZ * SIX   +
     &                O1XXZ * RJIX * RJIX * RJIZ * THREE +
     &                O1XZZ * RJIX * RJIZ * RJIZ * THREE
C
                  VOX1 = OXXX * RJIX * RLMX       + 
     &                   OXXY * RJIX * RLMY       +
     &                   OXXY * RJIY * RLMX       +
     &                   OXXZ * RJIX * RLMZ       +
     &                   OXXZ * RJIZ * RLMX       +
     &                   OXYY * RJIY * RLMY       +
     &                   OXYZ * RJIY * RLMZ       +
     &                   OXYZ * RJIZ * RLMY       +
     &                   OXZZ * RJIZ * RLMZ
C                                  
                  VOY1 = OXXY * RJIX * RLMX       + 
     &                   OXYY * RJIX * RLMY       +
     &                   OXYY * RJIY * RLMX       +
     &                   OXYZ * RJIX * RLMZ       +
     &                   OXYZ * RJIZ * RLMX       +
     &                   OYYY * RJIY * RLMY       +
     &                   OYYZ * RJIY * RLMZ       +
     &                   OYYZ * RJIZ * RLMY       +
     &                   OYZZ * RJIZ * RLMZ
C                                  
                  VOZ1 = OXXZ * RJIX * RLMX       + 
     &                   OXYZ * RJIX * RLMY       +
     &                   OXYZ * RJIY * RLMX       +
     &                   OXZZ * RJIX * RLMZ       +
     &                   OXZZ * RJIZ * RLMX       +
     &                   OYYZ * RJIY * RLMY       +
     &                   OYZZ * RJIY * RLMZ       +
     &                   OYZZ * RJIZ * RLMY       +
     &                   OZZZ * RJIZ * RLMZ
C
                  V1OX  =O1XXX * RJIX * RJIX       + 
     &                   O1XXY * RJIX * RJIY * TWO +
     &                   O1XXZ * RJIX * RJIZ * TWO +
     &                   O1XYY * RJIY * RJIY       +
     &                   O1XYZ * RJIY * RJIZ * TWO +
     &                   O1XZZ * RJIZ * RJIZ
C                                         
                  V1OY  =O1XXY * RJIX * RJIX       + 
     &                   O1XYY * RJIX * RJIY * TWO +
     &                   O1XYZ * RJIX * RJIZ * TWO +
     &                   O1YYY * RJIY * RJIY       +
     &                   O1YYZ * RJIY * RJIZ * TWO +
     &                   O1YZZ * RJIZ * RJIZ
C                                         
                  V1OZ  =O1XXZ * RJIX * RJIX       + 
     &                   O1XYZ * RJIX * RJIY * TWO +
     &                   O1XZZ * RJIX * RJIZ * TWO +
     &                   O1YYZ * RJIY * RJIY       +
     &                   O1YZZ * RJIY * RJIZ * TWO +
     &                   O1ZZZ * RJIZ * RJIZ
C
                  RRR  = THREE * RJIRI * RJI2
C
C                 CHARGE CONTRIBUTION
C
                  CHGMI = CHGI * RJI3
                  CHGM1 = CHG1I* RJI3 
                  FX = CHGMI * ( RLMX - RRR * RJIX ) - CHGM1 * RJIX
                  FY = CHGMI * ( RLMY - RRR * RJIY ) - CHGM1 * RJIY
                  FZ = CHGMI * ( RLMZ - RRR * RJIZ ) - CHGM1 * RJIZ
C
C                 DIPOLE CONTRIBUTION
C
                  FRDI = FIVE * RJI2 * DJRJI * RJIRI
C
                  FX = FX + THREE * RJI5 * ( DJRJI * RLMX +
     &                 DJRI * RJIX - FRDI * RJIX + DIPX * RJIRI )
C
                  FY = FY + THREE * RJI5 * ( DJRJI * RLMY +
     &                 DJRI * RJIY - FRDI * RJIY + DIPY * RJIRI )
C
                  FZ = FZ + THREE * RJI5 * ( DJRJI * RLMZ +
     &                 DJRI * RJIZ - FRDI * RJIZ + DIPZ * RJIRI )
C
                  FX = FX - THREE * RJI5 * RJID1 * RJIX + RJI3 * DIPX1
                  FY = FY - THREE * RJI5 * RJID1 * RJIY + RJI3 * DIPY1
                  FZ = FZ - THREE * RJI5 * RJID1 * RJIZ + RJI3 * DIPZ1
C
C                 QUADRUPOLE CONTRIBUTION
C
                  FVV = FIVE * RJI2
C
                  FX = FX + RJI5 * ( FVV * ( TWO   * RJIRI  * VQX + 
     &                                       TWO   * QJRJII * RJIX +
     &                                      QJRJI2 * RLMX ) - TWO * VQX1
     &                      - 35.00D+00 * RJIRI * QJRJI2 * RJI4 * RJIX )
C
                  FY = FY + RJI5 * ( FVV * ( TWO   * RJIRI  * VQY + 
     &                                       TWO   * QJRJII * RJIY +
     &                                      QJRJI2 * RLMY ) - TWO * VQY1
     &                      - 35.00D+00 * RJIRI * QJRJI2 * RJI4 * RJIY )
C
                  FZ = FZ + RJI5 * ( FVV * ( TWO   * RJIRI  * VQZ + 
     &                                       TWO   * QJRJII * RJIZ +
     &                                      QJRJI2 * RLMZ ) - TWO * VQZ1
     &                      - 35.00D+00 * RJIRI * QJRJI2 * RJI4 * RJIZ )
C
              FX = FX + TWO * RJI5 * Q1RJIX - FVV * RJI5 * Q1RJI2 * RJIX
              FY = FY + TWO * RJI5 * Q1RJIY - FVV * RJI5 * Q1RJI2 * RJIY
              FZ = FZ + TWO * RJI5 * Q1RJIZ - FVV * RJI5 * Q1RJI2 * RJIZ
C
C                 OCTUPOLE CONTRIBUTION
C
                  FX = FX + SEVEN*RJI9 * (SUM1*RLMX - 
     &                      QNINE*RJI2*SUM1*RJIRI*RJIX + 
     &                      THREE*SUM2*RJIX +
     &                      THREE*RJIRI*VOX -
     &                      SUM3*RJIX ) -
     &                 SIX*RJI7*VOX1 + THREE*RJI7*V1OX
C
                  FY = FY + SEVEN*RJI9 * (SUM1*RLMY - 
     &                      QNINE*RJI2*SUM1*RJIRI*RJIY + 
     &                      THREE*SUM2*RJIY +
     &                      THREE*RJIRI*VOY -
     &                      SUM3*RJIY ) -
     &                 SIX*RJI7*VOY1 + THREE*RJI7*V1OY
C
                  FZ = FZ + SEVEN*RJI9 * (SUM1*RLMZ - 
     &                      QNINE*RJI2*SUM1*RJIRI*RJIZ + 
     &                      THREE*SUM2*RJIZ +
     &                      THREE*RJIRI*VOZ -
     &                      SUM3*RJIZ ) -
     &                 SIX*RJI7*VOZ1 + THREE*RJI7*V1OZ
C
C                 WEIGHT EACH CONTRIBUTION BY MODE COEFFICIENT GRF
C
                  FLDX1 = FLDX1 - FX * GRF
                  FLDY1 = FLDY1 - FY * GRF
                  FLDZ1 = FLDZ1 - FZ * GRF
C
               ENDDO
            ENDDO
C
C           SAVE FIELD DERIVATIVES
C
            FIVEC(NJX3) = FLDX1
            FIVEC(NJY3) = FLDY1
            FIVEC(NJZ3) = FLDZ1
C
         ENDDO 
      ENDDO
C
C     ACCUMULATE -AVEC-
C
      CALL DGEMM('N','N',NDIM,NDIM,NDIM,ONE,
     &                   DIMAT,NDIM,DMAT,NDIM,
     &                   ZERO,MAT2,NDIM)
      CALL DGEMM('N','N',NDIM,NDIM,NDIM,ONE,
     &                   DMAT,NDIM,MAT2,NDIM,
     &                   ZERO,MAT3,NDIM)
      CALL DGMV('T',MAT3,FLDS,AVEC,NDIM)
C
      CALL DGMV('T',DMAT,FIVEC,VEC1,NDIM)
      CALL DGMV('N',DMAT,FIVEC,SDIPND,NDIM)
C
       DO IK=1,NDIM
          AVEC(IK) = AVEC(IK) - (VEC1(IK) + SDIPND(IK))
       ENDDO
C
c          CALL VECWRT(VEC1,NDIM,-1,"vec1.dat")
c          CALL VECWRT(SDIPND,NDIM,-1,"sdipnd.dat")
          CALL VECWRT(FIVEC,NDIM,-1,"fivec.dat")
c          CALL VECWRT(AVEC,NDIM,-1,"avec.dat")
c          CALL MATWRT(MAT1,NDIM,NDIM,-1,"rnew.dat")
C
C
19191 CONTINUE
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE SOLPOL(RDMA,CHG,DIP,QAD,OCT,RPOL,POL,EPOL,
     *                  DMAT,FLDS,DIPIND,
     *                  NMOLS,NDMA,NPOL,NDIM,NDMAS,
     *                  MDIP,MQAD,MOCT,MRPOL,MPOL,LWRITE)
C
C -----------------------------------------------------------------------------
C
C         ELECTROSTATIC POLARIZATION ENERGY FROM MULTIPOLE EXPANSION
C                  DISTRIBUTED DIPOLE POLARIZABILITY MODEL
C 
C              Bartosz Błasiak                        05.11.2013
C
C -----------------------------------------------------------------------------
C
C   Input variables:
C
C     ** Integer scalars:
C     NMOLS   - number of molecules
C     NDMAS   - total number of distributed sites
C     NPOLS   - total number of distributed polarizabilities (not used here for now)
C     NDIM    - dimension of D-matrix (NMOLS*NPOLS*3)
C     MDIP    - length of DIP vector = NDMAS * 3
C     MQAD    - length of QAD vector = NDMAS * 6
C     MOCT    - length of OCT vector = NDMAS * 10
C     MRPOL   - length of RPOL vector = NPOLS * 3
C     MPOL    - length of DPOL vector = NPOLS * 9
C
C     ** Integer arrays:
C     NDMA    - array of numbers of distributed electrostatic sites
C     NPOL    - array of numbers of distributed induction sites
C
C     ** Double precision scalar:
C     EPOL    - interaction energy
C
C     ** Double precision arrays:
C     RDMA    - array of DMA positions for each molecule
C     RPOL    - array of POL positions for each molecule
C     CHG,DIP,
C     QAD,OCT - distributed multipoles for each molecule
C     POL     - distributed polarizabilities for each molecule
C     
C   Returns:
C     EPOL    - polarization energy
C
C   External:
C     ILAENV  - evaluate optimum block size (lapack)
C     DGETRF  - LU decomposition (lapack) 
C     DGETRI  - rectanqular real matrix inversion (lapack)
C     DDOT    - dot product of two real vectors (lapack)
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
C            1   2   3   4   5   6   7   8   9   10
C     POL(i) XX  XY  XZ  YX  YY  YZ  ZX  ZY  ZZ 
C            1   2   3   4   5   6   7   8   9
C -----------------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RDMA(MDIP),CHG(NDMAS),DIP(MDIP),QAD(MQAD),OCT(MOCT),
     &          RPOL(MRPOL),POL(MPOL),DMAT(NDIM,NDIM),FLDS(NDIM),
     &          DIPIND(NDIM),NDMA(NMOLS),NPOL(NMOLS),
     &          IPIV(10000),WORK(5000000)
      PARAMETER (ZERO=0.0D+00,HALF=0.50D+00)
      DOUBLE PRECISION DDOT
      LOGICAL LWRITE
      EXTERNAL ILAENV,DGETRF,DGETRI,DDOT
Cf2py INTENT(OUT) EPOL
C     INTENT(IN,OUT) DMAT
      EPOL = ZERO
      DATA IPIV/10000*0/
      DATA WORK/5000000*0.0D+00/
C
C     CALCULATE FIELDS AT POLARIZABLE CENTERS AND D-MATRIX
C
      CALL FIELDS(RDMA,CHG,DIP,QAD,OCT,RPOL,POL,DMAT,FLDS,
     ^            NMOLS,NPOL,NDMA,NDIM,NDMAS,
     ^            MDIP,MQAD,MOCT,MRPOL,MPOL)
C
C     CALCULATE INDUCED DIPOLES AND POLARIZATION ENERGY
C
      NB = ILAENV(1,"DGETRI","NALIWKU",NDIM,-1,-1,-1)
      LWORK = NDIM*NB
      CALL DGETRF(NDIM,NDIM,DMAT,NDIM,IPIV,INFO)
C
      IF (LWRITE) WRITE(*,*) " LWORK= ", LWORK, "INFO= ",INFO
C
      CALL DGETRI(NDIM,DMAT,NDIM,IPIV,WORK,LWORK,INFO)
         IF (INFO.GT.0) THEN
            WRITE(*,*) "DMAT IS SINGULAR! QUITTING..."
            GOTO 1123
         ELSEIF (INFO.LT.0) THEN
            WRITE(*,*) "INDICES WRING! QUITTING..."
            GOTO 1123
         ENDIF
      CALL DGMV('N',DMAT,FLDS,DIPIND,NDIM)
      EPOL = - DDOT(NDIM,FLDS,1,DIPIND,1) * HALF
C
      IF (LWRITE) THEN
          CALL VECWRT(DIPIND,NDIM,-1,"dipind.dat")
          CALL VECWRT(FLDS,NDIM,-1,"fields.dat")
          CALL MATWRT(DMAT,NDIM,NDIM,-1,"dmat.dat")
      ENDIF
C
 1123 CONTINUE
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE FIELDS(RDMA,CHG,DIP,QAD,OCT,RPOL,POL,DMAT,FLDS,
     *                  NMOLS,NPOL,NDMA,NDIM,NDMAS,
     *                  MDIP,MQAD,MOCT,MRPOL,MPOL)
C
C     EVALUATE FIELDS DUE TO EFP FRAGMENTS. INITIAL FLDS VALUES ARE ZERO
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RDMA(MDIP),CHG(NDMAS),DIP(MDIP),
     &          QAD(MQAD),OCT(MOCT),
     &          RPOL(MRPOL),POL(MPOL),DMAT(NDIM,NDIM),FLDS(NDIM),
     &          NPOL(NMOLS),NDMA(NMOLS),
     &          WORKI(30),APOL(3,3),IPIVP(3)
      COMMON/SUMS  / VSUM1(3),VSUM2(3)
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,THREE=3.0D+00,
     &           FOUR=4.0D+00,FIVE=5.0D+00,SIX=6.0D+00,HALF=0.50D+00)
C
C     ITERATE OVER ALL MOLECULES AND THEIR POLARIZABLE CENTERS
C
      DATA WORKI/30*0.0D+00/
      DATA APOL /9*0.0D+00/ 
      DATA IPIVP/3*0/
C
      NPOLI = 0
      DO 9 IMOL=1,NMOLS
      NIM  = NPOL(IMOL)
      NPOLI = NPOLI + NIM
      DO 8 I=1,NIM
         NIX3 = 3*(NPOLI-NIM) + 3*(I-1) + 1
         NIX9 = 9*(NPOLI-NIM) + 9*(I-1) + 1
C
         NIY3 = NIX3 + 1
         NIZ3 = NIY3 + 1
C
         FVALX = ZERO
         FVALY = ZERO
         FVALZ = ZERO
C
C        ITERATE OVER PAIRS OF MOLECULES
C
         NDMAJ = 0
         DO JMOL=1,NMOLS
            NJM = NDMA(JMOL)
            NDMAJ = NDMAJ + NJM
C
C           CALCULATE FIELDS
C
            IF (JMOL.EQ.IMOL) GOTO 919
C
            DO J=1,NJM
               NJX0 =    (NDMAJ-NJM) +    J
               NJX3 =  3*(NDMAJ-NJM) + 3*(J-1) + 1
               NJX6 =  6*(NDMAJ-NJM) + 6*(J-1) + 1
               NJX10= 10*(NDMAJ-NJM) +10*(J-1) + 1
C                                                              
               NJY3 = NJX3 + 1
               NJZ3 = NJY3 + 1
C                                                              
               RIJX = RPOL(NIX3) - RDMA(NJX3) 
               RIJY = RPOL(NIY3) - RDMA(NJY3) 
               RIJZ = RPOL(NIZ3) - RDMA(NJZ3) 
C                                                              
               RIJ  = DSQRT(RIJX*RIJX+
     &                      RIJY*RIJY+
     &                      RIJZ*RIJZ)
C                                                              
               RIJ3 = ONE/(RIJ**3)
               RIJ5 = RIJ3/(RIJ*RIJ)
               RIJ7 = RIJ5/(RIJ*RIJ)
               RIJ9 = RIJ7/(RIJ*RIJ)
C                                                              
               DX = DIP(NJX3)
               DY = DIP(NJY3)
               DZ = DIP(NJZ3)
C                                                              
               QXX= QAD(NJX6  )
               QYY= QAD(NJX6+1)
               QZZ= QAD(NJX6+2)
               QXY= QAD(NJX6+3)
               QXZ= QAD(NJX6+4)
               QYZ= QAD(NJX6+5)
C                                                              
               OXXX=OCT(NJX10  ) 
               OYYY=OCT(NJX10+1) 
               OZZZ=OCT(NJX10+2) 
C                                                              
               OXXY=OCT(NJX10+3)
               OXXZ=OCT(NJX10+4)
               OXYY=OCT(NJX10+5)
C                                                              
               OYYZ=OCT(NJX10+6)
               OXZZ=OCT(NJX10+7)
               OYZZ=OCT(NJX10+8)
               OXYZ=OCT(NJX10+9)
C                                                              
C              AUXILIARY SUMS
C                                                              
               SUM1 = DX*RIJX+DY*RIJY+DZ*RIJZ
C                                                              
               SUM2 = QXX * RIJX * RIJX       +
     &                QXY * RIJX * RIJY * TWO +
     &                QXZ * RIJX * RIJZ * TWO +
     &                QYY * RIJY * RIJY       +
     &                QYZ * RIJY * RIJZ * TWO +
     &                QZZ * RIJZ * RIJZ
C                                                              
               SUM3 = OXXX * RIJX * RIJX * RIJX         +
     &                OXXY * RIJX * RIJX * RIJY * THREE +
     &                OXYY * RIJX * RIJY * RIJY * THREE +
     &                OYYY * RIJY * RIJY * RIJY         +
     &                OYYZ * RIJY * RIJY * RIJZ * THREE +
     &                OYZZ * RIJY * RIJZ * RIJZ * THREE +
     &                OZZZ * RIJZ * RIJZ * RIJZ         +
     &                OXYZ * RIJX * RIJY * RIJZ * SIX   +
     &                OXXZ * RIJX * RIJX * RIJZ * THREE +
     &                OXZZ * RIJX * RIJZ * RIJZ * THREE
C                                                              
               VSUM1(1) = QXX * RIJX + QXY * RIJY + QXZ * RIJZ
               VSUM1(2) = QXY * RIJX + QYY * RIJY + QYZ * RIJZ
               VSUM1(3) = QXZ * RIJX + QYZ * RIJY + QZZ * RIJZ
C                                                              
               VSUM2(1) = OXXX * RIJX * RIJX       +
     &                    OXXY * RIJX * RIJY * TWO +
     &                    OXXZ * RIJX * RIJZ * TWO +
     &                    OXYY * RIJY * RIJY       +
     &                    OXYZ * RIJY * RIJZ * TWO +
     &                    OXZZ * RIJZ * RIJZ
               VSUM2(2) = OXXY * RIJX * RIJX       +
     &                    OXYY * RIJX * RIJY * TWO +
     &                    OXYZ * RIJX * RIJZ * TWO +
     &                    OYYY * RIJY * RIJY       +
     &                    OYYZ * RIJY * RIJZ * TWO +
     &                    OYZZ * RIJZ * RIJZ
               VSUM2(3) = OXXZ * RIJX * RIJX       +
     &                    OXYZ * RIJX * RIJY * TWO +
     &                    OXZZ * RIJX * RIJZ * TWO +
     &                    OYYZ * RIJY * RIJY       +
     &                    OYZZ * RIJY * RIJZ * TWO +
     &                    OZZZ * RIJZ * RIJZ
C CHARGES   
               QC = CHG(NJX0) * RIJ3
               FVALX = FVALX + QC * RIJX
               FVALY = FVALY + QC * RIJY
               FVALZ = FVALZ + QC * RIJZ
C DIPOLES   
               FS5   = THREE * SUM1 * RIJ5
               FVALX = FVALX + RIJX * FS5 - DX * RIJ3
               FVALY = FVALY + RIJY * FS5 - DY * RIJ3
               FVALZ = FVALZ + RIJZ * FS5 - DZ * RIJ3
C QUADRUPOLES
               FS7   = FIVE * SUM2 * RIJ7
               T5    = TWO * RIJ5
               FVALX = FVALX + RIJX * FS7 - VSUM1(1) * T5
               FVALY = FVALY + RIJY * FS7 - VSUM1(2) * T5
               FVALZ = FVALZ + RIJZ * FS7 - VSUM1(3) * T5
C OCTUPOLES
               FS9 = FIVE * SUM3 * RIJ9
               T7  = THREE* RIJ7
               FVALX = FVALX + RIJX * FS9 - VSUM2(1) * T7
               FVALY = FVALY + RIJY * FS9 - VSUM2(2) * T7
               FVALZ = FVALZ + RIJZ * FS9 - VSUM2(3) * T7
            ENDDO
919         CONTINUE
         ENDDO
C
         FLDS(NIX3) = FVALX
         FLDS(NIY3) = FVALY
         FLDS(NIZ3) = FVALZ
C
C        CALCULATE D-MATRIX DIAGONALS
C
         APOL(1,1) = POL(NIX9  )
         APOL(1,2) = POL(NIX9+1)
         APOL(1,3) = POL(NIX9+2)
         APOL(2,1) = POL(NIX9+3)
         APOL(2,2) = POL(NIX9+4)
         APOL(2,3) = POL(NIX9+5)
         APOL(3,1) = POL(NIX9+6)
         APOL(3,2) = POL(NIX9+7)
         APOL(3,3) = POL(NIX9+8)
C
         CALL DGETRF(3,3,APOL,3,IPIVP,INFO)
         CALL DGETRI(3,APOL,3,IPIVP,WORKI,20,INFO)
         IF (INFO.GT.0) THEN
             WRITE(*,*) "HERE SINGULAR TOO!"
             GOTO 19191
         ELSEIF (INFO.LT.0) THEN
             WRITE(*,*) "HERE INDICES BAD TOO!"
             GOTO 19191
         ENDIF
         DMAT(NIX3,NIX3) = APOL(1,1)
         DMAT(NIY3,NIY3) = APOL(2,2)
         DMAT(NIZ3,NIZ3) = APOL(3,3)
C
         DMAT(NIX3,NIY3) = APOL(1,2)
         DMAT(NIX3,NIZ3) = APOL(1,3)
         DMAT(NIY3,NIZ3) = APOL(2,3)
C
         DMAT(NIY3,NIX3) = APOL(2,1)
         DMAT(NIZ3,NIX3) = APOL(3,1)
         DMAT(NIZ3,NIY3) = APOL(3,2)
C
C        CALCULATE D-MATRIX OFFDIAGONALS
C
         NPOLJ = 0
         DO JMOL=1,NMOLS
         NJM = NPOL(JMOL)
         NPOLJ = NPOLJ + NJM
         IF (IMOL.EQ.JMOL) GOTO 921
         DO J=1,NJM
            NJX3 = 3*(NPOLJ-NJM) + 3*(J-1) + 1
            NJY3 = NJX3 + 1
            NJZ3 = NJY3 + 1
C
            RIJX = RPOL(NIX3) - RPOL(NJX3)
            RIJY = RPOL(NIY3) - RPOL(NJY3)
            RIJZ = RPOL(NIZ3) - RPOL(NJZ3)
C
            RIJ  = DSQRT(RIJX*RIJX+
     &                   RIJY*RIJY+
     &                   RIJZ*RIJZ)
C
            RIJ3 = ONE/(RIJ*RIJ*RIJ)
            RIJ5 = RIJ3/(RIJ*RIJ)
C
            THR5 = THREE * RIJ5
            DMXX = RIJ3 - THR5 * RIJX*RIJX 
            DMXY =      - THR5 * RIJX*RIJY
            DMYY = RIJ3 - THR5 * RIJY*RIJY 
            DMXZ =      - THR5 * RIJX*RIJZ 
            DMYZ =      - THR5 * RIJY*RIJZ 
            DMZZ = RIJ3 - THR5 * RIJZ*RIJZ 
C
            DMAT(NIX3,NJX3) = DMXX
            DMAT(NIX3,NJY3) = DMXY
            DMAT(NIY3,NJX3) = DMXY
            DMAT(NIY3,NJY3) = DMYY
            DMAT(NIX3,NJZ3) = DMXZ
            DMAT(NIZ3,NJX3) = DMXZ
            DMAT(NIY3,NJZ3) = DMYZ
            DMAT(NIZ3,NJY3) = DMYZ
            DMAT(NIZ3,NJZ3) = DMZZ
         ENDDO
 921     CONTINUE
         ENDDO
C
 8       CONTINUE
 9       CONTINUE
19191 CONTINUE
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE MATWRT(A,M,N,L,PLIK)
C
C     WRITE THE SQUARE MATRIX OF DIMENSION M,N ON SCREEN (L>0) OR TO FILE (<0)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(M,N)
      CHARACTER(*) PLIK
C
      IF (L.GT.0) THEN
         GOTO 112
      ELSE
          IF (PLIK.EQ." ") THEN
             OPEN(8,FILE="mat.dat",ACCESS='sequential',FORM='formatted')
          ELSE
             OPEN(8,FILE=PLIK,ACCESS='sequential',FORM='formatted')
          ENDIF
          DO I=1,M
             DO J=1,N
                WRITE(8,*) A(I,J)
             ENDDO
          ENDDO
      ENDIF
C
 112  CONTINUE
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE VECWRT(V,N,L,PLIK)
C
C     WRITE THE VECTOR V OF LENGTH N ON SCREEN (L>0) OR TO FILE (<0)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION V(N)
      CHARACTER(*) PLIK
C
      IF (L.GT.0) THEN
          DO I=1,N/3
             K=3*(I-1)+1
             WRITE(*,*) K, V(K), V(K+1), V(K+2)
          ENDDO
      ELSE
          IF (PLIK.EQ." ") THEN
             OPEN(8,FILE="vec.dat",ACCESS='sequential',FORM='formatted')
          ELSE
             OPEN(8,FILE=PLIK,ACCESS='sequential',FORM='formatted')
          ENDIF
          DO I=1,N/3
             K=3*(I-1)+1
             VX = V(K)
             VY = V(K+1)
             VZ = V(K+2)
             VNORM = DSQRT(VX*VX+VY*VY+VZ*VZ)
             WRITE(8,*) I, VX, VY, VZ, VNORM
          ENDDO
      ENDIF
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE DGMV(TRANS,VM,VIN,VOUT,N)
C
C     MATRIX-VECTOR MULTIPLICATION M*VIN = VOUT
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*1 TRANS
      DIMENSION VM(N,N),VIN(N),VOUT(N)
      PARAMETER (ZERO=0.0D+00)
      IF (TRANS.EQ.'N') THEN
         DO I=1,N
            VAL = ZERO
            DO J=1,N
               VAL = VAL + VM(I,J) * VIN(J)
            ENDDO
            VOUT(I) = VAL        
         ENDDO
      ELSE
         DO I=1,N
            VAL = ZERO
            DO J=1,N
               VAL = VAL + VM(J,I) * VIN(J)
            ENDDO
            VOUT(I) = VAL        
         ENDDO
      ENDIF
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE CHECK(INFO)
C     CHECK THE MATRIX INVERSION
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IF (INFO.GT.0) THEN
         WRITE(*,*) "DMAT IS SINGULAR! QUITTING..."
      ELSEIF (INFO.LT.0) THEN
         WRITE(*,*) "INDICES WRONG! QUITTING..."
      ENDIF
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      BLOCK DATA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/SUMS  / VSUM1(3),VSUM2(3)
      DATA  VSUM1/3*0.0D+00/
      DATA  VSUM2/3*0.0D+00/
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|
