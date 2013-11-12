C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE MOLLST(RC,RM,IC,IP,NAT,CCUT,PCUT,NMOLS,NCOORD,NC)
C
C -----------------------------------------------------------------------------
C
C         DETERMINE THE MOLECULES LYING WITHIN COULOMB AND POLARIZATION RADII
C         FROM CENTRAL MOLECULE
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
C
C     ** Integer
C     IC         - array of condition numbers for Coulomb sphere
C     IP         - array of condition numbers for Polarization sphere
C     NMOLS      - number of other molecules (apart from central one)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RC(NC,3),RM(3*NCOORD),NAT(NMOLS),IC(NMOLS),IP(NMOLS)
      PARAMETER (ZERO=0.0D+00)
Cf2py INTENT(IN,OUT) IC,IP
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
         IF (RR.LT.CCUT) IC(I) = 1
         IF (RR.LT.PCUT) IP(I) = 1
C
99    CONTINUE
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
C              Bartosz Blasiak                        05.11.2013
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
     &          IPIV(10000),WORK(500000)
      PARAMETER (ZERO=0.0D+00,HALF=0.50D+00)
      DOUBLE PRECISION DDOT
      LOGICAL LWRITE
      EXTERNAL ILAENV,DGETRF,DGETRI,DDOT
Cf2py INTENT(OUT) EPOL
C     INTENT(IN,OUT) DMAT
      EPOL = ZERO
      DATA IPIV/10000*0/
      DATA WORK/500000*0.0D+00/
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
      CALL DGMV(DMAT,FLDS,DIPIND,NDIM)
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

      SUBROUTINE MATWRT(A,M,N,L,PLIK)
C
C     WRITE THE MATRIX OF DIMENSION M,N ON SCREEN (L>0) OR TO FILE (<0)
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
          DO I=1,N
             WRITE(8,*) V(I)
          ENDDO
      ENDIF
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE DGMV(VM,VIN,VOUT,N)
C
C     MATRIX-VECTOR MULTIPLICATION M*VIN = VOUT
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION VM(N,N),VIN(N),VOUT(N)
      PARAMETER (ZERO=0.0D+00)
      DO 1111 I=1,N
         VAL = ZERO
         DO 2222 J=1,N
            VAL = VAL + VM(I,J) * VIN(J)
 2222    CONTINUE
         VOUT(I) = VAL        
 1111 CONTINUE
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
         NIX0 =   (NPOLI-NIM) +    I
         NIX3 = 3*(NPOLI-NIM) + 3*(I-1) + 1
         NIX6 = 6*(NPOLI-NIM) + 6*(I-1) + 1
         NIX9 = 9*(NPOLI-NIM) + 9*(I-1) + 1
         NIX10=10*(NPOLI-NIM) +10*(I-1) + 1
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

      BLOCK DATA
C
C     MAX MOL - 100
C     MAX POL SITES PER MOL - 10
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/SUMS  / VSUM1(3),VSUM2(3)
      DATA  VSUM1/3*0.0D+00/
      DATA  VSUM2/3*0.0D+00/
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|
