C
      SUBROUTINE PREPOT
C
C   System:          ClH2
C   Functional form: Isaacson and Muckerman's diatomics-in-molecules version A
C   Common name:     ClH2 DIMA
C   Reference:       unpublished
C   Cross Reference: S. C. Tucker, D. G. Truhlar,
C                    G. C. Garrett, and A. D. Isaacson
C                    J. Chem. Phys. 82, 4102 (1985)
C
C   Calling Sequence: 
C      PREPOT - initializes the potential's variables and
C               must be called once before any calls to POT
C      POT    - driver for the evaluation of the energy and the derivatives 
C               of the energy with respect to the coordinates for a given 
C               geometric configuration
C
C   Units: 
C      energies    - hartrees
C      coordinates - bohr
C      derivatives - hartrees/bohr
C
C   Surfaces: 
C      ground electronic state
C
C   Zero of energy: 
C      The classical potential energy is set equal to zero for the Cl
C      infinitely far from the H2 diatomic and R(H2) set equal to the
C      H2 equilibrium diatomic value.
C
C   Parameters:
C      Set in the BLOCK DATA subprogram PTPACM
C
C   Coordinates:
C      Internal, Definition: R(1) = R(Cl-H)
C                            R(2) = R(H-H)
C                            R(3) = R(Cl-H)
C
C   Common Blocks (used between the calling program and this potential):
C      /PT31CM/ R(3), ENERGY, DEDR(3)
C        passes the coordinates, ground state electronic energy, and 
C        derivatives of the ground electronic state energy with respect 
C        to the coordinates.
C      /PT32CM/ NSURF, NDER, NFLAG(20), IDEBUG
C        passes the control flags where
C        NSURF   - not used
C        NDER   = 0 => no derivatives are computed
C               = 1 => derivatives of the energy for the ground electronic 
C                      state with respect to the coordinates are computed
C        NFLAG   - not used  
C        IDEBUG = 0 => do not print extra debug information
C        IDEBUG = 1 => print intermediate values to the file linked to
C                      FORTRAN unit IPRT
C      /PT34CM/ IPRT
C        passes the FORTRAN unit number used for potential output
C      /PT35CM/ EASYAB, EASYBC, EASYAC
C        passes the energy in the three asymptotic valleys for an A + BC system.
C        The energy in the AB valley, EASYAB, is equal to the energy of the 
C        C atom "infinitely" far from the AB diatomic and R(AB) set equal to 
C        Re(AB), the equilibrium bond length for the AB diatomic.  
C        In this potential the AB valley represents H infinitely far from
C        the ClH diatomic and R(ClH) equal to Re(ClH).  Similarly, the terms
C        EASYBC and EASYAC represent the energies in the H2 and the other ClH 
C        valleys, respectively.
C
C   Default Parameter Values:
C      Variable      Default value
C      NDER             1 
C      IDEBUG           0
C      IPRT             6
C
C*****
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT34CM/ IPRT
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20), IDEBUG
      COMMON /PT35CM/ EASYAB, EASYBC, EASYAC
      DIMENSION H1(2,2),H2(2,2),DH1(2,2),DH2(2,2)
      DIMENSION TH(6,6),VEC(6,6),W(6),EN(6),GVEC(6),DH(6,6),RR(3),DD(3)
      DIMENSION DH11(3),DH12(3),DH13(3),DH14(3),DH15(3),DH16(3),DH22(3)
      DIMENSION DH24(3),DH25(3),DH26(3),DH33(3),DH34(3),DH44(3),DH45(3)
      DIMENSION DH46(3),DH55(3),DH56(3),DH66(3)
      CHARACTER*80 TITLE(2)
      DIMENSION INDX(3)
      COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
      COMMON /POTCM2/ C(3)
      COMMON /POTCM3/ ZERO,ONE,TWO,TRE,FOR,SIX,EHT,HALF,FORTH,RT2I,
     $                SQRT3,DPUX,SSGX,DSGH,SSGH
      EQUIVALENCE (H1(1,1),SS11A),(H1(1,2),SS12A),(H1(2,2),SS22A)
      EQUIVALENCE (H2(1,1),SS11B),(H2(1,2),SS12B),(H2(2,2),SS22B)
      EQUIVALENCE (DH1(1,1),DSS11A),(DH1(1,2),DSS12A),(DH1(2,2),DSS22A)
      EQUIVALENCE (DH2(1,1),DSS11B),(DH2(1,2),DSS12B),(DH2(2,2),DSS22B)
      EQUIVALENCE (TH(1,1),DH(1,1))
      DATA NDIM/6/
      DATA NROOT/1/
      DATA TEST/0.999999999D0/
C
C   Initialize the array INDX which specifies the order of the 
C   input bond distances.  INDX=1,2,3 for Cl + H2 and INDX=132 for H + ClH
C
      DATA INDX / 1, 2, 3/
      DATA TITLE(1) /'Final BNL DIM surface for CL+H2 -- 5 and 7 kcal/mo
     *le barriers'/
      DATA TITLE(2) /'Spline fits replaced with analytic forms where pos
     *sible'/
C
C   Echo the name of the potential to the file linked to FORTRAN unit IPRT
C
       WRITE (IPRT, 600) C
       WRITE (IPRT, 610) TITLE
C
600   FORMAT(/,1X,'*****','Potential Energy Surface',1X,'*****',
     *      //,1X,T5,'ClH2 DIMA potential energy surface',
     *      //,1X,T5,'Parameters:',
     *        /,2X,T5,'Bond', T46, 'Cl-H', T58, 'H-H', T69, 'H-Cl',
     *        /,2X,T5,'Dissociation energies (atomic units):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5)
610   FORMAT(/,2X,T5,'Title card for the potential routine:',
     *       /,2X,T5,80A,/,2X,T5,80A,//,1X,'*****')
C
C     SET UP SPLINE FITS TO H2 CURVES
C
      CALL SETSPL
C
C   Initialize the energy in the three asymptotic valleys
C
      EASYAB = C(1)
      EASYBC = C(2)
      EASYAC = C(3)
C
      RETURN
C
      ENTRY POT
C
C***********************************************************************
C
C     ENTRY POT GIVE ENERGY AND PARTIALS W.R.T 3 DISTANCES
C
C     INTERNUCLEAR DISTANCES STORED IN R(3) AND RR(3)
C     NROOT-TH STATE ENERGY STORED IN ENERGY (OR TERM FOR DIATOMIC)
C     THREE PARTIALS STORED IN DD(3)
C
C***********************************************************************
C
C   Check the value of NDER
C
         IF (NDER .GT. 1) THEN
             WRITE (IPRT, 900) NDER
             STOP 'POT 1'
         ENDIF
C
C     IAPP=1
C
C     SET FLAG FOR EIGENVECTOR OUTPUT
C
      IFLAG=SIGN(1,NROOT)
      NROOT=ABS(NROOT)
C
C     FOLLOWING COORDINATE DEFINITIONS PROVIDE FOR PROPER SYSTEM
C
      DO 15 IX=1,3
   15 RR(IX)=R(INDX(IX))
      IF (IDEBUG .EQ. 1) WRITE(IPRT,6600) Q,R
6600  FORMAT(1X,'DIMIMA, Q=',1P3E13.5,/,9X,' R=',3E13.5)
C
C  CHECK FOR NEGATIVE R'S
C
      ENERGY = 1.D30                                                    23AUG83
      DEDR(1) = -ENERGY                                                 23AUG83
      DEDR(2) = -ENERGY                                                 23AUG83
      DEDR(3) = -ENERGY                                                 23AUG83
      IF(RR(1).LT.0.D0 .OR. RR(2).LT.0.D0 .OR. RR(3).LT.0.D0) RETURN    23AUG83
C
C     ZERO OUT UNUSED HAMILTONIAN ELEMENTS
C
      TH(2,3)=ZERO
      TH(3,5)=ZERO
      TH(3,6)=ZERO
C
C     COMPUTE SINE AND COSINE OF ANGLE AXB, AND POWERS
C
      CT=(RR(1)*RR(1)+RR(3)*RR(3)-RR(2)*RR(2))/(TWO*RR(1)*RR(3))
      IF(CT.GT.-TEST .AND. CT.LT.TEST) GO TO 10                         07AUG83
      CT = SIGN(TEST,CT)                                                07AUG83
   10 CONTINUE
      ARG=ONE-CT*CT
      ST=SQRT(ARG)
      CT2=CT*CT
      ST2=ST*ST
      CS=CT*ST
C
C     COMPUTE DIATOMIC ENERGIES FOR RR(3) ARRAY
C
      CALL HCL(RR(1),H1,TSXA,SPXA,TPXA,DH1,DTSXA,DSPXA,DTPXA)
      CALL HCL(RR(3),H2,TSXB,SPXB,TPXB,DH2,DTSXB,DSPXB,DTPXB)
      CALL HCLM(RR(1),DSXA,DDSXA)
      CALL HCLM(RR(3),DSXB,DDSXB)
      CALL H2CURV(RR(2),SSAB,TSAB,DGAB,DUAB,DSSAB,DTSAB,DDGAB,DDUAB)
      IF (IDEBUG .EQ. 1) WRITE(IPRT,6602) SS11A,SS11B,DSXA,DSXB,SSAB
6602  FORMAT(1X,'SS11A,SS11B,DSXA,DSXB,SSAB=',1P5E13.5)
C
C     CONSTRUCT DIM HAMILTONIAN IN ORTHOGONAL BASIS
C
      TH(1,1)=SS11A+FORTH*((SSAB+TRE*TSAB)+CT2*(SS11B+TRE*TSXB)
     1        +ST2*(SPXB+TRE*TPXB))-DPUX-TWO*DSGH
      TH(1,2)=FORTH*SQRT3*(SSAB-TSAB-CT2*(SS11B-TSXB)-ST2*(SPXB-TPXB))
      TH(1,3)=SS12A
      TH(1,4)=-HALF*CT*SS12B
      TH(1,5)=FORTH*CS*(SS11B+TRE*TSXB-SPXB-TRE*TPXB)
      TH(1,6)=-FORTH*SQRT3*CS*(SS11B-TSXB-SPXB+TPXB)
      TH(2,2)=TSXA+FORTH*((TRE*SSAB+TSAB)+CT2*(TRE*SS11B+TSXB)
     1        +ST2*(TRE*SPXB+TPXB))-DPUX-TWO*DSGH
      TH(2,4)=HALF*SQRT3*CT*SS12B
      TH(2,5)=TH(1,6)
      TH(2,6)=FORTH*CS*(TRE*SS11B+TSXB-TRE*SPXB-TPXB)
      TH(3,3)=SS22A+HALF*(DGAB+DUAB)+DSXB-SSGX-SSGH-DSGH
      TH(3,4)=HALF*(DGAB-DUAB)
      TH(4,4)=DSXA+HALF*(DGAB+DUAB)+SS22B-SSGX-SSGH-DSGH
      TH(4,5)=-HALF*ST*SS12B
      TH(4,6)=HALF*SQRT3*ST*SS12B
      TH(5,5)=SPXA+FORTH*((SSAB+TRE*TSAB)+ST2*(SS11B+TRE*TSXB)
     1        +CT2*(SPXB+TRE*TPXB))-DPUX-TWO*DSGH
      TH(5,6)=FORTH*SQRT3*(SSAB-TSAB-ST2*(SS11B-TSXB)-CT2*(SPXB-TPXB))
      TH(6,6)=TPXA+FORTH*((TRE*SSAB+TSAB)+ST2*(TRE*SS11B+TSXB)
     1        +CT2*(TRE*SPXB+TPXB))-DPUX-TWO*DSGH
C
C     FILL OUT MATRIX AND DIAGONALIZE
C
      DO 30 I=1,NDIM
      DO 30 J=1,I
      TH(I,J)=TH(J,I)
   30 CONTINUE
C
      CALL EIGN(NDIM,TH,VEC,EN,W,NDIM)
C
C     FOLLOWING STATEMENT ASSUMES ENERGY ZERO AT A + BC(R=RE)
C
      ENERGY=EN(NROOT)+C(2)
      IF (IDEBUG .EQ. 1) 
     *    WRITE(IPRT,6605) NROOT,EN(1),EN(NROOT),C(2),ENERGY
6605  FORMAT(1X,'NROOT,EN(1),EN(NROOT),C(2),ENERGY=',I5,1P4E13.5)
      IF(IFLAG.LT.ZERO) WRITE(IPRT,998) (VEC(I,NROOT),I=1,NDIM)
  998 FORMAT(1X,'VECTOR=',6D20.8)
C
C************************* END OF FPOT3N *******************************
C
C     IF(IAPP.EQ.0) RETURN
      IF(NDER.NE.1) RETURN
C
C     SET  VECTOR FOR STATE OF INTEREST
C
      DO 40 I=1,NDIM
      GVEC(I)=VEC(I,NROOT)
   40 CONTINUE
C
C     COMPUTE PARTIALS OF U=COS(THETA)
C
      DTDR1=(RR(1)*RR(1)+RR(2)*RR(2)-RR(3)*RR(3))/
     +      (TWO*RR(1)*RR(1)*RR(3))
      DTDR2=-RR(2)/(RR(1)*RR(3))
      DTDR3=(RR(2)*RR(2)+RR(3)*RR(3)-RR(1)*RR(1))/
     +      (TWO*RR(3)*RR(3)*RR(1))
      FACTR=-ONE/ST
      IF(ABS(ST).LT.1.D-6) WRITE(IPRT,905) ST
905   FORMAT(/,2X,T5,'Warning: Sine(THETA) = ',1PE20.10, 
     *               ' in the derivative calculation',
     *       /,2X,T14,'a serious error may result')
C
C     COMPUTE PARTIALS OF HAMILTONIAN ELEMENTS
C
      XTRA=+HALF*CT*(SS11B+TRE*TSXB-SPXB-TRE*TPXB)
      DH11(1)=DSS11A+XTRA*DTDR1
      DH11(2)=FORTH*(DSSAB+TRE*DTSAB)+XTRA*DTDR2
      DH11(3)=FORTH*(CT2*(DSS11B+TRE*DTSXB)+ST2*(DSPXB+TRE*DTPXB)) +
     1        XTRA*DTDR3
      XTRA=-HALF*SQRT3*CT*(SS11B-TSXB-SPXB+TPXB)
      DH12(1)=XTRA*DTDR1
      DH12(2)=FORTH*SQRT3*(DSSAB-DTSAB)+XTRA*DTDR2
      DH12(3)=-FORTH*SQRT3*(CT2*(DSS11B-DTSXB)+ST2*(DSPXB-DTPXB)) +
     1         XTRA*DTDR3
      DH13(1)=DSS12A
      DH13(2)=ZERO
      DH13(3)=ZERO
      XTRA=-HALF*SS12B
      DH14(1)=XTRA*DTDR1
      DH14(2)=XTRA*DTDR2
      DH14(3)=-HALF*CT*DSS12B+XTRA*DTDR3
      XTRA=FORTH*(CT2-ST2)*(SS11B+TRE*TSXB-SPXB-TRE*TPXB)*FACTR
      DH15(1)=XTRA*DTDR1
      DH15(2)=XTRA*DTDR2
      DH15(3)=FORTH*CS*(DSS11B+TRE*DTSXB-DSPXB-TRE*DTPXB)+XTRA*DTDR3
      XTRA=-FORTH*SQRT3*(CT2-ST2)*(SS11B-TSXB-SPXB+TPXB)*FACTR
      DH16(1)=XTRA*DTDR1
      DH16(2)=XTRA*DTDR2
      DH16(3)=-FORTH*SQRT3*CS*(DSS11B-DTSXB-DSPXB+DTPXB)+XTRA*DTDR3
      XTRA=+HALF*CT*(TRE*SS11B+TSXB-TRE*SPXB-TPXB)
      DH22(1)=DTSXA+XTRA*DTDR1
      DH22(2)=FORTH*(TRE*DSSAB+DTSAB)+XTRA*DTDR2
      DH22(3)=FORTH*(CT2*(TRE*DSS11B+DTSXB)+ST2*(TRE*DSPXB+DTPXB)) +
     1        XTRA*DTDR3
      XTRA=+HALF*SQRT3*SS12B
      DH24(1)=XTRA*DTDR1
      DH24(2)=XTRA*DTDR2
      DH24(3)=HALF*SQRT3*CT*DSS12B+XTRA*DTDR3
      DH25(1)=DH16(1)
      DH25(2)=DH16(2)
      DH25(3)=DH16(3)
      XTRA=FORTH*(CT2-ST2)*(TRE*SS11B+TSXB-TRE*SPXB-TPXB)*FACTR
      DH26(1)=XTRA*DTDR1
      DH26(2)=XTRA*DTDR2
      DH26(3)=FORTH*CS*(TRE*DSS11B+DTSXB-TRE*DSPXB-DTPXB)+XTRA*DTDR3
      DH33(1)=DSS22A
      DH33(2)=HALF*(DDGAB+DDUAB)
      DH33(3)=DDSXB
      DH34(1)=ZERO
      DH34(2)=HALF*(DDGAB-DDUAB)
      DH34(3)=ZERO
      DH44(1)=DDSXA
      DH44(2)=HALF*(DDGAB+DDUAB)
      DH44(3)=DSS22B
      XTRA=-HALF*CT*SS12B*FACTR
      DH45(1)=XTRA*DTDR1
      DH45(2)=XTRA*DTDR2
      DH45(3)=-HALF*ST*DSS12B+XTRA*DTDR3
      XTRA=HALF*SQRT3*CT*SS12B*FACTR
      DH46(1)=XTRA*DTDR1
      DH46(2)=XTRA*DTDR2
      DH46(3)=HALF*SQRT3*ST*DSS12B+XTRA*DTDR3
      XTRA=-HALF*CT*(SS11B+TRE*TSXB-SPXB-TRE*TPXB)
      DH55(1)=DSPXA+XTRA*DTDR1
      DH55(2)=FORTH*(DSSAB+TRE*DTSAB)+XTRA*DTDR2
      DH55(3)=FORTH*(ST2*(DSS11B+TRE*DTSXB)+CT2*(DSPXB+TRE*DTPXB)) +
     1        XTRA*DTDR3
      XTRA=+HALF*SQRT3*CT*(SS11B-TSXB-SPXB+TPXB)
      DH56(1)=XTRA*DTDR1
      DH56(2)=FORTH*SQRT3*(DSSAB-DTSAB)+XTRA*DTDR2
      DH56(3)=-FORTH*SQRT3*(ST2*(DSS11B-DTSXB)+CT2*(DSPXB-DTPXB)) +
     1         XTRA*DTDR3
      XTRA=-HALF*CT*(TRE*SS11B+TSXB-TRE*SPXB-TPXB)
      DH66(1)=DTPXA+XTRA*DTDR1
      DH66(2)=FORTH*(TRE*DSSAB+DTSAB)+XTRA*DTDR2
      DH66(3)=FORTH*(ST2*(TRE*DSS11B+DTSXB)+CT2*(TRE*DSPXB+DTPXB)) +
     1        XTRA*DTDR3
C
C     LOOP OVER RR(I) DERIVATIVES
C     FIRST ZERO OUT UNUSED ELEMENTS
C
      DH(2,3)=ZERO
      DH(3,5)=ZERO
      DH(3,6)=ZERO
      DO 50 I=1,3
C
C     SET UP DERIVATIVE MATRIX
C
      DH(1,1)=DH11(I)
      DH(1,2)=DH12(I)
      DH(1,3)=DH13(I)
      DH(1,4)=DH14(I)
      DH(1,5)=DH15(I)
      DH(1,6)=DH16(I)
      DH(2,2)=DH22(I)
      DH(2,4)=DH24(I)
      DH(2,5)=DH25(I)
      DH(2,6)=DH26(I)
      DH(3,3)=DH33(I)
      DH(3,4)=DH34(I)
      DH(4,4)=DH44(I)
      DH(4,5)=DH45(I)
      DH(4,6)=DH46(I)
      DH(5,5)=DH55(I)
      DH(5,6)=DH56(I)
      DH(6,6)=DH66(I)
C
C     FILL OUT DERIVATIVE MATRIX
C
      DO 70 J=1,NDIM
      DO 70 K=1,J
      DH(J,K)=DH(K,J)
   70 CONTINUE
C
      DNUMB=ZERO
      DO 80 J=1,NDIM
      DO 80 K=1,NDIM
      DNUMB=DNUMB+GVEC(J)*DH(J,K)*GVEC(K)
   80 CONTINUE
      DD(I)=DNUMB
C
C     END OF DERIVATIVE LOOP
C
   50 CONTINUE
C
C     PUT DERIVATIVES INTO CORRESPONDING SPOTS FOR CORRECT COORD CHOICE.
C
      DO 90 IX=1,3
   90 DEDR(INDX(IX))=DD(IX)
C
C************************ END OF DFPOT *********************************
C
900   FORMAT(/,1X,T5,'Error: POT has been called with NDER = ', I5,
     *       /,1X,T12,'only the first derivatives, NDER = 1, are ',
     *                'coded in this potential')
C
      RETURN
      END
C
      SUBROUTINE SETSPL
C
C     SETS UP SPLINE FITS FOR H2 CURVES
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /POTCM4/ V1H2(79),V2H2(87),V3H2(64),V3DUM(62),V4H2(64),
     $                V4DUM(62),R1H2(79),R2H2(87),R3H2(126)
      COMMON /POTCM3/ ZERO,ONE,TWO,TRE,FOR,SIX,EHT,HALF,QUAR,RT2I,RT3,
     $                DPUX,SSGX,DSGH,SSGH
      COMMON  /POTCM5/  S2SS(79),S3SS(79),DELSS(79),S2TS(87),S3TS(87),
     $                  DELTS(87),S2DG(126),S3DG(126),DELDG(126),
     $                  S2DU(126),S3DU(126),DELDU(126)
C
C     NOW COMPUTE SPLINE FIT DATA
C     H2 SINGLET STATE
C
      DO 10 I=1,67
   10 R1H2(I)=DBLE(I+3)*0.1D0
      R1H2(68)=7.2D-00
      R1H2(69)=7.4D-00
      R1H2(70)=7.6D-00
      R1H2(71)=7.8D-00
      R1H2(72)=8.0D-00
      R1H2(73)=8.25D-00
      R1H2(74)=8.5D-00
      R1H2(75)=9.0D-00
      R1H2(76)=9.5D-00
      R1H2(77)=10.D-00
      R1H2(78)=11.D-00
      R1H2(79)=12.D-00
      DO 30 I=1,79
 30   V1H2(I)=V1H2(I)+ONE
      S2SS(1)=3.2006875D+1
      S2SS(79)=-7.93478D-7                                              07AUG83
      CALL SPLNIN(79,R1H2,V1H2,T,SS,SS1,SS2,EPSLN,S2SS,S3SS,DELSS)
C
C     H2 TRIPLET STATE
C
      DO 70 I=1,81
 70   R2H2(I)=DBLE(I+9)*0.1D-00
      R2H2(82)=9.25D-00
      R2H2(83)=9.5D-00
      R2H2(84)=9.75D-00
      R2H2(85)=10.D-00
      R2H2(86)=11.D-00
      R2H2(87)=12.D-00
      DO 90 I=1,87
 90   V2H2(I)=V2H2(I)+ONE
      S2TS(1)=1.3039390D0
      S2TS(87)=-7.93478D-7                                              07AUG83
      CALL SPLNIN(87,R2H2,V2H2,T,SS,SS1,SS2,EPSLN,S2TS,S3TS,DELTS)
C
C     H2+ G STATE
C
      DO 130 I=1,100
 130  R3H2(I)=DBLE(I)*0.1D-00
      DO 140 I=1,26
 140  R3H2(I+100)=DBLE(I)*HALF+10.D-00
      DO 150 I=1,126
 150  V3H2(I)=V3H2(I)+ONE
      S2DG(1)=1.0742890D+3
      S2DG(126)=-3.03938D-7                                             07AUG83
      CALL SPLNIN(126,R3H2,V3H2,T,SS,SS1,SS2,EPSLN,S2DG,S3DG,DELDG)
C
C     H2+ U STATE
C
      DO 170 I=1,126
 170  V4H2(I)=V4H2(I)+ONE
      S2DU(1)=9.0876397D+2
      S2DU(126)=-3.0398D-7                                              07AUG83
      CALL SPLNIN(126,R3H2,V4H2,T,SS,SS1,SS2,EPSLN,S2DU,S3DU,DELDU)
      RETURN
      END
C
      SUBROUTINE HCLM(R,XHM2S,DXHM2S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     DOUBLET SIGMA STATE OF HCL-IN A.U.
C
      COMMON /POTCM3/ ZERO,ONE,TWO,TRE,FOR,SIX,EHT,HALF,QUAR,RT2I,RT3,
     $                DPUX,SSGX,DSGH,SSGH
      DATA A,B,C/2.52966D+1,-2.23D+0,-1.32774D-1/
      DATA A1,B1,C1/2.4913064D+4,7.1248778D+0,2.8830121/
      DATA RSMALL/0.94486D+0/
      IF(R.GE.RSMALL) GO TO 10
C
C     SMALL R FIT USED TO BE CONSISTENT WITH BNL SPLINE FIT VERSION
C
      ARG=A1*EXP(-B1*R)
      XHM2S=ARG+C1
      DXHM2S=-B1*ARG
      GO TO 20
   10 DX=ONE/R**7
      X=DX*R
      DY=X*R
      Y=DY*R
C
C     COMPUTE ENERGY AND DERIVATIVE
C
      XHM2S=A*X+B*Y+C
      DXHM2S=-SIX*A*DX-FOR*B*DY
   20 RETURN
      END
C
      SUBROUTINE HCL(R,H,XH3S,XH1P,XH3P,DH,DXH3S,DXH1P,DXH3P)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     CALCULATES 2X2 SINGLET SIGMA BLOCK, AND TRIPLET SIGMA, SINGLET
C     PI, AND TRIPLET PI CURVES FOR HCL IN A.U.
C
      DIMENSION H(2,2),DH(2,2)
      COMMON /POTCM3/ ZERO,ONE,TWO,TRE,FOR,SIX,EHT,HALF,QUAR,RT2I,RT3,
     $                DPUX,SSGX,DSGH,SSGH
      COMMON /PT34CM/ IPRT
C
      DATA DX,RX,BX/1.69552D-1,2.4082D+0,9.89104D-1/
      DATA DV,RV,BV,CV/1.84605D-1,4.746D+0,2.77985D-1,1.82621D-1/
      DATA C8,C6,C4,C0/-2.018460D+4,2.399677D+3,-5.351085D+1,3.67226D-1/
      DATA AA,BB,CC/25.20626238D-0,1.910495693D-0,1.828090476D-1/
C
C     TRIPLET SIGMA DE ALREADY CONTAINS SATO=0.06
C
      DATA RT,BT,DT/2.53D-0,1.047D-0,6.685927358D-2/
      DATA RA,BA,DA/2.274408D+0,1.069975D+0,5.3556D-2/
      DATA RMAX,GMAX,GB/2.5D+0,2.0212D-2,4.654D-1/
      DATA RV1/3.5D-0/
      DATA ANGCON/1.0D-0/
C
C     ANGCON = +1 FOR #2 PERP TO #1, = -1 FOR #1 PERP TO #2,
C        AND = 0 FOR SYMMETRIC ORTHOGONALIZATION
C
      DATA ZETA/2.0387D0/
      DATA Z1,Z2,Z3/7.222580119D+2,-2.749954467D+2,8.648352007D+2/
      DATA Z4,Z5,Z6/2.437196153D+2,3.144371653D+1,1.128252224D-1/
      DATA R0,B0,D0/2.53D-0,0.9897D-0,0.150789D-0/
C
C     SINGLET SIGMA BLOCK
C     GROUND STATE ADIABAT
C
      EX=EXP(-BX*(R-RX))
      OMEX=ONE-EX
      E1=DX*OMEX*OMEX-DX
      DE1=TWO*DX*BX*EX*OMEX
C
C     EXCITED STATE ADIABAT
C
      IF(R.GT.RV) GO TO 20
      IF(R.LT.RV1) GO TO 60
C
C     MORSE FIT
C
      EX=EXP(-BV*(R-RV))
      OMEX=ONE-EX
      E2=DV*OMEX*OMEX+CV
      DE2=TWO*DV*BV*EX*OMEX
      GO TO 30
C
C     LONG RANGE FIT
C
   20 RI=ONE/R
      RI2=RI*RI
      RI4=RI2*RI2
      RI6=RI4*RI2
      RI8=RI4*RI4
      E2=C8*RI8+C6*RI6+C4*RI4-RI+C0
      DE2=-RI*(EHT*C8*RI8+SIX*C6*RI6+FOR*C4*RI4-RI)
      GO TO 30
C
C     SHORT RANGE FIT
C
   60 ARG=AA*EXP(-BB*R)
      E2=ARG+CC
      DE2=-BB*ARG
   30 CONTINUE
C
C     LOWER COVALENT DIABAT
C
      EX=EXP(-B0*(R-R0))
      OMEX=ONE-EX
      H(1,1)=D0*OMEX*OMEX-D0
      DH(1,1)=TWO*D0*B0*EX*OMEX
C
C     OVERLAP
C
      R2=R*R
      T1=Z1*(ONE+R)+Z2*R2
      T2=Z1*(ONE+ZETA*R)+Z3*R2+Z4*R*R2+Z5*R2*R2
      EX1=EXP(-R)
      EX2=EXP(-ZETA*R)
      DT1=Z1+TWO*Z2*R
      DT2=Z1*ZETA+TWO*Z3*R+TRE*Z4*R2+FOR*Z5*R2*R
      FAC1=-T1*EX1+T2*EX2
      DFAC1=(-DT1+T1)*EX1+(DT2-ZETA*T2)*EX2
      OVLP=Z6*FAC1/R2
      DOVLP=Z6*(DFAC1-TWO*FAC1/R)/R2
      ROPS2=ONE/SQRT(ONE+OVLP*OVLP)
      S12=OVLP*ROPS2/RT2I
      IF(ABS(S12) .GT. 1.0D0) WRITE(IPRT,6600) R,S12 
 6600 FORMAT(1X,'In the DIM potential in the subprogram HCL: ',
     *          'R,S12=',F10.5,1PE15.7)
      DS12=DOVLP*ROPS2**3/RT2I
C
C     H12 AND H22 DIABATS IN NON-ORTHOGONAL REPRESENTATION
C
      E12=E1+E2
      DE12=DE1+DE2
      FAC1=E12*H(1,1)-E1*E2-H(1,1)*H(1,1)
      FAC2=ONE-S12*S12
      DFAC2=-TWO*S12*DS12
      TTT = FAC1*FAC2                                                   07AUG83
      ROOT= 0.D0                                                        07AUG83
      IF(TTT.GT.0.D0) ROOT=SQRT(TTT)                                    07AGU83
      H(1,2)=S12*H(1,1)+ROOT
      DH(1,2)=DS12*H(1,1)+S12*DH(1,1)                                   07AUG83
      TTT=HALF*(DFAC2*FAC1+FAC2*(DE12*H(1,1)                            07AUG83
     $+E12*DH(1,1)-E1*DE2-DE1*E2-TWO*H(1,1)*DH(1,1)))                   07AUG83
      IF(ROOT.NE.0.D0) DH(1,2) = DH(1,2) + TTT/ROOT                     07AUG83
      H(2,2)=E12*FAC2+TWO*S12*H(1,2)-H(1,1)
      DH(2,2)=DE12*FAC2+E12*DFAC2+TWO*(DS12*H(1,2)+S12*DH(1,2))
     $-DH(1,1)
C
C     TRANSFORM TO SYMMETRIC ORTHOGONAL REPRESENTATION
C
      T1=ONE/SQRT(ONE+S12)
      T2=ONE/SQRT(ONE-S12)
      DT1=-HALF*DS12*T1**3
      DT2=HALF*DS12*T2**3
      ALP=HALF*(T1+T2)
      BET=HALF*(T1-T2)
      DALP=HALF*(DT1+DT2)
      DBET=HALF*(DT1-DT2)
      A2=ALP*ALP
      B2=BET*BET
      AB=ALP*BET
      DA2=TWO*ALP*DALP
      DB2=TWO*BET*DBET
      DAB=ALP*DBET+DALP*BET
      G11=A2*H(1,1)+TWO*AB*H(1,2)+B2*H(2,2)
      G12=(A2+B2)*H(1,2)+AB*(H(1,1)+H(2,2))
      G22=B2*H(1,1)+TWO*AB*H(1,2)+A2*H(2,2)
      ARG=TWO*(AB*DH(1,2)+DAB*H(1,2))
      DG11=A2*DH(1,1)+DA2*H(1,1)+ARG+B2*DH(2,2)+DB2*H(2,2)
      DG12=(A2+B2)*DH(1,2)+(DA2+DB2)*H(1,2)+AB*(DH(1,1)+DH(2,2))
     $+DAB*(H(1,1)+H(2,2))
      DG22=B2*DH(1,1)+DB2*H(1,1)+ARG+A2*DH(2,2)+DA2*H(2,2)
      H(1,1)=G11
      H(1,2)=G12
      H(2,2)=G22
      DH(1,1)=DG11
      DH(1,2)=DG12
      DH(2,2)=DG22
C
C     PERFORM UNITARY TRANSFORMATION ON H
C
      Y=HALF*(ONE/T1+ONE/T2)
      IF(ABS(Y).GT.1.0D0) Y = SIGN(1.0D0,Y)                             07AUG83
      ANG=ANGCON*ACOS(Y)                                                05JUL83
      DANG=ZERO
      IF(ABS(ANG).LT.1.D-14) GO TO 90
      DANG=HALF*ANGCON*(DT1/(T1*T1)+DT2/(T2*T2))/SQRT(ONE-Y*Y)
      C=COS(ANG)
      S=SIN(ANG)
      DC=-S*DANG
      DS=C*DANG
      C2=C*C
      S2=S*S
      CS=C*S
      DC2=TWO*C*DC
      DS2=TWO*S*DS
      DCS=C*DS+S*DC
      H11=C2*H(1,1)+TWO*CS*H(1,2)+S2*H(2,2)
      DH11=C2*DH(1,1)+TWO*CS*DH(1,2)+S2*DH(2,2)
     $+DC2*H(1,1)+TWO*DCS*H(1,2)+DS2*H(2,2)
      H12=(C2-S2)*H(1,2)+CS*(H(2,2)-H(1,1))
      DH12=(C2-S2)*DH(1,2)+CS*(DH(2,2)-DH(1,1))
     $+(DC2-DS2)*H(1,2)+DCS*(H(2,2)-H(1,1))
      H22=S2*H(1,1)-TWO*CS*H(1,2)+C2*H(2,2)
      DH22=S2*DH(1,1)-TWO*CS*DH(1,2)+C2*DH(2,2)
     $+DS2*H(1,1)-TWO*DCS*H(1,2)+DC2*H(2,2)
      H(1,1)=H11
      H(2,2)=H22
      H(1,2)=H12
      H(2,1)=H12
      DH(1,1)=DH11
      DH(2,2)=DH22
      DH(1,2)=DH12
   90 DH(2,1)=DH12
C
C     TRIPLET SIGMA CURVE   (T STATE)
C
      EX=EXP(-BT*(R-RT))
      OPEX=ONE+EX
      XH3S=DT*OPEX*OPEX-DT
      DXH3S=-TWO*DT*BT*EX*OPEX
C
C     SINGLET PI CURVE   (A STATE)
C
      EX=EXP(-BA*(R-RA))
      OPEX=ONE+EX
      XH1P=DA*OPEX*OPEX-DA
      DXH1P=-TWO*DA*BA*EX*OPEX
C
C     TRIPLET PI CURVE   (SMALL A STATE)
C
      G=GMAX
      DG=ZERO
      IF(R.LE.RMAX)GO TO 100
      G=GMAX*EXP(-GB*(R-RMAX)**2)
      DG=-TWO*GB*(R-RMAX)*G
  100 XH3P=XH1P-G
      DXH3P=DXH1P-DG
      RETURN
      END
C
      SUBROUTINE H2CURV(R,AB1S,AB3S,AB2G,AB2U,DAB1S,DAB3S,DAB2G,DAB2U)
C
C     CALCULATES ENERGY AND DERIVATIVES FOR H2 SINGLET AND TRIPLET
C     SIGMA AND H2+ DOUBLET SIGMA G AND U STATES IN A.U.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /POTCM3/ ZERO,ONE,TWO,TRE,FOR,SIX,EHT,HALF,QUAR,RT2I,RT3,
     $                DPUX,SSGX,DSGH,SSGH
      COMMON /POTCM4/ V1H2(79),V2H2(87),V3H2(126),V4H2(126),
     $                R1H2(79),R2H2(87),R3H2(126)
      COMMON  /POTCM5/  S2SS(79),S3SS(79),DELSS(79),S2TS(87),S3TS(87),
     $                  DELTS(87),S2DG(126),S3DG(126),DELDG(126),
     $                  S2DU(126),S3DU(126),DELDU(126)
C
C     H2 EXPONENTIALS W.R.T. FINAL LIMITS, H2+ W.R.T. 0.0 LIMIT
C
      DATA B1,A1/6.0315697D-00,9.8214033D-00/
C**   DATA A2,B2,C2/2.41893D00,2.1329115D00,0.0918539D00/
      DATA C6,C8,C10/6.49903D-00,124.399D-00,1135.21D-00/
      DATA DE,BE,RE/3.21533D-2,8.3D-1,1.60D-0/
      DATA C4/2.25D0/
C
C     H2 CURVES
C     SINGLET SIGMA
C
      IF(R.GT.12.D0) GO TO 20
      IF(R.GE.0.4D0) GO TO 10
C
C     SHORT-RANGE EXPONENTIAL FIT
C
      AB1S=A1*EXP(-B1*R)
      DAB1S=-B1*AB1S
      GO TO 30
C
C     SPLINE FIT
C
   10 CALL SPLINE(79,R1H2,V1H2,R,AB1S,DAB1S,SS2,EPSLN,S2SS,S3SS,DELSS)
      GO TO 30
C
C     LONG-RANGE FORM
C
   20 R2=ONE/(R*R)
      R6=R2**3
      R8=R6*R2
      R10=R8*R2
      AB1S=-(C6*R6+C8*R8+C10*R10)
      DAB1S=(SIX*C6*R6+EHT*C8*R8+10.D0*C10*R10)/R
   30 CONTINUE
C
C     TRIPLET SIGMA
C     LINES STARTING WITH C** ARE FOR CORRECT TRIPLET CURVE
C
C**   IF(R.GT.12.D0) GO TO 50
C**   IF(R.GT.ONE) GO TO 40
C
C     SHORT-RANGE EXPONENTIAL FIT
C
C**   AB3S=A2*EXP(-B2*R)+C2
C**   DAB3S=-B2*(AB3S-C2)
C**   GO TO 60
C
C     SPLINE FIT
C
C**40 CALL SPLINE(87,R2H2,V2H2,R,AB3S,DAB3S,SS2,EPSLN,S2TS,S3TS,DELTS)
C**   GO TO 60
C
C     LONG-RANGE FORM -- SAME AS FOR SINGLET STATE
C
C**   AB3S=AB1S
C**   DAB3S=DAB1S
C
C     REPLACE TRIPLET CURVE WITH MODIFIED ANTI-MORSE
C
      BFAC=2.6D-2
      IF(R.LT.RE) BFAC=2.026D0
      ARG=R-RE
      EX=EXP(-BE*ARG)
      OPEX=ONE+EX
      F1=DE*OPEX*OPEX-DE
      DF1=-TWO*DE*BE*EX*OPEX
      ARG2=ARG*ARG
      ARG3=ARG2*ARG
      F2=EXP(-BFAC*ARG3)
      AB3S=F1*F2
      DAB3S=(DF1-F1*TRE*BFAC*ARG2)*F2
C
C     H2+ CURVES
C     DOUBLET G
C
   60 IF(R.GT.23.D0) GO TO 70
C
C     SPLINE FIT
C
      CALL SPLINE(126,R3H2,V3H2,R,AB2G,DAB2G,SS2,EPSLN,S2DG,S3DG,DELDG)
      GO TO 80
C
C     LONG-RANGE FORM
C
   70 R4=ONE/R**4
      AB2G=-C4*R4+HALF
      DAB2G=9.D0*R4/R
C
C     DOUBLET U
C
   80 IF(R.GT.23.D0) GO TO 90
C
C     SPLINE FIT
C
      CALL SPLINE(126,R3H2,V4H2,R,AB2U,DAB2U,SS2,EPSLN,S2DU,S3DU,DELDU)
      GO TO 100
C
C     LONG-RANGE FORM -- SAME AS FOR DOUBLET G STATE
C
   90 AB2U=AB2G
      DAB2U=DAB2G
  100 RETURN
      END
C
      SUBROUTINE SPLNIN(N,X,Y,T,SS,SS1,SS2,EPSLN,S2,S3,DELY)
C
C     JIM STINE'S INITIALIZATION ROUTINE FOR SPLINE FITTING
C     N=NO. POINTS (X,Y); T=ARGUMENT; SS=VALUE; SS1=1ST DERIV;
C     SS2=2ND DERIV; S2,S3, AND DELY ARE COEFFICIENT ARRAYS
C     SECOND DERIV'S AT ENDS STORED IN S2(1) AND S2(N) ON INPUT
C     EPSLN IS OVERALL FIT CRITERION
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT34CM/ IPRT
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20), IDEBUG
      DIMENSION C(126),H(126),H2(126),B(126)
      DIMENSION X(N),Y(N),S2(N),S3(N),DELY(N)
      DATA EX/1.1D0/,OMEGA/1.0717968D0/,HALF/0.5D0/
      EPSLN=1.D-12
      N1=N-1
      KNT=0
      DO 51 I=1,N1
      H(I)=X(I+1)-X(I)
   51 DELY(I)=(Y(I+1)-Y(I))/H(I)
      DO 52 I=2,N1
      H2(I)=H(I-1)+H(I)
      B(I)=HALF*H(I-1)/H2(I)
      DELSQY=(DELY(I)-DELY(I-1))/H2(I)
      S2(I)=DELSQY+DELSQY
   52 C(I)=DELSQY+S2(I)
    5 ETA=0.0D0
      DO 10 I=2,N1
      W=(C(I)-B(I)*S2(I-1)-(HALF-B(I))*S2(I+1)-S2(I))*OMEGA
      IF(ABS(W)-ETA) 10,10,9
    9 ETA=ABS(W)
   10 S2(I)=S2(I)+W
      KNT=KNT+1
      IF (KNT.GT.10) EPSLN=EPSLN*EX
      IF(ETA-EPSLN) 14,5,5
   14 DO 53 I=1,N1
   53 S3(I)=(S2(I+1)-S2(I))/H(I)
      IF (IDEBUG .EQ. 1) WRITE(IPRT,61) EPSLN
   61 FORMAT(5X,11HFinal EPS =,D15.6)
      RETURN
      END
C
      SUBROUTINE SPLINE(N,X,Y,T,SS,SS1,SS2,EPSLN,S2,S3,DELY)
C
C     JIM STINE'S SPLINE INTERPOLATOR
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT34CM/ IPRT
      DIMENSION X(N),Y(N),S2(N),S3(N),DELY(N)
      I=1
      IF(T.GT.X(1)) GO TO 10
      SS = Y(1)
      SS1 = DELY(1)
      WRITE(IPRT,600) T,X(1)
      RETURN
   10 IF(T.LT.X(N)) GO TO 20
      SS = Y(N)
      SS1 = DELY(N)
      WRITE(IPRT,600) T,X(1)
      RETURN
   20 CONTINUE
   56 IF(T-X(I)) 60,17,57
   57 I=I+1
      GO TO 56
   59 I=N
   60 I=I-1
   17 HT1=T-X(I)
      HT2=T-X(I+1)
      PROD=HT1*HT2
      SS2=S2(I)+HT1*S3(I)
      DELSQS=(S2(I)+S2(I+1)+SS2)/6.D-00
      SS=Y(I)+HT1*DELY(I)+PROD*DELSQS
      SS1=DELY(I)+(HT1+HT2)*DELSQS+PROD*S3(I)/6.D-00
   61 RETURN
  600 FORMAT(1X,'In the DIMIMA SPLINE routine ARG is out of range:',
     1       /,1X,' ARG =',1PE13.5,', X=',1PE13.5)
      END
C
      SUBROUTINE EIGN(NN,A,VEC,EIG,W,ND)
C
C     MATRIX DIAGONALIZATION ROUTINE FOR REAL SYMMETRIC CASE             VBDIM
C     HOUSEHOLDER METHOD                                                 VBDIM
C     RHO=UPPERLIMIT FOR OFF-DIAGONAL ELEMENT                            VBDIM
C     NN=SIZE OF MATRIX                                                  VBDIM
C     A=MATRIX (ONLY LOWER TRIANGLE IS USED ,THIS IS DESTROYED)          VBDIM
C     EIG=RETURNED EIGENVALUES IN ALGEBRAIC ASCENDING ORDER              VBDIM
C     VEC=RETURNED EIGENVECTORS IN COLUMNS                               VBDIM
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EIG(ND),W(ND)                                            VBDIM
      DIMENSION GAMMA(6),BETA(6),BETASQ(6)                               VBDIM
      DIMENSION P(6),Q(6)                                                VBDIM
      EQUIVALENCE (P(1),BETA(1)),(Q(1),BETA(1))
      DIMENSION IPOSV(6),IVPOS(6),IORD(6)                                VBDIM
      DIMENSION A(ND,ND),VEC(ND,ND)                                      VBDIM
      RHO=1.0D-14                                                        ED1
      RHOSQ=RHO*RHO                                                      VBDIM
      N=NN                                                               VBDIM
      IF(N.EQ.0) GO TO 560                                               VBDIM
    1 N1=N-1                                                             VBDIM
      N2=N-2                                                             VBDIM
      GAMMA(1)=A(1,1)                                                    VBDIM
      IF(N2) 280,270,120                                                 VBDIM
  120 DO 260 NR=1,N2                                                     VBDIM
      B=A(NR+1,NR)                                                       VBDIM
      S=0.0D0                                                            VBDIM
      DO 130 I=NR,N2                                                     VBDIM
C
C     PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION                      VBDIM
C
  130 S=S+A(I+2,NR)**2                                                   VBDIM
      A(NR+1,NR)=0.0D0                                                   VBDIM
      IF(S) 250,250,140                                                  VBDIM
  140 S=S+B*B                                                            VBDIM
      SGN=+1.0D0                                                         VBDIM
      IF(B) 150,160,160                                                  VBDIM
  150 SGN=-1.0D0                                                         VBDIM
  160 SQRTS=SQRT(S)                                                      VBDIM
      D=SGN/(SQRTS+SQRTS)                                                VBDIM
      TEMP=SQRT(.5D0+B*D)                                                VBDIM
      W(NR)=TEMP                                                         VBDIM
      A(NR+1,NR)=TEMP                                                    VBDIM
      D=D/TEMP                                                           VBDIM
      B=-SGN*SQRTS                                                       VBDIM
C
C     D IS FACTOR OF PROPORTIONALITY  NOW COMPUTE AND SAVE W VECTOR      VBDIM
C     EXTRA SINGLY SUBSCRIPTED W VECTOR USED FOR SPEED                   VBDIM
C
      DO 170 I=NR,N2                                                     VBDIM
      TEMP=D*A(I+2,NR)                                                   VBDIM
      W(I+1)=TEMP                                                        VBDIM
  170 A(I+2,NR)=TEMP                                                     VBDIM
C
C     PREMULTIPLY VECTOR W BY MATRIX A TO OBTAIN P VECTOR                VBDIM
C     SIMULTANEOUSLY ACCUMULATE DOT PRODUCT WP,(THE SCALAR K)            VBDIM
C
      WTAW=0.0D0                                                         VBDIM
      DO 220 I=NR,N1                                                     VBDIM
      SUM=0.0D0                                                          VBDIM
      DO 180 J=NR,I                                                      VBDIM
  180 SUM=SUM+A(I+1,J+1)*W(J)                                            VBDIM
      I1=I+1                                                             VBDIM
      IF(N1-I1) 210,190,190                                              VBDIM
  190 DO 200 J=I1,N1                                                     VBDIM
  200 SUM=SUM+A(J+1,I+1)*W(J)                                            VBDIM
  210 P(I)=SUM                                                           VBDIM
  220 WTAW=WTAW+SUM*W(I)                                                 VBDIM
C
C     P VECTOR AND SCALAR K NOW STORED, NEXT COMPUTE Q VECTOR.           VBDIM
C
      DO 230 I=NR,N1                                                     VBDIM
C
C     NOW FORM PAP MATRIX, REQUIRED PART                                 VBDIM
C
  230 Q(I)=P(I)-WTAW*W(I)                                                VBDIM
      DO 240 J=NR,N1                                                     VBDIM
      QJ=Q(J)                                                            VBDIM
      WJ=W(J)                                                            VBDIM
      DO 240 I=J,N1                                                      VBDIM
  240 A(I+1,J+1)=A(I+1,J+1)-2.0D0*(W(I)*QJ+WJ*Q(I))                      VBDIM
  250 BETA(NR)=B                                                         VBDIM
      BETASQ(NR)=B*B                                                     VBDIM
  260 GAMMA(NR+1)=A(NR+1,NR+1)                                           VBDIM
  270 B=A(N,N-1)                                                         VBDIM
      BETA(N-1)=B                                                        VBDIM
      BETASQ(N-1)=B*B                                                    VBDIM
      GAMMA(N)=A(N,N)                                                    VBDIM
  280 BETASQ(N)=0.0D0                                                    VBDIM
C
C     ADJOIN AN IDENTITY MATRIX TO BE POSTMULTIPLIED BY ROTATIONS        VBDIM
C
      DO 290 J=1,N                                                       VBDIM
      DO 290 I=1,N                                                       VBDIM
  290 VEC(I,J)=0.0D0                                                     VBDIM
      DO 300 I=1,N                                                       VBDIM
  300 VEC(I,I)=1.0D0                                                     VBDIM
      M=N                                                                VBDIM
      SUM=0.0D0                                                          VBDIM
      NPAS=1                                                             VBDIM
      GO TO 400                                                          VBDIM
  310 SUM=SUM+SHIFT                                                      VBDIM
      COSA=1.0D0                                                         VBDIM
      G=GAMMA(1)-SHIFT                                                   VBDIM
      PP=G                                                               VBDIM
      PPBS=PP*PP+BETASQ(1)                                               VBDIM
      PPBR=SQRT(PPBS)                                                    VBDIM
      DO 370 J=1,M                                                       VBDIM
      COSAP=COSA                                                         VBDIM
      IF(PPBS.NE.0.0D0) GO TO 320                                        VBDIM
  311 SINA=0.0D0                                                         VBDIM
      SINA2=0.0D0                                                        VBDIM
      COSA=1.0D0                                                         VBDIM
      GO TO 350                                                          VBDIM
  320 SINA=BETA(J)/PPBR                                                  VBDIM
      SINA2=BETASQ(J)/PPBS                                               VBDIM
      COSA=PP/PPBR                                                       VBDIM
C
C     POSTMULTIPLY IDAENTITY BY P-TRANSPOSE MATRIX                       VBDIM
C
      NT=J+NPAS                                                          VBDIM
      IF(NT.LT.N) GO TO 330                                              VBDIM
  321 NT=N                                                               VBDIM
  330 DO 340 I=1,NT                                                      VBDIM
      TEMP=COSA*VEC(I,J)+SINA*VEC(I,J+1)                                 VBDIM
      VEC(I,J+1)=-SINA*VEC(I,J)+COSA*VEC(I,J+1)                          VBDIM
  340 VEC(I,J)=TEMP                                                      VBDIM
  350 DIA=GAMMA(J+1)-SHIFT                                               VBDIM
      U=SINA2*(G+DIA)                                                    VBDIM
      GAMMA(J)=G+U                                                       VBDIM
      G=DIA-U                                                            VBDIM
      PP=DIA*COSA-SINA*COSAP*BETA(J)                                     VBDIM
      IF(J.NE.M) GO TO 360                                               VBDIM
  351 BETA(J)=SINA*PP                                                    VBDIM
      BETASQ(J)=SINA2*PP*PP                                              VBDIM
      GO TO 380                                                          VBDIM
  360 PPBS=PP*PP+BETASQ(J+1)                                             VBDIM
      PPBR=SQRT(PPBS)                                                    VBDIM
      BETA(J)=SINA*PPBR                                                  VBDIM
  370 BETASQ(J)=SINA2*PPBS                                               VBDIM
  380 GAMMA(M+1)=G                                                       VBDIM
C
C     TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT                      VBDIM
C
      NPAS=NPAS+1                                                        VBDIM
      IF(BETASQ(M).GT.RHOSQ) GO TO 410                                   VBDIM
  390 EIG(M+1)=GAMMA(M+1)+SUM                                            VBDIM
  400 BETA(M)=0.0D0                                                      VBDIM
      BETASQ(M)=0.0D0                                                    VBDIM
      M=M-1                                                              VBDIM
      IF(M.EQ.0) GO TO 430                                               VBDIM
  401 IF(BETASQ(M).LE.RHOSQ) GO TO 390                                   VBDIM
C
C     TAKE ROOT OF CORNER 2 BY 2 NEAREST TO LOWER DIAGONAL IN VALUE      VBDIM
C     AS ESTIMATE OF EIGENVALUE TO USE FOR SHIFT                         VBDIM
C
  410 A2=GAMMA(M+1)                                                      VBDIM
      R2=.5D0*A2                                                         VBDIM
      R1=.5D0*GAMMA(M)                                                   VBDIM
      R12=R1+R2                                                          VBDIM
      DIF=R1-R2                                                          VBDIM
      TEMP=SQRT(DIF*DIF+BETASQ(M))                                       VBDIM
      R1=R12+TEMP                                                        VBDIM
      R2=R12-TEMP                                                        VBDIM
      DIF=ABS(A2-R1)-ABS(A2-R2)                                          VBDIM
      IF(DIF.LT.0.0D0) GO TO 420                                         VBDIM
  411 SHIFT=R2                                                           VBDIM
      GO TO 310                                                          VBDIM
  420 SHIFT=R1                                                           VBDIM
      GO TO 310                                                          VBDIM
  430 EIG(1)=GAMMA(1)+SUM                                                VBDIM
C
C     INITIALIZE AUXILARY TABLES REQUIRED FOR REARRANGING THE VECTORS    VBDIM
C
      DO 440 J=1,N                                                       VBDIM
      IPOSV(J)=J                                                         VBDIM
      IVPOS(J)=J                                                         VBDIM
  440 IORD(J)=J                                                          VBDIM
C
C     USE A TRANSPOSITON SORT TO ORDER THE EIGENVALUES                   VBDIM
C
      M=N                                                                VBDIM
      GO TO 470                                                          VBDIM
  450 DO 460 J=1,M                                                       VBDIM
      IF(EIG(J).LE.EIG(J+1)) GO TO 460                                   VBDIM
  451 TEMP=EIG(J)                                                        VBDIM
      EIG(J)=EIG(J+1)                                                    VBDIM
      EIG(J+1)=TEMP                                                      VBDIM
      ITEMP=IORD(J)                                                      VBDIM
      IORD(J)=IORD(J+1)                                                  VBDIM
      IORD(J+1)=ITEMP                                                    VBDIM
  460 CONTINUE                                                           VBDIM
  470 M=M-1                                                              VBDIM
      IF(M.NE.0) GO TO 450                                               VBDIM
  471 IF(N1.EQ.0) GO TO 500                                              VBDIM
  472 DO 490 L=1,N1                                                      VBDIM
      NV=IORD(L)                                                         VBDIM
      NP=IPOSV(NV)                                                       VBDIM
      IF(NP.EQ.L) GO TO 490                                              VBDIM
  473 LV=IVPOS(L)                                                        VBDIM
      IVPOS(NP)=LV                                                       VBDIM
      IPOSV(LV)=NP                                                       VBDIM
      DO 480 I=1,N                                                       VBDIM
      TEMP=VEC(I,L)                                                      VBDIM
      VEC(I,L)=VEC(I,NP)                                                 VBDIM
  480 VEC(I,NP)=TEMP                                                     VBDIM
  490 CONTINUE                                                           VBDIM
C
C     BACK TRANSFORM THE VECTORS OF THE TRIPLE DIAGONAL MATRIX           VBDIM
C
  500 DO 550 NRR=1,N                                                     VBDIM
      K=N1                                                               VBDIM
  510 K=K-1                                                              VBDIM
      IF(K.LE.0) GO TO 550                                               VBDIM
  511 SUM=0.0D0                                                          VBDIM
      DO 520 I=K,N1                                                      VBDIM
  520 SUM=SUM+VEC(I+1,NRR)*A(I+1,K)                                      VBDIM
      SUM=SUM+SUM                                                        VBDIM
      DO 530 I=K,N1                                                      VBDIM
  530 VEC(I+1,NRR)=VEC(I+1,NRR)-SUM*A(I+1,K)                             VBDIM
      GO TO 510                                                          VBDIM
  550 CONTINUE                                                           VBDIM
  560 RETURN                                                             VBDIM
C
C  END OF EIGN
C
      END                                                                VBDIM
C
      BLOCK DATA PTPACM
C
C     CONTAINS FUNCTION VALUES FOR H2 SPLINE FITS AND OTHER CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT34CM/ IPRT
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20), IDEBUG
      COMMON /POTCM2/ C(3)
      COMMON /POTCM4/ V1H2(79),V2H2(87),V3H2(64),V3DUM(62),V4H2(64),
     $                V4DUM(62),R1H2(79),R2H2(87),R3H2(126)
      COMMON /POTCM3/ ZERO,ONE,TWO,TRE,FOR,SIX,EHT,HALF,QUAR,RT2I,RT3,
     $                DPUX,SSGX,DSGH,SSGH
C
C   Initialize the control flags for the potential
C
      DATA IPRT /6/
      DATA NDER /1/
      DATA IDEBUG /0/
      DATA NFLAG /20*0/ 
C
C   Initialize the parameters for the potential
C
      DATA ZERO,ONE,TWO,TRE,FOR,SIX/0.D0,1.D0,2.D0,3.D0,4.D0,6.D0/
      DATA EHT,HALF,QUAR/8.D0,0.5D0,0.25D0/
      DATA RT2I,RT3/7.0710 67811 86548D-1,1.7320 50807 56888D0/
      DATA DPUX,SSGX,DSGH,SSGH/0.D0,-0.132774D0,0.D0,0.5D0/
      DATA C/0.169552D0,0.1744734D0,0.169552D0/
      DATA V1H2/
     1 -.1202028D0, -.5266270D0, -.7696341D0, -.9220261D0,-1.0200556D0,
     2-1.0836422D0,-1.1245385D0,-1.1500562D0,-1.1649342D0,-1.1723459D0,
     3-1.1744744D0,-1.1728537D0,-1.1685799D0,-1.1624570D0,-1.1550670D0,
     4-1.1468496D0,-1.1381312D0,-1.1291562D0,-1.1201233D0,-1.1111725D0,
     5-1.1024127D0,-1.0939273D0,-1.0857810D0,-1.0780164D0,-1.0706700D0,
     6-1.0637641D0,-1.0573118D0,-1.0513185D0,-1.0457832D0,-1.0407003D0,
     7-1.0360578D0,-1.0318402D0,-1.0280272D0,-1.0245978D0,-1.0215297D0,
     8-1.0187967D0,-1.0163689D0,-1.0142247D0,-1.0123371D0,-1.0106810D0,
     9-1.0092303D0,-1.0079682D0,-1.0068703D0,-1.0059178D0,-1.0050923D0,
     1-1.0043782D0,-1.0037626D0,-1.0032309D0,-1.0027740D0,-1.0023800D0,
     2-1.0020423D0,-1.0017521D0,-1.0015030D0,-1.0012899D0,-1.0011069D0,
     3-1.0009498D0,-1.0008150D0,-1.0007002D0,-1.0006030D0,-1.0005162D0,
     4-1.0004466D0,-1.0003864D0,-1.0003328D0,-1.0002906D0,-1.0002466D0,
     5-1.0002154D0,-1.0001889D0,-1.0001434D0,-1.0001086D0,-1.0000868D0,
     6-1.0000682D0,-1.0000528D0,-1.0000404D0,-1.0000314D0,-1.0000185D0,
     7-1.0000121D0,-1.0000091D0,-1.0000045D0,-1.0000025D0/
      DATA V2H2/
     1 -.6215227D0, -.6757268D0, -.7189640D0, -.7544033D0, -.7841501D0,
     2 -.8096095D0, -.8317238D0, -.8511347D0, -.8682913D0, -.8835186D0,
     3 -.8970636D0, -.9091230D0, -.9198593D0, -.9294123D0, -.9379051D0,
     4 -.9454463D0, -.9521346D0, -.9580585D0, -.9632985D0, -.9679271D0,
     5 -.9720104D0, -.9756071D0, -.9787717D0, -.9815517D0, -.9839910D0,
     6 -.9861279D0, -.9879977D0, -.9896312D0, -.9910564D0, -.9922979D0,
     7 -.9933781D0, -.9943164D0, -.9951304D0, -.9958355D0, -.9964456D0,
     8 -.9969715D0, -.9974252D0, -.9978159D0, -.9981500D0, -.9984382D0,
     9 -.9986849D0, -.9988960D0, -.9990763D0, -.9992292D0, -.9993597D0,
     1 -.9994703D0, -.9995635D0, -.9996432D0, -.9997094D0, -.9997655D0,
     2 -.9998125D0, -.9998515D0, -.9998842D0, -.9999112D0, -.9999337D0,
     3 -.9999521D0, -.9999670D0, -.9999795D0, -.9999887D0, -.9999966D0,
     4-1.0000030D0,-1.0000076D0,-1.0000117D0,-1.0000143D0,-1.0000167D0,
     5-1.0000180D0,-1.0000191D0,-1.0000195D0,-1.0000197D0,-1.0000194D0,
     6-1.0000196D0,-1.0000192D0,-1.0000191D0,-1.0000180D0,-1.0000173D0,
     7-1.0000164D0,-1.0000158D0,-1.0000152D0,-1.0000143D0,-1.0000143D0,
     8-1.0000127D0,-1.0000109D0,-1.0000095D0,-1.0000076D0,-1.0000067D0,
     9-1.0000043D0,-1.0000025D0/
      DATA V3H2/
     1 8.0217579D0, 3.0713797D0, 1.4666292D0,  .6992459D0,
     2  .2650120D0, -.0048180D0, -.1826248D0, -.3044801D0, -.3902705D0,
     3 -.4517863D0, -.4964119D0, -.5289745D0, -.5527406D0, -.5699835D0,
     4 -.5823232D0, -.5909372D0, -.5966963D0, -.6002536D0, -.6021058D0,
     5 -.6026342D0, -.6021349D0, -.6008396D0, -.5989309D0, -.5965536D0,
     6 -.5938235D0, -.5908332D0, -.5876573D0, -.5843560D0, -.5809780D0,
     7 -.5775629D0, -.5741424D0, -.5707426D0, -.5673842D0, -.5640840D0,
     8 -.5608555D0, -.5577093D0, -.5546536D0, -.5516947D0, -.5488374D0,
     9 -.5460849D0, -.5434394D0, -.5409022D0, -.5384736D0, -.5361531D0,
     1 -.5339400D0, -.5318328D0, -.5298296D0, -.5279281D0, -.5261259D0,
     2 -.5244203D0, -.5228082D0, -.5212866D0, -.5198521D0, -.5185016D0,
     3 -.5172315D0, -.5160385D0, -.5149192D0, -.5138701D0, -.5128878D0,
     4 -.5119690D0, -.5111105D0, -.5103089D0, -.5095612D0, -.5088644D0/
      DATA V3DUM/
     5 -.5082155D0, -.5076116D0, -.5070501D0, -.5065283D0, -.5060437D0,
     6 -.5055940D0, -.5051769D0, -.5047902D0, -.5044319D0, -.5041001D0,
     7 -.5037930D0, -.5035087D0, -.5032458D0, -.5030027D0, -.5027780D0,
     8 -.5025704D0, -.5023785D0, -.5022014D0, -.5020377D0, -.5018866D0,
     9 -.5017472D0, -.5016184D0, -.5014996D0, -.5013900D0, -.5012888D0,
     1 -.5011954D0, -.5011093D0, -.5010298D0, -.5009564D0, -.5008887D0,
     2 -.5008262D0, -.5007685D0, -.5007153D0, -.5006661D0, -.5006207D0,
     3 -.5005787D0, -.5004122D0, -.5002992D0, -.5002219D0, -.5001683D0,
     4 -.5001306D0, -.5001035D0, -.5000837D0, -.5000690D0, -.5000577D0,
     5 -.5000489D0, -.5000420D0, -.5000364D0, -.5000318D0, -.5000279D0,
     6 -.5000247D0, -.5000220D0, -.5000196D0, -.5000176D0, -.5000158D0,
     7 -.5000142D0, -.5000129D0, -.5000117D0, -.5000106D0, -.5000097D0,
     8 -.5000089D0, -.5000081D0/
      DATA V4H2/
     1 9.4993326D0, 4.4973226D0, 2.8272874D0, 1.9892158D0,
     2 1.4831145D0, 1.1423563D0,  .8955959D0,  .7072541D0,  .5576760D0,
     3  .4351864D0,  .3324666D0,  .2447306D0,  .1687347D0,  .1022057D0,
     4  .0434985D0, -.0086173D0, -.0550906D0, -.0966755D0, -.1339863D0,
     5 -.1675344D0, -.1977530D0, -.2250137D0, -.2496386D0, -.2719091D0,
     6 -.2920720D0, -.3103453D0, -.3269215D0, -.3419716D0, -.3556477D0,
     7 -.3680850D0, -.3794043D0, -.3897134D0, -.3991089D0, -.4076771D0,
     8 -.4154957D0, -.4226343D0, -.4291554D0, -.4351155D0, -.4405653D0,
     9 -.4455506D0, -.4501129D0, -.4542896D0, -.4581145D0, -.4616184D0,
     1 -.4648290D0, -.4677716D0, -.4704692D0, -.4729427D0, -.4752111D0,
     2 -.4772916D0, -.4792001D0, -.4809509D0, -.4825572D0, -.4840309D0,
     3 -.4853831D0, -.4866238D0, -.4877622D0, -.4888066D0, -.4897648D0,
     4 -.4906439D0, -.4914502D0, -.4921898D0, -.4928681D0, -.4934901D0/
      DATA V4DUM/
     5 -.4940603D0, -.4945830D0, -.4950621D0, -.4955011D0, -.4959033D0,
     6 -.4962717D0, -.4966090D0, -.4969178D0, -.4972005D0, -.4974591D0,
     7 -.4976956D0, -.4979119D0, -.4981096D0, -.4982903D0, -.4984553D0,
     8 -.4986060D0, -.4987435D0, -.4988690D0, -.4989834D0, -.4990877D0,
     9 -.4991827D0, -.4992692D0, -.4993480D0, -.4994196D0, -.4994847D0,
     1 -.4995438D0, -.4995975D0, -.4996462D0, -.4996904D0, -.4997303D0,
     2 -.4997665D0, -.4997992D0, -.4998288D0, -.4998554D0, -.4998794D0,
     3 -.4999011D0, -.4999800D0, -.5000244D0, -.5000475D0, -.5000579D0,
     4 -.5000608D0, -.5000595D0, -.5000560D0, -.5000515D0, -.5000467D0,
     5 -.5000421D0, -.5000377D0, -.5000337D0, -.5000301D0, -.5000269D0,
     6 -.5000240D0, -.5000216D0, -.5000194D0, -.5000174D0, -.5000157D0,
     7 -.5000142D0, -.5000129D0, -.5000117D0, -.5000106D0, -.5000097D0,
     8 -.5000089D0, -.5000081D0/
      END

