C
      SUBROUTINE PREPOT
C
C   System:          ClH2
C   Functional form: Isaacson and Muckerman's diatomics-in-molecules version S
C   Common name:     ClH2 DIMS
C   Reference:       unpublished
C   Cross Reference: S. C. Tucker, D. G. Truhlar,
C                    B. C. Garrett, and A. D. Isaacson
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
C      Set in the BLOCK DATA subprogram PTPACM and read in from the file
C      potclh2dims.dat that is linked to FORTRAN unit IRD1
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
C      /PT34CM/ IPRT, IRD1
C        passes the FORTRAN unit number used for potential output
C      /PT5CM/ EASYAB, EASYBC, EASYAC
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
C      IRD1             4
C
C*****
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
      COMMON /PT34CM/ IPRT, IRD1
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20), IDEBUG
      COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
      COMMON /POTCM2/ C(3)
      COMMON/POTCM3/ZERO,ONE,TWO,TRE,HALF,FORTH,SQRT3,
     $              DPUX,SSGX,DSGH,SSGH
     $              SABFIT(121),S2SAB(121),S3SAB(121),DELSAB(121),
     $              TABFIT(121),S2TAB(121),S3TAB(121),DELTAB(121),
     $              GABFIT(121),S2GAB(121),S3GAB(121),DELGAB(121),
     $              UABFIT(121),S2UAB(121),S3UAB(121),DELUAB(121),
     $              S11FIT(121),S2S11(121),S3S11(121),DELS11(121),
     $              S12FIT(121),S2S12(121),S3S12(121),DELS12(121),
     $              S22FIT(121),S2S22(121),S3S22(121),DELS22(121),
     $              TSXFIT(121),S2TSX(121),S3TSX(121),DELTSX(121),
     $              SPXFIT(121),S2SPX(121),S3SPX(121),DELSPX(121),
     $              TPXFIT(121),S2TPX(121),S3TPX(121),DELTPX(121),
     $              DSXFIT(121),S2DSX(121),S3DSX(121),DELDSX(121),NS
      COMMON/POTCM4/X(121)
      DIMENSION TH(6,6),VEC(6,6),W(6),EN(6),GVEC(6),DH(6,6),RR(3),DD(3)
      DIMENSION DH11(3),DH12(3),DH13(3),DH14(3),DH15(3),DH16(3),DH22(3)
      DIMENSION DH24(3),DH25(3),DH26(3),DH33(3),DH34(3),DH44(3),DH45(3)
      DIMENSION DH46(3),DH55(3),DH56(3),DH66(3)
      DIMENSION INDX(3)
      CHARACTER*80 TITLE(2) 
      EQUIVALENCE (TH(1,1),DH(1,1))
C
      DATA NDIM/6/,NROOT/1/
      DATA TEST/0.999999999D0/
      DATA INDX /1, 2, 3/
      DATA TITLE(1) /'    Final BNL DIM surface for Cl+H2 -- 5 and 7 kca
     *l/mole barriers'/
      DATA TITLE(2) /'    Spline fits used in atomic units'/
C
C   Echo the name of the potential to the file linked to FORTRAN unit IPRT
C
      WRITE (IPRT, 600)
      WRITE (IPRT, 610) TITLE
C
C   Open the file with the potential input data
C
      OPEN (UNIT=IRD1, FILE='potclh2dims.dat', FORM='FORMATTED', 
     *      STATUS='OLD', ERR=100)
C
C     SET UP SPLINE FITS TO H2 CURVES
C
      CALL SETSPL
C
C   Close the input data file
C
      CLOSE (UNIT = IRD1)
C
C   Initialize the energies in the three asymptotic valleys
C
      EASYAB = C(1)
      EASYBC = C(2)
      EASYAC = C(3)
C
600   FORMAT(/,1X,'*****','Potential Energy Surface',1X,'*****',
     *      //,1X,T5,'ClH2 DIMS potential energy surface')
610   FORMAT(/,2X,T5,'Title cards for the potential routine:',
     *       /,80A,/,80A,//,1X,'*****')
6000  FORMAT(/,2X,T5,'Error opening input data file')
      RETURN
100   WRITE (IPRT, 6000)
      STOP 'PREPOT 1'
C
      ENTRY POT
C
C***********************************************************************
C
C     ENTRY POT GIVES NROOT-TH POTENTIAL ENERGY FOR THE CLH2 SYSTEM
C     ENTRY POT ALSO GIVES PARTIALS W.R.T. THREE DISTANCES
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
C     SET FLAG FOR EIGENVECTOR OUTPUT
C
      IFLAG=SIGN(1,NROOT)
      NROOT=ABS(NROOT)
C
C     FOLLOWING COORDINATE DEFINITIONS PROVIDE FOR PROPER SYSTEM
C
      DO 15 IX=1,3
   15 RR(IX)=R(INDX(IX))
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
      IF(ARG.LT.ZERO) ARG=ZERO
      ST=SQRT(ARG)
      CT2=CT*CT
      ST2=ST*ST
      CS=CT*ST
C
C     COMPUTE DIATOMIC ENERGIES FOR RR(3) ARRAY
C
      CALL SPLINE(NS,X,S11FIT,RR(1),SS11A,DSS11A,S3,EPS,S2S11,S3S11,
     $DELS11)
      CALL SPLINE(NS,X,S11FIT,RR(3),SS11B,DSS11B,S3,EPS,S2S11,S3S11,
     $DELS11)
      CALL SPLINE(NS,X,S12FIT,RR(1),SS12A,DSS12A,S3,EPS,S2S12,S3S12,
     $DELS12)
      CALL SPLINE(NS,X,S12FIT,RR(3),SS12B,DSS12B,S3,EPS,S2S12,S3S12,
     $DELS12)
      CALL SPLINE(NS,X,S22FIT,RR(1),SS22A,DSS22A,S3,EPS,S2S22,S3S22,
     $DELS22)
      CALL SPLINE(NS,X,S22FIT,RR(3),SS22B,DSS22B,S3,EPS,S2S22,S3S22,
     $DELS22)
      CALL SPLINE(NS,X,TSXFIT,RR(1),TSXA,DTSXA,S3,EPS,S2TSX,S3TSX,
     +            DELTSX)
      CALL SPLINE(NS,X,TSXFIT,RR(3),TSXB,DTSXB,S3,EPS,S2TSX,S3TSX,
     +            DELTSX)
      CALL SPLINE(NS,X,SPXFIT,RR(1),SPXA,DSPXA,S3,EPS,S2SPX,S3SPX,
     +            DELSPX)
      CALL SPLINE(NS,X,SPXFIT,RR(3),SPXB,DSPXB,S3,EPS,S2SPX,S3SPX,
     +            DELSPX)
      CALL SPLINE(NS,X,TPXFIT,RR(1),TPXA,DTPXA,S3,EPS,S2TPX,S3TPX,
     +            DELTPX)
      CALL SPLINE(NS,X,TPXFIT,RR(3),TPXB,DTPXB,S3,EPS,S2TPX,S3TPX,
     +            DELTPX)
      CALL SPLINE(NS,X,DSXFIT,RR(1),DSXA,DDSXA,S3,EPS,S2DSX,S3DSX,
     +            DELDSX)
      CALL SPLINE(NS,X,DSXFIT,RR(3),DSXB,DDSXB,S3,EPS,S2DSX,S3DSX,
     +            DELDSX)
      CALL SPLINE(NS,X,SABFIT,RR(2),SSAB,DSSAB,S3,EPS,S2SAB,S3SAB,
     +            DELSAB)
      CALL SPLINE(NS,X,TABFIT,RR(2),TSAB,DTSAB,S3,EPS,S2TAB,S3TAB,
     +            DELTAB)
      CALL SPLINE(NS,X,GABFIT,RR(2),DGAB,DDGAB,S3,EPS,S2GAB,S3GAB,
     +            DELGAB)
      CALL SPLINE(NS,X,UABFIT,RR(2),DUAB,DDUAB,S3,EPS,S2UAB,S3UAB,
     +            DELUAB)
C
C     CONSTRUCT DIM HAMILTONIAN IN ORTHOGONAL BASIS
C
      TH(1,1)=SS11A+FORTH*((SSAB+TRE*TSAB)+CT2*(SS11B+TRE*TSXB)
     1 +ST2*(SPXB+TRE*TPXB))-DPUX-TWO*DSGH
      TH(1,2)=FORTH*SQRT3*(SSAB-TSAB-CT2*(SS11B-TSXB)-ST2*(SPXB-TPXB))
      TH(1,3)=SS12A
      TH(1,4)=-HALF*CT*SS12B
      TH(1,5)=FORTH*CS*(SS11B+TRE*TSXB-SPXB-TRE*TPXB)
      TH(1,6)=-FORTH*SQRT3*CS*(SS11B-TSXB-SPXB+TPXB)
      TH(2,2)=TSXA+FORTH*((TRE*SSAB+TSAB)+CT2*(TRE*SS11B+TSXB)
     1 +ST2*(TRE*SPXB+TPXB))-DPUX-TWO*DSGH
      TH(2,4)=HALF*SQRT3*CT*SS12B
      TH(2,5)=TH(1,6)
      TH(2,6)=FORTH*CS*(TRE*SS11B+TSXB-TRE*SPXB-TPXB)
      TH(3,3)=SS22A+HALF*(DGAB+DUAB)+DSXB-SSGX-SSGH-DSGH
      TH(3,4)=HALF*(DGAB-DUAB)
      TH(4,4)=DSXA+HALF*(DGAB+DUAB)+SS22B-SSGX-SSGH-DSGH
      TH(4,5)=-HALF*ST*SS12B
      TH(4,6)=HALF*SQRT3*ST*SS12B
      TH(5,5)=SPXA+FORTH*((SSAB+TRE*TSAB)+ST2*(SS11B+TRE*TSXB)
     1 +CT2*(SPXB+TRE*TPXB))-DPUX-TWO*DSGH
      TH(5,6)=FORTH*SQRT3*(SSAB-TSAB-ST2*(SS11B-TSXB)-CT2*(SPXB-TPXB))
      TH(6,6)=TPXA+FORTH*((TRE*SSAB+TSAB)+ST2*(TRE*SS11B+TSXB)
     1 +CT2*(TRE*SPXB+TPXB))-DPUX-TWO*DSGH
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
      IF(IFLAG.LT.ZERO) WRITE(IPRT,998) (VEC(I,NROOT),I=1,NDIM)
998   FORMAT(/,2X,T5,'Vector:',(T15,3E20.8))
C
C************************* END OF FPOT3N *******************************
C
      IF (NDER .NE. 1) RETURN
C     IF(IAPP.EQ.0) RETURN
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
 905  FORMAT(/,2X,T5,'Warning: A serious error may result because ',
     *               'sine(THETA) = ',1PE20.10)
C
C     COMPUTE PARTIALS OF HAMILTONIAN ELEMENTS
C
      XTRA=+HALF*CT*(SS11B+TRE*TSXB-SPXB-TRE*TPXB)
      DH11(1)=DSS11A+XTRA*DTDR1
      DH11(2)=FORTH*(DSSAB+TRE*DTSAB)+XTRA*DTDR2
      DH11(3)=FORTH*(CT2*(DSS11B+TRE*DTSXB)+ST2*(DSPXB+TRE*DTPXB))+XTRA*
     1 DTDR3
      XTRA=-HALF*SQRT3*CT*(SS11B-TSXB-SPXB+TPXB)
      DH12(1)=XTRA*DTDR1
      DH12(2)=FORTH*SQRT3*(DSSAB-DTSAB)+XTRA*DTDR2
      DH12(3)=-FORTH*SQRT3*(CT2*(DSS11B-DTSXB)+ST2*(DSPXB-DTPXB))+XTRA*
     1 DTDR3
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
      DH22(3)=FORTH*(CT2*(TRE*DSS11B+DTSXB)+ST2*(TRE*DSPXB+DTPXB))+XTRA*
     1 DTDR3
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
      DH55(3)=FORTH*(ST2*(DSS11B+TRE*DTSXB)+CT2*(DSPXB+TRE*DTPXB))+XTRA*
     1 DTDR3
      XTRA=+HALF*SQRT3*CT*(SS11B-TSXB-SPXB+TPXB)
      DH56(1)=XTRA*DTDR1
      DH56(2)=FORTH*SQRT3*(DSSAB-DTSAB)+XTRA*DTDR2
      DH56(3)=-FORTH*SQRT3*(ST2*(DSS11B-DTSXB)+CT2*(DSPXB-DTPXB))+XTRA*
     1 DTDR3
      XTRA=-HALF*CT*(TRE*SS11B+TSXB-TRE*SPXB-TPXB)
      DH66(1)=DTPXA+XTRA*DTDR1
      DH66(2)=FORTH*(TRE*DSSAB+DTSAB)+XTRA*DTDR2
      DH66(3)=FORTH*(ST2*(TRE*DSS11B+DTSXB)+CT2*(TRE*DSPXB+DTPXB))+XTRA*
     1 DTDR3
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
900   FORMAT(/,1X,T5,'Error: POT has been called with NDER = ', I5,
     *       /,1X,T12,'only the first derivatives, NDER = 1, are ',
     *                'coded in this potential')
C
      RETURN
      END
C
      SUBROUTINE SETSPL
C
C     SETS UP SPLINE FITS FOR DIATOMIC CURVES
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT34CM/ IPRT, IRD1
      COMMON/POTCM4/X(121),
     $              SABFIT(121),S2SAB(121),S3SAB(121),DELSAB(121),
     $              TABFIT(121),S2TAB(121),S3TAB(121),DELTAB(121),
     $              GABFIT(121),S2GAB(121),S3GAB(121),DELGAB(121),
     $              UABFIT(121),S2UAB(121),S3UAB(121),DELUAB(121),
     $              S11FIT(121),S2S11(121),S3S11(121),DELS11(121),
     $              S12FIT(121),S2S12(121),S3S12(121),DELS12(121),
     $              S22FIT(121),S2S22(121),S3S22(121),DELS22(121),
     $              TSXFIT(121),S2TSX(121),S3TSX(121),DELTSX(121),
     $              SPXFIT(121),S2SPX(121),S3SPX(121),DELSPX(121),
     $              TPXFIT(121),S2TPX(121),S3TPX(121),DELTPX(121),
     $              DSXFIT(121),S2DSX(121),S3DSX(121),DELDSX(121),NS
      COMMON/POTCM3/ZERO,ONE,TWO,TRE,HALF,QUAR,RT3,
     $              DPUX,SSGX,DSGH,SSGH
      DATA AUTERG/4.3598283D-11/,ANGTA0/1.88972634D0/
      NS=121
C
C     SET UP X ARRAY IN BOHR
C
      FAC=0.1D0*ANGTA0
      DO 20 I=1,NS
   20 X(I)=DBLE(I-1)*FAC
C
C     READ IN FIT DATA
C
  999 FORMAT(2D20.13)
      READ(IRD1,999) (SABFIT(I),TABFIT(I),I=1,NS)
      READ(IRD1,999) (GABFIT(I),UABFIT(I),I=1,NS)
      READ(IRD1,998) (S11FIT(I),S12FIT(I),S22FIT(I),TSXFIT(I),I=1,NS)
  998 FORMAT(4D20.13)
      READ(IRD1,997) (SPXFIT(I),TPXFIT(I),DSXFIT(I),I=1,NS)
  997 FORMAT(3D20.13)
C
C     CONVERT TO ATOMIC UNITS
C
      FAC=ONE/AUTERG
      DO 10 I=1,NS
      SABFIT(I)=SABFIT(I)*FAC
      TABFIT(I)=TABFIT(I)*FAC
      GABFIT(I)=GABFIT(I)*FAC
      UABFIT(I)=UABFIT(I)*FAC
      S11FIT(I)=S11FIT(I)*FAC
      S12FIT(I)=S12FIT(I)*FAC
      S22FIT(I)=S22FIT(I)*FAC
      TSXFIT(I)=TSXFIT(I)*FAC
      SPXFIT(I)=SPXFIT(I)*FAC
      TPXFIT(I)=TPXFIT(I)*FAC
      DSXFIT(I)=DSXFIT(I)*FAC
   10 CONTINUE
C
C     READ IN SECOND DERIVATIVES AT ENDS OF FIT CURVES
C     AND CONVERT TO ATOMIC UNITS
C
      FAC=FAC*1.D-16/(ANGTA0*ANGTA0)
      READ(IRD1,999) S2SAB(1),S2TAB(1),S2SAB(NS),S2TAB(NS)
      READ(IRD1,999) S2GAB(1),S2UAB(1),S2GAB(NS),S2UAB(NS)
      READ(IRD1,998) S2S11(1),S2S12(1),S2S22(1),S2TSX(1)
      READ(IRD1,998) S2S11(NS),S2S12(NS),S2S22(NS),S2TSX(NS)
      READ(IRD1,997) S2SPX(1),S2TPX(1),S2DSX(1)
      READ(IRD1,997) S2SPX(NS),S2TPX(NS),S2DSX(NS)
      S2SAB(1)=S2SAB(1)*FAC
      S2TAB(1)=S2TAB(1)*FAC
      S2GAB(1)=S2GAB(1)*FAC
      S2UAB(1)=S2UAB(1)*FAC
      S2S11(1)=S2S11(1)*FAC
      S2S12(1)=S2S12(1)*FAC
      S2S22(1)=S2S22(1)*FAC
      S2TSX(1)=S2TSX(1)*FAC
      S2SPX(1)=S2SPX(1)*FAC
      S2TPX(1)=S2TPX(1)*FAC
      S2DSX(1)=S2DSX(1)*FAC
C
C     GENERATE SPLINE FITS
C
      CALL SPLNIN(NS,X,SABFIT,R,R1,R2,R3,EPS,S2SAB,S3SAB,DELSAB)
      CALL SPLNIN(NS,X,TABFIT,R,R1,R2,R3,EPS,S2TAB,S3TAB,DELTAB)
      CALL SPLNIN(NS,X,GABFIT,R,R1,R2,R3,EPS,S2GAB,S3GAB,DELGAB)
      CALL SPLNIN(NS,X,UABFIT,R,R1,R2,R3,EPS,S2UAB,S3UAB,DELUAB)
      CALL SPLNIN(NS,X,S11FIT,R,R1,R2,R3,EPS,S2S11,S3S11,DELS11)
      CALL SPLNIN(NS,X,S12FIT,R,R1,R2,R3,EPS,S2S12,S3S12,DELS12)
      CALL SPLNIN(NS,X,S22FIT,R,R1,R2,R3,EPS,S2S22,S3S22,DELS22)
      CALL SPLNIN(NS,X,TSXFIT,R,R1,R2,R3,EPS,S2TSX,S3TSX,DELTSX)
      CALL SPLNIN(NS,X,SPXFIT,R,R1,R2,R3,EPS,S2SPX,S3SPX,DELSPX)
      CALL SPLNIN(NS,X,TPXFIT,R,R1,R2,R3,EPS,S2TPX,S3TPX,DELTPX)
      CALL SPLNIN(NS,X,DSXFIT,R,R1,R2,R3,EPS,S2DSX,S3DSX,DELDSX)
      RETURN
      END
C
      SUBROUTINE SPLINE(N,X,Y,T,SS,SS1,SS2,EPSLN,S2,S3,DELY)
C
C     JIM STINE'S SPLINE INTERPOLATOR
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),Y(N),S2(N),S3(N),DELY(N)
      I=1
      IF(T.GT.X(1)) GO TO 10
      SS = Y(1)
      SS1 = DELY(1)
      RETURN
   10 IF(T.LT.X(N)) GO TO 20
      SS = Y(N)
      SS1 = DELY(N)
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
      COMMON /PT34CM/ IPRT, IRD1
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20), IDEBUG
      DIMENSION C(121),H(121),H2(121),B(121)
      DIMENSION X(N),Y(N),S2(N),S3(N),DELY(N)
      DATA EX/1.1D0/,OMEGA/1.0717968D0/,HALF/0.5D0/
      EPSLN=6.4229D-13
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
61    FORMAT(2X,T5,'Final EPS = ', 1PE15.6)
      RETURN
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
C     POST-MULTIPLY IDENTITY BY P-TRANSPOSE MATRIX                       VBDIM
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT34CM/IPRT, IRD1
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20), IDEBUG
      COMMON /POTCM2/ C(3)
      COMMON/POTCM3/ZERO,ONE,TWO,TRE,HALF,QUAR,RT3,
     $              DPUX,SSGX,DSGH,SSGH
C   Initialize the control flags for the potential
      DATA IPRT, IRD1 /6, 4/
      DATA NDER /1/
      DATA IDEBUG /0/
      DATA NFLAG /20*0/ 
C   Initialize the parameters for the potential
      DATA ZERO,ONE,TWO,TRE/0.D0,1.D0,2.D0,3.D0/
      DATA HALF,QUAR/0.5D0,0.25D0/
      DATA RT3/1.7320 50807 56888D0/
      DATA DPUX,SSGX,DSGH,SSGH/0.D0,-0.132774D0,0.D0,0.5D0/
      DATA C/0.169552D0,0.1744734D0,0.169552D0/
      END
