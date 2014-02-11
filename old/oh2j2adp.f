C
      SUBROUTINE PREPOT
C
C   System:          OH2 A'' surface
C   Common name:     J2
C   Functional form: Rotated-Morse-oscillator-spline (RMOS) potential
C                    for collinear geometeries plus anti-Morse 
C                    bending potential. Gamma is a function of the
C                    collinear OH distance.
C   References:      T. Joseph, D. G. Truhlar, and B. C. Garrett
C                    J. Chem. Phys. 88, 6982-6990 (1988).
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in the block data subprogram PTPACM.
C   Coordinates, potential energy, and derivatives are passed 
C   through the common block PT31CM:
C                  COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
C   The potential energy in the three asymptotic valleys are 
C   stored in the common block PT35CM:
C                  COMMON /PT35CM/ EASYAB, EASYBC, EASYAC
C   The potential energy in the AB valley, EASYAB, is equal to the potential 
C   energy of the H "infinitely" far from the OH diatomic, with the 
C   OH diatomic at its equilibrium configuration.  Similarly, the terms 
C   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys, 
C   respectively.
C   All the information passed through the common blocks PT31CM and PT35CM 
C   is in hartree atomic units.  
C
C        This potential is written such that:
C                       R(1) = R(O-H)
C                       R(2) = R(H-H)
C                       R(3) = R(H-O)
C   The classical potential energy is set equal to zero for the O
C   infinitely far from the equilibrium H2 diatomic.
C
C   The flags that indicate what calculations should be carried out in 
C   the potential routine are passed through the common block PT32CM:
C                  /PT32CM/ NSURF, NDER, NFLAG(20)
C   where:
C        NSURF - which electronic state should be used.
C                This option is not used for this potential as only the 
C                ground electronic state is available.
C        NDER  = 0 => no derivatives should be calculated
C        NDER  = 1 => calculate first derivatives
C        NFLAG - these integer values can be used to flag options
C                within the potential; in this potential these options 
C                are not used.
C   The common block PT34CM contains the FORTRAN unit numbers for the 
C   potential output.  In this potential PT34CM contains one variable, IPRT,
C                      /PT34CM/ IPRT
C
C   Potential parameters' default settings
C                  Variable            Default value
C                  NSURF               0
C                  NDER                1
C                  IPRT                6
C
C
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
         COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
         COMMON /PT34CM/ IPRT
         COMMON /PT35CM/ EASYAB, EASYBC, EASYAC
         COMMON /J2SWCM/ ASW, IOPSW
C
      DIMENSION DE1(3), DE2(3), DCOSDR(3), DFDR(3)
      DATA ONE,TOL/1.0D0,1.D-8/
C
      HPI = 2.0D0*ATAN(ONE)
      WRITE (IPRT, 600)
      IOPSW = MAX(0, IOPSW)
      ASW = MAX(AMN, ASW)
      WRITE (IPRT, 601) IOPSW, ASW
      CALL PREM2
      ASQ = ASW*ASW
      RETURN
C
600   FORMAT (/,2X,T5,'PREPOT has been called for the A'''' OH2 ',
     *                'potential energy surface J2',
     *       //,2X,T5,'Potential energy surface parameters:')
601   FORMAT (/,2X,T5,'Switching function parameter, n=', T50,I5,
     *        /,2X,T5,'Switching function, Fn(chi) = ',
     *                'sin**2((pi/2)*Fn-1(chi))',
     *        /,2X,T5,'where F0(chi) = sin**2(chi/2)', 
     *        /,2X,T5,'cos(chi) = 0.5*(R(1)**2 - R(3)**2)/R(2)*rho', 
     *        /,2X,T5,'rho**2 = R**2 + a**2,  a=', T50,1PE15.7)
C
      ENTRY POT
C
C   Check the values of NSURF and NDER for validity.
C
         IF (NSURF .NE. 0) THEN
             WRITE (IPRT, 900) NSURF
             STOP 'POT 1'
         ENDIF
         IF (NDER .GT. 1) THEN
             WRITE (IPRT, 910) NDER
             STOP 'POT 2'
         ENDIF
C
C  CHECK FOR SMALL VALUES OF R
C
      IF (R(1) .LT. 0.001D0 .OR. R(2). LT. 0.001D0 .OR. 
     *    R(3) .LT. 0.001D0) THEN
         ENERGY = 1.D10
         IF (NDER .EQ. 1) THEN
             DEDR(1) = -ENERGY
             DEDR(2) = -ENERGY
             DEDR(3) = -ENERGY
             RETURN
         ENDIF
      END IF
C
C    EVALUATE SWITCHING FUNCTION AND ITS DERIVATIVE
C
      R12 = R(1)*R(1)
      R22 = R(2)*R(2)
      R32 = R(3)*R(3)
      BIGR2 = 0.5D0*(R12 - 0.5D0*R22 + R32)
      BIGR2 = MAX(ZERO, BIGR2)
      RHO2 = BIGR2 + ASQ
      FRHO2 = 4.0D0*RHO2
      RHO = SQRT(RHO2)
      DR13 = R12 - R32
      COSCHI = DR13/(2.0D0*R(2)*RHO)
C
C    DERIVATIVES OF COS
C
      IF (NDER .EQ. 1) THEN
          T = 0.25D0/(R(2)*RHO2)
          DCOSDR(2) = -COSCHI*(FRHO2 - R22)*T
          T = T/RHO
          DCOSDR(1) =   R(1)*(FRHO2 - DR13)*T
          DCOSDR(3) = - R(3)*(FRHO2 + DR13)*T
      ENDIF
C
C    EVALUATE SWITCHING FUNCTION
C
      G = COSCHI
      DGDX = 1.0D0
      IF (IOPSW .GT. 0) THEN
         DO 20 I = 1,IOPSW
            T = HPI*G
            G = SIN(T)
            IF (NDER .EQ. 1) DGDX = HPI*COS(T)*DGDX
   20    CONTINUE
      END IF
      F = 0.5D0*(1.0D0 - G)
      IF (NDER .EQ. 1) THEN
          DFDX = -0.5D0*DGDX
          DO 30 I = 1,3
                DFDR(I) = DFDX*DCOSDR(I)
30        CONTINUE
      ENDIF
      OMF = 1.0D0 - F
      IF (F .GT. TOL) THEN
C
C    EVALUATE POTENTIAL FOR O + H2 DIRECTION
C
         CALL POTM3
         V1 = ENERGY
         IF (NDER .EQ. 1) THEN
             DO 50 I = 1,3
50                 DE1(I) = DEDR(I)
         ENDIF
      ELSE
         V1 = 0.0D0
         IF (NDER .EQ. 1) THEN
             DO 55 I = 1,3
55                 DE1(I) = 0.0D0
         ENDIF
      END IF
      IF (F-0.5D0 .LT. TOL) THEN
         V2 = ENERGY
         IF (NDER .EQ. 1) THEN
             DE2(1) = DEDR(3)
             DE2(2) = DEDR(2)
             DE2(3) = DEDR(1)
         ENDIF
      ELSE IF (OMF .GT. TOL) THEN
C
C    INTERCHANGE R(1) AND R(3) AND RECOMPUTE POTENTIAL
C
         T = R(1)
         R(1) = R(3)
         R(3) = T
         CALL POTM3
         V2 = ENERGY
         IF (NDER .EQ. 1) THEN
             DE2(1) = DEDR(3)
             DE2(2) = DEDR(2)
             DE2(3) = DEDR(1)
         ENDIF
         T = R(1)
         R(1) = R(3)
         R(3) = T
      ELSE
         V2 = 0.0D0
         IF (NDER .EQ. 1) THEN
         DO 60 I = 1,3
60             DE2(I) = 0.0D0
         ENDIF
      END IF
C
      ENERGY = F*V1 + OMF*V2
      IF (NDER .EQ. 1) THEN
          DELV = V1 - V2
          DO 300 I = 1,3
                 DEDR(I) = F*DE1(I) + OMF*DE2(I) + DELV*DFDR(I)
300       CONTINUE
      ENDIF
C
900   FORMAT(/,2X,T5,'NSURF has been set equal to ',I5,
     *       /,2X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,2X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/,2X,T5,'POT has been called with NDER = ',I5,
     *       /,2X,T5,'This value of NDER is not allowed in this ',
     *               'version of the potential.')
C
      RETURN
      END
C
      SUBROUTINE PREM2
C
C   The subprogram PREM2 determines the bending correction to the 
C   collinear O + H2 potential subroutine.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION NAB,NAB1,NBC,NBC1
      DIMENSION DBCOOR(2)
      COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
      COMMON /PT34CM/ IPRT
      COMMON /PT35CM/ EASYAB, EASYBC, EASYAC
      COMMON /J2BDCM/ D(3), BETA, REQ(3), A, GAMP(3), RBB
C
C   Conversion for Angstroms to Bohr (Angstroms/Bohr)
C
      PARAMETER (BOHR = .52917706D0)
C
C   Conversion for kcal to Hartree (kcal/Hartree)
C
      PARAMETER (HARTRE = 627.5095D0)
      DATA TOL /1.D-5/
C
C   Initialize the potential parameters for the 
C   collinear portion of the potential energy
C 
      CALL RMOS
C
C   Echo the potential parameters for the bending correction.
C
      WRITE (IPRT,600) REQ, D
      WRITE (IPRT,601) BETA, A, RBB, GAMP
C
C       CONVERT TO ATOMIC UNITS
C
      DO 50 IT = 1,3
         REQ(IT) = REQ(IT)/BOHR
         D(IT)   = D(IT)/HARTRE
   50 CONTINUE
      BETA = BETA * BOHR
      DE  = D(1)
      A = A/BOHR
      RHH = 0.74144D0 / BOHR
      RBB = RBB /BOHR
C
C   Set the values of the energy in the asymptotic valleys
C
      EASYAB = D(1)
      EASYBC = D(2)
      EASYAC = D(3)
C
600   FORMAT(/,2X,T5,'Parameters for the bending correction',
     *       /,2X,T5,'Bond', T47, 'O-H', T58, 'H-H', T69, 'H-O',
     *       /,2X,T5,'Equilibrium bond lengths (Angstroms):', 
     *       T44, F10.5, T55, F10.5, T66, F10.5,
     *       /,2X,T5,'Dissociation energy (kcal/mol):', 
     *       T44, F10.5, T55, F10.5, T66, F10.5)
601   FORMAT(2X,T5,'Morse beta parameter (Angstroms**-1):', T44,1PE15.8,
     *     /,2X,T5,'Pauling parameter (Angstroms):',T44,1PE15.8,
     *     /,2X,T5,'RBB (Angstroms):',T44,1PE15.8,
     *     /,2X,T5,'Gamma for the bending correction:',T44,1PE15.8,
     *          T60,1PE15.8,/,2X,T44,1PE15.8)
C
      RETURN
C
      ENTRY POTM3
C
      IF (R(1) .GT. TOL .AND. R(2) .GT. TOL .AND. R(3) .GT. TOL) THEN
         RACCOL = R(1) + R(2)
         BCOOR = 0.0D0
         DBCOOR(1) = 0.0D0
         DBCOOR(2) = 0.0D0
         DAC = 0.0D0
         ROB = 1.0D0
         IF(ABS(RHH-RBB).GT.1.D-10) THEN
C
C     FIRST CALCULATE SOME USEFUL NUMBERS
C        NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS
C        AB AND BC RESPECTIVELY. (EQUATION 2.2)
C
            NAB = EXP((REQ(1) - R(1))/A)
            IF (NAB.LT.1.D-15) NAB = 1.D-15
            NBC = EXP((REQ(2) - R(2))/A)
            IF (NBC.LT.1.D-15) NBC = 1.D-15
C
C       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE
C       SEE EQUATION 2.7A, 2.7B
C
            C = 1.D0 - NAB/NBC
            IF (ABS(C) .LT. 1.D-14) THEN
               NAB1=.5D0
               NBC1=.5D0
            ELSE
               C = C/NAB
               NAB1 = (2.D0 + C - SQRT(4.D0+C*C))/(2.D0*C)
               NBC1 = 1.D0 - NAB1
            END IF
C
C          ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS
C          FIRST TRAP OUT ANY ZERO ARGUMENTS
C
            IF (NAB1*NBC1 .LE. 0.0D0) THEN
               ROB = 1.D0
            ELSE
               STUFF = A * LOG(NAB1*NBC1)
               T = RHH - STUFF
               ROB = 1.0D0
               IF(ABS(T).GT.1.D-15) ROB = (RBB - STUFF)/T
            END IF
         END IF
C
C        CALCULATE BENDING CORRECTION
C        EVALUATE GAMMA AND DERIVATIVE
C
         EX = EXP(RACCOL*(GAMP(2) + GAMP(3)*RACCOL))
         GAMMA = GAMP(1) + EX
         GAMMA = GAMMA*DE
         IF (NDER .EQ. 1) THEN
             DGAM = (GAMP(2) + 2.D0*GAMP(3)*RACCOL)*EX
             DGAM = DGAM*DE
         ENDIF
C
C           CALCULATE V(R(3)) TERM
C
         EX = EXP(-BETA*(R(3)-REQ(3))/ROB)
         EX1 = EX*(1.D0 + 0.5D0*EX)
         IF (NDER .EQ. 1) DEX1 = EX*(1.D0 + EX)
C 
C           CALCULATE V(R(1)+R(2)) TERM
C
         EX = EXP(BETA*(REQ(3)-R(2)-R(1))/ROB)
         EX2 = EX*(1.D0 + 0.5D0*EX)
         IF (NDER .EQ. 1) DEX2 = EX*(1.D0 + EX)
C
C        HERE IS THE BENDING CORRECTION
C
         BCOOR = GAMMA*(EX1 - EX2)
C
C        WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE
C        OF BCOOR WITH RESPECT TO R(1),  R(2), AND R(3).
C        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO R(3))
C
         IF (NDER .EQ. 1) THEN
             GAMMA = GAMMA*BETA/ROB
             DAC = -GAMMA*DEX1
C
C          CALCULATE CORRECTIONS TO DAB, DBC
C
             DBCOOR(1) = GAMMA*DEX2 + DGAM*(EX1-EX2)
             DBCOOR(2)= DBCOOR(1)
         ENDIF
C
C       NOW CALCULATE THE COLLINEAR ENERGY
C          STORE R(3) IN RACTMP THEN SET R(3) TO THE COLLINEAR
C          GEOMETRY AND CALL POT2D
C
         RACTMP = R(3)
         R(3) = R(2) + R(1)
         CALL POT2D
C
         R(3) = RACTMP
C
         ENERGY = ENERGY + BCOOR
         IF (NDER .NE. 1) RETURN
C
         DO 300 IT = 1,2
              DEDR(IT) = DBCOOR(IT) + DEDR(IT) + DEDR(3)
  300    CONTINUE
         DEDR(3) = DAC
      ELSE
        ENERGY = 1.D30
        IF (NDER .NE. 1) RETURN
        DEDR(1) = 0.0D0
        IF (R(1) .LE. TOL) DEDR(1) = -ENERGY
        DEDR(2) = 0.0D0
        IF (R(2) .LE. TOL) DEDR(2) = -ENERGY
        DEDR(3) = 0.0D0
        IF (R(3) .LE. TOL) DEDR(3) = -ENERGY
      END IF
      RETURN
      END
C
      SUBROUTINE RMOS
C
C     RMO-SPLINE  FIT for collinear O + H2
C     the calling program uses atomic units
C     Appropriate unit conversions are carried out first
C     THE ORIGINAL FIT WAS DONE IN KCAL/MOL AND ANG. UNITS
C     R(1) = R(O-H) ; R(2) = R(H-H)
C
C     CALL RMOS ONLY ONCE. FOR SUBSEQUENT CALLS USE ENTRY POT2D
C
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION THETA(11),DSP(11),BSP(11),XLSP(11)
      DIMENSION AD(11),BD(11),CD(11),DD(11),AB(11),BB(11),
     1          CB(11),DB(11),AL(11),BL(11),CL(11),DL(11)
      DIMENSION TAB1(3),TAB2(3),TAB3(3),W(11),IOP(2)
C
      COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
      COMMON /PT34CM/ IPRT
      COMMON /MCM/ NR
      DATA DETRAD,CKCAL,BHOR/0.01745329D0,627.510D0,0.529177D0/
      DATA IOP,W(1),W(11)/3,3,0.0D+0,0.0D+0/
C
      DATA R1S,R2S,DOH,DHH,BOH,BHH,REOH,REHH/2.21208D+0,2.21531D+0,
     1 1.0656D+2,1.0947D+2,2.2957D+0,1.9444D+0,0.96966D+0,0.74144D+0/
      DATA THETA/0.0D+0,0.10D+2,0.20D+2,0.30D+2,0.3733D+2,0.40D+2,
     10.50D+2,0.60D+2,0.70D+2,0.80D+2,0.90D+2/
      DATA DSP/0.10947D+3,0.1076490D+3,0.1046819D+3,0.993946D+2,
     10.968872D+2,0.971891D+2,0.1005656D+3,0.1041727D+3,
     20.1055045D+3,0.1061649D+3,0.10656D+3/
      DATA BSP/1.9444D+0,1.9065D+0,1.7702D+0,1.5485D+0,1.4360D+0,
     1 1.4230D+0,1.6218D+0,1.8728D+0,2.0717D+0,2.2289D+0,2.2957D+0/
      DATA XLSP/1.47387D+0,1.4755D+0,1.5305D+0,1.6064D+0,1.6349D+0,
     1 1.6323D+0,1.5346D+0,1.3966D+0,1.2998D+0,1.2438D+0,1.24242D+0/
      DATA NSP/11/
C
C     CONVERT TO ATOMIC UNITS
C
      R1S = R1S/BHOR
      R2S = R2S/BHOR
      REOH = REOH/BHOR
      REHH = REHH/BHOR
      BOH = BOH*BHOR
      BHH = BHH*BHOR
      DOH = DOH/CKCAL
      DHH = DHH/CKCAL
      DO 10 I=1,NSP
      THETA(I) = THETA(I)*DETRAD
      DSP(I) = DSP(I)/CKCAL
      BSP(I) = BSP(I)*BHOR
      XLSP(I) = XLSP(I)/BHOR
 10   CONTINUE
C
      WRITE (IPRT, 600)
      WRITE (IPRT, 610) THETA
      WRITE (IPRT, 620)DSP
      WRITE (IPRT, 630)XLSP
      WRITE (IPRT, 640)BSP
C  
C     DETERMINE THE SPLINE COEFFICIENTS
C
      CALL SPL1D1(NSP,THETA,DSP,W,IOP,1,AD,BD,CD)
      CALL SPL1B1(NSP,THETA,DSP,W,1,AD,BD,CD,DD)
      W(1) = 0.0D+0
      W(11) = 0.0D+0
      CALL SPL1D1(NSP,THETA,BSP,W,IOP,1,AB,BB,CB)
      CALL SPL1B1(NSP,THETA,BSP,W,1,AB,BB,CB,DB)
      W(1) = 0.0D+0
      W(11) = 0.0D+0
      CALL SPL1D1(NSP,THETA,XLSP,W,IOP,1,AL,BL,CL)
      CALL SPL1B1(NSP,THETA,XLSP,W,1,AL,BL,CL,DL)
C
600   FORMAT(/,2X,T5,'Parameters for the rotated-Morse-',
     *               'oscillator-spline fit',
     *       /,2X,T5,'Accurate Morse parameters are used for ',
     *               'both asymptotes',
     *       /,2X,T5,'Potential data in atomic units, THETA ',
     *               'in Radians')
610   FORMAT(/,2X,T5,'THETA:',(T15,4(F11.5,1X)))
620   FORMAT(/,2X,T5,'DE:',(T15,4(F11.5,1X)))
630   FORMAT(/,2X,T5,'LE:',(T15,4(F11.5,1X)))
640   FORMAT(/,2X,T5,'BETA:',(T15,4(F11.5,1X)))
C  
      RETURN
C
C**
C
      ENTRY POT2D
C
      IF (R(1) .GE. R1S) THEN
C
C     REACTANT ASYMPTOTIC REGION
C
       NR = 1
       DIF = R(2) - REHH
       IF (DIF*BHH .LT. -70.0D0) WRITE (IPRT,60)R(1),R(2),R(3),DIF
60    FORMAT(/,2X,T5,'Warning: In POT2D R(1), R(2), R(3), DIF = ',
     +       4(F12.4,1X))
      EXP1 = EXP(-BHH*DIF)
      TERM = 1.D+0 - EXP1
      ENERGY = DHH*TERM**2
      IF (NDER .EQ. 1) THEN
          DEDR(1) = 0.0D+0
          DEDR(2) = 2.D+0*DHH*TERM*BHH*EXP1
      ENDIF
        ELSE IF (R(2) .GE. R2S) THEN
C
C     PRODUCT ASYMPTOTIC REGION
C
      NR = 3
      DIF = R(1) - REOH
      IF (DIF*BOH.LT.-70.0D+0) WRITE(IPRT,60)R(1),R(2),R(3),DIF
      EXP2 = EXP(-BOH*DIF)
      TERM = 1.D+0 - EXP2
      ENERGY = DOH*TERM**2 - DOH + DHH
           IF (NDER .EQ. 1) THEN
               DEDR(1) = 2.D+0*DOH*TERM*BOH*EXP2
               DEDR(2) = 0.0D+0
           ENDIF
      ELSE
C
C     SPLINE INTERPOLATION REGION
C
      NR = 2
      R1DIF = R1S - R(1)
      R2DIF = R2S - R(2)
      ZZ = ATAN(R1DIF/R2DIF)
      XL = SQRT(R1DIF**2 + R2DIF**2)
      CALL SPL1B2(NSP,THETA,AD,BD,CD,DD,ZZ,TAB1,1)
      CALL SPL1B2(NSP,THETA,AB,BB,CB,DB,ZZ,TAB2,1)
      CALL SPL1B2(NSP,THETA,AL,BL,CL,DL,ZZ,TAB3,1)
      DIF = XL - TAB3(1)
      EXP3 = EXP(TAB2(1)*DIF)
      TERM = 1.0D+0 - EXP3
      ENERGY = TAB1(1)*TERM**2 - TAB1(1) + DHH
C
C     DERIVATIVES
C     DTH1,DTH2,DL1,DL2 ARE THE DERIVATIVES OF THETA AND L WRT. R(1),R(2)
C     THE DERIVATIVES OF D,B&LE WRT. THETA ARE OBTAINED FROM THE SPLINE
C     FIT. THEY ARE STORED IN TAB1(2),TAB2(2) AND TAB3(2).
C
      IF (NDER .EQ. 1) THEN
          DTH1 = -R2DIF/XL**2
          DTH2 = R1DIF/XL**2
          DL1 = -R1DIF/XL
          DL2 = -R2DIF/XL
C
          PROD = DIF*TAB2(2) - TAB2(1)*TAB3(2)
          DVDTH =TAB1(2)*TERM**2 - 2.0D+0*TAB1(1)*TERM*EXP3*PROD-TAB1(2)
          DVDL = -2.0D+0*TAB1(1)*TAB2(1)*TERM*EXP3
          DEDR(1) = DVDL*DL1 + DVDTH*DTH1
          DEDR(2) = DVDL*DL2 + DVDTH*DTH2
      ENDIF
C
      END IF
      IF(R(1).GT.R1S.AND.R(2).GT.R2S)THEN
      ENERGY = DHH
      IF (NDER .EQ. 1) THEN
          DEDR(1) = 0.0D+0
          DEDR(2) = 0.0D+0
      ENDIF
C
      ENDIF
      DEDR(3) = 0.0D0
      RETURN
      END
C
         BLOCK DATA PTPACM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
         COMMON /PT34CM/ IPRT
         COMMON /J2SWCM/ ASW, IOPSW
         COMMON /J2BDCM/ D(3), BETA, REQ(3), A, GAMP(3), RBB
C
C   Initialize the flags and the I/O unit numbers for the potential
C
         DATA IPRT /6/
         DATA NDER /1/
         DATA NFLAG /20*0/
C
C   Initialize the potential parameters; the energy parameters 
C   are in kcal/mol, and the lengths are in Angstroms.
C   These potential parameters are for the A'' surface of O + H2
C
         DATA IOPSW / 2/
         DATA ASW / 0.4D0/
         DATA D / 106.56D0, 109.472D0, 106.56D0/
         DATA BETA / 2.07942D0/
         DATA REQ / 0.96966D0, 0.74144D0, 0.96966D0/
         DATA A / 0.26D0/
         DATA GAMP / 0.210820D0, 5.42027D0, -1.44235D0/
         DATA RBB / 0.74144D0/
C
         END
C
C Spline utility subprograms
C
C***********************************************************************
C  SPL1D1
C***********************************************************************
C
      SUBROUTINE SPL1D1 (N,X,F,W,IOP,IJ,A,B,C)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C
C     WHERE N = NUMBER OF POINTS IN THE INTERPOLATION
C           X = ORIGIN OF TABLE OF INDEPENDENT VARIALE
C           F = ORIGIN OF TABLE OF DEPENDENT VARIABLE
C           W = AN ARRAY OF DIMENSION N WHICH CONTAINS THE CALCULATED
C               SECOND DERIVATIVES UPON RETURN
C         IOP = AN ARRAY OF DIMENSION 2 WHICH CONTAINS COMBINATIONS OF
C               THE INTEGERS 1 THRU 5 USED TO SPECIFY THE BOUNDARY
C               CONDITIONS
C     IF IOP(1)=1 AND IOP(2)=1, W(1) AND W(K) ARE THE VALUES OF THE
C     SECOND DERIVATIVE AT X(1) AND X(K), RESPECTIVELY.
C     IF IOP(1)=2 AND IOP(2)=2, W(1) DETERMINES F""(X(1)) BY THE
C     RELATION F""(X(1))=W(1)*F""(X(2)), AND W(K) DOES F""(X(K)) BY THE
C     RELATION F""(X(N))=W(K)*F""(X(N-1)).
C     IF IOP(1)=3 AND IOP(2)=3, W(1) AND W(K) ARE THE VALUES OF THE
C     FIRST DERIVATIVE AT X(1) AND X(K), RESPECTIVELY.
C     IF IOP(1)=4 AND IOP(2)=4, THE RERIODIC BOUNDARY CONDITIONS ARE
C     ALLOWED, F""(1)=F""(N), F(1)=F(N).
C     IF IOP(1)=5 AND IOP(2)=5, THE FIRST DERIVATIVES AT X(1) AND X(N)
C     ARE CALCULATED BY USING A DIFFERENTIATED FOUR-POINT LAGRANGIAN
C     INTERPOLATION FORMULA EVALUATED AT X(1) AND X(N), RESPECTIVELYR.
C          IJ = SPACING IN THE F AND W TABLES
C       A,B,C = ARRAYS OF DIMENSION N USED FOR TEMPORARY STORAGE
C
C   Output from this subprogram is written to the file linked to
C   FORTRAN unit IPRT
C
      DIMENSION IOP(2),X(N),F(N),W(N),A(N),B(N),C(N)
      COMMON /PT34CM/ IPRT
C     PARAMETER (IUNIT = 6)
      DATA BOB / 0.0D0 /, BILL / 0.0D0 /, EPS/ 1.0D-11 /
      K = N-1
      SXIV = (1.0D0/6.0D0)-EPS
      THIV = 2.0D0*SXIV
      A(2) = -(X(2)-X(1))*SXIV
      B(2) = (X(3)-X(1))*THIV
      W(IJ+1) = (F(2*IJ+1)-F(IJ+1))/(X(3)-X(2))-(F(IJ+1)-F(1))/(X(2)-X(1
     *            ))
      IF (N-3) 10, 30, 10
   10 DO 20 I = 3, K
         M = (I-1)*IJ+1
         J1 = M+IJ
         J2 = M-IJ
         CON = (X(I+1)-X(I-1))*THIV
         DIFX = X(I)-X(I-1)
         DON = DIFX*SXIV
         BIMI = 1.0D0/B(I-1)
         B(I) = CON-(DON**2)*BIMI
         E = (F(J1)-F(M))/(X(I+1)-X(I))-(F(M)-F(J2))/DIFX
         W(M) = E-(DON*W(J2))*BIMI
         A(I) = -(DON*A(I-1))*BIMI
   20 CONTINUE
   30 K1 = (N-2)*IJ+1
      BNMI = 1.0D0/B(N-1)
      C(N-1) = -(X(N)-X(N-1))*SXIV*BNMI
      W(K1) = W(K1)*BNMI
      A(N-1) = A(N-1)*BNMI
      K2 = K-1
      IF (N-3) 40, 60, 40
   40 DO 50 I = 2, K2
         J = N-I
         CON = (X(J+1)-X(J))*SXIV
         BJI = 1.0D0/B(J)
         A(J) = (A(J)-CON*A(J+1))*BJI
         C(J) = -(CON*C(J+1))*BJI
         K3 = (J-1)*IJ+1
         M = K3+IJ
         W(K3) = (W(K3)-CON*W(M))*BJI
   50 CONTINUE
   60 K4 = (N-1)*IJ+1
      IF (IOP(1)-5) 70, 90, 70
   70 C1 = W(1)
      IF (IOP(2)-5) 80, 110, 80
   80 C2 = W(K4)
      GO TO 130
   90 IF (N-4) 570, 100, 100
  100 A1 = X(1)-X(2)
      A2 = X(1)-X(3)
      A3 = X(1)-X(4)
      A4 = X(2)-X(3)
      A5 = X(2)-X(4)
      A6 = X(3)-X(4)
      W(1) = F(1)*(1.0D0/A1+1.0D0/A2+1.0D0/A3)-A2*A3*F(IJ+1)/(A1*A4*A5)+
     *       A1*A3*F(2*IJ+1)/(A2*A4*A6)-A1*A2*F(3*IJ+1)/(A3*A5*A6)
      GO TO 70
  110 IF (N-4) 570, 120, 120
  120 B1 = X(N)-X(N-3)
      B2 = X(N)-X(N-2)
      B3 = X(N)-X(N-1)
      B4 = X(N-1)-X(N-3)
      B5 = X(N-1)-X(N-2)
      B6 = X(N-2)-X(N-3)
      L1 = K4-IJ
      L2 = L1-IJ
      L3 = L2-IJ
      W(K4) = -B2*B3*F(L3)/(B6*B4*B1)+B1*B3*F(L2)/(B6*B5*B2)-B1*B2*F(L1)
     *        /(B4*B5*B3)+F(K4)*(1.0D0/B1+1.0D0/B2+1.0D0/B3)
      GO TO 80
  130 DO 160 I = 1, K
         M = (I-1)*IJ+1
  170 MK = IOP(1)
      GO TO (180,210,260,310,260), MK
  180 IF (I-1) 200, 190, 200
  190 A(1) = -1.0D0
      C(1) = 0.0D0
      GO TO 340
  200 BOB = 0.0D0
      GO TO 340
  210 IF (I-1) 230, 220, 230
  220 A(1) = -1.0D0
      C(1) = 0.0D0
      W(1) = 0.0D0
      GO TO 340
  230 IF (I-2) 240, 240, 250
  240 BOB = -C1
      GO TO 340
  250 BOB = 0.0D0
      GO TO 340
  260 IF (I-1) 280, 270, 280
  270 XDTO = X(2)-X(1)
      A(1) = -XDTO*THIV
      C(1) = 0.0D0
      W(1) = -C1+(F(IJ+1)-F(1))/XDTO
      GO TO 340
  280 IF (I-2) 290, 290, 300
  290 BOB = (X(2)-X(1))*SXIV
      GO TO 340
  300 BOB = 0.0D0
      GO TO 340
  310 IF (I-1) 330, 320, 330
  320 A(1) = -1.0D0
      C(1) = 1.0D0
      W(1) = 0.0D0
      GO TO 340
  330 BOB = 0.0D0
  340 ML = IOP(2)
      GO TO (350,380,430,480,430), ML
  350 IF (I-1) 370, 360, 370
  360 A(N) = 0.0D0
      C(N) = -1.0D0
      GO TO 140
  370 BILL = 0.0D0
      GO TO 140
  380 IF (I-1) 400, 390, 400
  390 A(N) = 0.0D0
      C(N) = -1.0D0
      W(K4) = 0.0D0
      GO TO 140
  400 IF (I-K) 420, 410, 420
  410 BILL = -C2
      GO TO 140
  420 BILL = 0.0D0
      GO TO 140
  430 IF (I-1) 450, 440, 450
  440 A(N) = 0.0D0
      C(N) = (X(N-1)-X(N))*THIV
      W(K4) = C2-(F(K4)-F(K1))/(X(N)-X(N-1))
      GO TO 140
  450 IF (I-K) 470, 460, 470
  460 BILL = (X(N)-X(N-1))*SXIV
      GO TO 140
  470 BILL = 0.0D0
      GO TO 140
  480 IF (I-1) 500, 490, 500
  490 A(N) = 0.0D0
      C(N) = (X(N-1)+X(1)-X(N)-X(2))*THIV
      W(K4) = (F(IJ+1)-F(1))/(X(2)-X(1))-(F(K4)-F(K1))/(X(N)-X(N-1))
      GO TO 140
  500 IF (I-2) 520, 510, 520
  510 BILL = (X(2)-X(1))*SXIV
      GO TO 140
  520 IF (I-K) 540, 530, 540
  530 BILL = (X(N)-X(N-1))*SXIV
      GO TO 140
  540 BILL = 0.0D0
      GO TO 140
  140    IF (I-1) 150, 160, 150
  150    W(1) = W(1)-BOB*W(M)
         W(K4) = W(K4)-BILL*W(M)
         A(1) = A(1)-BOB*A(I)
         A(N) = A(N)-BILL*A(I)
         C(1) = C(1)-BOB*C(I)
         C(N) = C(N)-BILL*C(I)
  160 CONTINUE
  550 CON = A(1)*C(N)-C(1)*A(N)
      D1 = -W(1)
      D2 = -W(K4)
      W(1) = (D1*C(N)-C(1)*D2)/CON
      W(K4) = (A(1)*D2-D1*A(N))/CON
      DO 560 I = 2, K
         M = (I-1)*IJ+1
         W(M) = W(M)+A(I)*W(1)+C(I)*W(K4)
  560 CONTINUE
      GO TO 580
  570 WRITE (IPRT,1000)
  580 RETURN
C
 1000 FORMAT(1H ,39HSPL1D1 N LESS THAN 4, RESULTS INCORRECT  )
C
      END
C
      SUBROUTINE SPL1B1(N,X,F,W,IJ,A,B,C,D)
C
C     WHERE N = NUMBER OF POINTS IN THE INTERPOLATION
C           X = ORIGIN OF TABLE OF THE INDEPENDENT VARIABLE
C           F = ORIGIN OF TABLE OF THE DEPENDENT VARIABLE
C           W = ORIGIN OF TABLE OF SECOND DERIVATIVES AS CALCULATED
C               BY SPL1D1
C          IJ = SPACING IN THE TABLES F AND W
C           A = ORIGIN OF TABLE OF THE CUBIC COEFFICIENT
C           B = ORIGIN OF TABLE OF THE QUADRATIC COEFFICIENT
C           C = ORIGIN OF TABLE OF THE LINEAR COEFFICIENT
C           D = ORIGIN OF TABLE OF THE CONSTANT COEFFICIENT
C
C
C     SPL1B1 CONVERTS THE SPLINE FIT DATA SUPPLIED BY SPL1D1
C     FROM THE FUNCTION AND ITS SECOND DERIVATIVE AT THE KNOTS,
C     TO THE FOUR COEFFICINTS OF THE CUBIC - A,B,C,D, WHERE
C      F IS NOW APPROXIMATED BY A*X**3 + B*X**2 + C*X + D
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(3),F(3),W(3),A(3),B(3),C(3),D(3)
      DATA C6 /.1666 6666 6666D0/ , C3 /.3333 3333 3333D0/
C
C      GIVEN X,F,W COMPUTE A,B,C,D
C
      NM1 = N - 1
      DO 100 I = 1,NM1
      M = (I-1)*IJ + 1
      MPIJ = M + IJ
      IP1 = I + 1
      XIP1 = X(IP1)
      XI = X(I)
      H = XIP1 - XI
       FI = F(I)
      FIP1 = F(MPIJ)
      WIP1 = W(MPIJ)
      WI = W(M)
      AI = (WIP1 - WI)*C6
      T1 = XIP1*WI
      T2 = XI*WIP1
      BI = .5D0*(T1-T2)
      T1 = XIP1*T1
      T2 = XI*T2
      C(I) = (.5D0*(T2-T1) + FIP1 - FI)/H - AI*H
      A(I) = AI/H
      T1 = XIP1*T1
      T2 = XI*T2
      D(I) = (C6*(T1-T2) + FI*XIP1 - FIP1*XI)/H - C3*H*BI
      B(I) = BI/H
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE SPL1B2(N,X,A,B,C,D,Y,TAB,IOP)
C
C     WHERE N = NUMBER OF POINTS IN THE INTERPOLATION
C           X = ORIGIN OF TABLE OF THE INDEPENDENT VARIABLE
C           A = ORIGIN OF TABLE OF CUBIC COEFFICIENTS
C           B = ORIGIN OF TABLE OF QUADRATIC COEFFICIENTS
C           C = ORIGIN OF TABLE OF LINEAR COEFFICIENTS
C           D = ORIGIN OF TABLE OF CONSTANT COEFFICIENTS
C           Y = THE POINT AT WHICH INTERPOLATION IS DESIRED
C         TAB = AN ARRAY OF DIMENSION 3 WHICH CONTAINS THE FUNCTION
C               VALUE, FIRST DERIVATIVE, AND SECOND DERIVATIVE AT Y
C         IOP = INTEGER SPECIFYING WHETHER DERIVATIVES ARE TO BE COMPUTED
C            .LE. 0 , F ONLY COMPUTED
C             = 1 , F AND FIRST DERIVATIVE ARE COMPUTED
C            .GE. 2 , F AND FIRST AND SECOND DERIVATIVES COMPUTED
C
C       THE FUNCTION F IS NOW APPROXIMATED BY THE CUBIC EQUATION
C          A*X**3 + B*X**2 + C*X + D
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(3),A(3),B(3),C(3),D(3),TAB(3)
C
C     LOCATE Y IN THE X TABLE
C
      NM1 = N - 1
      IF(Y-X(1)) 10,10,20
   10 I = 1
      GO TO 30
   20 IF(Y-X(N)) 15,40,40
   40 I = NM1
      GO TO 30
   15 I = 1
      DO 25 K = 2,NM1
      IF(X(K).GT.Y) GO TO 30
   25 I = I + 1
   30 CONTINUE
C
C      CALCULATE F(Y)
C
      T = A(I)*Y
      TAB(1) = ((T + B(I))*Y + C(I))*Y + D(I)
C
C      CALCULATE THE FIRST DERIVATIVE OF F(Y)
C
      IF(IOP.LE.0) RETURN
      T = 3.D0*T
      TAB(2) = (T+2.D0*B(I))*Y + C(I)
C
C      CALCULATE THE SECOND DERIVATIVE OF F(Y)
C
      IF(IOP.EQ.1) RETURN
      TAB(3) = 2.D0*(T+B(I))
      RETURN
      END
