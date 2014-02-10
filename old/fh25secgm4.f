C
      SUBROUTINE PREPOT
C
C   System:             FH2
C   Functional form:    Double many-body expansion
C   Common name:        5SECGM4
C   Reference:          unpublished
C   Cross reference:    S. L. Mielke, G. C. Lynch, 
C                       D. G. Truhlar, and D. W. Schwenke
C                       Chem. Phys. Lett. 213, 10-16 (1993)
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in the block data subprogram PTPACM.
C   Coordinates, potential energy, and derivatives are passed 
C   through the common block PT31CM:
C                  COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
C   The potential energy in the three asymptotic valleys are 
C   stored in the common block PT35CM:
C                  COMMON /PT5COM/ EASYAB, EASYBC, EASYAC
C   The potential energy in the AB valley, EASYAB, is equal to the potential 
C   energy of the H "infinitely" far from the FH diatomic, with the 
C   FH diatomic at its equilibrium configuration.  Similarly, the terms 
C   EASYBC and EASYAC represent the H2 and the HF asymptotic valleys, 
C   respectively.
C   All the information passed through the common blocks PT31CM and PT35CM 
C   are in hartree atomic units.  
C
C        This potential is written such that:
C                       R(1) = R(F-H)
C                       R(2) = R(H-H)
C                       R(3) = R(H-F)
C   The classical potential energy is set equal to zero for the F
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
C        NFLAG - these 20 integer values can be used to flag options
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
C   ENERGY  = E(EHF) + E(CORR) + De(H2)
C   E(EHF)  = E(London) + E(3-body)
C   E(CORR) = ECORR(2-body) + ECORR(3-body)
C
C   Note:  The zero of energy for this surface is F infinitely far
C          from H2 at R(H2) = Re(H2).
C          To get the zero of energy R1 must be at least 50 a.u.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
      COMMON /PT34CM/ IPRT
      COMMON /PT35CM/ EASYAB, EASYBC, EASYAC
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
      COMMON /POTDCM/ VDUM(10)
      COMMON /DECOM/ DISSE(3)
C
      COMMON /DIST/RM(3),R0(3)                                          
      COMMON /REF/RO(3)        
      COMMON /TPAR/XHH(3),XHF(3),P(4)      
      COMMON /ATAT/CHH(10),CHF(10)                 
      COMMON /COEFF/CFHH,CHHF  
      COMMON /TRIAL/ETA1(10),ETA2(10),ETA3(10),
     1              ETAP1(10),ETAP2(10),ETAP3(10)  
      COMMON /KKP/CK1(10),CK2(10),CK3(10),CKP1(10),CKP2(10),CKP3(10)      
      COMMON /DIAT/D(3),A1(3),A2(3),A3(3),GAMMA(3)
      COMMON /SWPAR/ DEC(3),RB(3),SW(3),DSW(3,3)
      COMMON /CSI/FCSI,FCSI1,FCSI2
      COMMON /LOCAL/C3C,ALF,ALFP,ALFT,P3C,Q3C,BETP
      COMMON /DISPC/ALPH0,ALPH1,BET0,BET1                           
C
      COMMON /DAMPC/AD1,AD2,AD3,BD1,BD2,BD3
      COMMON /DPC/DP1(3),DP2(3),DP3(3),DDP1(3),DDP2(3),DDP3(3)
      COMMON /DEDRCM/DCORR2(3),DCORR3(3),DLDR(3),DE3C(3)
      COMMON /ENGYCM/ECOR2B, ECOR3B, ELOND, E3C
      COMMON /PAR/BHH(5),BHF(5)
      COMMON /VTANH/RTANH(3,10),DRTANH(3,10)
C
      WRITE(IPRT, 1000) 
      WRITE(IPRT, 1100) (XHH(I),I=1,3)
      WRITE(IPRT, 1200) (XHF(I),I=1,3)
      WRITE(IPRT, 1300) (P(I),I=1,4)
      WRITE(IPRT, 1400) C3C, P3C, ALF, ALFT, ALFP, BETP
C
C     CONSTANTS NEEDED FOR THE DAMPING FUNCTION
C     A=ALPH0/FLOAT(N)**ALPH1        N=6,8,10
C     B=BET0*EXP(-BET1*FLOAT(N))     N=6,8,10
C
      AD1= ALPH0/6.0D0**ALPH1
      AD2= ALPH0/8.0D0**ALPH1
      AD3= ALPH0/10.0D0**ALPH1
      BD1= BET0*EXP(-BET1*6.0D0)
      BD2= BET0*EXP(-BET1*8.0D0)
      BD3= BET0*EXP(-BET1*10.0D0)
C
C    Initialize the energy in the three asymptotic valleys
C
               EASYAB = DISSE(1)
               EASYBC = DISSE(2)
               EASYAC = DISSE(3)
C
1000  FORMAT(/,2X,T5,'PREPOT has been called for the FH2 surface ',
     1               '5SECGM4',
     2       /,2X,T5,'SCALE = 1.0, a = 0.103, b = 1/sqrt(.38), ',
     3               'R0 = 2.0')
1100  FORMAT(/,2X,T5,'Potential parameters',/,
     1       /,2X,T5,'Parameters for the H2 diatomic:',
     2       /,2X,T8,'a4',T13,'=',T15,F10.6,T29,'a5',T33,'=',T35,F10.6,
     3            T48,'a6',T52,'=',T54,F10.6)
1200  FORMAT(/,2X,T5,'Parameters for the FH diatomic:',
     1       /,2X,T8,'Deff',T13,'=',T15,F10.6,T29,'a7',T33,'=',
     2            T35,F10.6,T48,'a8',T52,'=',T54,F10.6)
1300  FORMAT(2X,T8,'c11',T13,'=',T15,F10.6,T29,'c12',T33,'=',
     1            T35,F10.6,/,2X,T8,'c21',T13,'=',T15,F10.6,T29,
     2            'c22',T33,'=',T35,F10.6)
1400  FORMAT(/,2X,T5,'Parameters for the 3-body term:',
     1       /,2X,T8,'J',T13,'=',T15,F10.6,T29,'p',T33,'=',T35,F10.6,
     2            T48,'c1',T52,'=',T54,F10.6,
     3       /,2X,T8,'c2',T13,'=',T15,F10.6,T29,'c3',T33,'=',T35,F10.6,
     4            T48,'c4',T52,'=',T54,F10.6,/)
C
      RETURN
C
C==============================================================================
C   POT: driver for computing the energy and the derivatives.
C==============================================================================
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
C   I: Compute the damping function for the N th. dispersion coefficient
C      and, if NDER = 1, compute the derivative of the damping function
C      w.r.t. R(i)
C
      CALL DAMP(R,RM)
C
C   II: Compute the two-body dynamical correlation term and, if NDER = 1, 
C       compute the derivative of the two-body dynamical correlation term
C       w.r.t. R(i)
C
      CALL CORR2
C
C   III: Compute the auxiliary function, hn(R), and, if NDER = 1,
C            compute the derivative of hn(R) w.r.t. R
C
      CALL XTANH(1,ETA1)
      CALL XTANH(2,ETA2)
      CALL XTANH(3,ETA3)
C
C   IV: Compute the three-body dynamical correlation term and, if NDER = 1, 
C       compute the derivative of the three-body dynamical correlation term
C       w.r.t. R(i)
C
      CALL CORR3
C
C   V: Compute the switching function, S(R), which is used in the computation 
C      of the H2 triplet state curves.  If NDER = 1, compute the derivatives
C      of S w.r.t. R
C
      CALL SWITCH
C
C   VI: Compute the lower eigenvalue of the London eq. and, if NDER = 1, 
C       compute the derivatives of the London energy w.r.t. R.
C
      CALL VLOND
C
C   VII: Compute the three-body term and, if NDER = 1, compute the 
C        the derivatives of the three-body energy w.r.t. R.
C
      CALL THCENT
C
C   IX: Compute the total energy which is a sum of the extended-Hartree-Fock 
C       and the dynamical correlation terms.
C
       ENERGY = ECOR2B + ECOR3B + ELOND + E3C + DISSE(2)
C
C   Compute the derivatives of the energy with respect to the coordinates.
C
      IF (NDER .EQ. 1) THEN
          DO 10 I = 1, 3
                DEDR(I) = DCORR2(I) + DCORR3(I) + DLDR(I) + DE3C(I)
10        CONTINUE
      ENDIF
C
900   FORMAT(/,1X,T5,'NSURF has been set equal to ',I5,
     *       /,1X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,1X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/,1X,T5, 'POT has been called with NDER = ',I5,
     *       /,1X,T5, 'This value of NDER is not allowed in this ',
     *                'version of the potential.')
C
      RETURN
      END
C
C*****
C
      SUBROUTINE CORR2
C==============================================================================
C   Compute the two-body dynamical energies and derivatives.
C==============================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
      COMMON /POTDCM/ VDUM(10)
      COMMON /ATAT/CHH(10),CHF(10)
      COMMON /DEDRCM/DCORR2(3),DCORR3(3),DLDR(3),DE3C(3)
      COMMON /DPC/DP1(3),DP2(3),DP3(3),DDP1(3),DDP2(3),DDP3(3)
      COMMON /ENGYCM/ECOR2B, ECOR3B, ELOND, E3C
C
      R1P6  = R(1)**6
      R1P8  = R(1)**8
      R1P10 = R(1)**10
      R2P6  = R(2)**6
      R2P8  = R(2)**8
      R2P10 = R(2)**10
      R3P6  = R(3)**6
      R3P8  = R(3)**8
      R3P10 = R(3)**10
C
      COR2R1 = -DP1(1)*CHF(6)/R1P6
     1         -DP1(2)*CHF(8)/R1P8
     2         -DP1(3)*CHF(10)/R1P10
C
      COR2R2 = -DP2(1)*CHH(6)/R2P6
     1         -DP2(2)*CHH(8)/R2P8
     2         -DP2(3)*CHH(10)/R2P10
C
      COR2R3 = -DP3(1)*CHF(6)/R3P6
     1         -DP3(2)*CHF(8)/R3P8
     2         -DP3(3)*CHF(10)/R3P10
C
      ECOR2B = COR2R1 + COR2R2 + COR2R3
C
C   Compute the derivatives w.r.t. R(i) if NDER = 1
C
      IF (NDER .NE. 1) GO TO 100
C
      DCORR2(1) = -CHF(6)*(DDP1(1)-6.0D0*DP1(1)/R(1))/R1P6
     1            -CHF(8)*(DDP1(2)-8.0D0*DP1(2)/R(1))/R1P8
     2            -CHF(10)*(DDP1(3)-10.0D0*DP1(3)/R(1))/R1P10
C
      DCORR2(2) = -CHH(6)*(DDP2(1)-6.0D0*DP2(1)/R(2))/R2P6
     1            -CHH(8)*(DDP2(2)-8.0D0*DP2(2)/R(2))/R2P8
     2            -CHH(10)*(DDP2(3)-10.0D0*DP2(3)/R(2))/R2P10
C
      DCORR2(3) = -CHF(6)*(DDP3(1)-6.0D0*DP3(1)/R(3))/R3P6
     1            -CHF(8)*(DDP3(2)-8.0D0*DP3(2)/R(3))/R3P8
     2            -CHF(10)*(DDP3(3)-10.0D0*DP3(3)/R(3))/R3P10
C
100   CONTINUE
C
      RETURN
      END
C
C****
C      
      SUBROUTINE DAMP(RT1,RT2)
C==============================================================================
C   Compute the damping function for the Nth. dispersion coefficient and the
C   derivative w.r.t. R(i)
C==============================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
      COMMON /POTDCM/ VDUM(10)
      COMMON /DAMPC/ AD1,AD2,AD3,BD1,BD2,BD3
      COMMON /DPC/ DP1(3),DP2(3),DP3(3),DDP1(3),DDP2(3),DDP3(3)
      COMMON /DIST/RM(3),R0(3)
      DIMENSION RT1(3),RT2(3),X(3),DX(3),POL1(3),POL2(3),POL3(3)
      DIMENSION DPOL1(3),DPOL2(3),DPOL3(3)
C
            X(1) = 2.0D0*RT1(1)/(RT2(1)+2.5D0*R0(1))
            X(2) = 2.0D0*RT1(2)/(RT2(2)+2.5D0*R0(2))
            X(3) = 2.0D0*RT1(3)/(RT2(3)+2.5D0*R0(3))
C
            POL1(1) = AD1*X(1)+BD1*X(1)**2
            POL2(1) = AD2*X(1)+BD2*X(1)**2
            POL3(1) = AD3*X(1)+BD3*X(1)**2
C
            POL1(2) = AD1*X(2)+BD1*X(2)**2
            POL2(2) = AD2*X(2)+BD2*X(2)**2
            POL3(2) = AD3*X(2)+BD3*X(2)**2
C
            POL1(3) = AD1*X(3)+BD1*X(3)**2
            POL2(3) = AD2*X(3)+BD2*X(3)**2
            POL3(3) = AD3*X(3)+BD3*X(3)**2
C
            DP1(1) = (1.0D0-EXP(-POL1(1)))**6
            DP1(2) = (1.0D0-EXP(-POL2(1)))**8
            DP1(3) = (1.0D0-EXP(-POL3(1)))**10
C
            DP2(1) = (1.0D0-EXP(-POL1(2)))**6
            DP2(2) = (1.0D0-EXP(-POL2(2)))**8
            DP2(3) = (1.0D0-EXP(-POL3(2)))**10
C
            DP3(1) = (1.0D0-EXP(-POL1(3)))**6
            DP3(2) = (1.0D0-EXP(-POL2(3)))**8
            DP3(3) = (1.0D0-EXP(-POL3(3)))**10
C
      IF (NDER .EQ. 1) THEN
          DX(1) = 2.0D0/(RT2(1)+2.5D0*R0(1))
          DX(2) = 2.0D0/(RT2(2)+2.5D0*R0(2))
          DX(3) = 2.0D0/(RT2(3)+2.5D0*R0(3))
C
          DPOL1(1) = AD1+2.0D0*BD1*X(1)
          DPOL2(1) = AD2+2.0D0*BD2*X(1)
          DPOL3(1) = AD3+2.0D0*BD3*X(1)
C
          DPOL1(2) = AD1+2.0D0*BD1*X(2)
          DPOL2(2) = AD2+2.0D0*BD2*X(2)
          DPOL3(2) = AD3+2.0D0*BD3*X(2)
C
          DPOL1(3) = AD1+2.0D0*BD1*X(3)
          DPOL2(3) = AD2+2.0D0*BD2*X(3)
          DPOL3(3) = AD3+2.0D0*BD3*X(3)
C
          DDP1(1) = 6.0D0*DP1(1)*EXP(-POL1(1))*DPOL1(1)*DX(1)/             
     1              (1.0D0-EXP(-POL1(1)))
          DDP1(2) = 8.0D0*DP1(2)*EXP(-POL2(1))*DPOL2(1)*DX(1)/
     1              (1.0D0-EXP(-POL2(1)))
          DDP1(3) = 10.0D0*DP1(3)*EXP(-POL3(1))*DPOL3(1)*DX(1)/
     1              (1.0D0-EXP(-POL3(1)))
C
          DDP2(1) = 6.0D0*DP2(1)*EXP(-POL1(2))*DPOL1(2)*DX(2)/
     1              (1.0D0-EXP(-POL1(2)))
          DDP2(2) = 8.0D0*DP2(2)*EXP(-POL2(2))*DPOL2(2)*DX(2)/
     1              (1.0D0-EXP(-POL2(2)))
          DDP2(3) = 10.0D0*DP2(3)*EXP(-POL3(2))*DPOL3(2)*DX(2)/
     1              (1.0D0-EXP(-POL3(2)))
C
          DDP3(1) = 6.0D0*DP3(1)*EXP(-POL1(3))*DPOL1(3)*DX(3)/
     1              (1.0D0-EXP(-POL1(3)))
          DDP3(2) = 8.0D0*DP3(2)*EXP(-POL2(3))*DPOL2(3)*DX(3)/
     1              (1.0D0-EXP(-POL2(3)))
          DDP3(3) = 10.0D0*DP3(3)*EXP(-POL3(3))*DPOL3(3)*DX(3)/
     1              (1.0D0-EXP(-POL3(3)))
      ENDIF
C
      RETURN
      END
C
C*****
C
      SUBROUTINE CORR3
C==============================================================================
C   Compute the three-body dynamical energies and derivatives.
C==============================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
      COMMON /POTDCM/ VDUM(10)
      COMMON /ATAT/CHH(10),CHF(10)
      COMMON /DEDRCM/DCORR2(3),DCORR3(3),DLDR(3),DE3C(3)
      COMMON /ENGYCM/ECOR2B, ECOR3B, ELOND, E3C
      COMMON /DPC/DP1(3),DP2(3),DP3(3),DDP1(3),DDP2(3),DDP3(3)
      COMMON /DIST/RM(3),R0(3)
      COMMON /KKP/CK1(10),CK2(10),CK3(10),CKP1(10),CKP2(10),CKP3(10)  
      COMMON /REF/RO(3)
      COMMON /TRIAL/ETA1(10),ETA2(10),ETA3(10),
     1              ETAP1(10),ETAP2(10),ETAP3(10)
      COMMON /VTANH/ RTANH(3,10),DRTANH(3,10)
      DIMENSION CO3(3),DC3DR1(3),DC3DR2(3),DC3DR3(3)
C
      JJ=0
C
      DO 10 I=6,10,2
            JJ = JJ+1
            FI = DBLE(I)
            R1PI = R(1)**I
            R2PI = R(2)**I
            R3PI = R(3)**I
C
C   Compute the auxiliary function gn(R)
C
            G1 = 1.0D0+CK1(I)*EXP(-CKP1(I)*(R(1)-RO(1)))
            G2 = 1.0D0+CK2(I)*EXP(-CKP2(I)*(R(2)-RO(2)))
            G3 = 1.0D0+CK3(I)*EXP(-CKP3(I)*(R(3)-RO(3)))
C
C   Compute the auxiliary function hn(R)
C
            H1 = RTANH(1,I)**ETAP1(I)
            H2 = RTANH(2,I)**ETAP2(I)
            H3 = RTANH(3,I)**ETAP3(I)
C
            T1 = CHF(I)*(1.0D0-0.5D0*(G2*H3+G3*H2))
     1           *DP1(JJ)/R1PI
            T2 = CHH(I)*(1.0D0-0.5D0*(G3*H1+G1*H3))
     1           *DP2(JJ)/R2PI
            T3 = CHF(I)*(1.0D0-0.5D0*(G1*H2+G2*H1))
     1           *DP3(JJ)/R3PI
C
            CO3(JJ) = T1+T2+T3
C
            IF (NDER .EQ. 1) THEN
C
C   Compute the partial derivatives w.r.t. R(i)
C
                DG1 = -CKP1(I)*(G1-1.0D0)
                DG2 = -CKP2(I)*(G2-1.0D0)
                DG3 = -CKP3(I)*(G3-1.0D0)
C
                DH1 = ETAP1(I)*(RTANH(1,I)**(ETAP1(I)-1))*DRTANH(1,I)
                DH2 = ETAP2(I)*(RTANH(2,I)**(ETAP2(I)-1))*DRTANH(2,I)
                DH3 = ETAP3(I)*(RTANH(3,I)**(ETAP3(I)-1))*DRTANH(3,I)
C
                DT1DR1 = T1*(DDP1(JJ)/DP1(JJ)-FI/R(1))
                DT1DR2 = -0.5D0*CHF(I)*DP1(JJ)*(DG2*H3+G3*DH2)/R1PI
                DT1DR3 = -0.5D0*CHF(I)*DP1(JJ)*(G2*DH3+DG3*H2)/R1PI
C
                DT2DR1 = -0.5D0*CHH(I)*DP2(JJ)*(G3*DH1+DG1*H3)/R2PI
                DT2DR2 = T2*DDP2(JJ)/DP2(JJ)-T2*FI/R(2)
                DT2DR3 = -0.5D0*CHH(I)*DP2(JJ)*(DG3*H1+G1*DH3)/R2PI
C
                DT3DR1 = -0.5D0*CHF(I)*DP3(JJ)*(DG1*H2+G2*DH1)/R3PI
                DT3DR2 = -0.5D0*CHF(I)*DP3(JJ)*(G1*DH2+DG2*H1)/R3PI
                DT3DR3 = T3*(DDP3(JJ)/DP3(JJ)-FI/R(3))
C
                DC3DR1(JJ) = DT1DR1+DT2DR1+DT3DR1
                DC3DR2(JJ) = DT1DR2+DT2DR2+DT3DR2
                DC3DR3(JJ) = DT1DR3+DT2DR3+DT3DR3
            ENDIF
10    CONTINUE
C
C   Compute the three-body dynamical correlation term by summing over all
C   C6, C8, and C10 terms. 
C
       ECOR3B = CO3(1)+CO3(2)+CO3(3)
C
      IF (NDER .EQ. 1) THEN
C
C   Sum all the partial derivatives.
C
          DCORR3(1) = DC3DR1(1)+DC3DR1(2)+DC3DR1(3)
          DCORR3(2) = DC3DR2(1)+DC3DR2(2)+DC3DR2(3)
          DCORR3(3) = DC3DR3(1)+DC3DR3(2)+DC3DR3(3)
      ENDIF
C
      RETURN
      END 
C   
C*****
C
      SUBROUTINE VLOND
C==============================================================================
C   Compute the energy expressed by the lower eignenvalue of the London eq.
C==============================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
      COMMON /PT2CCM/ NSURF, NDER, NFLAG(20)
      COMMON /POTDCM/ VDUM(10)
C                     
      COMMON /ATAT/CHH(10),CHF(10)       
      COMMON /CSI/FCSI,FCSI1,FCSI2
      COMMON /DAMPC/AD1,AD2,AD3,BD1,BD2,BD3
      COMMON /DEDRCM/DCORR2(3),DCORR3(3),DLDR(3),DE3C(3)
      COMMON /DIAT/D,A1,A2,A3,GAMMA 
      COMMON /DIST/RM,R0                                                 
      COMMON /DPC/DP1(3),DP2(3),DP3(3),DDP1(3),DDP2(3),DDP3(3)
      COMMON /ENGYCM/ECOR2B, ECOR3B, ELOND, E3C
      COMMON /LOCAL/C3C,ALF,ALFP,ALFT,P3C,Q3C,BETP
      COMMON /PAR/BHH(5),BHF(5)
      COMMON /SWPAR/ DEC(3),RB(3),SW(3),DSW(3,3)
      COMMON /TPAR/XHH(3),XHF(3),P(4)
C
      DIMENSION RM(3),R0(3)
      DIMENSION D(3),A1(3),A2(3),A3(3),GAMMA(3)
      DIMENSION DCSI(3),DFCSI1(3),DFCSI2(3),EHF(3),DEHF(3)
      DIMENSION EHFT(3),DEHFT(3),DSUMQ(3),TI(3),DTI(3,3)
      DIMENSION ETRIP(3),DETRIP(3,3),DX(3),DRI(3),X(3)
      DIMENSION Q(3),DQ(3,3),EX(3),DEX(3,3)
      DIMENSION TSCALE(3), DTSCLE(3,3)
C
C   Initialize the arrays that will be used to store the partial 
C   derivatives w.r.t. R(i)
C
      DO 11 I = 1, 3
      DO 12 J = 1, 3
            DETRIP(I,J) = 0.0D0
            DTSCLE(I,J) = 0.0D0
12    CONTINUE
11    CONTINUE
C
C   Compute the pseudo-angular coordinate
C
      CSI = (R(1)-R(3))/R(2)
      CSIP2 = CSI**2
      CSIP3 = CSI**3
      CSIP4 = CSI**4
C 
C    Compute the F1 and F2 terms
C
      FCSI1T = 1.0D0+P(1)*(1.0D0-CSIP2)+P(2)*(1.0D0-CSIP4)
      FCSI2T = 1.0D0+P(3)*(1.0D0-CSIP2)+P(4)*(1.0D0-CSIP4)
C
      FCSI1 = FCSI1T**2
      FCSI2 = FCSI2T**2
C
      IF (NDER .EQ. 1) THEN
C
C   Compute the derivatives w.r.t. R of the terms computed to this point.
C
          DCSI(1) = 1.0D0/R(2)
          DCSI(2) = -CSI/R(2)
          DCSI(3) = -1.0D0/R(2)
C
          TMP1 = 2.0D0*FCSI1T*(-2.0D0*P(1)*CSI-4.0D0*P(2)*CSIP3)
          TMP2 = 2.0D0*FCSI2T*(-2.0D0*P(3)*CSI-4.0D0*P(4)*CSIP3)
C
          DO 10 I=1,3
                DFCSI1(I) = TMP1*DCSI(I)
                DFCSI2(I) = TMP2*DCSI(I)
10        CONTINUE
      ENDIF
C
      R1P2 = R(1)**2
      R2P2 = R(2)**2
      R2P3 = R(2)**3
      R2P4 = R(2)**4
      R3P2 = R(3)**2
C
C   Compute the extended-Hartree-Fock curves and the derivatives w.r.t. R(i) 
C   for the singlet state of H2 and HF.
C   Ref: A. J. C. Varandas and J. D. Silva, J. Chem. Soc. Faraday Trans. II,
C        82, 593 (1986).
C
      DO 20 I=1,3
            DR = R(I)-RM(I)
            DRP2 = DR**2
            EHF(I) = -D(I)*(1.0D0+A1(I)*DR+A2(I)*DRP2+A3(I)*DRP2*DR)
     1               *EXP(-GAMMA(I)*DR)
            IF (NDER .EQ. 1) 
     1          DEHF(I) = -D(I)*(A1(I) + 2.0D0*A2(I)*DR +
     2                    3.0D0*A3(I)*DRP2)*EXP(-GAMMA(I)*DR) -
     3                    GAMMA(I)*EHF(I)
C
20    CONTINUE
C
C   Compute the triplet state curve and the derivative of the triplet state
C   curve w.r.t. R(i) for the H2 diatomic
C   Ref.: A. J. C. Varandas and J. Brando
C         Mol. Phys. 45, 857 (1982)
C
      SCALE = 1.0D0+XHH(1)*R(2)+XHH(2)*R2P2+XHH(3)*R2P3
C
      ARG = BHH(2)*R(2)+BHH(3)*R2P2+BHH(4)*R2P3+BHH(5)*R2P4
C
      EHFT(2) = BHH(1)*EXP(-ARG)/R(2)
C
      ETRIP(2) = SCALE*EHFT(2)
C
      IF (NDER .EQ. 1) THEN
C
C   Compute the derivatives of the H2 triplet curve w.r.t. R(i)
C
          DSCALE = XHH(1)+2.0D0*XHH(2)*R(2)+3.0D0*XHH(3)*R2P2
          DARG = BHH(2)+2.0D0*BHH(3)*R(2)+3.0D0*BHH(4)*R2P2+
     1           4.0D0*BHH(5)*R2P3
          DEHFT(2) = EHFT(2)*(-DARG-1.0D0/R(2))
C
          DETRIP(2,2) = DSCALE*EHFT(2)+SCALE*DEHFT(2)
          DETRIP(2,1) = 0.0D0
          DETRIP(2,3) = 0.0D0
      ENDIF
C
C   Compute the triplet state curve and the derivatives of the triplet 
C   state curve w.r.t. R(i) for the HF diatomic.
C      
      DO 30 I=1, 3, 2
            ARG      = BHF(2)*R(I)+BHF(3)*R(I)**2
            EHFT(I)  = BHF(1)*EXP(-ARG)/R(I)
            X(I)     = EXP(-XHF(2)*R(I))
            ETRIP(I) = XHF(1)*(FCSI1*X(I)*X(I)+FCSI2*XHF(3)*X(I))/R(I)
30    CONTINUE
C
      IF (NDER .EQ. 1) THEN
C
C   Compute the derivatives of the HF triplet state curve w.r.t. R(i)
C
          DO 50 I = 1, 3, 2
                DARG     = BHF(2)+2.0D0*BHF(3)*R(I)
                DEHFT(I) = EHFT(I)*(-1.0D0/R(I)-DARG)
                DX(I)    = -XHF(2)*X(I)
                DRI(I)   = 1.0D0
                DO 40 J = 1, 3
                      IF (J .NE. I) THEN
                          DX(J) = 0.0D0
                          DRI(J) = 0.0D0
                      ENDIF
                      DETRIP(I,J) = XHF(1) * (DFCSI1(J)*X(I)*X(I) +
     1                                        2.0D0*FCSI1*X(I)*DX(J) + 
     2                                        DFCSI2(J)*XHF(3)*X(I) +
     3                                        FCSI2*XHF(3)*DX(J))/R(I) -
     4                              ETRIP(I)*DRI(J)/R(I)
40              CONTINUE
50    CONTINUE
      ENDIF
C
C   Compute the scaling factor for the triplet state curves
C
      DO 51 I = 1, 3
            TSCALE(I) = 1.0D0
51    CONTINUE
      ASC1   = 1.0D0
      ASC2   = 0.103D0
      TMP1   = 1.0D0/0.38D0
      ASC3   = TMP1*TMP1*TMP1*TMP1
      TMP1   = (R(2) - 2.D0)
      TMP1P2 = TMP1 * TMP1
      TMP2   = EXP(-ASC3*TMP1P2*TMP1P2)
      TSCALE(2) = ASC1 + ASC2*TMP2
      IF (NDER .EQ. 1) 
     1    DTSCLE(2,2) = ASC2*TMP2*(-ASC3*4.0D0*TMP1P2*TMP1)
C
C   Compute the triplet state curves and combine the triplet and the
C   singlet state curves to determine the eigenvalue of the 
C   London equation.
C
      DO 60 I = 1, 3
C
C   Compute the triplet state curves
C
            TI(I) = (1.0D0-SW(I))*ETRIP(I)+SW(I)*EHFT(I)
C
C   Scale the triplet state curves
C
            TI(I) = TI(I)*TSCALE(I)
C
C   Compute the coulomb integrals
C
            Q(I) = 0.5D0*(EHF(I)+TI(I))
C
C   Compute the exchange integrals
C
            EX(I) = 0.5D0*(EHF(I)-TI(I))
60    CONTINUE
C
C  Diagnostics
C
      VDUM(1) = TSCALE(2)
      VDUM(2) = TI(2)
C
      FF1 = (EX(1) - EX(2)) ** 2
      FF2 = (EX(2) - EX(3)) ** 2
      FF3 = (EX(3) - EX(1)) ** 2
C
      EXCH = SQRT(0.5D0*(FF1+FF2+FF3))
C
      ELOND = Q(1) + Q(2) + Q(3) - EXCH
C
      IF (NDER .EQ. 1) THEN
C
C   Compute the derivatives of the lowest eignenvalue of the London equation
C   w.r.t. R(i)
C
          DO 70 I = 1, 3
          DO 71 J = 1, 3
                DTI(I,J) = DSW(I,J)*(EHFT(I)-ETRIP(I)) + 
     1                     DETRIP(I,J)*(1.0D0 - SW(I)) 
                IF (J .EQ. I) DTI(I,J) = DTI(I,J) + SW(I)*DEHFT(I)
                DTI(I,J) = DTI(I,J)*TSCALE(I) + 
     1                     TI(I)*DTSCLE(I,J)/TSCALE(I)
71        CONTINUE
70        CONTINUE
C
          DO 80 I = 1, 3
          DO 81 J = 1, 3
                IF(J .EQ. I)THEN
                   DQ(I,J) = 0.5D0*(DEHF(J)+DTI(I,J))
                   DEX(I,J) = 0.5D0*(DEHF(J)-DTI(I,J))
                ELSE
                   DQ(I,J) = 0.5D0*DTI(I,J)
                   DEX(I,J) = -DQ(I,J)
                ENDIF
81        CONTINUE
80        CONTINUE
C
          DSUMQ(1) = DQ(1,1)+DQ(2,1)+DQ(3,1)
          DSUMQ(2) = DQ(1,2)+DQ(2,2)+DQ(3,2)
          DSUMQ(3) = DQ(1,3)+DQ(2,3)+DQ(3,3)
C
          XTERM1 = 2.0D0*EX(1)-EX(2)-EX(3)
          XTERM2 = 2.0D0*EX(2)-EX(1)-EX(3)
          XTERM3 = 2.0D0*EX(3)-EX(1)-EX(2)
C
          DO 90 I = 1, 3
                TEXCH = EXCH
                IF (TEXCH .EQ. 0.0D0)TEXCH = 1.0D0
                XTERM = 0.5D0*(XTERM1*DEX(1,I) + XTERM2*DEX(2,I)+
     1                         XTERM3*DEX(3,I))/TEXCH
                DLDR(I) = DSUMQ(I)-XTERM
90        CONTINUE
      ENDIF
C
      RETURN
      END
C
C****
C
      SUBROUTINE THCENT
C==============================================================================
C   Compute the localized three-body term and the derivatives w.r.t. R(i)
C==============================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
      COMMON /POTDCM/ VDUM(10)
      COMMON /DEDRCM/DCORR2(3),DCORR3(3),DLDR(3),DE3C(3)
      COMMON /ENGYCM/ECOR2B, ECOR3B, ELOND, E3C
      COMMON /LOCAL/C3C,ALF,ALFP,ALFT,P3C,Q3C,BETP
C
      OMP3C = 1.D0 - P3C
      OMQ3C = 1.D0 - Q3C   
      R12 = R(1) * R(1)
      R22 = R(2) * R(2)
      R32 = R(3) * R(3)
      RSUM = R(1) + R(3)
      RDIF = R(1) - R(3)
      RDIF2 = RDIF * RDIF
      T1 = -ALFP * RSUM
      EX1 = C3C * EXP(T1 * RSUM)
      T2 = -ALF * RDIF
      EX2 = OMP3C * EXP(T2 * RDIF)
      T3 = -ALFT * RDIF2
      EX3 = P3C * EXP(T3 * RDIF2)
      R1I = 1.D0 / R(1)
      R3I = 1.D0 / R(3)
      T4 = R1I*R3I
C
C   COSTH IS THE ANGLE BETWEEN R1 AND R3
C
      COSTH = 0.5D0 * (R12 + R32 - R22) * T4
      T5 = OMQ3C * COSTH
C
C   SET MODIFICATION ...
C
      T = Q3C + T5 * COSTH
C
C   E3C CONTAINS THE THREE CENTER CORRECTION TO V
C
      E3C = (EX2 + EX3) * T * EX1
      F3=EXP(-BETP*(R(1)+R(3)-R(2))**2)
      E3C=E3C*F3
      IF (NDER .NE. 1) GO TO 280
C
C   COMPUTE DERIVATIVE OF THE 3C TERM
C   FIRST, DERIVATIVE OF COSTH
C
         DCOS1 = R3I - COSTH*R1I
         DCOS2 = - R(2)*T4
         DCOS3 = R1I - COSTH*R3I
         T4 = EX2 + EX3
         T5 = T4 * T5
C
C   NOW, DERIVATIVES OF E3C
C
         DE3C(1) = 2.D0*((T2 * EX2 + 2.D0 * T3 * RDIF * EX3 + T1 * T4) 
     *      * T + T5 * DCOS1) * EX1*F3+
     1      2.0D0*E3C*(-BETP*(R(1)+R(3)-R(2)))
         DE3C(2) = 2.D0 * T5 * DCOS2 * EX1*F3+
     1           2.0D0*E3C*(BETP*(R(1)+R(3)-R(2)))
         DE3C(3) = 2.D0*((-T2 * EX2 - 2.D0 * T3 * RDIF * EX3 + T1 * T4)
     *      * T + T5 * DCOS3) * EX1*F3+
     1      2.0D0*E3C*(-BETP*(R(1)+R(3)-R(2)))
  280 CONTINUE
      RETURN
      END
C
C*****
C
      SUBROUTINE XTANH(I,A)
C==============================================================================
C   Compute the auxiliary function, hn, and the derivatives of hn w.r.t. R(i)
C==============================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
      COMMON /POTDCM/ VDUM(10)
      COMMON /VTANH/ RTANH(3,10),DRTANH(3,10)
      DIMENSION A(10)
C
      DO 10 J = 6, 10, 2
            RVAL = R(I) * A(J)
            IF(RVAL .GT. 0.0D0) THEN
               EX = EXP( -2.0D0 * RVAL)
               SGN = 1.0D0
            ELSE
               EX = EXP( 2.0D0 * RVAL)
               SGN = -1.0D0
            ENDIF
            T = 1.0D0 + EX
            TH = 1.0D0/T
            IF (NDER .EQ. 1) DRTANH(I,J) = 4.0D0 * A(J) * EX * TH * TH
            RTANH(I,J) = SGN * TH * (1.0D0 - EX)
10    CONTINUE
      RETURN
      END
C
C*****
C
      SUBROUTINE SWITCH
C==============================================================================
C   Compute the switching function, S, and the derivatives of S w.r.t. R(i)
C==============================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
      COMMON /POTDCM/ VDUM(10)
      COMMON /SWPAR/ DEC(3),RB(3),SW(3),DSW(3,3)
C
C   Initialize the array that will be used to store the derivatives of 
C   the switching function w.r.t. R(i)
C
      DO 10 I = 1, 3
      DO 10 J = 1, 3
            DSW(I,J) = 0.0D0
10    CONTINUE
C
C   Compute the switching function.
      DO 20 I=1,3
            RMRB = R(I)-RB(I)
            RVAL = RMRB*DEC(I)
            IF (RVAL .GT. 0.0D0)THEN
                EX2U = EXP(-2.0D0*RVAL)
                EXU  = EXP(-RVAL)
                SGN  = 1.0D0
            ELSE 
                EX2U = EXP(2.0D0*RVAL)
                EXU  = EXP(RVAL)
                SGN = -1.0D0
            ENDIF
            TEMP1 = 1.0D0 + EX2U
            HYTAN = SGN*(1.0D0 - EX2U)/TEMP1
            SW(I) = 0.5D0*(1.0D0 + HYTAN)
   20 CONTINUE
      IF (NDER .NE. 1) GO TO 100
      DO 30 I = 1, 3
            RMRB = R(I)-RB(I)
            RVAL = RMRB*DEC(I)
            IF (RVAL .GT. 0.0D0)THEN
                EXU  = EXP(-RVAL)
                SGN  = 1.0D0
            ELSE 
                EXU  = EXP(RVAL)
                SGN = -1.0D0
            ENDIF
            TEMP2 = 1.0D0 + EXU*EXU
            HYSEC = SGN*2.0D0*EXU/TEMP2
            DSW(I,I) = 0.5D0*DEC(I)*HYSEC*HYSEC
30    CONTINUE
C
100   CONTINUE
      RETURN
      END
C
C*****
C==============================================================================
      BLOCK DATA PTPACM
C==============================================================================
C     INPUT DATA FOR FH2
C
C==============================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT34CM/ IPRT
      COMMON /PT2CCM/ NSURF, NDER, NFLAG(20)
      COMMON /DECOM/ DISSE(3)
      COMMON /ATAT/CHH(10),CHF(10)
      COMMON /DIAT/D(3),A1(3),A2(3),A3(3),GAMMA(3)
      COMMON /DIST/RM(3),R0(3)
      COMMON /DISPC/ALPH0,ALPH1,BET0,BET1                           
      COMMON /KKP/CK1(10),CK2(10),CK3(10),
     1            CKP1(10),CKP2(10),CKP3(10)      
      COMMON /LOCAL/C3C,ALF,ALFP,ALFT,P3C,Q3C,BETP
      COMMON /PAR/BHH(5),BHF(5)
      COMMON /REF/RO(3)       
      COMMON /SWPAR/ DEC(3),RB(3),SW(3),DSW(3,3)
      COMMON /TPAR/XHH(3),XHF(3),P(4)
      COMMON /TRIAL/ETA1(10),ETA2(10),ETA3(10),
     1              ETAP1(10),ETAP2(10),ETAP3(10)
C
      DATA IPRT /6/
      DATA NSURF /0/
      DATA NDER /1/
      DATA NFLAG /20*0/
      DATA DISSE /0.224989335D0, 0.174472692D0, 0.224989355D0/
C
      DATA XHF /2.633449D0, 0.417523D0, -0.090289D0/
      DATA XHH /-0.393828D0, 0.167810D0, -0.017081D0/
      DATA P  /-0.071547D0, -0.108550D0, 0.050461D0, -0.370255D0/
      DATA C3C, ALFP /0.296668D0, 0.078606D0/
C
      DATA ALF,ALFT/0.18D0,2.14D0/
      DATA P3C,Q3C,BETP/0.95D0,1.0D0,0.18D0/
      DATA DEC/2.0D0,2.0D0,2.0D0/
      DATA RB/6.80D0,5.00D0,6.80D0/
      DATA BHH/0.448467D0,-0.056687D0,0.246831D0,-0.018419D0,0.000598D0/
      DATA BHF/0.9149555D0,0.3222197D0,0.1273001D0,0.0D0,0.0D0/
      DATA CHH(6),CHH(8),CHH(10)/6.499027D0,1.243991D2,3.2858D3/
      DATA CHF(6),CHF(8),CHF(10)/6.68799D0,1.0289812D2,2.07391451D3/
      DATA RM(1),RM(2),RM(3)/1.7329D0,1.4010D0,1.7329D0/
      DATA R0(1),R0(2),R0(3)/5.9D0,6.9283D0,5.9D0/
      DATA RO(1),RO(2),RO(3)/1.7329D0,1.449D0,1.7329D0/
      DATA ALPH0,ALPH1,BET0,BET1/2.59528D1,1.1868D0,1.57381D1,9.729D-2/
      DATA ETA1(6),ETA1(8),ETA1(10)/0.9438421608D0,0.9774440533D0,
     1                              0.9452240321D0/              
      DATA ETA2(6),ETA2(8),ETA2(10)/0.1085550150D1,0.1117269260D1,
     1                              0.1138536961D1/  
      DATA ETA3(6),ETA3(8),ETA3(10)/0.9438421608D0,0.9774440533D0,
     1                              0.9452240321D0/  
      DATA ETAP1(6),ETAP1(8),ETAP1(10)/6.0D0,6.0D0,6.0D0/  
      DATA ETAP2(6),ETAP2(8),ETAP2(10)/6.0D0,6.0D0,6.0D0/  
      DATA ETAP3(6),ETAP3(8),ETAP3(10)/6.0D0,6.0D0,6.0D0/
      DATA CK1(6),CK1(8),CK1(10)/-0.2610767389D-1,-0.3428701964D-1,
     1                           -0.4858536634D-1/
      DATA CK2(6),CK2(8),CK2(10)/-0.9709847270D-1,-0.1248186655D 0,
     1                           -0.1426680807D 0/  
      DATA CK3(6),CK3(8),CK3(10)/-0.2610767389D-1,-0.3428701964D-1,
     1                           -0.4858536634D-1/
      DATA CKP1(6),CKP1(8),CKP1(10)/0.9438421608D 0,0.9774440533D 0,
     1                            0.9452240321D 0/     
      DATA CKP2(6),CKP2(8),CKP2(10)/0.1085550150D 1,0.1117269260D 1,
     1                            0.1138536961D 1/  
      DATA CKP3(6),CKP3(8),CKP3(10)/0.9438421608D 0,0.9774440533D 0,
     1                            0.9452240321D 0/   
      DATA D(1),A1(1),A2(1),A3(1),GAMMA(1)/0.19383609D 0,2.4822787D 0,
     1                        1.5435337D 0,0.83093855D 0,2.3992999D 0/
      DATA D(2),A1(2),A2(2),A3(2),GAMMA(2)/0.15796326D0,2.1977034D0,
     1                        1.2932502D0,0.64375666D0,2.1835071D0/
      DATA D(3),A1(3),A2(3),A3(3),GAMMA(3)/0.19383609D 0,2.4822787D 0,
     1                        1.5435337D 0,0.83093855D 0,2.3992999D 0/
      END
C==============================================================================
