
      SUBROUTINE PREPOT
C
C   System:          H3 collinear
C   Functional form: Diatomics-in-molecules plus three-center term
C   Common name:     LSTH
C   Reference:       D. G. Truhlar and C. J. Horowitz
C                    J. Chem. Phys. 68, 2466 (1978); 71, 1514(E) (1979)
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in the block data subprogram PTPACM.
C   Coordinates, potential energy, and derivatives are passed 
C   through the common block PT31CM:
C                  COMMON /PT31CM/ RSV(3), ENERGY, DEDR(3)
C   The potential energy in the three asymptotic valleys are 
C   stored in the common block PT35CM:
C                  COMMON /PT35CM/ EASYAB, EASYBC, EASYAC
C   The potential energy in the AB valley, EASYAB, is equal to the potential 
C   energy of the H "infinitely" far from the H2 diatomic, with the 
C   H2 diatomic at its equilibrium configuration.  For this potential energy
C   surface the potential energy in the three asymptotic valleys are equivalent.
C   All the information passed through the common blocks PT31CM and PT35CM
C   is in Hartree atomic units.  
C
C   This potential is written such that:
C                  RSV(1) = R(H-H)
C                  RSV(2) = R(H-H)
C                  RSV(3) = R(H-H)
C   The zero of energy is defined at H "infinitely" far from the H2 diatomic.
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
C
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
C   Units: 
C         Energy: 1 Hartree = 27.21161 ev = 2625.5 kJ/mol
C         Length: 1 bohr = 0.5291771E-10 m
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT31CM/ RSV(3), ENERGY, DEDR(3)
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
      COMMON /PT34CM/ IPRT
      COMMON /PT35CM/ EASYAB, EASYBC, EASYAC
      DIMENSION S1(3),S2(3),S3(3)
      COMMON/POTCOM/C6,C8,RKW(87),EKW(87),WKW(87)
      COMMON/VCOM/C,A,A1,F,FNS,AN1,AN2,AN3,AN4
C
C   Echo the name of the potential to the file linked to FORTRAN unit IPRT
C
      WRITE (IPRT, 600)
C
C   Initialize potential variables.
C
      NDER = 1
      NSURF = 0
      DO 20 I = 1, 20
            NFLAG(I) = 0
20    CONTINUE
C
      R = 1.40105D0
      IC = 0
    5 CALL SPL1D2(87,RKW,EKW,WKW,1,R,S1)
      IF(ABS(S1(2)) .LT. 1.D-10) GO TO 10
      IC = IC + 1
      IF(IC.LT.30) GO TO 7
      WRITE(6,6000)
      STOP 'PREPOT 1'
    7 R = R - S1(2)/S1(3)
      GO TO 5
   10 EZ = S1(1)
C
C   Initialize the energy in the asymptotic valleys
C   The energy in the asymptotic valley is set equal to -EZ because the
C   potential energy is defined as ENERGY  = E - EZ in this routine,
C   but the convention is ENERGY = E + EZ
C
      EASYAB = -EZ
      EASYBC = -EZ
      EASYAC = -EZ
C
600   FORMAT (/,2X,T5,'PREPOT has been called for the H3 ',
     *                'potential energy surface LSTH')
6000  FORMAT (/,1X,T5,'Error: In PREPOT, cannot find minimum of VH2')
      RETURN
C
      ENTRY POT
C
C   Initialize the energy and derivatives
C
      ENERGY  = 0.0D0
      DEDR(1) = 0.0D0
      DEDR(2) = 0.0D0
      DEDR(3) = 0.0D0
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
C E LONDON
C
      EF1=EXP(-F*RSV(1))
      EF2=EXP(-F*RSV(2))
      EF3=EXP(-F*RSV(3))
      X21=RSV(1)*RSV(1)
      X22=RSV(2)*RSV(2)
      X23=RSV(3)*RSV(3)
      T1=C*(A+RSV(1)+A1*X21)*EF1
      T2=C*(A+RSV(2)+A1*X22)*EF2
      T3=C*(A+RSV(3)+A1*X23)*EF3
      CALL VH2(RSV,S1,S2,S3)
      XQ1=S1(1)+T1
      XQ2=S2(1)+T2
      XQ3=S3(1)+T3
      XJ1=S1(1)-T1
      XJ2=S2(1)-T2
      XJ3=S3(1)-T3
      XQ=(XQ1+XQ2+XQ3)/2.D0
      XJ=SQRT(((XJ1-XJ2)**2+(XJ2-XJ3)**2+(XJ3-XJ1)**2)/8.D0)
      ELOND=XQ-XJ
C
C ENS
C
      R=RSV(1)+RSV(2)+RSV(3)
      R2=R*R
      R3=R2*R
      EXNS=EXP(-FNS*R3)
      ENS = 0.D0
      IF(EXNS.EQ.0.D0) GO TO 201
      WNT=(RSV(1)-RSV(2))*(RSV(2)-RSV(3))*(RSV(3)-RSV(1))
      WN=ABS(WNT)
      WN2=WN*WN
      WN3=WN2*WN
      WN4=WN3*WN
      WN5=WN4*WN
      ENS=(AN1*WN2+AN2*WN3+AN3*WN4+AN4*WN5)*EXNS
 201  CONTINUE

      ENERGY = ELOND+ENS-EZ
C
C DERIVATIVES
C
C E LONDON DERIVATIVES
C
      IF (NDER .EQ. 1) THEN
          IF (XJ.EQ.0.0D0) GO TO 1
          XJS=(XJ1+XJ2+XJ3)/8.D0
          T1P=C*(1.D0+2.D0*A1*RSV(1))*EF1-F*T1
          T2P=C*(1.D0+2.D0*A1*RSV(2))*EF2-F*T2
          T3P=C*(1.D0+2.D0*A1*RSV(3))*EF3-F*T3
          ELON1P=(S1(2)+T1P)/2.D0-(S1(2)-T1P)*(.375D0*XJ1-XJS)/XJ
          ELON2P=(S2(2)+T2P)/2.D0-(S2(2)-T2P)*(.375D0*XJ2-XJS)/XJ
          ELON3P=(S3(2)+T3P)/2.D0-(S3(2)-T3P)*(.375D0*XJ3-XJS)/XJ
C
C ENS DERIVATIVES
C
          ENSP1 = 0.D0
          ENSP2 = 0.D0
          ENSP3 = 0.D0
          IF(EXNS.EQ.0.D0) GO TO 30
          ENSPWN=(AN1*WN*2.D0+AN2*3.D0*WN2+AN3*4.D0*WN3+
     *            AN4*5.D0*WN4)*EXNS
          ENSPR=(-3.D0)*FNS*R2*ENS
C
C WN DERIVATIVES
C
          DELTA =-1.D0
          IF (WN.EQ.WNT) DELTA =1.D0
          WNP1=(2.D0*RSV(1)*(RSV(3)-RSV(2))+X22-X23)*DELTA
          WNP2=(2.D0*RSV(2)*(RSV(1)-RSV(3))+X23-X21)*DELTA
          WNP3=(2.D0*RSV(3)*(RSV(2)-RSV(1))+X21-X22)*DELTA
C
C DENS/DXI=(DENS/DWN)(DWN/DXI)+(DENS/DR)(DR/DXI)
C
          ENSP1=ENSPWN*WNP1+ENSPR
          ENSP2=ENSPWN*WNP2+ENSPR
          ENSP3=ENSPWN*WNP3+ENSPR
30        CONTINUE
          DEDR(1)=ELON1P+ENSP1
          DEDR(2)=ELON2P+ENSP2
          DEDR(3)=ELON3P+ENSP3
      ENDIF
      RETURN
1     WRITE (6,2)
2     FORMAT(/,2X,T5,'Warning: The coordinates form an equilateral ',
     *               'triangle and ',
     *       /,2X,T14,'the derivatives will be infinite.',
     *       /,2X,T14,'The derivatives at this geometry have been set ',
     *                'equal to zero.')
      DEDR(1) = 0.0D0
      DEDR(2) = 0.0D0
      DEDR(3) = 0.0D0
C
900   FORMAT(/,2X,T5,'NSURF has been set equal to ',I5,
     *       /,2X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,2X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,
     *       /, 2X, 'This value of NDER is not allowed in this ',
     *              'version of the potential.')
 RETURN
      END
C
      SUBROUTINE VH2(R,S1,S2,S3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/POTCOM/C6,C8,RKW(87),EKW(87),WKW(87)
      DIMENSION R(3),S1(3),S2(3),S3(3)
1     IF (X(1).GT.10.D0) CALL VBIGR(R(1),S1)
      IF (X(1).GT.10.D0) GO TO 2
      CALL SPL1D2(87,RKW,EKW,WKW,1,R(1),S1)
    2 IF (X(2).GT.10.D0) CALL VBIGR(R(2),S2)
      IF (X(2).GT.10.D0) GO TO 3
      CALL SPL1D2(87,RKW,EKW,WKW,1,R(2),S2)
    3 IF (X(3).GT.10.D0) CALL VBIGR(R(3),S3)
      IF (X(3).GT.10.D0) RETURN
      CALL SPL1D2(87,RKW,EKW,WKW,1,R(3),S3)
      RETURN
      END
C
      SUBROUTINE SPL1D2(N,X,F,W,IJ,Y,TAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),F(N),W(N),TAB(3)
      IF (Y-X(1)) 10,10,20
   10 I=1
      GO TO 30
   20 IF (Y-X(N)) 15,40,40
   40 I=N-1
      GO TO 30
   15 I=0
      DO 25 K=1,N
      IF (X(K).GT.Y) GO TO 30
   25 I=I+1
   30 MI =(I-1)*IJ+1
      KI=MI+IJ
      FLK=X(I+1)-X(I)
      A=(W(MI)*(X(I+1)-Y)**3+W(KI)*(Y-X(I))**3)/(6.D0*FLK)
      B=(F(KI)/FLK-W(KI)*FLK/6.D0)*(Y-X(I))
      C=(F(MI)/FLK-FLK*W(MI)/6.D0)*(X(I+1)-Y)
      TAB(1)=A+B+C
      A=(W(KI)*(Y-X(I))**2-W(MI)*(X(I+1)-Y)**2)/(2.D0*FLK )
      B=(F(KI)-F(MI))/FLK
      C=FLK*(W(MI)-W(KI))/6.D0
      TAB(2)=A+B+C
      TAB(3)=(W(MI)*(X(I+1)-Y)+W(KI)*(Y-X(I)))/FLK
      RETURN
      END
C
      SUBROUTINE VBIGR(X,S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/POTCOM/C6,C8,RKW(87),EKW(87),WKW(87)
      DIMENSION S(3)
      X2=X*X
      X3=X2*X
      X6=X3*X3
      C8A=C8/X2
      S(1) = -(C6+C8A)/X6
      S(2)=(C6*6.D0+C8A*8.D0)/X6/X
      RETURN
      END
C
      BLOCK DATA PTPACM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /POTCOM/ C6,C8,RKW(87),EKW(87),WKW(87)
      COMMON/VCOM/C,A,A1,F,FN,XN1,XN2,XN3,XN4
      DATA C6,C8/6.89992032D0,219.9997304D0/
      DATA C,A,A1,F/-1.2148730613D0,-1.514663474D0,-1.46D0,2.088442D0/
      DATA FN,XN1,XN2,XN3,XN4/.0035D0,.0012646477D0,-.0001585792D0,
     1 .0000079707D0,-.0000001151D0/
      DATA (RKW(I),I=1,40)/
     1 .400000000D+00, .450000000D+00, .500000000D+00, .550000000D+00,
     2 .600000000D+00, .650000000D+00, .700000000D+00, .750000000D+00,
     3 .800000000D+00, .900000000D+00, .100000000D+01, .110000000D+01,
     4 .120000000D+01, .130000000D+01, .135000000D+01, .139000000D+01,
     5 .140000000D+01, .140100001D+01, .140109999D+01, .141000000D+01,
     6 .145000000D+01, .150000000D+01, .160000000D+01, .170000000D+01,
     7 .180000000D+01, .190000000D+01, .200000000D+01, .210000000D+01,
     8 .220000000D+01, .230000000D+01, .240000000D+01, .250000000D+01,
     9 .260000000D+01, .270000000D+01, .280000000D+01, .290000000D+01,
     9 .300000000D+01, .310000000D+01, .320000000D+01, .330000000D+01/
      DATA (RKW(I),I=41,87)/
     1 .340000000D+01, .350000000D+01, .360000000D+01, .370000000D+01,
     2 .380000000D+01, .390000000D+01, .400000000D+01, .410000000D+01,
     3 .420000000D+01, .430000000D+01, .440000000D+01, .450000000D+01,
     4 .460000000D+01, .470000000D+01, .480000000D+01, .490000000D+01,
     5 .500000000D+01, .510000000D+01, .520000000D+01, .530000000D+01,
     6 .540000000D+01, .550000000D+01, .560000000D+01, .570000000D+01,
     7 .580000000D+01, .590000000D+01, .600000000D+01, .610000000D+01,
     8 .620000000D+01, .630000000D+01, .640000000D+01, .650000000D+01,
     9 .660000000D+01, .670000000D+01, .680000000D+01, .690000000D+01,
     9 .700000000D+01, .720000000D+01, .740000000D+01, .760000000D+01,
     1 .780000000D+01, .800000000D+01, .824999991D+01, .850000000D+01,
     2 .900000000D+01, .950000000D+01, .100000000D+02/
      DATA (EKW(I),I=1,40)/
     1 .879796188D+00, .649071056D+00, .473372447D+00, .337228924D+00,
     2 .230365628D+00, .145638432D+00, .779738117D-01, .236642733D-01,
     3-.200555771D-01,-.836421044D-01,-.124538356D+00,-.150056027D+00,
     4-.164934012D+00,-.172345701D+00,-.173962500D+00,-.174451499D+00,
     5-.174474200D+00,-.174474400D+00,-.174474400D+00,-.174459699D+00,
     6-.174055600D+00,-.172853502D+00,-.168579707D+00,-.162456813D+00,
     7-.155066822D+00,-.146849432D+00,-.138131041D+00,-.129156051D+00,
     8-.120123163D+00,-.111172372D+00,-.102412583D+00,-.939271927D-01,
     9-.857809026D-01,-.780163108D-01,-.706699181D-01,-.637640270D-01,
     9-.573117349D-01,-.513184414D-01,-.457831464D-01,-.407002530D-01/
      DATA (EKW(I),I=41,87)/
     1-.360577581D-01,-.318401624D-01,-.280271683D-01,-.245977718D-01,
     2-.215296753D-01,-.187966785D-01,-.163688812D-01,-.142246837D-01,
     3-.123370858D-01,-.106809878D-01,-.923028934D-02,-.796819096D-02,
     4-.687029215D-02,-.591779314D-02,-.509229414D-02,-.437819496D-02,
     5-.376259562D-02,-.323089623D-02,-.277399691D-02,-.237999732D-02,
     6-.204229767D-02,-.175209799D-02,-.150299828D-02,-.128989853D-02,
     7-.110689874D-02,-.949798920D-03,-.814999069D-03,-.700199190D-03,
     8-.602999302D-03,-.516199400D-03,-.446599479D-03,-.386399548D-03,
     9-.332799617D-03,-.290599668D-03,-.246599722D-03,-.215399753D-03,
     9-.188899784D-03,-.143399836D-03,-.108599875D-03,-.867998994D-04,
     1-.681999214D-04,-.527999393D-04,-.403999540D-04,-.313999636D-04,
     2-.184999787D-04,-.120999861D-04,-.909998949D-05/
      DATA (WKW(I),I=1,40)/
     1 .308019605D+02, .214419954D+02, .154937452D+02, .115151545D+02,
     2 .871827707D+01, .673831756D+01, .527864661D+01, .419929947D+01,
     3 .333940643D+01, .219403463D+01, .149861953D+01, .103863661D+01,
     4 .730647471D+00, .518552387D+00, .441110777D+00, .383461006D+00,
     5 .373946396D+00, .358559402D+00, .372215569D+00, .356670198D+00,
     6 .312744133D+00, .261523038D+00, .180817537D+00, .124665543D+00,
     7 .807794104D-01, .486562494D-01, .251952492D-01, .452257820D-02,
     8-.854560161D-02,-.196001146D-01,-.276538076D-01,-.344244662D-01,
     9-.381080935D-01,-.421628973D-01,-.441600287D-01,-.454966841D-01,
     9-.460129217D-01,-.458513118D-01,-.453815149D-01,-.440623159D-01/
      DATA (WKW(I),I=41,87)/
     1-.426089183D-01,-.404417185D-01,-.383839285D-01,-.361823035D-01,
     2-.336666088D-01,-.302110314D-01,-.286090554D-01,-.255125522D-01,
     3-.233005599D-01,-.201850499D-01,-.191990995D-01,-.161784216D-01,
     4-.146071006D-01,-.126330766D-01,-.110605069D-01,-.996481997D-02,
     5-.818014482D-02,-.765454189D-02,-.608163613D-02,-.575887028D-02,
     6-.466284400D-02,-.408972107D-02,-.363824334D-02,-.295728079D-02,
     7-.259261281D-02,-.221225014D-02,-.193837141D-02,-.203425060D-02,
     8-.484614204D-03,-.226728547D-02,-.766232140D-03,-.307779418D-03,
     9-.196264565D-02, .131836977D-02,-.223083472D-02,-.750220030D-04,
     9-.289074004D-03,-.220265690D-03,-.434861384D-03, .971346041D-05,
     1-.839919101D-04,-.153745275D-03,-.369227366D-04,-.249634065D-04,
     2-.290482724D-04,-.148433244D-04, .682166282D-05/
      END
