C
      SUBROUTINE PREPOT
C
C   System:          Modified Kauppi-Halonen surface for H2O
C   Common Name:     h20mkh
C   Reference:       E. Kauppi and L. Halonen, J.Phys.Chem. 94, 5779 (1990)
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in DATA statements.
C   Coordinates (TWO BOND LENGTHS, AND A BOND ANGLE), 
C   potential energy, and derivatives are passed through the common block
C   PT31CM:
C
C                  COMMON /PT31CM/ R1, R2, THETA, ENERGY, DEDR(3)
C
C   All the information passed through the common block PT31CM
C   is in Hartree atomic units.
C
C   The flags that indicate what calculations should be carried out in 
C   the potential routine are passed through the common block PT32CM:
C                  /PT32CM/ NSURF, NDER, NFLAG(20), IDBUG
C   where:
C        NSURF - which electronic state should be used.
C                This option is not used for this potential as only the 
C                ground electronic state is available.
C        NDER  = 0 => no derivatives should be calculated
C        NDER  = 1 => calculate first derivatives
C        IDBUG = 0 => do not print extra information
C        NFLAG  - these 7 integer values can be used to flag options
C                within the potential; in this potential these options 
C                are not used.
C        IDBUG > 0 => print details of the energy calculation 
C        IDBUG > 1 => print details of the derivative calculation
C                     All extra output is written to FORTRAN unit 7
C
C   The common block PT34CM contains the FORTRAN unit numbers for the 
C   potential output.  In this potential PT34CM contains one variable, IPRT,
C                      /PT34CM/ IPRT
C
C   Potential parameters' default settings
C                  Variable            Default value
C                  NSURF               0
C                  NDER                1
C                  IDBUG               0
C                  IPRT                6
C
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT31CM/ R1, R2, THETA, ENERGY, DEDR(3)
         COMMON /PT32CM/ NSURF, NDER, NFLAG(20), IDBUG
         COMMON /PT34CM/ IPRT
C
C
         NSURF = 0
         NDER  = 1
         IDBUG = 0
         IPRT  = 6
C
         RETURN
         END
c
	SUBROUTINE POT
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT31CM/ R1, R2, THETA, ENERGY, DEDR(3)
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20), IDBUG
      COMMON /PT34CM/ IPRT
C
C	Coded by R.Q. Topper.
C	Calculates the potential energy of the H2S molecule
C	according to the Kauppi-Halonen surface.
C	Input are the two O-H bond lengths and the H-O-H angle; R1, R2, THETA.
C	Output is ENERGY, the potential energy.
C	References:
C	(1) E. Kauppi and L. Halonen, J.Phys.Chem. 94, 5779(1990)
C	(2) Zhao, Gonzales-Lafont, Truhlar, and Steckler, JCP 94, 5544(1991).
C	All calculations in atomic units; angles in RADIANS. 
C
C	These are the potential parameters. T=theta,R=r,RP=r' (see Kauppi).
C
	DATA      DE,A,T1,T2,
     +            ALPHA,BETA,GAMMA,DELTA,EPSILON,RKAPPA,
     +            RLAMBDA,RMU,RNU,RHO, 
     +            RE,PHIE/
     +    0.159286720904555423     ,       !DE
     +    0.929431016596844017     ,       !A
     +    0.000000000000000000E+00 ,       !T1
     +    0.118745389259935549E-01 ,       !T2
     +   -0.149749727501814604E-02 ,       !ALPHA
     +   -0.748748637509073019E-03 ,       !BETA
     +    0.874694987069420393E-01 ,       !GAMMA
     +   -0.849055136658370853E-02 ,       !DELTA
     +   -0.735898723129194826E-02 ,       !EPSILON
     +    0.157692024067037119E-01 ,       !RKAPPA
     +   -0.174212128271397614E-01 ,       !RLAMBDA
     +   -0.319649409478022645E-01 ,       !RMU
     +    0.788460120335185595E-02 ,       !RNU
     +   -0.154210988516505863E-01 ,       !RHO
     +    2.52391806779835504     ,        !RE
     +    1.60779735487407605     /        !THETAE
C
C  The potential is expanded in terms of Y1, Y2, THETA in 
C  the paper by Kauppi and Halonen, as below.
C
	Y1=1.0D00-DEXP(-A*(R1-RE))
	Y2=1.0D00-DEXP(-A*(R2-RE))
	THETA=PHI-PHIE
C
C  For efficiency, next we use the Greek letter parameters 
C  as "force constants" and calculate the potential energy. 
C
	Y1SQ=Y1*Y1
	Y2SQ=Y2*Y2
	Y1Y2=Y1*Y2
C	
	Y1CU=Y1SQ*Y1
	Y2CU=Y2SQ*Y2
	Y1QU=Y1CU*Y1
	Y2QU=Y2CU*Y2
C
	TSQ=THETA*THETA
	TCU=TSQ*THETA
	TQU=TCU*THETA
C
	ENERGY = DE*(Y1SQ + Y2SQ) + ALPHA*Y1Y2 
     +       + T1*(Y1CU + Y2CU) + T2*(Y1QU + Y2QU)
     +       + BETA*(Y1SQ*Y2 + Y1*Y2SQ)
     +       + GAMMA*TSQ + DELTA*TCU + EPSILON*TQU
     +       + RKAPPA*(Y1 + Y2)*THETA + RLAMBDA*(Y1 + Y2)*TSQ
     +       + RMU*(Y1SQ + Y2SQ)*TSQ + RNU*(Y1SQ+Y2SQ)*THETA
     +       + RHO*Y1Y2*THETA
C
	RETURN 
	END
