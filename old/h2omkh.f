C
      SUBROUTINE PREPOT
C
C   System:          Modified Kauppi-Halonen surface for H2O
C   Common Name:     h20mkh
C   Reference:       E. Kauppi and L. Halonen, J.Phys.Chem. 94, 5779 (1990)
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in DATA statements.
C   Coordinates (Here, Cartesian Coordinates are used by the routine), 
C   potential energy, and derivatives are passed through the common block
C   PT31CM:
C
C                  COMMON /PT31CM/ CARTX(9), ENERGY, DEDR(3)
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
         COMMON /PT31CM/ CARTX(9), ENERGY, DEDR(3)
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
      subroutine pot
c
c      Modified Kauppi-Halonen surface for H2O
c
c      Input are the two O-H bond lengths and the H-O-H angle: r1,r2,theta.
c      All calculations in atomic units; angles in radians.
c 
      implicit double precision(a-h,o-z)
      COMMON /PT31CM/ CARTX(9), ENERGY, DEDR(3)
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20), IDBUG
      COMMON /PT34CM/ IPRT
c
      data      de,a,t1,t2,alpha,beta,gamma,delta,epsilon,rkappa,
     +          rlambda,rmu,rnu,rho,re,phie/
     +     0.126154072681395890d0   , 1.46555639874357269d0   , 
     +     0.000000000000000000d+00 , 0.127759668062261463d-01,
     +    -0.255114862220033774d-02 ,-0.341515428195341688d-02,
     +     0.822960364451488097d-01 ,-0.249937976602641748d-01,
     +    -0.104937253831824621d-01 , 0.301797454253811617d-01,
     +    -0.120503652480248528d-01 ,-0.761011850276373681d-02,
     +     0.150898727126905809d-01 ,-0.945877795944889514d-02,
     +     1.80941259980574709d0   ,  1.82404363854352902d0/ 
c
       c3mc9 = cartx(3) - cartx(9)
       c2mc8 = cartx(2) - cartx(8)  
       c1mc7 = cartx(1) - cartx(7)
       c6mc9 = cartx(6) - cartx(9)
       c5mc8 = cartx(5) - cartx(8)
       c4mc7 = cartx(4) - cartx(7)
c
       r1 = dsqrt(c3mc9*c3mc9+c2mc8*c2mc8+c1mc7*c1mc7)
       r2 = dsqrt(c6mc9*c6mc9+c5mc8*c5mc8+c4mc7*c4mc7)
       r1dotr2=(c1mc7*c4mc7+c2mc8*c5mc8+c3mc9*c6mc9)
       cosphi=r1dotr2/(r1*r2)
       phi=acos(cosphi)
c
c  The potential is expanded in terms of y1, y2, theta in 
c  the paper by Kauppi and Halonen, as below.
c
        y1=1.d0-dexp(-a*(r1-re))
        y2=1.d0-dexp(-a*(r2-re))
        theta=phi-phie
c 
        y1sq=y1*y1
        y2sq=y2*y2
        y1y2=y1*y2
        tsq=theta*theta
c
      ENERGY = (y1sq+y2sq)*(rmu*tsq+rnu*theta+de)
     +       + (y1+y2)*(rkappa*theta+rlambda*tsq+beta*y1y2)
     +       + y1y2*(alpha+rho*theta)
     +       + tsq*(gamma+delta*theta+epsilon*tsq)
     +       + t2*(y1sq*y1sq+y2sq*y2sq)
      return
      end
