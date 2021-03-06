c****************************************************         
c System: H2Se
c Common Name: h2sekh,
c Functional Form: Expanded into 2 bond length and 1 bond angle               
c Internal Coords: 2 bond length and 1 bond angle
c Features: Calibrated for the equilibrium state          
c****************************************************         
      SUBROUTINE PREPOT
      IMPLICIT REAL*8 (A-H,O-Z)
      character text1,text2,text3,text4,text5
      CHARACTER*75 REF(5)
      PARAMETER (NATOM=25,N3ATOM=3*NATOM)
      PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+1)/2)
      COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),     
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER  
      COMMON /INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),IRCTNT,  
     X                NATOMS,ICARTR,MDER,MSURF,REF
      COMMON /DATACM/ de,a,t1,t2,                                                          
     +            alpha,beta,gamma,delta,epsilon,rkappa,                        
     +            rlambda,rmu,rnu,rho,                                          
     +            re,phie               
c****************************************************         
c any data or labeled commons for PTPACM or POT               
c****************************************************         
      REF(1)= 'Potential Energy for H2Se'
      REF(2)= ' '
      REF(3)= ' '
      REF(4)= ' '
      REF(5)= ' '
      IRCTNT = 3
      INDEXES(1) = 1
      INDEXES(2) = 34                                   
      INDEXES(3) = 1
      CALL POTINFO
      CALL ANCVRT
c****************************************************         
c any assignments for labeled common share with pot           
c****************************************************         
      RETURN
      END
      SUBROUTINE POT
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*75 REF(5)
      PARAMETER (NATOM=25,N3ATOM=3*NATOM)
      PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+1)/2)
      COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),     
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER  
      COMMON /INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),IRCTNT,  
     X                NATOMS,ICARTR,MDER,MSURF,REF            
      COMMON /USROCM/ PENGYGS,PENGYES(ISURF),PENGYIJ(JSURF),  
     X                DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     X                DIJCART(NATOM,3,JSURF)                  
      COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)          
      COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)       
      COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
      COMMON /DATACM/ de,a,t1,t2,                                                          
     +            alpha,beta,gamma,delta,epsilon,rkappa,                        
     +            rlambda,rmu,rnu,rho,                                          
     +            re,phie          
c****************************************************         
c any labeled common statements for PTPACM or POT             
c****************************************************         
      CALL CARTOU
      CALL CARTTOR
      
      cosphi=((R(1)**2 + R(2)**2 - R(3)**2) / (2.d0 * R(1) * R(2)))                                                 
      phi=acos(cosphi)  
c****************************************************         
c code for V,dVdr as functions of int coord R                 
c or for polyatomics as functions of cartesian coord          
c****************************************************
c START OF ORIGINAL CODE
                           
                                                       
                                                                                
c                                                                               
c  The potential is expanded in terms of y1, y2, theta in                       
c  the paper by Kauppi and Halonen, as below.                                   
c                                                                               
	y1=1.0d00-dexp(-a*(R(1)-re))                                                     
	y2=1.0d00-dexp(-a*(R(2)-re))                                                     
	theta=phi-phie                                                                 
c                                                                               
c  for efficiency, next we use the greek letter parameters                      
C  as "force constants" and calculate the potential energy.                     
c                                                                               
	y1sq=y1*y1                                                                     
	y2sq=y2*y2                                                                     
	y1y2=y1*y2                                                                     
	tsq=theta*theta                                                                
c                                                                               
	ENGYGS = (y1sq+y2sq)*(rmu*tsq+rnu*theta+de)                                    
     +       + (y1+y2)*(rkappa*theta+rlambda*tsq+beta*y1y2)                     
     +       + y1y2*(alpha+rho*theta)                                           
     +       + tsq*(gamma+delta*theta+epsilon*tsq)                              
     +       + t2*(y1sq*y1sq+y2sq*y2sq)                                         
c END OF ORIGINAL CODE
         
      CALL EUNITZERO
      IF(NDER.NE.0) THEN
         CALL RTOCART
         IF(NFLAG(1)+NFLAG(2).NE.0) CALL DEDCOU
      ENDIF
      RETURN
      END
C START OF ORIGINAL CODE                                                                               
      SUBROUTINE NPREPOT                                                        
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
      subroutine npot                                                           
c                                                                               
	implicit double precision (a-h,o-z)                                            
      COMMON /PT31CM/ CARTX(9), ENERGY, DEDR(3)                                 
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20), IDBUG                             
      COMMON /PT34CM/ IPRT                                                      
c                                                                               
c	coded by R.Q. Topper.                                                         
c	calculates the potential energy of the H2S molecule                           
c	according to the Kauppi-Halonen surface.                                      
c	input are the two O-H bond lengths and the H-O-H angle; r1, r2, theta.        
c	output is ENERGY, the potential energy.                                       
c	references:                                                                   
c	(1) E. Kauppi and L. Halonen, J.Phys.Chem. 94, 5779 (1990)                    
c	(2) Zhao, Gonzales-Lafont, Truhlar, and Steckler, JCP 94, 5544 (1991).        
c	all calculations in atomic units; angles in degrees.                          
c	parameters created by program module convert.                                 
c                                                                               
c	These are the potential parameters: t=theta, r=r, rp=r' (see Kauppi).         
c                                                                               
	data      de,a,t1,t2,                                                          
     +            alpha,beta,gamma,delta,epsilon,rkappa,                        
     +            rlambda,rmu,rnu,rho,                                          
     +            re,phie/                                                      
     +    0.144760654398338406     ,       !de                                  
     +    0.881450520988266817     ,       !a                                   
     +    0.000000000000000000e+00 ,       !t1                                  
     +    0.110900893492383024e-01 ,       !t2                                  
     +   -0.179227352954911988e-02 ,       !alpha                               
     +   -0.896136764774559938e-03 ,       !beta                                
     +    0.802328447657812416e-01 ,       !gamma                               
     +   -0.615096649861561571e-02 ,       !delta                               
     +    0.000000000000000000e+00 ,       !epsilon                             
     +    0.110148223195307097e-01 ,       !rkappa                              
     +   -0.188996679643331431e-01 ,       !rlambda                             
     +   -0.250743588903988915e-01 ,       !rmu                                 
     +    0.270905047360925745e-02 ,       !rnu                                 
     +   -0.274545233066495094e-01 ,       !rho                                 
     +    2.75900001813446805     ,        !re                                  
     +    1.58074469820494201     /        !thetae                              
c                                                                               
       c3mc9=cartx(3)-cartx(9)                                                  
       c2mc8=cartx(2)-cartx(8)                                                  
       c1mc7=cartx(1)-cartx(7)                                                  
       c6mc9=cartx(6)-cartx(9)                                                  
       c5mc8=cartx(5)-cartx(8)                                                  
       c4mc7=cartx(4)-cartx(7)                                                  
c                                                                               
       r1=sqrt(c3mc9*c3mc9+c2mc8*c2mc8+c1mc7*c1mc7)                             
       r2=sqrt(c6mc9*c6mc9+c5mc8*c5mc8+c4mc7*c4mc7)                             
       r1dotr2=(c1mc7*c4mc7+c2mc8*c5mc8+c3mc9*c6mc9)                            
       cosphi=r1dotr2/(r1*r2)                                                   
       phi=acos(cosphi)                                                         
                                                                                
c                                                                               
c  The potential is expanded in terms of y1, y2, theta in                       
c  the paper by Kauppi and Halonen, as below.                                   
c                                                                               
	y1=1.0d00-dexp(-a*(r1-re))                                                     
	y2=1.0d00-dexp(-a*(r2-re))                                                     
	theta=phi-phie                                                                 
c                                                                               
c  for efficiency, next we use the greek letter parameters                      
C  as "force constants" and calculate the potential energy.                     
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
c                                                                               
	return                                                                         
	end        
c END OR ORIGINAL CODE                                                                    
      BLOCK DATA PTPACM
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*75 REF(5)
      PARAMETER (NATOM=25,N3ATOM=3*NATOM)
      PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+1)/2)
c     PARAMETER (ia=1,ib=1,ic=1,id=1)
      COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),     
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER  
      COMMON /INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),IRCTNT,  
     X                NATOMS,ICARTR,MDER,MSURF,REF
      COMMON /DATACM/ de,a,t1,t2,                                                          
     +            alpha,beta,gamma,delta,epsilon,rkappa,                        
     +            rlambda,rmu,rnu,rho,                                          
     +            re,phie          
c****************************************************         
c all other labeled commons for PREPOT or POT                 
c****************************************************         
      DATA NDER,ANUZERO,NFLAG,NULBL/0,0.0,1,1,15*0,6,0,0,25*0/
      DATA NATOMS,ICARTR,MDER,MSURF/3,3,0,0/
      DATA NASURF/1,35*0/
c****************************************************         
c all other data statements                                   
c****************************************************
      data      de,a,t1,t2,                                                          
     +            alpha,beta,gamma,delta,epsilon,rkappa,                        
     +            rlambda,rmu,rnu,rho,                                          
     +            re,phie/                                                      
     +    0.144760654398338406     ,       !de                                  
     +    0.881450520988266817     ,       !a                                   
     +    0.000000000000000000e+00 ,       !t1                                  
     +    0.110900893492383024e-01 ,       !t2                                  
     +   -0.179227352954911988e-02 ,       !alpha                               
     +   -0.896136764774559938e-03 ,       !beta                                
     +    0.802328447657812416e-01 ,       !gamma                               
     +   -0.615096649861561571e-02 ,       !delta                               
     +    0.000000000000000000e+00 ,       !epsilon                             
     +    0.110148223195307097e-01 ,       !rkappa                              
     +   -0.188996679643331431e-01 ,       !rlambda                             
     +   -0.250743588903988915e-01 ,       !rmu                                 
     +    0.270905047360925745e-02 ,       !rnu                                 
     +   -0.274545233066495094e-01 ,       !rho                                 
     +    2.75900001813446805     ,        !re                                  
     +    1.58074469820494201     /        !thetae          
      END
