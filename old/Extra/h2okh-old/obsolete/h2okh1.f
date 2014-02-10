      subroutine surf(ndof,cartx,v,n)                                           
c     Kauppi-Halonen surface for H2O                                            
      implicit double precision(a-h,o-z)                                        
C                                                 
      CHARACTER*75 REF(5)                         
C                                                 
      PARAMETER (I3ATOMS = 3)                     
      PARAMETER (ISURF = 1)                       
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)       
      PARAMETER (PI = 3.141592653589793D0)        
      PARAMETER (ANUZERO=0.0D0)                   
      PARAMETER (NATOMS = 3)                      
      PARAMETER (IREORDER = 0)                    
C                                                 
      COMMON /PT1CM/ R(I3ATOMS),ENGYGS,DEGSDR(I3ATOMS)        
      COMMON /PT2CM/ NSURF, NDER, NFLAG(20)       
      COMMON /PT3CM/ EZERO                        
      COMMON /PT4CM/ ENGYES(ISURF),DEESDR(I3ATOMS,ISURF)      
      COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(I3ATOMS,JSURF)      
      COMMON /ASYCM/ EASYAB,EASYBC,EASYAC         
C                                                 
      COMMON /INFOCM/ CART(NATOMS,3),DGSCART(NATOMS,3),       
     +                DESCART(NATOMS,3,ISURF),    
     +                DIJCART(NATOMS,3,JSURF),    
     +                NULBL(NATOMS),INDEXES(NATOMS),          
     +                IAMAN(NATOMS,2),IASYM(NATOMS,2),REF     
C                                                 
      dimension cartx(ndof+3,*),v(*)                                            
      do 8888 i=1,n                                                             
       c3mc9=cartx(3,i)-cartx(9,i)                                              
       c2mc8=cartx(2,i)-cartx(8,i)                                              
       c1mc7=cartx(1,i)-cartx(7,i)                                              
       c6mc9=cartx(6,i)-cartx(9,i)                                              
       c5mc8=cartx(5,i)-cartx(8,i)                                              
       c4mc7=cartx(4,i)-cartx(7,i)                                              
       r1=dsqrt(c3mc9*c3mc9+c2mc8*c2mc8+c1mc7*c1mc7)                            
       r2=dsqrt(c6mc9*c6mc9+c5mc8*c5mc8+c4mc7*c4mc7)                            
       r1dotr2=(c1mc7*c4mc7+c2mc8*c5mc8+c3mc9*c6mc9)                            
       cosphi=r1dotr2/(r1*r2)                                                   
       phi=acos(cosphi)                                                         
c      input are the two o-h bond lengths and the h-o-h angle; r1,r2,theta.     
c      all calculations in atomic units; angles in radians.                     
       data      de,a,t1,t2,alpha,beta,gamma,delta,epsilon,rkappa,              
     &            rlambda,rmu,rnu,rho,re,phie/                                  
     &    0.195286504346741885d0   , 1.17789560286810513d0   ,                  
     &    0.000000000000000000d+00 , 0.127759668062261463d-01 ,                 
     &   -0.394936461650084985d-02 ,-0.609580676207029798d-02 ,                 
     &    0.822960364451488097d-01 ,-0.249937976602641748d-01 ,                 
     &   -0.104937253831824621d-01 , 0.375501181199095618d-01 ,                 
     &   -0.149932556445900381d-01 ,-0.995022445892074572d-02 ,                 
     &    0.187750590599547809d-01 ,-0.146428799417284827d-01 ,                 
     &    1.80941259980574709d0   ,  1.82404363854352902d0   /                  
c  the potential is expanded in terms of y1,y2,theta in the paper by kauppi     
c  and halonen, as below.                                                       
        y1=1.d0-dexp(-a*(r1-re))                                                
        y2=1.d0-dexp(-a*(r2-re))                                                
        theta=phi-phie                                                          
c                                                                               
        y1sq=y1*y1                                                              
        y2sq=y2*y2                                                              
        y1y2=y1*y2                                                              
        tsq=theta*theta                                                         
c                                                                               
      v(i) = (y1sq+y2sq)*(rmu*tsq+rnu*theta+de)                                 
     &       + (y1+y2)*(rkappa*theta+rlambda*tsq+beta*y1y2)                     
     &       + y1y2*(alpha+rho*theta)                                           
     &       + tsq*(gamma+delta*theta+epsilon*tsq)                              
     &       + t2*(y1sq*y1sq+y2sq*y2sq)                                         
8888  continue                                                                  
      return                                                                    
      end                                                                       
