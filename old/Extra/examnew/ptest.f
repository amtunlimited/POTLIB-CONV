      Program test
      Implicit Real*8 (A-H,O-Z)
C
      CHARACTER*75 REF(5)
      CHARACTER*5 PERIODIC_1(7,32,5)
C
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
      PARAMETER (PI = 3.141592653589793D0)
C
      CHARACTER*2 NAME1(NATOM)
      CHARACTER*2 NAME2(NATOM)
      CHARACTER*1 IBLANK
C
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON /PT3CM/  EZERO(ISURF+1)
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IAMAN(NATOM,2),IRCTNT,
     +               REF
C     
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NSURF,NDER,
     +               I3ATOMS,NATOMS
C
         nflag(1)=1
         nflag(2)=1
C         nflag(19)=1
C         Call prepot
         cart(1,1)=-0.91666666667d0
         cart(1,2)=1.77756075064d0
         cart(1,3)=0.0d0
         cart(2,1)=0.0d0
         cart(2,2)=0.0d0
         cart(2,3)=0.0d0
         cart(3,1)=1.5d0
         cart(3,2)=0.0d0
         cart(3,3)=0.0d0
C         r(2)=1.5d0
C         r(1)=2.0d0
C         r(3)=3.0d0         
C         Call pot
         write(6,100) engygs           
100      format(2x,'Energy',5x,g25.15,/)
         do i=1,3
         write(6,200) degsdr(i)
200      format(2x, 'Derivative',5x,g25.15)
         enddo
         write(6,250)
250      format(2x)
         write(6,300) easyab
300      format(2x, 'asy-AB',5x,g25.15) 
         write(6,325) easybc
325      format(2x, 'asy-BC',5x,g25.15) 
         write(6,350) easyac
350      format(2x, 'asy-AC',5x,g25.15) 
         stop
         end
