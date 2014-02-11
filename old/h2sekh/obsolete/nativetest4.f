       PROGRAM test
       
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (NATOM=25,N3ATOM=3*NATOM)
       PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+1)/2)
       COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),     
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER
       COMMON /USROCM/ PENGYGS,PENGYES(ISURF),PENGYIJ(JSURF),  
     X                DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     X                DIJCART(NATOM,3,JSURF)
       COMMON /PT31CM/ CARTX(9), ENERGY, DEDR(3)
       
       do i=1,9
         read (5,*) CARTX(i)
       enddo
       
       CART(1,1) = CARTX(1)
       CART(1,2) = CARTX(2)
       CART(1,3) = CARTX(3)
       CART(2,1) = CARTX(7)
       CART(2,2) = CARTX(8)
       CART(2,3) = CARTX(9)
       CART(3,1) = CARTX(4)
       CART(3,2) = CARTX(5)
       CART(3,3) = CARTX(6)
       
       call PREPOT
       call POT
       write (6,*) 'ENERGY: ', PENGYGS
       
       call NPREPOT
       call NPOT
       write (6,*) 'ENERGY: ', ENERGY
       
       end
