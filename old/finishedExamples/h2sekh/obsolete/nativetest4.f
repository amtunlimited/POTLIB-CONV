c to test h2sekh.f 
c
      implicit real*8 (a-h,o-z)
      PARAMETER (NATOM=25,ISURF=5,JSURF=ISURF*(ISURF+1)/2)
      COMMON /PT31CM/ CARTX(9), ENERGY, DEDR(3)
      COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),     
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER  
      COMMON /USROCM/ PENGYGS,PENGYES(ISURF),PENGYIJ(JSURF),  
     X                DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     X                DIJCART(NATOM,3,JSURF)                  
c
c read in input
c
      write(6,*) 'Input Geometry Raw (angstroms)'
      do i = 1,7,3
         read(5,*) (cartx(j),j=i,i+2)
         write(6,101) (cartx(j),j=i,i+2)
101      format(8f10.5)
      enddo
c
c contract and expand by 10% two H's from equilibrium value in raw
c convert to bohr
c
      write(6,*) 'Input Geometry final (bohr)'
      cartx(2) = 1.1d0*cartx(2)
      cartx(4) =  .9d0*cartx(4)
      d = .5291771d0
      do i = 1,7,3
         do j = i,i+2
            cartx(j) = cartx(j)/d
         enddo
         write(6,101) (cartx(j),j=i,i+2)
      enddo
c
c convert cartx to cart and in doing so put Se as the middle atom
c this will allow using CARTTOR options
c
      cart(1,1) = cartx(1)
      cart(1,2) = cartx(2)
      cart(1,3) = cartx(3)
      cart(2,1) = cartx(7)
      cart(2,2) = cartx(8)
      cart(2,3) = cartx(9)
      cart(3,1) = cartx(4)
      cart(3,2) = cartx(5)
      cart(3,3) = cartx(6)
c
c get energy
c
      call prepot
      call pot
      write(6,*) 'Output Energy in au, kcal from pot'
      write(6,101) pengygs,pengygs*627.5096d0
      call ppot
      write(6,*) 'Output Energy in au, kcal'
      write(6,101) energy,energy*627.5096d0
      stop
      end
