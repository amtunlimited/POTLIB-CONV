c to test the native form of h2okh
c
      implicit real*8 (a-h,o-z)
      dimension cartx(9,10),v(10)
c
c info required to be consistent with PREPOT
c
      PARAMETER (NATOM=25)
      PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+1)/2)
      COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER
      COMMON /USROCM/ PENGYGS,PENGYES(ISURF),PENGYIJ(JSURF),
     X                DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     X                DIJCART(NATOM,3,JSURF)
c
c set up test geometry
c
      ndof = 6
      n = 1
      read(5,*) (cartx(i,1),i=1,9)
      write(6,*) 'test geometry in atomic units'
      write(6,101) (i,cartx(i,1),cartx(i+1,1),cartx(i+2,1),i=1,7,3)
101   format(i5,3f10.5)
      write(6,*) 'test geometry in angstroms'
      d = .5291771d0
      write(6,101) (i,cartx(i,1)*d,cartx(i+1,1)*d,
     x                             cartx(i+2,1)*d,i=1,7,3)
c
c fill in CART in an H,O,H order (not the H,H,O order of cartx)
c
      cart(1,1) = cartx(1,1)
      cart(1,2) = cartx(2,1)
      cart(1,3) = cartx(3,1)
      cart(2,1) = cartx(7,1)
      cart(2,2) = cartx(8,1)
      cart(2,3) = cartx(9,1)
      cart(3,1) = cartx(4,1)
      cart(3,2) = cartx(5,1)
      cart(3,3) = cartx(6,1)
c
c call the potential
c
      call prepot
      call pot
      write(6,103) pengygs,pengygs*627.5096d0
      call surf(ndof,cartx,v,n)
      write(6,103) v(1),v(1)*627.5096d0
103   format('V(au or kcal):',f15.7,f15.5)
      stop
      end
