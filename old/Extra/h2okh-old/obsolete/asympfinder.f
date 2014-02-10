c to test the native form of h2okh
c
      implicit real*8 (a-h,o-z)
      dimension cartx(9,10),v(10)
c
c set up test geometry
c
      ndof = 6
      n = 1
      read(5,*) (cartx(i,1),i=1,9)
c      cartx(9,1) = 0.d0
c      cartx(8,1) = 0.d0
c      cartx(7,1) = 0.d0
c      cartx(6,1) = 0.d0
c      cartx(5,1) = 0.d0
c      cartx(4,1) = 1.889d0
c      cartx(3,1) = 0.d0
cc      cartx(2,1) = 1.d0
c      cartx(2,1) = .9510d0
c      cartx(1,1) = -sqrt(cartx(4,1)**2 - cartx(2,1)**2)
      write(6,*) 'test geometry in atomic units'
      write(6,101) (i,cartx(i,1),cartx(i+1,1),cartx(i+2,1),i=1,7,3)
101   format(i5,3f10.5)
      write(6,*) 'test geometry in angstroms'
      d = .5291771d0
      write(6,101) (i,cartx(i,1)*d,cartx(i+1,1)*d,
     x                             cartx(i+2,1)*d,i=1,7,3)
c
c call the potential
c
      call surf(ndof,cartx,v,n)
      write(6,103) v(1),v(1)*627.5096d0
103   format('V(au or kcal):',f15.7,f15.5)
c
c loop over separation of one H
c
      carte = cartx(4,1)
      cartinc = .1d0
      do i = -5,151
         cartx(4,1) = carte + (i-1)*cartinc
         call surf(ndof,cartx,v,n)
         write(6,105) cartx(4,1)*d,v(1)*627.5096d0
105      format(f10.5,',',f10.5)
      enddo
      stop
      end
