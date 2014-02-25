      implicit double precision (a-h,o-z)
      dimension r(3), e(3)
      
      data r /2.0d0,1.5d0,3.0d0/
      data e /0.0d0,0.0d0,0.0d0/
c nativetest for ClH2g3ta4

c 2, 1.5, 3

      call prepot
      call pot(r,e)
      
      write (6,*) "e(1) is ", e(1)
      write (6,*) "e(2) is ", e(2)
      write (6,*) "e(3) is ", e(3)
      
      end
