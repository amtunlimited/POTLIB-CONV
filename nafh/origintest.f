      double precision r(1,3), e(1), De(3,1)

      r(1,1) = 1.80278
      r(1,2) = 1
      r(1,3) = 1.5

      call prepot()
      call apot(r,e,De,1,1)
      
      write(6,*) "The energy is ", e(1)
      end
