      COMMON /POTCM/ R(3),ENERGY,DER(3)
      
      call prepot
      
      R(1)=1
      R(2)=1.4
      R(3)=2
      
      call pot
      
      
