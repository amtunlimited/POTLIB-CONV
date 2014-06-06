      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /POTCM/ R(3),ENERGY,DER(3)
      

      
      call prepot
      R(1)=1.0d0
      R(2)=2.0d0
      R(3)=3.0d0
      call pot
      
      write(6, *) "ENERGY IS ", ENERGY
      write(6, *) "DER(1) IS ", DER(1)
      write(6, *) "DER(2) IS ", DER(2)
      write(6, *) "DER(3) IS ", DER(3)
      
      END
