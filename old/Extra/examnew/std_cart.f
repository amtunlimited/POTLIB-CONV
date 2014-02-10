      Program Standard
C**********************************
C                                 *
C   Written by:  R. J. Duchovic   *
C   Date:        2006-06-24       *
C                2006-06-30       *
C                                 *
C   Program used to calculate the *
C   energy from 2  unmodified PES *
C   functions                     *
C                                 *
C      h2omkh.f                   *
C      h2sekh.f                   *
C                                 *
C**********************************
      Implicit Real*8 (A-N, O-Z)
      COMMON /PT31CM/ CARTX(9), ENERGY, DEDR(3)
      COMMON /PT61CM/ e(3)
C**********************************
C                                 *
C   INITIALIZE VARIABLES          *
C                                 *
C**********************************
      ENERGY=0.0d0
      DEDR(1)=0.0D0
      DEDR(2)=0.0D0
      DEDR(3)=0.0D0
      e(1)=0.0d0
      e(2)=0.0d0
      e(3)=0.0d0
C**********************************
C                                 *
C   Set Geometry                  *
C                                 *
C**********************************
         cartx(1)=-0.91666666667d0
         cartx(2)=1.77756075064d0
         cartx(3)=0.0d0
         cartx(4)=0.0d0
         cartx(5)=0.0d0
         cartx(6)=0.0d0
         cartx(7)=1.5d0
         cartx(8)=0.0d0
         cartx(9)=0.0d0

      Call Prepot
      Call Pot
      Write(6,1000) cartx(1), cartx(2), cartx(3),
     +              cartx(4), cartx(5), cartx(6), 
     +              cartx(7), cartx(8), cartx(9), ENERGY
 1000 Format(2x, "Geometry",/,
     +       2x,"cc1 = ", g20.10,/,
     +       2x,"cc2 = ", g20.10,/,
     +       2x,"cc3 = ", g20.10,/,
     +       2x,"cc4 = ", g20.10,/,
     +       2x,"cc5 = ", g20.10,/,
     +       2x,"cc6 = ", g20.10,/,
     +       2x,"cc7 = ", g20.10,/,
     +       2x,"cc8 = ", g20.10,/,
     +       2x,"cc9 = ", g20.10,//,
     +       2x, "Energy",//,
     +       2x, "Energy = ", g20.10,//)
      IF(ENERGY.eq.0.0d0) then
         Write(6,1500) e(1),e(2),e(3)
      ENDIF
 1500    Format(2x, "Multiple Energies",/,
     +          2x, "e1 = ", g20.10,/,
     +          2x, "e2 = ", g20.10,/,
     +          2x, "e3 = ",g20.10,//)
      stop
      end