      Program Standard
C**********************************
C                                 *
C   Written by:  R. J. Duchovic   *
C   Date:        2006-06-24       *
C                2006-06-30       *
C                                 *
C   Program used to calculate the *
C   energy from 19 unmodified PES *
C   functions                     *
C                                 *
C      clh2dima.f                 *
C      clh2dims.f                 *
C      clh2g1.f                   *
C      clh2g3ta1.f                *
C      clh2g3ta2.f                *
C      clh2g3ta3.f                *
C      clh2g3ta4.f                *
C      clhbr.f                    *
C      fh25secgm4.f               *
C      h2skh.f                    *
C      h3bkmp.f                   *
C      h3.f                       *
C      h3lsth1d.f                 *
C      h3lstha.f                  *
C      h3lsthb.f                  *
C      h3lsth.f                   *
C      hcl2sfpolci.f              *
C      oh2dimrmos.f               *
C      oh2j2adp.f                 *
C      oh2j3ap.f                  *
C      oh2modpolci.f              *
C                                 *
C**********************************
      Implicit Real*8 (A-N, O-Z)
      COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
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
      R(1)=1.5d0
      R(2)=2.0d0
      R(3)=3.0d0
      Call Prepot
      Call Pot
      Write(6,1000) R(1), R(2), R(3), ENERGY
 1000 Format(2x, "Geometry",/,
     +       2x,"R1 = ", g20.10,/,
     +       2x,"R2 = ", g20.10,/,
     +       2x,"R3 = ", g20.10,//,
     +       2x, "Energy",/,
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