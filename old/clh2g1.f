C
      SUBROUTINE PREPOT
C
C   System:           ClH2
C   Functional form:  Extended LEPS (London-Erying-Polyani-Sato)
C   Common name:      G1
C
C   References for the potential parameters and the potential functional form:
C   Morse Parameters: 
C   Sato Parameters:  
C   Functional form:  P. J. Kuntz, E. M. Nemth, J. C. Polanyi, 
C                     S. D. Rosner, and C. E. Young 
C                     J. Chem. Phys. 44, 1168 (1966)
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in the block data subprogram PTPACM.
C   Coordinates, potential energy, and derivatives are passed 
C   through the common block PT31CM:
C                  COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
C   The potential energy in the three asymptotic valleys are 
C   stored in the common block PT5COM:
C                  COMMON /PT35CM/ EASYAB, EASYBC, EASYAC
C   The potential energy in the AB valley, EASYAB, is equal to the potential 
C   energy of the H "infinitely" far from the ClH diatomic, with the 
C   ClH diatomic at its equilibrium configuration.  Similarly, the terms 
C   EASYBC and EASYAC represent the H2 and the HCl asymptotic valleys, 
C   respectively.
C   The other potential parameters are passed through the common blocks
C   SATOCM and LEPSCM.
C   All the information passed through the common blocks PT31CM, PT35CM, 
C   SATOCM, and LEPSCM are in hartree atomic units.  
C
C        This potential is written such that:
C                       R(1) = R(Cl-H)
C                       R(2) = R(H-H)
C                       R(3) = R(H-Cl)
C   The classical potential energy is set equal to zero for the Cl
C   infinitely far from the equilibrium H2 diatomic.
C
C   The flags that indicate what calculations should be carried out in 
C   the potential routine are passed through the common block PT32CM:
C                  /PT32CM/ NSURF, NDER, NFLAG(20)
C   where:
C        NSURF - which electronic state should be used.
C                This option is not used for this potential as only the 
C                ground electronic state is available.
C        NDER  = 0 => no derivatives should be calculated
C        NDER  = 1 => calculate first derivatives
C        NFLAG - these integer values can be used to flag options
C                 within the potential; in this potential these options 
C                 are not used.
C   The common block IOCOM contains the FORTRAN unit numbers for the 
C   potential output.  In this potential IOCOM contains one variable, IPRT,
C                      /PT34CM/ IPRT
C
C   Potential parameters' default settings
C                  Variable            Default value
C                  NSURF               0
C                  NDER                1
C                  IPRT                6
C
C
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
         COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
         COMMON /PT34CM/ IPRT
         COMMON /PT35CM/ EASYAB, EASYBC, EASYAC
         COMMON /SATOCM/ DE(3), RE(3), BETA(3), Z(3) 
         COMMON /LEPSCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3), 
     *                   DO4Z(3), B(3)
         PARAMETER (CKCAL = 627.5095D0)
         PARAMETER (CANGS =   0.529177106D0)
C
C   Echo the potential parameters
C
         WRITE (IPRT, 100) DE, RE, BETA, Z
C
100   FORMAT (/, 2X, T5, 'PREPOT has been called for the ClH2 ',
     *                   'G1 potential',
     *        //, 2X, T5, 'Potential energy surface parameters:',
     *        /, 2X, T5, 'Bond', T46, 'Cl-H', T58, 'H-H', T69, 'H-Cl',
     *        /, 2X, T5, 'Dissociation energies (kcal/mol):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /, 2X, T5, 'Equilibrium bond lengths (Angstroms):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /, 2X, T5, 'Morse beta parameters (Angstroms**-1):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /, 2X, T5, 'Sato parameters:', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,/)
C
      DO  10 I = 1,3
C
C   Convert to atomic units
C
             DE(I)  = DE(I)/CKCAL
             RE(I)   = RE(I)/CANGS
             BETA(I) = BETA(I)*CANGS
C
C   Compute useful constants
C
             ZPO(I)   = 1.0D0 + Z(I)
             OP3Z(I)  = 1.0D0 + 3.0D0*Z(I)
             TOP3Z(I) = 2.0D0*OP3Z(I)
             ZP3(I)   = Z(I) + 3.0D0
             TZP3(I)  = 2.0D0*ZP3(I)
             DO4Z(I)  = DE(I)/4.0D0/ZPO(I)
             B(I)     = BETA(I)*DO4Z(I)*2.0D0
10    CONTINUE
C
C    Set the values of the classical energy in the three asymptotic valleys
C
             EASYAB = DE(1)
             EASYBC = DE(2)
             EASYAC = DE(3)
C
      RETURN
      END
C
C*****
C
         SUBROUTINE POT
C
C   System:          ABC
C   Functional form: Extended LEPS (London-Erying-Polyani-Sato)
C   References:      P. J. Kuntz, E. M. Nemth, J. C. Polanyi, S. D. Rosner,
C                    and C. E. Young 
C                    J. Chem. Phys. 44, 1168 (1966)
C
C   The potential parameters must be passed through the common blocks
C   PT31CM, PT35CM, PT32CM, PT34CM, SATOCM, and LEPSCM.  
C   All information passed through the common blocks PT31CM, PT35CM, 
C   SATOCM, and LEPSCM must be in Hartree atomic units.
C
C        For the reaction: A + BC -> AB + C we write:
C                          R(1) = R(A-B)
C                          R(2) = R(B-C)
C                          R(3) = R(C-A)
C
C   NOTE: The potential energy at the reactant asympotote, that is at 
C         A infinitely far from the BC diatomic, BC diatomic at its
C         equilibrium configuration, is set equal to zero.
C
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT31CM/ R(3), ENERGY, DEDR(3)
         COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
         COMMON /PT34CM/ IPRT
         COMMON /PT35CM/ EASYAB, EASYBC, EASYAC
         COMMON /SATOCM/ DE(3), RE(3), BETA(3), Z(3) 
         COMMON /LEPSCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3), 
     *                   DO4Z(3), B(3)
         DIMENSION X(3), COUL(3), EXCH(3)
         PARAMETER (R2 = 1.41421356D0)
C
C   Initialize the variable used in the calculation of the energy.
C
         ENERGY = 0.D0
C
C   Check the values of NSURF and NDER for validity.
C
         IF (NSURF .NE. 0) THEN
             WRITE (IPRT, 900) NSURF
             STOP 'POT 1'
         ENDIF
         IF (NDER .GT. 1) THEN
             WRITE (IPRT, 910) NDER
             STOP 'POT 2'
         ENDIF
C
C   Compute the energy.
C
         DO 10 I = 1,3
               X(I)    = EXP(-BETA(I)*(R(I)-RE(I)))
               COUL(I) = DO4Z(I)*(ZP3(I)*X(I)-TOP3Z(I))*X(I)
               EXCH(I) = DO4Z(I)*(OP3Z(I)*X(I)-TZP3(I))*X(I)
               ENERGY  = ENERGY + COUL(I)
10       CONTINUE
C
         RAD = SQRT((EXCH(1)-EXCH(2))**2 + (EXCH(2)-EXCH(3))**2 +
     *              (EXCH(3)-EXCH(1))**2)
C
         ENERGY = ENERGY - RAD/R2 + EASYBC
C
C   Compute the derivatives of the energy with respect to the internal
C   coordinates.
C
         IF (NDER .EQ. 1) THEN
             S = EXCH(1) + EXCH(2) + EXCH(3)
             DO 20 I = 1,3
                   DEDR(I) = B(I)*X(I)*((3.0D0*EXCH(I)-S)/R2*
     *                       (OP3Z(I)*X(I)-ZP3(I))/RAD-
     *                       ZP3(I)*X(I)+OP3Z(I))
20           CONTINUE
         ENDIF
C
900   FORMAT(/,2X,T5,'NSURF has been set equal to ',I5,
     *       /,2X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,2X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,
     *       /, 2X, 'This value of NDER is not allowed in this ',
     *              'version of the potential.')
C
      RETURN
      END
C
C*****
C
         BLOCK DATA PTPACM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
         COMMON /PT34CM/ IPRT
         COMMON /SATOCM/ DE(3), RE(3), BETA(3), Z(3) 
C
C   Initialize the flags and the I/O unit numbers for the potential
C
         DATA IPRT /6/
         DATA NDER /1/
         DATA NFLAG /20*0/ 
C
C   Initialize the potential parameters; the energy 
C   parameters are in kcal/mol, and the lengths are in Angstroms.
C
         DATA DE/ 106.447D0, 109.458D0, 106.447D0/
         DATA RE/ 1.2732D0, 0.74127D0, 1.2732D0/
         DATA BETA/ 1.8674D0, 1.94130D0, 1.8674D0/
         DATA Z/ 0.277D0, 0.077D0, 0.277D0/
C
         END
C
C*****
