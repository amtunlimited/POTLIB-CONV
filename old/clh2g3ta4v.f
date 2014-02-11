c This is the scalar version of the clh2g3 multiple electronic
c potential energy surface code.  TCA, 16 August 1994.
c
      subroutine prepot

      implicit double precision (a-h,o-z)

      logical usesep
c
c If the parameter usesep = .true. the adiabatic surfaces will be separated
c by (wvnsep) cm**-1 at the saddle point.  If usesep = .false. the surfaces
c will be separated by their normal value.
c
      parameter (usesep=.true.)
      parameter (wvnsep=100.d0)

      dimension de(3),re(3),b(3),s(3),pf(3),spl(3),r(nt,3),e(nt)

      data de     / 106.447d0, 109.458d0, 106.447d0 /
      data re     / 1.2732d0, 0.74127d0, 1.2732d0 /
      data b      / 1.8674d0, 1.9413d0, 1.8674d0 /
      data s      / 0.1835d0, 0.167d0, 0.1835d0 /
      data rj     / 0.0758016022d0 /
      data ap     / 0.0008627355d0 /
      data at     / 0.2981969860d0 /
      data bp     / 0.1439106543d0 /
      data q2     / 0.6940323070d0 /
      data q4     / 1.6963092005d0 /
      data ckcau  / 627.5095d0 /
      data cangau / 0.529177249d0 /
      data c1     / 3.1053877618397071175758280546d+0 /
      data c2     / -1.8942608491155350536449828878d+0 /
      data spl    / 2.6477d0, 1.8701d0, 4.5178d0 /
      save

      do 10 i=1,3
        de(i)=de(i)/ckcau
        re(i)=re(i)/cangau
        b(i)=b(i)*cangau
        pf(i)=0.5d0*de(i)*((1.d0-s(i))/(1.d0+s(i)))
10    continue

      esep=136.748116495694973d0/ckcau
      if (usesep) then
        zl=(esep-wvnsep*2.859144d-3/ckcau)/esep
      else
        zl=0.d0
      end if

      sqrt3=dsqrt(3.d0)

      write(6,*) 'Potential code called for ClH2 surface TA4.'
      write(6,*) 'Last modified 16 August 1994 - Version 1.1'
      if (usesep) write(6,*) 'Surfaces separated by ',wvnsep,' cm**-1.'
      write(6,*)

      return

      entry pot(r,e,nt,nsurf)

      do 20 i=1,nt

      f1=dexp(-ap*(r(i,1)+r(i,3))**4)
      f2=dexp(-at*(r(i,1)-r(i,3))**2)
      f3=dexp(-bp*(r(i,1)+r(i,3)-r(i,2))**2)
      ctheta=(r(i,1)**2+r(i,3)**2-r(i,2)**2)/(2.d0*r(i,1)*r(i,3))
      if (ctheta .lt. -1.d0) then
        ctheta=-1.d0
      else if (ctheta .gt. 1.d0) then
        ctheta=1.d0
      end if
      stheta=dsin(dacos(ctheta))
      f4=1.d0+q2*stheta**2+q4*stheta**4
      v3b=f4*(rj*f1*f2*f3)

      x=(r(i,1)-r(i,3))/r(i,2)
      g=(0.5d0*(1.d0+ctheta))**6
      fmod=1.d0+c1*g*(1.d0-x**2)+c2*g*(1.d0-x**4)

      y1=dexp(-b(1)*(r(i,1)-re(1)))
      vs1=de(1)*(y1-2.d0)*y1
      vt1=pf(1)*(fmod*y1+2.d0)*y1
      y2=dexp(-b(2)*(r(i,2)-re(2)))
      vs2=de(2)*(y2-2.d0)*y2
      vt2=pf(2)*(y2+2.d0)*y2
      y3=dexp(-b(3)*(r(i,3)-re(3)))
      vs3=de(3)*(y3-2.d0)*y3
      vt3=pf(3)*(fmod*y3+2.d0)*y3
c
c Calculate original Hamiltonian matrix and obtain eigenvalues and
c eigenvectors
c
      h11=de(2)+vs2+0.25d0*(vs1+vs3)+0.75d0*(vt1+vt3)+v3b
      h22=de(2)+0.75d0*(vs1+vs3)+vt2+0.25d0*(vt1+vt3)+v3b
      h12=sqrt3*0.25d0*(vs1-vs3+vt3-vt1)
      h12p=-dsqrt(h12**4/((1.227d-2*dexp(-2.d0*b(2)*r(i,2))**2)+h12**2))
      call dlaev2(h11,h12p,h22,rt1,rt2,cs1,sn1)
c
c Calculate new value of e2 and backtransform to new Hamiltonian matrix
c
      z=rt1-rt2
      cfd=2.d0*r(i,1)**2+2.d0*r(i,3)**2-r(i,2)**2
      if (cfd .lt. 0.d0) cfd=0.d0
      coschi=(r(i,3)**2-r(i,1)**2)/(r(i,2)*dsqrt(cfd))
      if (coschi .gt. 1.d0) coschi=1.d0
      if (coschi .lt. -1.d0) coschi=-1.d0
      sinchi=dsqrt(1.d0-coschi**2)
      gmsp=esep*dexp(1.d0-esep)
      gm=z*dexp(1.d0-z)
      e2new=rt1-(1.d0-sinchi**4)*esep*zl*(2.d0*gmsp*gm/(gmsp**2+gm**2))
c    &      *dexp(-(b(1)*(r(i,1)+r(i,3))+b(2)*r(i,2)))

      tp11=rt2*sn1*sn1+e2new*cs1*cs1
      tp12=(e2new-rt2)*sn1*cs1
      tp22=rt2*cs1*cs1+e2new*sn1*sn1
c
c Calculate scaling factor for coupling and transform matrix
c
      chi=0.5d0*(1.d0-dtanh(b(2)*(r(i,2)-5.3d0)))
      if (chi .lt. 0.d0) chi=0.d0
      if (chi .gt. 1.d0) chi=1.d0
      ab=tp11-tp22
      ac=ab*ab+4.d0*(1.d0+chi)*tp12*tp12
      ad=ab*dsqrt(ab*ab+4.d0*(1.d0-chi*chi)*tp12*tp12)
      ae1=ac+ad
      if (ae1 .lt. 0.d0) ae1=0.d0
      ae2=ac-ad
      if (ae2 .lt. 0.d0) ae2=0.d0
      af=8.d0*tp12*tp12+2.d0*ab*ab
      ax1=dsqrt(ae1/af)
      if (ax1 .gt. 1.d0) ax1=1.d0
      ax2=dsqrt(ae2/af)
      if (ax2 .gt. 1.d0) ax2=1.d0
      ax=max(ax1,ax2)
      if (tp12 .ge. 0.d0) then
        ay=-dsqrt(1.d0-ax*ax)
      else
        ay=dsqrt(1.d0-ax*ax)
      end if
c
c Return value of appropriate surface
c
      if (nsurf .eq. 1) then
        e(i)=tp11*ax*ax-2.d0*tp12*ax*ay+tp22*ay*ay
      else if (nsurf .eq. 2) then
        e(i)=(tp11-tp22)*ax*ay+tp12*(ax*ax-ay*ay)
      else
        e(i)=tp11*ay*ay+2.d0*tp12*ax*ay+tp22*ax*ax
      end if

20    continue

      return

      end

      SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
*     ..
*
*  Purpose
*  =======
*
*  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
*     [  A   B  ]
*     [  B   C  ].
*  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
*  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
*  eigenvector for RT1, giving the decomposition
*
*     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
*     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*          The (1,1) element of the 2-by-2 matrix.
*
*  B       (input) DOUBLE PRECISION
*          The (1,2) element and the conjugate of the (2,1) element of
*          the 2-by-2 matrix.
*
*  C       (input) DOUBLE PRECISION
*          The (2,2) element of the 2-by-2 matrix.
*
*  RT1     (output) DOUBLE PRECISION
*          The eigenvalue of larger absolute value.
*
*  RT2     (output) DOUBLE PRECISION
*          The eigenvalue of smaller absolute value.
*
*  CS1     (output) DOUBLE PRECISION
*  SN1     (output) DOUBLE PRECISION
*          The vector (CS1, SN1) is a unit right eigenvector for RT1.
*
*  Further Details
*  ===============
*
*  RT1 is accurate to a few ulps barring over/underflow.
*
*  RT2 may be inaccurate if there is massive cancellation in the
*  determinant A*C-B*B; higher precision or correctly rounded or
*  correctly truncated arithmetic would be needed to compute RT2
*  accurately in all cases.
*
*  CS1 and SN1 are accurate to a few ulps barring over/underflow.
*
*  Overflow is possible only if RT1 is within a factor of 5 of overflow.
*  Underflow is harmless if the input data is 0 or exceeds
*     underflow_threshold / macheps.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            SGN1, SGN2
      DOUBLE PRECISION   AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,
     $                   TB, TN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Compute the eigenvalues
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        Includes case AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
         SGN1 = -1
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        Includes case RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF
*
*     Compute the eigenvector
*
      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
*
*     End of DLAEV2
*
      END


