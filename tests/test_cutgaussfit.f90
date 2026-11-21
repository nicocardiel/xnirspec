!------------------------------------------------------------------------------
! Unit tests for CUTCONTFIT (linear continuum fitted by least squares to two
! continuum regions) and CUTGAUSSFIT (gaussian fitted to the profile minus
! the continuum in the line region).
!------------------------------------------------------------------------------
        PROGRAM TEST_CUTGAUSSFIT
        USE TEST_UTILS
        IMPLICIT NONE
!
        INTEGER N
        PARAMETER (N=100)
        REAL YRMSTOL
        PARAMETER (YRMSTOL=1.E-5)
!
        INTEGER I
        INTEGER NCONT,NLINE,IERR
        REAL X(N),Y(N)
        REAL XLIM(6)
        REAL A0,A1,RMSCONT
        REAL X0,SIGMA,AMP,EX0,ESIGMA,EAMP
        LOGICAL LERRSIM
!------------------------------------------------------------------------------
! emission line on a linear continuum, without noise:
! y = 10 + 0.5*x + 50*exp[-(x-40)^2/(2*3^2)]
        DO I=1,N
          X(I)=REAL(I)
          Y(I)=10.0+0.5*X(I)+50.0*EXP(-(X(I)-40.0)**2/(2.0*9.0))
        END DO
! the limits of each region can be given in any order
        XLIM=(/25.0,5.0,55.0,90.0,50.0,30.0/)
        CALL CUTGAUSSFIT(N,1,N,X,Y,XLIM,YRMSTOL,0,A0,A1,NCONT,RMSCONT, &
         X0,SIGMA,AMP,EX0,ESIGMA,EAMP,LERRSIM,NLINE,IERR)
        CALL CHECK_INT('emission: IERR',IERR,0)
        CALL CHECK_INT('emission: continuum points',NCONT,21+36)
        CALL CHECK_INT('emission: line points',NLINE,21)
        CALL CHECK_REAL('emission: a0',A0,10.0,1.E-4,1.E-4)
        CALL CHECK_REAL('emission: a1',A1,0.5,1.E-4,1.E-5)
        CALL CHECK_REAL('emission: rms continuum',RMSCONT,0.0,0.0,1.E-4)
        CALL CHECK_REAL('emission: x0',X0,40.0,0.0,1.E-2)
        CALL CHECK_REAL('emission: sigma',SIGMA,3.0,1.E-2,0.0)
        CALL CHECK_REAL('emission: amp',AMP,50.0,1.E-2,0.0)
        CALL CHECK_TRUE('emission: errors from DOWNHILL (NSIMUL=0)', &
         .NOT.LERRSIM)
!------------------------------------------------------------------------------
! absorption line on a decreasing continuum:
! y = 200 - 1.2*x - 30*exp[-(x-70.3)^2/(2*2.5^2)]
        DO I=1,N
          Y(I)=200.0-1.2*X(I)-30.0*EXP(-(X(I)-70.3)**2/(2.0*6.25))
        END DO
        XLIM=(/45.0,58.0,83.0,98.0,62.0,79.0/)
        CALL CUTGAUSSFIT(N,1,N,X,Y,XLIM,YRMSTOL,0,A0,A1,NCONT,RMSCONT, &
         X0,SIGMA,AMP,EX0,ESIGMA,EAMP,LERRSIM,NLINE,IERR)
        CALL CHECK_INT('absorption: IERR',IERR,0)
        CALL CHECK_REAL('absorption: a0',A0,200.0,1.E-4,0.0)
        CALL CHECK_REAL('absorption: a1',A1,-1.2,1.E-4,0.0)
        CALL CHECK_REAL('absorption: x0',X0,70.3,0.0,1.E-2)
        CALL CHECK_REAL('absorption: sigma',SIGMA,2.5,1.E-2,0.0)
        CALL CHECK_REAL('absorption: amp (negative)',AMP,-30.0,1.E-2,0.0)
!------------------------------------------------------------------------------
! with (deterministic) noise of amplitude 0.5: the continuum rms is not zero,
! so the errors are computed with simulations and must be positive
        DO I=1,N
          Y(I)=10.0+0.5*X(I)+50.0*EXP(-(X(I)-40.0)**2/(2.0*9.0))+ &
           0.5*SIN(7.3*REAL(I))
        END DO
        XLIM=(/5.0,25.0,55.0,90.0,30.0,50.0/)
        CALL CUTGAUSSFIT(N,1,N,X,Y,XLIM,YRMSTOL,50,A0,A1,NCONT,RMSCONT, &
         X0,SIGMA,AMP,EX0,ESIGMA,EAMP,LERRSIM,NLINE,IERR)
        CALL CHECK_INT('noise: IERR',IERR,0)
        CALL CHECK_TRUE('noise: rms continuum > 0',RMSCONT.GT.0.0)
        CALL CHECK_TRUE('noise: errors from simulations',LERRSIM)
        CALL CHECK_TRUE('noise: errors > 0', &
         (EX0.GT.0.0).AND.(ESIGMA.GT.0.0).AND.(EAMP.GT.0.0))
        CALL CHECK_REAL('noise: x0',X0,40.0,0.0,0.2)
        CALL CHECK_REAL('noise: sigma',SIGMA,3.0,0.05,0.0)
        CALL CHECK_REAL('noise: amp',AMP,50.0,0.05,0.0)
!------------------------------------------------------------------------------
! using only part of the profile (N1:N2): the right continuum region lies
! outside, so only the left one is used
        DO I=1,N
          Y(I)=10.0+0.5*X(I)+50.0*EXP(-(X(I)-40.0)**2/(2.0*9.0))
        END DO
        XLIM=(/5.0,25.0,55.0,90.0,30.0,50.0/)
        CALL CUTGAUSSFIT(N,1,54,X,Y,XLIM,YRMSTOL,0,A0,A1,NCONT,RMSCONT, &
         X0,SIGMA,AMP,EX0,ESIGMA,EAMP,LERRSIM,NLINE,IERR)
        CALL CHECK_INT('N1:N2: IERR',IERR,0)
        CALL CHECK_INT('N1:N2: continuum points (left region only)', &
         NCONT,21)
        CALL CHECK_REAL('N1:N2: x0',X0,40.0,0.0,1.E-2)
!------------------------------------------------------------------------------
! CUTCONTFIT alone: exact straight line plus a line outside the continuum
! regions, which must not affect the fit
        CALL CUTCONTFIT(N,1,N,X,Y,(/90.0,55.0,25.0,5.0/),A0,A1,NCONT, &
         RMSCONT,IERR)
        CALL CHECK_INT('CUTCONTFIT: IERR',IERR,0)
        CALL CHECK_INT('CUTCONTFIT: points',NCONT,21+36)
        CALL CHECK_REAL('CUTCONTFIT: a0',A0,10.0,1.E-4,1.E-4)
        CALL CHECK_REAL('CUTCONTFIT: a1',A1,0.5,1.E-4,1.E-5)
! two points: exact line through them, rms set to zero
        CALL CUTCONTFIT(N,1,N,X,Y,(/10.0,10.0,20.0,20.0/),A0,A1,NCONT, &
         RMSCONT,IERR)
        CALL CHECK_INT('CUTCONTFIT 2 points: IERR',IERR,0)
        CALL CHECK_INT('CUTCONTFIT 2 points: points',NCONT,2)
        CALL CHECK_REAL('CUTCONTFIT 2 points: a1',A1,0.5,1.E-4,1.E-5)
        CALL CHECK_REAL('CUTCONTFIT 2 points: rms',RMSCONT,0.0,0.0,0.0)
! a single point (twice the same X): error
        CALL CUTCONTFIT(N,1,N,X,Y,(/10.0,10.0,10.0,10.0/),A0,A1,NCONT, &
         RMSCONT,IERR)
        CALL CHECK_INT('CUTCONTFIT 1 point: IERR=1',IERR,1)
!------------------------------------------------------------------------------
! errors: no continuum points, and fewer than 3 points in the line region
        XLIM=(/200.0,210.0,300.0,310.0,30.0,50.0/)
        CALL CUTGAUSSFIT(N,1,N,X,Y,XLIM,YRMSTOL,0,A0,A1,NCONT,RMSCONT, &
         X0,SIGMA,AMP,EX0,ESIGMA,EAMP,LERRSIM,NLINE,IERR)
        CALL CHECK_INT('no continuum points: IERR=1',IERR,1)
        XLIM=(/5.0,25.0,55.0,90.0,39.5,41.5/)
        CALL CUTGAUSSFIT(N,1,N,X,Y,XLIM,YRMSTOL,0,A0,A1,NCONT,RMSCONT, &
         X0,SIGMA,AMP,EX0,ESIGMA,EAMP,LERRSIM,NLINE,IERR)
        CALL CHECK_INT('2 points in line region: IERR=2',IERR,2)
        CALL CHECK_REAL('IERR=2: amp set to zero',AMP,0.0,0.0,0.0)
!
        CALL TEST_SUMMARY()
        END PROGRAM TEST_CUTGAUSSFIT
