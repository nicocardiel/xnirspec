!------------------------------------------------------------------------------
! SUBROUTINE CUTGAUSSFIT(N,N1,N2,X,Y,XLIM,YRMSTOL,NSIMUL,
!                        A0,A1,NCONT,RMSCONT,
!                        X0,SIGMA,AMP,EX0,ESIGMA,EAMP,LERRSIM,NLINE,IERR)
!
! Fit a gaussian to a line in a 1-D profile after subtracting a linear
! continuum. The continuum is a straight line, Y = A0 + A1*X, fitted by least
! squares (CUTCONTFIT) to the points of the profile within two continuum
! regions (left and right of the line). The gaussian, Y - (A0 + A1*X) =
! AMP*EXP[-(X-X0)^2/(2*SIGMA^2)], is fitted (GAUSSFIT) to the continuum-
! subtracted points within the line region.
!
! Parameters
! ----------
! N : INTEGER
!     Dimension of X and Y.
! N1, N2 : INTEGER
!     Range of points of the profile to be used (X(N1:N2), Y(N1:N2)).
! X, Y : REAL(N)
!     Profile.
! XLIM : REAL(6)
!     Limits (in X units, in any order within each pair) of the left
!     continuum region (XLIM(1:2)), the right continuum region (XLIM(3:4))
!     and the line region (XLIM(5:6)).
! YRMSTOL : REAL
!     Stopping criterion for DOWNHILL (used by GAUSSFIT).
! NSIMUL : INTEGER
!     Number of simulations to compute the errors of the gaussian
!     parameters, adding to the line points gaussian noise with standard
!     deviation RMSCONT. If NSIMUL<2 or RMSCONT=0, the errors are the rms of
!     the final simplex of DOWNHILL instead.
!
! Returns
! -------
! A0, A1 : REAL
!     Coefficients of the continuum, Y = A0 + A1*X.
! NCONT : INTEGER
!     Number of points used to fit the continuum.
! RMSCONT : REAL
!     Rms of the residuals of the continuum fit (with NCONT-2 degrees of
!     freedom; zero if NCONT=2).
! X0, SIGMA, AMP : REAL
!     Center, sigma (positive) and amplitude (negative for absorption lines)
!     of the gaussian.
! EX0, ESIGMA, EAMP : REAL
!     Errors of X0, SIGMA and AMP.
! LERRSIM : LOGICAL
!     .TRUE. if the errors come from simulations, .FALSE. if they are the
!     rms of the final simplex of DOWNHILL.
! NLINE : INTEGER
!     Number of points used to fit the gaussian.
! IERR : INTEGER
!     0: no error; 1: fewer than 2 continuum points or all of them with the
!     same X; 2: fewer than 3 points in the line region; 3: more than NXYMAX
!     points in the line region. If IERR/=0 the gaussian parameters are
!     not computed (set to zero).
!------------------------------------------------------------------------------
        SUBROUTINE CUTGAUSSFIT(N,N1,N2,X,Y,XLIM,YRMSTOL,NSIMUL, &
         A0,A1,NCONT,RMSCONT, &
         X0,SIGMA,AMP,EX0,ESIGMA,EAMP,LERRSIM,NLINE,IERR)
        IMPLICIT NONE
        INTEGER N,N1,N2
        REAL X(N),Y(N)
        REAL XLIM(6)
        REAL YRMSTOL
        INTEGER NSIMUL
        REAL A0,A1
        INTEGER NCONT
        REAL RMSCONT
        REAL X0,SIGMA,AMP,EX0,ESIGMA,EAMP
        LOGICAL LERRSIM
        INTEGER NLINE
        INTEGER IERR
!
        INCLUDE 'largest.inc'
!
        INTEGER I,NSIM
        REAL XL1,XL2
        REAL EEX0,EESIGMA,EEAMP
        REAL, SAVE :: XF(NXYMAX),YF(NXYMAX),EYF(NXYMAX)
!------------------------------------------------------------------------------
        A0=0.0
        A1=0.0
        RMSCONT=0.0
        X0=0.0
        SIGMA=0.0
        AMP=0.0
        EX0=0.0
        ESIGMA=0.0
        EAMP=0.0
        LERRSIM=.FALSE.
        NLINE=0
        IERR=0
! limits of the line region in increasing order
        XL1=MIN(XLIM(5),XLIM(6))
        XL2=MAX(XLIM(5),XLIM(6))
!------------------------------------------------------------------------------
! continuum: least-squares straight line
        CALL CUTCONTFIT(N,N1,N2,X,Y,XLIM(1:4),A0,A1,NCONT,RMSCONT,IERR)
        IF(IERR.NE.0) RETURN
!------------------------------------------------------------------------------
! line: continuum-subtracted points within the line region
        DO I=N1,N2
          IF((X(I).GE.XL1).AND.(X(I).LE.XL2))THEN
            NLINE=NLINE+1
            IF(NLINE.GT.NXYMAX)THEN
              IERR=3
              RETURN
            END IF
            XF(NLINE)=X(I)
            YF(NLINE)=Y(I)-(A0+A1*X(I))
            EYF(NLINE)=RMSCONT
          END IF
        END DO
        IF(NLINE.LT.3)THEN
          IERR=2
          RETURN
        END IF
        IF((NSIMUL.GE.2).AND.(RMSCONT.GT.0.0))THEN
          NSIM=NSIMUL
        ELSE
          NSIM=0
        END IF
        CALL GAUSSFIT(NLINE,XF,YF,EYF,X0,SIGMA,AMP, &
         EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP,YRMSTOL,NSIM)
        SIGMA=ABS(SIGMA)           !the fitted function only depends on SIGMA^2
        LERRSIM=(NSIM.GE.2)
        IF(.NOT.LERRSIM)THEN
          EX0=EEX0
          ESIGMA=EESIGMA
          EAMP=EEAMP
        END IF
        END
