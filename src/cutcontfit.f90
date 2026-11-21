!------------------------------------------------------------------------------
! SUBROUTINE CUTCONTFIT(N,N1,N2,X,Y,XLIM,A0,A1,NCONT,RMSCONT,IERR)
!
! Fit a straight line, Y = A0 + A1*X, by least squares to the points of a
! 1-D profile within two continuum regions (used by CUTGAUSSFIT and
! PCUTGAUSS).
!
! Parameters
! ----------
! N : INTEGER
!     Dimension of X and Y.
! N1, N2 : INTEGER
!     Range of points of the profile to be used (X(N1:N2), Y(N1:N2)).
! X, Y : REAL(N)
!     Profile.
! XLIM : REAL(4)
!     Limits (in X units, in any order within each pair) of the first
!     (XLIM(1:2)) and second (XLIM(3:4)) continuum regions.
!
! Returns
! -------
! A0, A1 : REAL
!     Coefficients of the continuum, Y = A0 + A1*X (zero if IERR/=0).
! NCONT : INTEGER
!     Number of points used in the fit.
! RMSCONT : REAL
!     Rms of the residuals of the fit (with NCONT-2 degrees of freedom;
!     zero if NCONT=2 or IERR/=0).
! IERR : INTEGER
!     0: no error; 1: fewer than 2 points or all of them with the same X.
!------------------------------------------------------------------------------
        SUBROUTINE CUTCONTFIT(N,N1,N2,X,Y,XLIM,A0,A1,NCONT,RMSCONT,IERR)
        IMPLICIT NONE
        INTEGER N,N1,N2
        REAL X(N),Y(N)
        REAL XLIM(4)
        REAL A0,A1
        INTEGER NCONT
        REAL RMSCONT
        INTEGER IERR
!
        INTEGER I
        REAL XC1,XC2,XC3,XC4
        DOUBLE PRECISION SX,SY,SXX,SXY,XMEAN,DX,RES,SRES2
        LOGICAL LCONT
!------------------------------------------------------------------------------
        A0=0.0
        A1=0.0
        RMSCONT=0.0
        IERR=0
! limits of each region in increasing order
        XC1=MIN(XLIM(1),XLIM(2))
        XC2=MAX(XLIM(1),XLIM(2))
        XC3=MIN(XLIM(3),XLIM(4))
        XC4=MAX(XLIM(3),XLIM(4))
! least-squares straight line (X referred to its mean value to avoid
! rounding errors)
        NCONT=0
        SX=0.D0
        DO I=N1,N2
          LCONT=((X(I).GE.XC1).AND.(X(I).LE.XC2)).OR. &
                ((X(I).GE.XC3).AND.(X(I).LE.XC4))
          IF(LCONT)THEN
            NCONT=NCONT+1
            SX=SX+DBLE(X(I))
          END IF
        END DO
        IF(NCONT.LT.2)THEN
          IERR=1
          RETURN
        END IF
        XMEAN=SX/DBLE(NCONT)
        SY=0.D0
        SXX=0.D0
        SXY=0.D0
        DO I=N1,N2
          LCONT=((X(I).GE.XC1).AND.(X(I).LE.XC2)).OR. &
                ((X(I).GE.XC3).AND.(X(I).LE.XC4))
          IF(LCONT)THEN
            DX=DBLE(X(I))-XMEAN
            SY=SY+DBLE(Y(I))
            SXX=SXX+DX*DX
            SXY=SXY+DX*DBLE(Y(I))
          END IF
        END DO
        IF(SXX.EQ.0.D0)THEN
          IERR=1
          RETURN
        END IF
        A1=REAL(SXY/SXX)
        A0=REAL(SY/DBLE(NCONT)-(SXY/SXX)*XMEAN)
        IF(NCONT.GT.2)THEN
          SRES2=0.D0
          DO I=N1,N2
            LCONT=((X(I).GE.XC1).AND.(X(I).LE.XC2)).OR. &
                  ((X(I).GE.XC3).AND.(X(I).LE.XC4))
            IF(LCONT)THEN
              RES=DBLE(Y(I))-(DBLE(A0)+DBLE(A1)*DBLE(X(I)))
              SRES2=SRES2+RES*RES
            END IF
          END DO
          RMSCONT=REAL(DSQRT(SRES2/DBLE(NCONT-2)))
        END IF
        END
