!------------------------------------------------------------------------------
! SUBROUTINE FINDMMLX(N,N1,N2,X,Y,X1,X2,YMIN,YMAX,LFOUND)
!
! Return the minimum and maximum of Y(I), with I in the range [N1,N2], for
! the points whose abscissa X(I) lies within the interval [X1,X2] (limits
! included).
!
! Parameters
! ----------
! N : INTEGER
!     Number of elements of X and Y.
! N1, N2 : INTEGER
!     First and last elements (1-based) of X and Y to be considered.
! X, Y : REAL arrays (N)
!     Abscissae and ordinates of the points.
! X1, X2 : REAL
!     Interval of X (X1 <= X2) where the points are considered.
!
! Returns
! -------
! YMIN, YMAX : REAL
!     Minimum and maximum of the selected Y values (unchanged if LFOUND is
!     .FALSE.).
! LFOUND : LOGICAL
!     .TRUE. if at least one point has been found within [X1,X2].
!------------------------------------------------------------------------------
        SUBROUTINE FINDMMLX(N,N1,N2,X,Y,X1,X2,YMIN,YMAX,LFOUND)
        IMPLICIT NONE
        INTEGER N,N1,N2
        REAL X(N),Y(N)
        REAL X1,X2
        REAL YMIN,YMAX
        LOGICAL LFOUND
!
        INTEGER I
!------------------------------------------------------------------------------
        LFOUND=.FALSE.
        DO I=MAX(N1,1),MIN(N2,N)
          IF((X(I).GE.X1).AND.(X(I).LE.X2))THEN
            IF(LFOUND)THEN
              IF(Y(I).LT.YMIN) YMIN=Y(I)
              IF(Y(I).GT.YMAX) YMAX=Y(I)
            ELSE
              YMIN=Y(I)
              YMAX=Y(I)
              LFOUND=.TRUE.
            END IF
          END IF
        END DO
        END
