!------------------------------------------------------------------------------
!omment
!
! REAL FUNCTION FPOLY(NDEG,COEFF,X)
!
! Input: NDEG,COEFF,X
! Output: FPOLY (function)
!
! Evaluate the polynomial of degree NDEG and coefficients COEFF at X.
!
! INTEGER NDEG -> polynomial degree
! REAL    COEFF(NDEG+1) -> polynomial coefficients
! REAL    X -> abscissa at which the polynomial is going to be evaluated
!
!omment
!------------------------------------------------------------------------------
        REAL FUNCTION FPOLY(NDEG,COEFF,X)
        IMPLICIT NONE
        INTEGER NDEG
        REAL COEFF(NDEG+1)
        REAL X
!
        INTEGER K
        DOUBLE PRECISION DSUM
!------------------------------------------------------------------------------
        DSUM=DBLE(COEFF(NDEG+1))
        IF(NDEG.GT.0)THEN
          DO K=NDEG,1,-1
            DSUM=DSUM*DBLE(X)+DBLE(COEFF(K))
          END DO
        END IF
!
        FPOLY=REAL(DSUM)
        END
