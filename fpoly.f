C------------------------------------------------------------------------------
Comment
C
C REAL FUNCTION FPOLY(NDEG,COEFF,X)
C
C Input: NDEG,COEFF,X
C Output: FPOLY (function)
C
C Evaluate the polynomial of degree NDEG and coefficients COEFF at X.
C
C INTEGER NDEG -> polynomial degree
C REAL    COEFF(NDEG+1) -> polynomial coefficients
C REAL    X -> abscissa at which the polynomial is going to be evaluated
C
Comment
C------------------------------------------------------------------------------
        REAL FUNCTION FPOLY(NDEG,COEFF,X)
        IMPLICIT NONE
        INTEGER NDEG
        REAL COEFF(NDEG+1)
        REAL X
C
        INTEGER K
        DOUBLE PRECISION DSUM
C------------------------------------------------------------------------------
        DSUM=DBLE(COEFF(NDEG+1))
        IF(NDEG.GT.0)THEN
          DO K=NDEG,1,-1
            DSUM=DSUM*DBLE(X)+DBLE(COEFF(K))
          END DO
        END IF
C
        FPOLY=REAL(DSUM)
        END
