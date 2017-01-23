C Calcula la derivada del polinomio en el valor de X
        REAL FUNCTION FDERI(NDEG,COEFF,X)
        IMPLICIT NONE
        INTEGER NDEG
        REAL COEFF(NDEG+1)
        REAL X
C
        INTEGER K
C------------------------------------------------------------------------------
        FDERI=COEFF(NDEG+1)*REAL(NDEG)
        IF(NDEG.GT.1)THEN
          DO K=NDEG,2,-1
            FDERI=FDERI*X+COEFF(K)*REAL(K-1)
          END DO
        END IF
C
        END
