! Calcula la longitud de arco de un polinomio. Para una constante y una recta
! no hace falta hacer nada. Para polinomios de segundo grado o mayor, es
! necesario hacer una integral (aqui se resuelve por el metodo de Simpson).
        REAL FUNCTION ARCLENGTH(NDEG,COEFF,XMIN,XMAX,NINTERVAL)
        IMPLICIT NONE
        INTEGER NDEG
        REAL COEFF(NDEG+1)
        REAL XMIN,XMAX
        INTEGER NINTERVAL     !no. of intervals for the extended Simpson's rule
!
        INTEGER I,K
        REAL S
        REAL SUM1,SUM2
        REAL FDER
        REAL X,Y1,Y2
!------------------------------------------------------------------------------
        IF(NDEG.EQ.0)THEN
          S=(XMAX-XMIN)
!..............................................................................
        ELSEIF(NDEG.EQ.1)THEN
          S=(XMAX-XMIN)*SQRT(1.+COEFF(2)*COEFF(2))
!..............................................................................
        ELSE
          SUM1=0.
          DO I=1,NINTERVAL-1,2
            X=XMIN+(XMAX-XMIN)*REAL(I)/REAL(NINTERVAL)
            FDER=COEFF(NDEG+1)*REAL(NDEG)
            DO K=NDEG,2,-1
              FDER=FDER*X+COEFF(K)*REAL(K-1)
            END DO
            SUM1=SUM1+SQRT(1.+FDER*FDER)
          END DO
          SUM2=0.
          DO I=2,NINTERVAL-2,2
            X=XMIN+(XMAX-XMIN)*REAL(I)/REAL(NINTERVAL)
            FDER=COEFF(NDEG+1)*REAL(NDEG)
            DO K=NDEG,2,-1
              FDER=FDER*X+COEFF(K)*REAL(K-1)
            END DO
            SUM2=SUM2+SQRT(1.+FDER*FDER)
          END DO
          X=XMIN
          FDER=COEFF(NDEG+1)*REAL(NDEG)
          DO K=NDEG,2,-1
            FDER=FDER*X+COEFF(K)*REAL(K-1)
          END DO
          Y1=SQRT(1.+FDER*FDER)
          X=XMAX
          FDER=COEFF(NDEG+1)*REAL(NDEG)
          DO K=NDEG,2,-1
            FDER=FDER*X+COEFF(K)*REAL(K-1)
          END DO
          Y2=SQRT(1.+FDER*FDER)
          S=(Y1+Y2+4.*SUM1+2.*SUM2)*(XMAX-XMIN)/REAL(3*NINTERVAL)
!..............................................................................
        END IF
!------------------------------------------------------------------------------
        ARCLENGTH=S
!
        END
