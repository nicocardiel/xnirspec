! Calcula la interseccion de dos polinomios y=f(x) y x=g(y). Para ello: 
! - parte del valor medio del recorrido de x en f(x), cacula y=f(x)
! - con el valor anterior calcular x=g(y)
! - con el valor anterior calcular y=f(x)
! y asi sucesivamente hasta alcanzar la precision deseada
!------------------------------------------------------------------------------
        SUBROUTINE INTERSEC(NDEGBA,COEFFBA,NDEGBL,COEFFBL,X0,XFIN,YFIN)
        IMPLICIT NONE
        INTEGER NDEGBA
        REAL COEFFBA(NDEGBA+1)                                !polinomio x=g(y)
        INTEGER NDEGBL
        REAL COEFFBL(NDEGBL+1)                                !polinomio y=f(x)
        REAL X0                      !valor inicial de comienzo de la iteracion
        REAL XFIN,YFIN                                                !solucion
!
        INTEGER NITERMAX
        PARAMETER (NITERMAX=500)
        REAL THRESHOLD
        PARAMETER (THRESHOLD=0.0001)
!
        REAL FPOLY
!
        INTEGER ILOOP
        REAL X1,Y1,X2,Y2
        LOGICAL LOOP
!------------------------------------------------------------------------------
        Y1=FPOLY(NDEGBL,COEFFBL,X0)
        X1=FPOLY(NDEGBA,COEFFBA,Y1)
        LOOP=.TRUE.
        ILOOP=0
!
        DO WHILE(LOOP)
          Y2=FPOLY(NDEGBL,COEFFBL,X1)
          X2=FPOLY(NDEGBA,COEFFBA,Y2)
          LOOP=(ABS(X2-X1).GT.THRESHOLD).OR.(ABS(Y2-Y1).GT.THRESHOLD)
          IF(LOOP)THEN
            X1=X2
            Y1=Y2
            ILOOP=ILOOP+1
            LOOP=(ILOOP.LT.NITERMAX)
          END IF
        END DO
!
        XFIN=X2
        YFIN=Y2
!
        END
