C Calcula el valor X0 que hace que el valor del polinomio sea Y0 (utiliza el
C metodo de Newton, por lo que un valor inicial XINI es necesario)
        SUBROUTINE FPOLYINV(NDEG,COEFF,Y0,XINI,X0)
        IMPLICIT NONE
        INTEGER NDEG
        REAL COEFF(NDEG+1)
        REAL Y0
        REAL XINI
        REAL X0
C
        INTEGER NITERMAX
        PARAMETER (NITERMAX=500)
        REAL THRESHOLD
        PARAMETER (THRESHOLD=0.0001)
C
        REAL FPOLY
        REAL FDERI
C
        INTEGER ILOOP
        REAL X1,X2
        REAL FDER,FPOL
        LOGICAL LOOP
C------------------------------------------------------------------------------
        X1=XINI
        ILOOP=0
C
        LOOP=.TRUE.
        DO WHILE(LOOP)
          FPOL=FPOLY(NDEG,COEFF,X1)-Y0
          FDER=FDERI(NDEG,COEFF,X1)
          X2=X1-FPOL/FDER
          LOOP=(ABS(X2-X1).GT.THRESHOLD)
          IF(LOOP)THEN
            X1=X2
            ILOOP=ILOOP+1
            LOOP=(ILOOP.LT.NITERMAX)
          END IF
        END DO
C
        X0=X2
C
        END
