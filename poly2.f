C Calcula los coeficientes del polinomio de segundo grado que pasa por
C los puntos (X1,Y1), (X2,Y2) y (X3,Y3).
        SUBROUTINE POLY2(X1,Y1,X2,Y2,X3,Y3,COEFF)
        IMPLICIT NONE
        REAL X1,Y1
        REAL X2,Y2
        REAL X3,Y3
        REAL COEFF(3)
C------------------------------------------------------------------------------
        IF(X1.EQ.X2)THEN
          WRITE(*,101) '***FATAL ERROR***'
          WRITE(*,100) '=> X1,Y1:'
          WRITE(*,*) X1,Y1
          WRITE(*,100) '=> X2,Y2:'
          WRITE(*,*) X2,Y2
          WRITE(*,100) '=> X3,Y3:'
          WRITE(*,*) X3,Y3
          WRITE(*,101) '=> X1=X2 in subroutine POLY2'
          STOP
        END IF
C
        COEFF(3)=((X2-X1)*(Y3-Y1)-(Y2-Y1)*(X3-X1))/
     +   ((X3*X3-X1*X1)*(X2-X1)-(X2*X2-X1*X1)*(X3-X1))
        COEFF(2)=((Y2-Y1)-COEFF(3)*(X2*X2-X1*X1))/(X2-X1)
        COEFF(1)=Y1-COEFF(2)*X1-COEFF(3)*X1*X1
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
