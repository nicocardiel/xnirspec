! Calcula la distancia de un punto a una recta, devolviendo ademas las
! coordenadas del punto mas cercano sobre la recta.
! XC,YC: coordenadas de un punto generico C
! XA,YA, XB,YB: coordenadas de dos puntos A y B que definen la recta
! D: distancia del punto C a la recta que pasa por A y B
! X0,Y0: coordenadas del punto de la recta mas cercano al punto C
!------------------------------------------------------------------------------
        SUBROUTINE DISTPR(XC,YC,XA,YA,XB,YB,D,X0,Y0)
        IMPLICIT NONE
        REAL XC,YC
        REAL XA,YA
        REAL XB,YB
        REAL D
        REAL X0,Y0
!
        REAL L,R,S
!------------------------------------------------------------------------------
        L=((XB-XA)*(XB-XA)+(YB-YA)*(YB-YA)) !distancia AB al cuadrado
        IF(L.EQ.0.0)THEN
          WRITE(*,101) '***FATAL ERROR***'
          WRITE(*,101) '=> A=B in subroutine DISTPR'
          INCLUDE 'deallocate_arrays.inc'
          STOP
        END IF
        R=((YA-YC)*(YA-YB)-(XA-XC)*(XB-XA))/L
        S=((YA-YC)*(XB-XA)-(XA-XC)*(YB-YA))/L
        X0=XA+R*(XB-XA)
        Y0=YA+R*(YB-YA)
        D=S*SQRT(L) !distance from C to the line
!------------------------------------------------------------------------------
101     FORMAT(A)
        END
