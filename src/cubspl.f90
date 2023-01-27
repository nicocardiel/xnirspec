!------------------------------------------------------------------------------
! Version 18-June-1998                                           File: cubspl.f
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This routine is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!omment
!
! SUBROUTINE CUBSPL(X,Y,N,IMODE,S,A,B,C)
!
! Input: X,Y,N,IMODE,S
! Output: A,B,C,S
!
! This subroutine computes the coefficients of a cubic spline. See C.F. Gerald
! and P. O. Wheatley, in Applied Numerical Analysis, 4th edition, pag. 207.
! The subroutine returns the spline coefficients, where the spline defined
! in the interval between X(I),Y(I) and X(I+1),Y(I+1) is given by:
!
!      Y = A(I)*(X-X(I))**3 + B(I)*(X-X(I))**2 + C(I)*(X-X(I)) + D(I)
!
! REAL X(N) -> X-values to be fitted
! REAL Y(N) -> Y-values to be fitted
! INTEGER N -> number of data points
! INTEGER IMODE -> End conditions mode: if S(I) represent the second derivative
!                  at the point X(I),Y(I), the following four possibilites 
!                  are available:
!                  1) IMODE=1: S(0)=0, S(N)=0. This is called natural cubic
!                     spline. It is equivalent to assuming that the end cubics
!                     aproach linearity at their extremities.
!                  2) IMODE=2: S(0)=S(1), S(N)=S(N-1). This is equivalent to
!                     assuming that the cubics approach parabolas at their
!                     extremities.
!                  3) IMODE=3: S(0) is a linear extrapolation from S(1) and
!                     S(2), and S(N) is a linear extrapolation from S(N-2) 
!                     and S(N-1).
!                  4) IMODE=4: Force the slopes at each end to assume certain
!                     values.
! REAL S(N) -> if IMODE=4, in input S(1) and S(N) contain the first derivatives
!              at X(1) and X(N). In output, this matrix contains the second
!              derivatives
! REAL A(N) -> spline coefficients
! REAL B(N) -> spline coefficients
! REAL C(N) -> spline coefficients
!
!omment
!------------------------------------------------------------------------------
        SUBROUTINE CUBSPL(X,Y,N,IMODE,S,A,B,C)
        IMPLICIT NONE
!
        INTEGER N
        REAL X(N),Y(N)
        INTEGER IMODE
        REAL S(N)
        REAL A(N),B(N),C(N)
! local parameters
        INTEGER NMAX
        PARAMETER (NMAX=100)
! local variables
        INTEGER I,I1,I2
        REAL DX1,DX2,DY1,DY2
        REAL SS(0:NMAX,4)
        REAL H
!------------------------------------------------------------------------------
! evitamos algun warning de compilacion
        SS(1,1)=0.0
! verificamos que al menos tenemos dos puntos
        IF(N.LT.2)THEN
          WRITE(*,100)'FATAL ERROR in subroutine CUBSPL: '
          WRITE(*,100)' no. of data points too small: '
          WRITE(*,*)N
          STOP
        END IF
! verificamos que el numero de puntos no sea demasiado grande
        IF(N.GT.NMAX)THEN        !en realidad vale con N-1, pero lo dejamos asi
          WRITE(*,100)'FATAL ERROR in subroutine CUBSPL: '
          WRITE(*,100)' no. of data points too large: '
          WRITE(*,*)N
          STOP
        END IF
! verificamos que IMODE toma un valor posible
        IF((IMODE.LT.1).OR.(IMODE.GT.4))THEN
          WRITE(*,100)'FATAL ERROR in subroutine CUBSPL: '
          WRITE(*,100)' invalid IMODE value:'
          WRITE(*,*)IMODE
          STOP
        END IF
!------------------------------------------------------------------------------
! calculamos la matriz reducida, que contiene N-2 filas x 4 columnas
        DX1=X(2)-X(1)
        DY1=6.*(Y(2)-Y(1))/DX1
        DO I=1,N-2
          DX2=X(I+2)-X(I+1)
          DY2=6.*(Y(I+2)-Y(I+1))/DX2
          SS(I,1)=DX1
          SS(I,2)=2.*(DX1+DX2)
          SS(I,3)=DX2
          SS(I,4)=DY2-DY1
          DX1=DX2
          DY1=DY2
        END DO
!------------------------------------------------------------------------------
! modificamos la primera y ultima fila segun el valor de IMODE
        IF(IMODE.EQ.1)THEN              !no hay que modificar nada en este caso
!..............................................................................
        ELSEIF(IMODE.EQ.2)THEN
          SS(1,2)=SS(1,2)+(X(2)-X(1))
          SS(N-2,2)=SS(N-2,2)+(X(N)-X(N-1))
!..............................................................................
        ELSEIF(IMODE.EQ.3)THEN
          DX1=X(2)-X(1)
          DX2=X(3)-X(2)
          SS(1,2)=(DX1+DX2)*(DX1+2.*DX2)/DX2
          SS(1,3)=(DX2*DX2-DX1*DX1)/DX2
          DX1=X(N-1)-X(N-2)
          DX2=X(N)-X(N-1)
          SS(N-2,1)=(DX1*DX1-DX2*DX2)/DX1
          SS(N-2,2)=(DX2+DX1)*(DX2+2.*DX1)/DX1
!..............................................................................
        ELSEIF(IMODE.EQ.4)THEN
! notar que en este caso tambien cambia el taman~o de la matriz
          DX1=X(2)-X(1)
          DY1=(Y(2)-Y(1))/DX1
          SS(0,1)=1.
          SS(0,2)=2.*DX1
          SS(0,3)=DX1
          SS(0,4)=6.*(DY1-S(1))
          DX1=X(N)-X(N-1)
          DY1=(Y(N)-Y(N-1))/DX1
          SS(N-1,1)=DX1
          SS(N-1,2)=2.*DX1
          SS(N-1,3)=0.
          SS(N-1,4)=6.*(S(N)-DY1)
!..............................................................................
        END IF
!------------------------------------------------------------------------------
! dimensiones de la matriz a resolver
        IF(IMODE.EQ.4)THEN
          I1=1
          I2=N-1
        ELSE
          I1=2
          I2=N-2
        END IF
! resolvemos el sistema tridiagonal, eliminando en primer lugar los elementos
! que se encuentran por debajo de la diagonal
        DO I=I1,I2
          SS(I,1)=SS(I,1)/SS(I-1,2)
          SS(I,2)=SS(I,2)-SS(I,1)*SS(I-1,3)
          SS(I,4)=SS(I,4)-SS(I,1)*SS(I-1,4)
        END DO
! y ahora hacemos la sustitucion desde atras
        SS(I2,4)=SS(I2,4)/SS(I2,2)
        DO I=I2-1,I1-1,-1
          SS(I,4)=(SS(I,4)-SS(I,3)*SS(I+1,4))/SS(I,2)
        END DO
!------------------------------------------------------------------------------
! los valores de las derivadas segundas se almacenan en la matriz S
        DO I=I1-1,I2
          S(I+1)=SS(I,4)
        END DO
!
        IF(IMODE.EQ.1)THEN
          S(1)=0.
          S(N)=0.
        ELSEIF(IMODE.EQ.2)THEN
          S(1)=S(2)
          S(N)=S(N-1)
        ELSEIF(IMODE.EQ.3)THEN
          DX1=X(2)-X(1)
          DX2=X(3)-X(2)
          S(1)=((DX1+DX2)*S(2)-DX1*S(3))/DX2
          DX1=X(N-1)-X(N-2)
          DX2=X(N)-X(N-1)
          S(N)=((DX1+DX2)*S(N-1)-DX2*S(N-2))/DX1
        ELSEIF(IMODE.EQ.4)THEN                         !no hay nada que cambiar
        END IF
!------------------------------------------------------------------------------
! finalmente calculamos los coeficientes
        DO I=1,N-1
          H=X(I+1)-X(I)
          A(I)=(S(I+1)-S(I))/(6.*H)
          B(I)=S(I)/2.
          C(I)=(Y(I+1)-Y(I))/H-H*(2.*S(I)+S(I+1))/6.
        END DO
!------------------------------------------------------------------------------
100     FORMAT(A,$)
        END
