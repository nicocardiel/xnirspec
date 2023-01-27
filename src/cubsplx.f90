!------------------------------------------------------------------------------
! Version 18-June-1998                                          File: cubsplx.f
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
! SUBROUTINE CUBSPLX(X,Y,A,B,C,N,I0,X0,Y0)
!
! Input: X,Y,A,B,C,N,I0,X0
! Output: Y0
!
! The subroutine returns the cubic spline evaluated at X0, using the
! coefficients determined in a previous call to CUBSPL. The spline defined in 
! the interval between X(I),Y(I) and X(I+1),Y(I+1) is given by:
!
!      Y = A(I)*(X-X(I))**3 + B(I)*(X-X(I))**2 + C(I)*(X-X(I)) + D(I)
!
! If X0.LT.X(1), I=1 is employed (first computed spline)
! If X0.GT.X(N), I=N-1 is employed (last computed spline)
!
! REAL X(N) -> X-values fitted with CUBSPL
! REAL Y(N) -> Y-values fitted with CUBSPL
! REAL A(N) -> spline coefficients
! REAL B(N) -> spline coefficients
! REAL C(N) -> spline coefficients
! INTEGER N -> number of data points
! INTEGER I0 -> initial location to start the search of the place of X0 in
!               the X array
! REAL X0 -> X-value where the spline function will be evaluated
! REAL Y0 -> spline value at X0
!
!omment
!------------------------------------------------------------------------------
        SUBROUTINE CUBSPLX(X,Y,A,B,C,N,I0,X0,Y0)
        IMPLICIT NONE
!
        INTEGER N
        REAL X(N),Y(N)
        REAL A(N),B(N),C(N)
        INTEGER I0
        REAL X0,Y0
! local variables
        REAL DX
!------------------------------------------------------------------------------
! buscamos el lugar en la tabla en la que se encuentra X0, para lo cual
! empleamos la subrutina BINSEARCH, la cual permite emplear un valor de prueba
! para iniciar la busqueda, lo cual acelera el proceso cuando se realizan
! llamadadas sucesivas a esta funcion, con valores de X0 consecutivos.
        CALL BINSEARCH(X,N,X0,I0)
        IF(I0.EQ.0) I0=1
        IF(I0.EQ.N) I0=N-1
!------------------------------------------------------------------------------
! evaluate the spline
        DX=X0-X(I0)
        Y0=Y(I0)+DX*(C(I0)+DX*(B(I0)+DX*A(I0)))
!------------------------------------------------------------------------------
        END
