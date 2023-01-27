!------------------------------------------------------------------------------
! Version 18-June-1998                                           File: lusolv.f
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
! SUBROUTINE LUSOLV(A,N,NDIM,ORDER,SCALEROW,B,X)
!
! Input: A,N,NDIM,ORDER,SCALEROW,B
! Output: X
!
! This subroutine solves the set of N linear equations A X = B, where the
! A matrix corresponds to the LU decomposition of the initial coefficient
! matrix. See C.F. Gerald and P. O. Wheatley, in Applied Numerical
! Analysis, 4th edition, pag. 110. The matrix A remains unchanged (also ORDER
! and SCALEROW), so subsequent calls to this subroutine, variying the B matrix,
! can be performed.
!
! REAL A(NDIM,NDIM) -> matrix of coefficients (LU in compact scheme)
! INTEGER N -> logical dimension of A
! INTEGER NDIM -> physical dimension of A in the calling program
! INTEGER ORDER(N) -> vector holding row order after pivoting in LUDCMP
! REAL SCALEROW(N) -> vector holding scaling factors applied to each row
! REAL B(N) -> right-hand side vector B
! REAL X(N) -> solution vector X
!
!omment
!------------------------------------------------------------------------------
        SUBROUTINE LUSOLV(A,N,NDIM,ORDER,SCALEROW,B,X)
        IMPLICIT NONE
!
        INTEGER N,NDIM
        REAL A(NDIM,NDIM)
        INTEGER ORDER(N)
        REAL SCALEROW(N)
        REAL B(N),X(N)
! local variables
        INTEGER I,I0,J
        DOUBLE PRECISION DSUM
!------------------------------------------------------------------------------
! reordenamos y reescalamos los elementos del vector B siguiendo la informacion
! contenida en los vectores ORDER y SCALEROW
        DO I=1,N
          I0=ORDER(I)
          X(I)=B(I0)*SCALEROW(I0)
        END DO
!------------------------------------------------------------------------------
! calculamos el vector Bprima
        X(1)=X(1)/A(1,1)
        IF(N.EQ.1) RETURN
        DO I=2,N
          DSUM=0.D0
          DO J=1,I-1
            DSUM=DSUM+DBLE(A(I,J))*DBLE(X(J))
          END DO
          X(I)=(X(I)-REAL(DSUM))/A(I,I)
        END DO
!------------------------------------------------------------------------------
! finalmente obtenemos las soluciones para X
        DO I=2,N
          DSUM=0.D0
          DO J=N-I+2,N
            DSUM=DSUM+DBLE(A(N-I+1,J))*DBLE(X(J))
          END DO
          X(N-I+1)=X(N-I+1)-REAL(DSUM)
        END DO
!------------------------------------------------------------------------------
        END
