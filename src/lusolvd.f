C------------------------------------------------------------------------------
C Version 18-June-1998                                           File: lusolv.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This routine is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C SUBROUTINE LUSOLVD(A,N,NDIM,ORDER,SCALEROW,B,X)
C
C Input: A,N,NDIM,ORDER,SCALEROW,B
C Output: X
C
C This subroutine solves the set of N linear equations A X = B, where the
C A matrix corresponds to the LU decomposition of the initial coefficient
C matrix. See C.F. Gerald and P. O. Wheatley, in Applied Numerical
C Analysis, 4th edition, pag. 110. The matrix A remains unchanged (also ORDER
C and SCALEROW), so subsequent calls to this subroutine, variying the B matrix,
C can be performed.
C
C DOUBLE PRECISION A(NDIM,NDIM) -> matrix of coefficients (LU in compact scheme)
C INTEGER N -> logical dimension of A
C INTEGER NDIM -> physical dimension of A in the calling program
C INTEGER ORDER(N) -> vector holding row order after pivoting in LUDCMP
C DOUBLE PRECISION SCALEROW(N) -> vector holding scaling factors 
C         applied to each row
C DOUBLE PRECISION B(N) -> right-hand side vector B
C DOUBLE PRECISION X(N) -> solution vector X
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE LUSOLVD(A,N,NDIM,ORDER,SCALEROW,B,X)
        IMPLICIT NONE
C
        INTEGER N,NDIM
        DOUBLE PRECISION A(NDIM,NDIM)
        INTEGER ORDER(N)
        DOUBLE PRECISION SCALEROW(N)
        DOUBLE PRECISION B(N),X(N)
C local variables
        INTEGER I,I0,J
        DOUBLE PRECISION DSUM
C------------------------------------------------------------------------------
C reordenamos y reescalamos los elementos del vector B siguiendo la informacion
C contenida en los vectores ORDER y SCALEROW
        DO I=1,N
          I0=ORDER(I)
          X(I)=B(I0)*SCALEROW(I0)
        END DO
C------------------------------------------------------------------------------
C calculamos el vector Bprima
        X(1)=X(1)/A(1,1)
        IF(N.EQ.1) RETURN
        DO I=2,N
          DSUM=0.D0
          DO J=1,I-1
            DSUM=DSUM+A(I,J)*X(J)
          END DO
          X(I)=(X(I)-DSUM)/A(I,I)
        END DO
C------------------------------------------------------------------------------
C finalmente obtenemos las soluciones para X
        DO I=2,N
          DSUM=0.D0
          DO J=N-I+2,N
            DSUM=DSUM+A(N-I+1,J)*X(J)
          END DO
          X(N-I+1)=X(N-I+1)-DSUM
        END DO
C------------------------------------------------------------------------------
        END
