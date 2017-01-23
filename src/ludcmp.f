C------------------------------------------------------------------------------
C Version 13-September-1999                                      File: ludcmp.f
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
C SUBROUTINE LUDCMP(A,N,NDIM,ORDER,SCALEROW,IOK,IPAR)
C
C Input: A,N,NDIM
C Output: A,ORDER,SCALEROW,IOK,IPAR
C
C This subroutine computes the L and U triangular matrices equivalent to the
C A matrix, such that LU = A.  These matrices are returned in the space of A,
C in compact form. See C.F. Gerald and P. O. Wheatley, in Applied Numerical
C Analysis, 4th edition, pag. 106.
C
C REAL A(NDIM,NDIM) -> matrix of coefficients
C INTEGER N -> logical dimension of A
C INTEGER NDIM -> physical dimension of A in the calling program
C INTEGER ORDER(N) -> vector holding row order after pivoting
C REAL SCALEROW(N) -> vector holding scaling factors applied to each row
C INTEGER IOK -> returns 0 if everything works properly, +(the row number)
C                if all elements in a row are zero, or -(the row number) if 
C                the pivot value is zero.
C INTEGER IPAR -> returns as +1 or -1 depending on whether the number of row
C                 interchanges was even or odd, respectively
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE LUDCMP(A,N,NDIM,ORDER,SCALEROW,IOK,IPAR)
        IMPLICIT NONE
C
        INTEGER N,NDIM
        REAL A(NDIM,NDIM)
        INTEGER ORDER(N)
        REAL SCALEROW(N)
        INTEGER IOK,IPAR
C local variables
        INTEGER I,J,K
        REAL AMIN,AMAX,ATEST
        DOUBLE PRECISION DSUM
        LOGICAL LSCALE,LOK
C------------------------------------------------------------------------------
C el orden inicial es el de partida
        DO I=1,N
          ORDER(I)=I
        END DO
        IOK=0
        IPAR=1
C------------------------------------------------------------------------------
C si N es igual a 1, no hay nada que hacer
        IF(N.EQ.1)THEN
          SCALEROW(1)=1.0
          RETURN
        END IF
C------------------------------------------------------------------------------
C Antes de realizar el proceso, es interesante verificar si existen numeros
C con diferentes ordenes de magnitud en una misma fila. Si es asi, realizamos
C un proceso de reescalado de forma que el valor maximo en cada fila sea igual
C a uno. Si la diferencia no es muy grande no reescalamos, para evitar asi 
C introducir errores de redondeo al realizar precisamente dicho reescalado.
        LSCALE=.FALSE.
C
        DO I=1,N
          AMIN=ABS(A(I,1))
          AMAX=ABS(A(I,1))
          DO J=2,N
            ATEST=ABS(A(I,J))
            IF(ATEST.LT.AMIN) AMIN=ATEST
            IF(ATEST.GT.AMAX) AMAX=ATEST
          END DO
          IF(AMAX.EQ.0.0)THEN
            WRITE(*,100)'ERROR: all coefficientes are zero in row '
            WRITE(*,*)I
            WRITE(*,100)'Press <CR> to continue...'
            READ(*,*)
            IOK=I
            RETURN
          ELSEIF(AMIN.EQ.0.0)THEN
            LSCALE=.TRUE.
          ELSE
            IF(ALOG10(AMAX)-ALOG10(AMIN).GT.1.5) LSCALE=.TRUE.
          END IF
          SCALEROW(I)=1./AMAX
        END DO
C
        IF(LSCALE)THEN
          DO I=1,N
            DO J=1,N
              A(I,J)=A(I,J)*SCALEROW(I)
            END DO
          END DO
        ELSE
          DO I=1,N
            SCALEROW(I)=1.
          END DO
        END IF
C------------------------------------------------------------------------------
C hacemos pivote sobre la primera columna, siempre que N sea mayor que 1
!!!10      CALL PIVOTE(A,N,NDIM,ORDER,1,IPAR)
        CALL PIVOTE(A,N,NDIM,ORDER,1,IPAR)
C------------------------------------------------------------------------------
C si el pivote es muy pequen~o, podemos tener una matriz singular o casi 
C singular
        CALL LITTLEPIVOTE(A(1,1),LOK)
        IF(.NOT.LOK)THEN
          IOK=1
          RETURN
        END IF
C------------------------------------------------------------------------------
C calculamos los coeficientes de la primera fila de la matriz U
        DO J=2,N
          A(1,J)=A(1,J)/A(1,1)
        END DO
C------------------------------------------------------------------------------
C continuamos el calculo de los coeficientes de L y U (salvo el ultimo
C coeficiente de L, que lo dejamos para despues). Notar como el sumatorio
C puede realizarse en doble precision (Ver Gerald and Wheatley, pag. 109).
        DO J=2,N-1
C calculamos una columna de la matriz L
          DO I=J,N
            DSUM=0.D0
            DO K=1,J-1
              DSUM=DSUM+DBLE(A(I,K))*DBLE(A(K,J))
            END DO
            A(I,J)=A(I,J)-REAL(DSUM)
          END DO
C en caso necesario podemos intercambiar filas
          CALL PIVOTE(A,N,NDIM,ORDER,J,IPAR)
C comprobamos si el pivote es demasiado pequen~o
          CALL LITTLEPIVOTE(A(J,J),LOK)
          IF(.NOT.LOK)THEN
            IOK=J
            RETURN
          END IF
C calculamos una fila de la matriz U
          DO K=J+1,N
            DSUM=0.D0
            DO I=1,J-1
              DSUM=DSUM+DBLE(A(J,I))*DBLE(A(I,K))
            END DO
            A(J,K)=(A(J,K)-REAL(DSUM))/A(J,J)
          END DO
        END DO
C finalmente, calculamos el ultimo elemento de la matriz L
        DSUM=0.D0
        DO K=1,N-1
          DSUM=DSUM+DBLE(A(N,K))*DBLE(A(K,N))
        END DO
        A(N,N)=A(N,N)-REAL(DSUM)
        CALL LITTLEPIVOTE(A(N,N),LOK)
        IF(.NOT.LOK)THEN
          IOK=N
          RETURN
        END IF
C------------------------------------------------------------------------------
100     FORMAT(A,$)
        END
C
C******************************************************************************
C Esta subrutina encuentra el elemento mayor para realizar pivote en la columna
C J0, intercambia los elementos de la matriz A y guarda un registro de los
C cambios realizados en la matriz ORDER. La variable IPAR retorna +1 o -1
C dependiendo del numero de cambios realizados en las filas.
        SUBROUTINE PIVOTE(A,N,NDIM,ORDER,J0,IPAR)
        IMPLICIT NONE
C
        INTEGER N,NDIM
        REAL A(NDIM,NDIM)
        INTEGER ORDER(N),J0,IPAR
C local variables
        INTEGER I,IPIVOTE,J,I0
        REAL AMAX,ATEST
C------------------------------------------------------------------------------
        IF(J0.EQ.N)THEN
          WRITE(*,100)'ERROR: no pivot row can be found for the last '
          WRITE(*,101)'element on diagonal.'
          WRITE(*,100)'Press <CR> to continue...'
          READ(*,*)
        END IF
C------------------------------------------------------------------------------
C buscamos la fila para realizar pivote, considerando solamente los elementos
C que se encuentran sobre y bajo la diagonal
        AMAX=ABS(A(J0,J0))
        IPIVOTE=J0
        DO I=J0+1,N
          ATEST=ABS(A(I,J0))
          IF(ATEST.GT.AMAX)THEN
            AMAX=ATEST
            IPIVOTE=I
          END IF
        END DO
C------------------------------------------------------------------------------
C Si la fila con el elemento mayor es J0, no hay nada que hacer. En caso
C contrario, realizamos el intercambio, tanto en la matriz A como en el
C vector ORDER.
        IF(IPIVOTE.EQ.J0) RETURN
C
        DO J=1,N
          ATEST=A(J0,J)
          A(J0,J)=A(IPIVOTE,J)
          A(IPIVOTE,J)=ATEST
        END DO
        I0=ORDER(J0)
        ORDER(J0)=ORDER(IPIVOTE)
        ORDER(IPIVOTE)=I0
        IPAR=-IPAR
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C Comprueba si el valor de A es demasiado pequen~o como para poder tener
C una matriz problematica.
        SUBROUTINE LITTLEPIVOTE(A,LOK)
        IMPLICIT NONE
C
        REAL A
        LOGICAL LOK
C------------------------------------------------------------------------------
        LOK=.TRUE.
C
        IF(ABS(A).EQ.0.0)THEN
          WRITE(*,101)'ERROR: singular matrix.'
          WRITE(*,100)'Press <CR> to continue...'
          READ(*,*)
          LOK=.FALSE.
        ELSEIF(ABS(A).LT.1.E-7)THEN
          WRITE(*,101)'WARNING: nearly singular matrix.'
          WRITE(*,100)'Press <CR> to continue...'
          READ(*,*)
          LOK=.FALSE.
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
