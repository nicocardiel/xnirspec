C------------------------------------------------------------------------------
C Version 16-June-1998                                         File:binsearch.f
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
C SUBROUTINE BINSEARCH(X,N,X0,N0)
C
C Input: X,N,X0,N0
C Output: N0
C
C Given the array X(N), and the test value X0, this subroutine returns an
C integer N0, such that X0 is between X(N0) and X(N0+1). As input N0 is
C employed to start the searching. If X0.LT.X(1) then N0=0 on output, whereas 
C if X0.GT.X(N) then N0=N.
C
C REAL    X(N) -> ordered input array (not necesarilly equally-spaced)
C INTEGER N -> no. of points in input array
C REAL    X0 -> argument to be searched for
C INTEGER N0 -> location of X0 in the input array
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE BINSEARCH(X,N,X0,N0)
        IMPLICIT NONE
C
        INTEGER N,N0
        REAL X(N),X0
C local variables
        INTEGER L,U,I
        INTEGER STEP
        LOGICAL LOOP
C------------------------------------------------------------------------------
        IF(N0.LT.1)THEN
ccc          WRITE(*,100)'* WARNING: in subroutine BINSEARCH: '
ccc          WRITE(*,101)'N0.LT.1'
          N0=1
        END IF
        IF(N0.GT.N)THEN
ccc          WRITE(*,100)'* WARNING: in subroutine BINSEARCH: '
ccc          WRITE(*,101)'N0.GT.N'
          N0=N
        END IF
C------------------------------------------------------------------------------
C Buscamos el intervalo inicial duplicando el paso de busqueda
        STEP=1
        L=N0
        LOOP=.TRUE.
c..............................................................................
        IF((X(1).LT.X(N)).EQV.(X0.GE.X(L)))THEN
          DO WHILE(LOOP)
            U=L+STEP
            IF(U.GT.N)THEN
              U=N+1
              LOOP=.FALSE.
            ELSE
              IF((X(1).LT.X(N)).EQV.(X0.GE.X(U)))THEN
                L=U
                STEP=2*STEP
              ELSE
                LOOP=.FALSE.
              END IF
            END IF
          END DO
c..............................................................................
        ELSE
          U=L
          DO WHILE(LOOP)
            L=U-STEP
            IF(L.LT.1)THEN
              L=0
              LOOP=.FALSE.
            ELSE
              IF((X(1).LT.X(N)).EQV.(X0.LT.X(L)))THEN
                U=L
                STEP=2*STEP
              ELSE
                LOOP=.FALSE.
              END IF
            END IF
          END DO
c..............................................................................
        END IF
C------------------------------------------------------------------------------
C Ahora buscamos el valor de N0 dividiendo el paso de busqueda
        DO WHILE(U-L.GT.1)
          I=(U+L)/2
          IF(X0.EQ.X(I))THEN
            N0=I
            RETURN
          END IF
          IF((X(1).LT.X(N)).EQV.(X0.GT.X(I)))THEN
            L=I
          ELSE
            U=I
          END IF
        END DO
C
        IF(U.LT.L)THEN
          WRITE(*,101)'FATAL ERROR: in subroutine BINSEARCH.'
          STOP
        END IF
C
        N0=L
C
ccc100     FORMAT(A,$)
101     FORMAT(A)
        END
