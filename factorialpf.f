C------------------------------------------------------------------------------
C Version 28-February-1997                                  File: factorialpf.f
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
C DOUBLE PRECISION FUNCTION FACTORIALPF(N)
C
C Input: N
C Output: FACTORIALPF (function)
C
C Calculate N factorial
C
C INTEGER N
C
Comment
C------------------------------------------------------------------------------
        DOUBLE PRECISION FUNCTION FACTORIALPF(N)
        IMPLICIT NONE
        INTEGER N
C
        INTEGER I
C------------------------------------------------------------------------------
        IF(N.LT.0)THEN
          WRITE(*,101)'FATAL ERROR: factorial(n<0)!'
          STOP
        END IF
        IF(N.GT.30)THEN
          WRITE(*,101)'FATAL ERROR: factorial(n>30)!'
          STOP
        END IF
        FACTORIALPF=1.D0
        IF(N.EQ.0) RETURN
        DO I=1,N
          FACTORIALPF=DBLE(I)*FACTORIALPF 
        END DO
C
101     FORMAT(A)
        END
