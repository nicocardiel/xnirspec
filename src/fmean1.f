C------------------------------------------------------------------------------
C Version 28-February-1997                                       File: fmean1.f
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
C REAL FUNCTION FMEAN1(N,X)
C
C Input: N,X
C Output: FMEAN1 (function)
C
C Calculate the mean value of X(N)
C
C INTEGER N -> no. of elements
C REAL    X(N) -> input matrix
C
Comment
C------------------------------------------------------------------------------
        REAL FUNCTION FMEAN1(N,X)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
C
        INTEGER I
        REAL SUM
C------------------------------------------------------------------------------
        IF(N.EQ.0) STOP 'FATAL ERROR: in function FMEAN1: N=0.'
        SUM=0.
        DO I=1,N
          SUM=SUM+X(I)
        END DO
        FMEAN1=SUM/REAL(N)
C
        END
