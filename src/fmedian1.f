C------------------------------------------------------------------------------
C Version 16-June-1998                                         File: fmedian1.f
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
C REAL FUNCTION FMEDIAN1(N,X)
C
C Input: N,X
C Output: FMEDIAN (function), X (sorted)
C
C Calculate the median value of X(N). It is important to note that this 
C subroutine rearranges the matrix X which is returned sorted.
C
C INTEGER N -> no. of elements
C REAL    X(N) -> input matrix
C
Comment
C------------------------------------------------------------------------------
        REAL FUNCTION FMEDIAN1(N,X)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
C variables locales
        INTEGER NN
C------------------------------------------------------------------------------
        IF(N.EQ.0) STOP 'FATAL ERROR: in function FMEDIAN: N=0.'
        CALL ORDENA1F(N,X)
        NN=N/2
        IF(MOD(N,2).EQ.0)THEN
          FMEDIAN1=(X(NN)+X(NN+1))/2.
        ELSE
          FMEDIAN1=X(NN+1)
        END IF
C
        END
