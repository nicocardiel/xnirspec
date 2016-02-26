C------------------------------------------------------------------------------
C Version 28-February-1997                                       File: combpf.f
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
C DOUBLE PRECISION FUNCTION COMBPF(N,K)
C
C Input: N,K
C Output: COMBPF (function)
C
C Calculate the binomial coefficient N over K
C
C INTEGER N
C INTEGER K
C
Comment
C------------------------------------------------------------------------------
        DOUBLE PRECISION FUNCTION COMBPF(N,K)
        IMPLICIT NONE
C        
        DOUBLE PRECISION FACTORIALPF
        INTEGER N,K
C------------------------------------------------------------------------------
        IF(K.GT.N)THEN
          WRITE(*,101)'FATAL ERROR: in function COMBPF(N,K), K>N'
          STOP
        END IF
        COMBPF=FACTORIALPF(N)/(FACTORIALPF(K)*FACTORIALPF(N-K))
101     FORMAT(A)
        END
