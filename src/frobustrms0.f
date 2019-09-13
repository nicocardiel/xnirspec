C------------------------------------------------------------------------------
C Version 13-September-2019                                  File: frobustrms.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel, Departamento de Física de la Tierra y Astrofísica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: cardiel@ucm.es
C------------------------------------------------------------------------------
C This routine is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C REAL FUNCTION FROBUSTRMS0(N,X)
C
C Input: N,X
C Output: FROBUSTRMS0 (function)
C
C Calculate the robust r.m.s. of X(N) using Eq.3.36 from 'Statistics, Data
C Mining, and Machine Learning in Astronomy', by Ivezic, Connolly, VanderPlas &
C Gray (page 84).
C
C INTEGER N -> no. of elements
C REAL    X(N) -> input matrix
C
Comment
C------------------------------------------------------------------------------
        REAL FUNCTION FROBUSTRMS0(N,X)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
C
        INTEGER I25,I75
C------------------------------------------------------------------------------
        IF(N.LE.0) STOP 'FATAL ERROR: in function FROBUSTRMS: N.LE.0.'
C
        IF(N.EQ.1)THEN
          FROBUSTRMS0=0.0
          RETURN
        END IF
C 
        CALL ORDENA1F(N,X)
C compute percentiles using the nearest-rank method (note that if X has fewer 
C than 100 distinct values, this methos assign the same value to more than one
C percentile; see https://en.wikipedia.org/wiki/Percentile)
        I25=INT(0.25*FLOAT(N)+0.5)
        I75=INT(0.75*FLOAT(N)+0.5)
        FROBUSTRMS0=0.7413*(X(I75)-X(I25))
C
        END
