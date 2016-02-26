C------------------------------------------------------------------------------
C Version 8-October-1998                                         File: fmean0.f
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
C REAL FUNCTION FMEAN0E(N,X,EX,SIGMA,EFMEAN)
C
C Input: N,X,EX
C Output: FMEAN0E (function),SIGMA,EFMEAN
C
C Calculate the mean value of X(N) and its r.m.s. Using the errors EX(N), the
C routine also computes the error of FMEAN0E.
C
C INTEGER N -> no. of elements
C REAL    X(N) -> input matrix
C REAL    EX(N) -> input error matrix
C REAL SIGMA -> r.m.s. around the mean value
C REAL EFMEAN -> error of FMEAN0E
C
Comment
C------------------------------------------------------------------------------
        REAL FUNCTION FMEAN0E(N,X,EX,SIGMA,EFMEAN)
        IMPLICIT NONE
        INTEGER N
        REAL X(N),EX(N)
        REAL SIGMA
        REAL EFMEAN
C
        INTEGER I
        DOUBLE PRECISION SUM,SUME
C------------------------------------------------------------------------------
        IF(N.LE.0) STOP 'FATAL ERROR: in function FMEAN0E: N.LE.0'
        SUM=0.D0
        SUME=0.D0
        DO I=1,N
          SUM=SUM+DBLE(X(I))
          SUME=SUME+DBLE(EX(I))*DBLE(EX(I))
        END DO
        FMEAN0E=REAL(SUM/DBLE(N))
        EFMEAN=REAL(DSQRT(SUME)/DBLE(N))
C
        IF(N.EQ.1)THEN
          SIGMA=0.
        ELSE
          SUM=0.
          DO I=1,N
            SUM=SUM+(X(I)-FMEAN0E)*(X(I)-FMEAN0E)
          END DO
          SIGMA=SQRT(REAL(SUM)/REAL(N-1))
        END IF
C
        END
