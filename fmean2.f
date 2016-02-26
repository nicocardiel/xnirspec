C------------------------------------------------------------------------------
C Version 28-February-1997                                       File: fmean2.f
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
C REAL FUNCTION FMEAN2(N,X,TIMES,SIGMA)
C
C Input: N,X,TIMES
C Output: FMEAN2 (function), SIGMA
C
C Calculate the mean value of X(N) rejecting points at TIMES sigma. The
C function can recover points rejected in previous iterations.
C
C INTEGER N -> no. of elements
C REAL    X(N) -> input matrix
C REAL    TIMES -> times sigma to reject points before calculating the mean
C REAL    SIGMA -> resulting sigma after removing outliers
C
Comment
C------------------------------------------------------------------------------
        REAL FUNCTION FMEAN2(N,X,TIMES,SIGMA)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL TIMES
C
        INCLUDE 'dimensions.inc'
C
        INTEGER NMAX
        PARAMETER(NMAX=NXMAX*NYMAX)
C
        INTEGER I,NN
        REAL SUM,SIGMA
        LOGICAL IFX(NMAX),IFXX(NMAX)
        LOGICAL LREPEAT
C------------------------------------------------------------------------------
        IF(N.EQ.0) STOP 'FATAL ERROR in function FMEAN2: N=0.'
        IF(N.GT.NMAX)THEN
          WRITE(*,101)'FATAL ERROR in function FMEAN2: '//
     +      'N too large.'
          STOP
        END IF
C
        DO I=1,N
          IFX(I)=.TRUE.
        END DO
C
10        NN=0
        SUM=0.
        DO I=1,N
          IF(IFX(I))THEN
            NN=NN+1
            SUM=SUM+X(I)
          END IF
        END DO
        FMEAN2=SUM/REAL(NN)
        IF(N.EQ.1) RETURN
C
        SIGMA=0.
        IF(NN.GT.1)THEN
          DO I=1,N
            IF(IFX(I)) SIGMA=SIGMA+(X(I)-FMEAN2)*(X(I)-FMEAN2)
          END DO
          SIGMA=SQRT(SIGMA/REAL(NN-1))
        END IF
C
        DO I=1,N
          IFXX(I)=(ABS(X(I)-FMEAN2).LE.TIMES*SIGMA)
        END DO
C
        LREPEAT=.FALSE.
        DO I=1,N
          IF(IFX(I).NEQV.IFXX(I)) LREPEAT=.TRUE.
        END DO
        IF(.NOT.LREPEAT) RETURN
C
        DO I=1,N
          IFX(I)=IFXX(I)
        END DO
        GOTO 10
C
101     FORMAT(A)
        END
