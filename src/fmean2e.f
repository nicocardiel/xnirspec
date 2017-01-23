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
C REAL FUNCTION FMEAN2E(N,X,EX,TIMES,SIGMA,EFMEAN)
C
C Input: N,X,TIMES
C Output: FMEAN2E (function), SIGMA, EFMEAN
C
C Calculate the mean value of X(N) rejecting points at TIMES sigma. The
C function can recover points rejected in previous iterations. This routine
C also computes the error of FMEAN2E using the errors EX(N).
C
C INTEGER N -> no. of elements
C REAL    X(N) -> input matrix
C REAL    EX(N) -> input error matrix
C REAL    TIMES -> times sigma to reject points before calculating the mean
C REAL    SIGMA -> resulting sigma after removing outliers
C REAL    EFMEAN -> error of FMEAN2E
C
Comment
C------------------------------------------------------------------------------
        REAL FUNCTION FMEAN2E(N,X,EX,TIMES,SIGMA,EFMEAN)
        IMPLICIT NONE
        INTEGER N
        REAL X(N),EX(N)
        REAL TIMES
        REAL EFMEAN
C
        INCLUDE 'dimensions.inc'
C
        INTEGER NMAX
        PARAMETER(NMAX=NXMAX*NYMAX)
C
        INTEGER I,NN
        REAL SIGMA
        DOUBLE PRECISION SUM,ESUM
        LOGICAL IFX(NMAX),IFXX(NMAX)
        LOGICAL LREPEAT
C------------------------------------------------------------------------------
        IF(N.EQ.0) STOP 'FATAL ERROR in function FMEAN2E: N=0.'
        IF(N.GT.NMAX)THEN
          WRITE(*,101)'FATAL ERROR in function FMEAN2E: '//
     +      'N too large.'
          STOP
        END IF
C
        DO I=1,N
          IFX(I)=.TRUE.
        END DO
C
10        NN=0
        SUM=0.D0
        ESUM=0.D0
        DO I=1,N
          IF(IFX(I))THEN
            NN=NN+1
            SUM=SUM+DBLE(X(I))
            ESUM=ESUM+DBLE(EX(I))*DBLE(EX(I))
          END IF
        END DO
        FMEAN2E=REAL(SUM/DBLE(NN))
        EFMEAN=REAL(DSQRT(ESUM)/DBLE(NN))
        IF(N.EQ.1) RETURN
C
        SIGMA=0.
        IF(NN.GT.1)THEN
          DO I=1,N
            IF(IFX(I)) SIGMA=SIGMA+(X(I)-FMEAN2E)*(X(I)-FMEAN2E)
          END DO
          SIGMA=SQRT(SIGMA/REAL(NN-1))
        END IF
C
        DO I=1,N
          IFXX(I)=(ABS(X(I)-FMEAN2E).LE.TIMES*SIGMA)
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
