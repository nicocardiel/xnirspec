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
C REAL FUNCTION FROBUSTRMS1(N,X,TSIGMA)
C
C Input: N,X,TSIGMA
C Output: FROBUSTRMS1 (function)
C
C Calculate the robust r.m.s. of X(N) using Eq.3.36 from 'Statistics, Data
C Mining, and Machine Learning in Astronomy', by Ivezic, Connolly, VanderPlas &
C Gray (page 84). After the initial computation, points are rejected and the
C r.m.s. is computed again (note that this process is iterated).
C
C INTEGER N -> no. of elements
C REAL    X(N) -> input matrix
C REAL    TSIGMA -> times sigma to reject points
C
Comment
C------------------------------------------------------------------------------
        REAL FUNCTION FROBUSTRMS1(N,X,TSIGMA)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL TSIGMA
C
        INCLUDE 'dimensions.inc'
C
        INTEGER NMAX
        PARAMETER(NMAX=NXMAX*NYMAX)
C
        REAL FMEAN2
C
        INTEGER I
        INTEGER I25,I50,I75
        INTEGER NEFF
        INTEGER NITER,NITERMAX
        REAL SIGMA,FMEDIAN
        REAL XEFF(NMAX)
        LOGICAL IFX(NMAX),IFXX(NMAX)
        LOGICAL LREPEAT
C------------------------------------------------------------------------------
        IF(N.LE.0) STOP 'FATAL ERROR: in function FROBUSTRMS1: N.LE.0'
C
        IF(N.EQ.1)THEN
          FROBUSTRMS1=0.0
          RETURN
        END IF
C 
        CALL ORDENA1F(N,X)
C
        DO I=1,N
          IFX(I)=.TRUE.
        END DO
C
        NITERMAX=10  !avoid infinite loops
        NITER=0
10      NEFF=0
        DO I=1,N
          IF(IFX(I))THEN
            NEFF=NEFF+1
            XEFF(NEFF)=X(I)
          END IF
        END DO
C compute percentiles using the nearest-rank method (note that if X has fewer 
C than 100 distinct values, this methos assign the same value to more than one
C percentile; see https://en.wikipedia.org/wiki/Percentile)
        I25=INT(0.25*FLOAT(NEFF)+0.5)
        I75=INT(0.75*FLOAT(NEFF)+0.5)
        FROBUSTRMS1=0.7413*(XEFF(I75)-XEFF(I25))
        IF(NEFF.EQ.1) RETURN
C
        I50=INT(0.50*FLOAT(NEFF)+0.5)
        FMEDIAN=XEFF(I50)
C
        DO I=1,N
          IFXX(I)=(ABS(X(I)-FMEDIAN).LE.TSIGMA*FROBUSTRMS1)
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
        NITER=NITER+1
        IF(NITER.LT.NITERMAX) GOTO 10
C
        END
