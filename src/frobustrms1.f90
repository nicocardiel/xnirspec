!------------------------------------------------------------------------------
! Version 13-September-2019                                  File: frobustrms.f
!------------------------------------------------------------------------------
! Copyright N. Cardiel, Departamento de Física de la Tierra y Astrofísica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: cardiel@ucm.es
!------------------------------------------------------------------------------
! This routine is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!omment
!
! REAL FUNCTION FROBUSTRMS1(N,X,TSIGMA)
!
! Input: N,X,TSIGMA
! Output: FROBUSTRMS1 (function)
!
! Calculate the robust r.m.s. of X(N) using Eq.3.36 from 'Statistics, Data
! Mining, and Machine Learning in Astronomy', by Ivezic, Connolly, VanderPlas &
! Gray (page 84). After the initial computation, points are rejected and the
! r.m.s. is computed again (note that this process is iterated).
!
! INTEGER N -> no. of elements
! REAL    X(N) -> input matrix
! REAL    TSIGMA -> times sigma to reject points
!
!omment
!------------------------------------------------------------------------------
        REAL FUNCTION FROBUSTRMS1(N,X,TSIGMA)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL TSIGMA
!
        INCLUDE 'dimensions.inc'
!
        INTEGER NMAX
        PARAMETER(NMAX=NXMAX*NYMAX)
!
        INTEGER I
        INTEGER I25,I50,I75
        INTEGER NEFF
        INTEGER NITER,NITERMAX
        REAL FMEDIAN
        REAL XEFF(NMAX)
        LOGICAL IFX(NMAX),IFXX(NMAX)
        LOGICAL LREPEAT
!------------------------------------------------------------------------------
        IF(N.LE.0) STOP 'FATAL ERROR: in function FROBUSTRMS1: N.LE.0'
!
        IF(N.EQ.1)THEN
          FROBUSTRMS1=0.0
          RETURN
        END IF
! 
        CALL ORDENA1F(N,X)
!
        DO I=1,N
          IFX(I)=.TRUE.
        END DO
!
        NITERMAX=10  !avoid infinite loops
        NITER=0
10      NEFF=0
        DO I=1,N
          IF(IFX(I))THEN
            NEFF=NEFF+1
            XEFF(NEFF)=X(I)
          END IF
        END DO
! compute percentiles using the nearest-rank method (note that if X has fewer 
! than 100 distinct values, this methos assign the same value to more than one
! percentile; see https://en.wikipedia.org/wiki/Percentile)
        I25=INT(0.25*FLOAT(NEFF)+0.5)
        I75=INT(0.75*FLOAT(NEFF)+0.5)
        FROBUSTRMS1=0.7413*(XEFF(I75)-XEFF(I25))
        IF(NEFF.EQ.1) RETURN
!
        I50=INT(0.50*FLOAT(NEFF)+0.5)
        FMEDIAN=XEFF(I50)
!
        DO I=1,N
          IFXX(I)=(ABS(X(I)-FMEDIAN).LE.TSIGMA*FROBUSTRMS1)
        END DO
!
        LREPEAT=.FALSE.
        DO I=1,N
          IF(IFX(I).NEQV.IFXX(I)) LREPEAT=.TRUE.
        END DO
        IF(.NOT.LREPEAT) RETURN
!
        DO I=1,N
          IFX(I)=IFXX(I)
        END DO
        NITER=NITER+1
        IF(NITER.LT.NITERMAX) GOTO 10
!
        END
