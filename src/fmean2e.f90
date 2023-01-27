!------------------------------------------------------------------------------
! Version 28-February-1997                                       File: fmean2.f
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This routine is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!omment
!
! REAL FUNCTION FMEAN2E(N,X,EX,TIMES,SIGMA,EFMEAN)
!
! Input: N,X,TIMES
! Output: FMEAN2E (function), SIGMA, EFMEAN
!
! Calculate the mean value of X(N) rejecting points at TIMES sigma. The
! function can recover points rejected in previous iterations. This routine
! also computes the error of FMEAN2E using the errors EX(N).
!
! INTEGER N -> no. of elements
! REAL    X(N) -> input matrix
! REAL    EX(N) -> input error matrix
! REAL    TIMES -> times sigma to reject points before calculating the mean
! REAL    SIGMA -> resulting sigma after removing outliers
! REAL    EFMEAN -> error of FMEAN2E
!
!omment
!------------------------------------------------------------------------------
        REAL FUNCTION FMEAN2E(N,X,EX,TIMES,SIGMA,EFMEAN)
        IMPLICIT NONE
        INTEGER N
        REAL X(N),EX(N)
        REAL TIMES
        REAL EFMEAN
!
        INCLUDE 'dimensions.inc'
!
        INTEGER NMAX
        PARAMETER(NMAX=NXMAX*NYMAX)
!
        INTEGER I,NN
        REAL SIGMA
        DOUBLE PRECISION SUM,ESUM
        LOGICAL IFX(NMAX),IFXX(NMAX)
        LOGICAL LREPEAT
!------------------------------------------------------------------------------
        IF(N.EQ.0) STOP 'FATAL ERROR in function FMEAN2E: N=0.'
        IF(N.GT.NMAX)THEN
          WRITE(*,101)'FATAL ERROR in function FMEAN2E: N too large.'
          STOP
        END IF
!
        DO I=1,N
          IFX(I)=.TRUE.
        END DO
!
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
!
        SIGMA=0.
        IF(NN.GT.1)THEN
          DO I=1,N
            IF(IFX(I)) SIGMA=SIGMA+(X(I)-FMEAN2E)*(X(I)-FMEAN2E)
          END DO
          SIGMA=SQRT(SIGMA/REAL(NN-1))
        END IF
!
        DO I=1,N
          IFXX(I)=(ABS(X(I)-FMEAN2E).LE.TIMES*SIGMA)
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
        GOTO 10
!
101     FORMAT(A)
        END
