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
! REAL FUNCTION FMEAN2(N,X,TIMES,SIGMA)
!
! Input: N,X,TIMES
! Output: FMEAN2 (function), SIGMA
!
! Calculate the mean value of X(N) rejecting points at TIMES sigma. The
! function can recover points rejected in previous iterations.
!
! INTEGER N -> no. of elements
! REAL    X(N) -> input matrix
! REAL    TIMES -> times sigma to reject points before calculating the mean
! REAL    SIGMA -> resulting sigma after removing outliers
!
!omment
!------------------------------------------------------------------------------
        REAL FUNCTION FMEAN2(N,X,TIMES,SIGMA)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL TIMES
!
        INCLUDE 'dimensions.inc'
!
        INTEGER NMAX
        PARAMETER(NMAX=NXMAX*NYMAX)
!
        INTEGER I,NN
        REAL SUM,SIGMA
        LOGICAL IFX(NMAX),IFXX(NMAX)
        LOGICAL LREPEAT
!------------------------------------------------------------------------------
        IF(N.EQ.0) STOP 'FATAL ERROR in function FMEAN2: N=0.'
        IF(N.GT.NMAX)THEN
          WRITE(*,101)'FATAL ERROR in function FMEAN2: N too large.'
          STOP
        END IF
!
        DO I=1,N
          IFX(I)=.TRUE.
        END DO
!
10      NN=0
        SUM=0.
        DO I=1,N
          IF(IFX(I))THEN
            NN=NN+1
            SUM=SUM+X(I)
          END IF
        END DO
        FMEAN2=SUM/REAL(NN)
        IF(N.EQ.1) RETURN
!
        SIGMA=0.
        IF(NN.GT.1)THEN
          DO I=1,N
            IF(IFX(I)) SIGMA=SIGMA+(X(I)-FMEAN2)*(X(I)-FMEAN2)
          END DO
          SIGMA=SQRT(SIGMA/REAL(NN-1))
        END IF
!
        DO I=1,N
          IFXX(I)=(ABS(X(I)-FMEAN2).LE.TIMES*SIGMA)
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
