!------------------------------------------------------------------------------
! Version 8-October-1998                                         File: fmean0.f
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
! REAL FUNCTION FMEAN0E(N,X,EX,SIGMA,EFMEAN)
!
! Input: N,X,EX
! Output: FMEAN0E (function),SIGMA,EFMEAN
!
! Calculate the mean value of X(N) and its r.m.s. Using the errors EX(N), the
! routine also computes the error of FMEAN0E.
!
! INTEGER N -> no. of elements
! REAL    X(N) -> input matrix
! REAL    EX(N) -> input error matrix
! REAL SIGMA -> r.m.s. around the mean value
! REAL EFMEAN -> error of FMEAN0E
!
!omment
!------------------------------------------------------------------------------
        REAL FUNCTION FMEAN0E(N,X,EX,SIGMA,EFMEAN)
        IMPLICIT NONE
        INTEGER N
        REAL X(N),EX(N)
        REAL SIGMA
        REAL EFMEAN
!
        INTEGER I
        DOUBLE PRECISION SUM,SUME
!------------------------------------------------------------------------------
        IF(N.LE.0) STOP 'FATAL ERROR: in function FMEAN0E: N.LE.0'
        SUM=0.D0
        SUME=0.D0
        DO I=1,N
          SUM=SUM+DBLE(X(I))
          SUME=SUME+DBLE(EX(I))*DBLE(EX(I))
        END DO
        FMEAN0E=REAL(SUM/DBLE(N))
        EFMEAN=REAL(DSQRT(SUME)/DBLE(N))
!
        IF(N.EQ.1)THEN
          SIGMA=0.
        ELSE
          SUM=0.
          DO I=1,N
            SUM=SUM+(X(I)-FMEAN0E)*(X(I)-FMEAN0E)
          END DO
          SIGMA=SQRT(REAL(SUM)/REAL(N-1))
        END IF
!
        END
