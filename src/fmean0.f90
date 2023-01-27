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
! REAL FUNCTION FMEAN0(N,X,SIGMA)
!
! Input: N,X
! Output: FMEAN0 (function),SIGMA
!
! Calculate the mean value of X(N) and its r.m.s.
!
! INTEGER N -> no. of elements
! REAL    X(N) -> input matrix
! REAL SIGMA -> r.m.s. around the mean value
!
!omment
!------------------------------------------------------------------------------
        REAL FUNCTION FMEAN0(N,X,SIGMA)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL SIGMA
!
        INTEGER I
        REAL SUM
!------------------------------------------------------------------------------
        IF(N.LE.0) STOP 'FATAL ERROR: in function FMEAN0: N.LE.0'
        SUM=0.
        DO I=1,N
          SUM=SUM+X(I)
        END DO
        FMEAN0=SUM/REAL(N)
!
        IF(N.EQ.1)THEN
          SIGMA=0.
        ELSE
          SUM=0.
          DO I=1,N
            SUM=SUM+(X(I)-FMEAN0)*(X(I)-FMEAN0)
          END DO
          SIGMA=SQRT(SUM/REAL(N-1))
        END IF
!
        END
