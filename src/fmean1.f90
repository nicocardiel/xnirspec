!------------------------------------------------------------------------------
! Version 28-February-1997                                       File: fmean1.f
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
! REAL FUNCTION FMEAN1(N,X)
!
! Input: N,X
! Output: FMEAN1 (function)
!
! Calculate the mean value of X(N)
!
! INTEGER N -> no. of elements
! REAL    X(N) -> input matrix
!
!omment
!------------------------------------------------------------------------------
        REAL FUNCTION FMEAN1(N,X)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
!
        INTEGER I
        REAL SUM
!------------------------------------------------------------------------------
        IF(N.EQ.0)THEN
          INCLUDE 'deallocate_arrays.inc'
          STOP 'FATAL ERROR: in function FMEAN1: N=0.'
        END IF
        SUM=0.
        DO I=1,N
          SUM=SUM+X(I)
        END DO
        FMEAN1=SUM/REAL(N)
!
        END
