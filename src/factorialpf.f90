!------------------------------------------------------------------------------
! Version 28-February-1997                                  File: factorialpf.f
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
! DOUBLE PRECISION FUNCTION FACTORIALPF(N)
!
! Input: N
! Output: FACTORIALPF (function)
!
! Calculate N factorial
!
! INTEGER N
!
!omment
!------------------------------------------------------------------------------
        DOUBLE PRECISION FUNCTION FACTORIALPF(N)
        IMPLICIT NONE
        INTEGER N
!
        INTEGER I
!------------------------------------------------------------------------------
        IF(N.LT.0)THEN
          WRITE(*,101)'FATAL ERROR: factorial(n<0)!'
          STOP
        END IF
        IF(N.GT.30)THEN
          WRITE(*,101)'FATAL ERROR: factorial(n>30)!'
          STOP
        END IF
        FACTORIALPF=1.D0
        IF(N.EQ.0) RETURN
        DO I=1,N
          FACTORIALPF=DBLE(I)*FACTORIALPF 
        END DO
!
101     FORMAT(A)
        END
