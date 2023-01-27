!------------------------------------------------------------------------------
! Version 28-February-1997                                       File: combpf.f
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
! DOUBLE PRECISION FUNCTION COMBPF(N,K)
!
! Input: N,K
! Output: COMBPF (function)
!
! Calculate the binomial coefficient N over K
!
! INTEGER N
! INTEGER K
!
!omment
!------------------------------------------------------------------------------
        DOUBLE PRECISION FUNCTION COMBPF(N,K)
        IMPLICIT NONE
!        
        DOUBLE PRECISION FACTORIALPF
        INTEGER N,K
!------------------------------------------------------------------------------
        IF(K.GT.N)THEN
          WRITE(*,101)'FATAL ERROR: in function COMBPF(N,K), K>N'
          STOP
        END IF
        COMBPF=FACTORIALPF(N)/(FACTORIALPF(K)*FACTORIALPF(N-K))
101     FORMAT(A)
        END
