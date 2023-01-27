!------------------------------------------------------------------------------
! Version 16-June-1998                                         File: fmedian1.f
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
! REAL FUNCTION FMEDIAN1(N,X)
!
! Input: N,X
! Output: FMEDIAN (function), X (sorted)
!
! Calculate the median value of X(N). It is important to note that this 
! subroutine rearranges the matrix X which is returned sorted.
!
! INTEGER N -> no. of elements
! REAL    X(N) -> input matrix
!
!omment
!------------------------------------------------------------------------------
        REAL FUNCTION FMEDIAN1(N,X)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
! variables locales
        INTEGER NN
!------------------------------------------------------------------------------
        IF(N.EQ.0) STOP 'FATAL ERROR: in function FMEDIAN: N=0.'
        CALL ORDENA1F(N,X)
        NN=N/2
        IF(MOD(N,2).EQ.0)THEN
          FMEDIAN1=(X(NN)+X(NN+1))/2.
        ELSE
          FMEDIAN1=X(NN+1)
        END IF
!
        END
