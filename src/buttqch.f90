!------------------------------------------------------------------------------
! Version 26-November-1996                                      File: buttqch.f
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
! SUBROUTINE BUTTQCH(SIZE)
!
! Output: SIZE
!
! Return the current character font size in buttons.
!
! REAL SIZE -> the current font size (dimensionless multiple of the default
!      size)
!
!omment
!------------------------------------------------------------------------------

! Return the current character height in buttons
!
        SUBROUTINE BUTTQCH(SIZE)
        IMPLICIT NONE
        REAL SIZE
        INCLUDE 'button.inc'
!------------------------------------------------------------------------------
        SIZE=PGSCH_BUTT
        END
