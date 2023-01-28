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
! REAL FUNCTION FROBUSTRMS0(N,X)
!
! Input: N,X
! Output: FROBUSTRMS0 (function)
!
! Calculate the robust r.m.s. of X(N) using Eq.3.36 from 'Statistics, Data
! Mining, and Machine Learning in Astronomy', by Ivezic, Connolly, VanderPlas &
! Gray (page 84).
!
! INTEGER N -> no. of elements
! REAL    X(N) -> input matrix
!
!omment
!------------------------------------------------------------------------------
        REAL FUNCTION FROBUSTRMS0(N,X)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
!
        INTEGER I25,I75
!------------------------------------------------------------------------------
        IF(N.LE.0)THEN
          INCLUDE 'deallocate_arrays.inc'
          STOP 'FATAL ERROR: in function FROBUSTRMS: N.LE.0.'
        END IF
!
        IF(N.EQ.1)THEN
          FROBUSTRMS0=0.0
          RETURN
        END IF
! 
        CALL ORDENA1F(N,X)
! compute percentiles using the nearest-rank method (note that if X has fewer 
! than 100 distinct values, this methos assign the same value to more than one
! percentile; see https://en.wikipedia.org/wiki/Percentile)
        I25=INT(0.25*FLOAT(N)+0.5)
        I75=INT(0.75*FLOAT(N)+0.5)
        FROBUSTRMS0=0.7413*(X(I75)-X(I25))
!
        END
