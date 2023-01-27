!------------------------------------------------------------------------------
! Version 25-November-1996                                      File: findmml.f
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
! SUBROUTINE FINDMML(N,N1,N2,X,XMIN,XMAX)
!
! Input: N,N1,N2,X
! Output: XMIN,XMAX
!
! Return the maximum and minimum value of matrix X of N elements (in the
! range from N1 to N2 exclusively)
!
! INTEGER N -> no. of elements of matrix X
! INTEGER N1 -> first element of X() to search minimum/maximum
! INTEGER N2 -> last element of X() to search minimum/maximum
! REAL    X(N) -> data matrix
! REAL    XMIN -> minimum value of X()
! REAL    XMAX -> maximum value of X()
!
!omment
!------------------------------------------------------------------------------
        SUBROUTINE FINDMML(N,N1,N2,X,XMIN,XMAX)
        IMPLICIT NONE
        INTEGER N,N1,N2
        REAL X(N)
        REAL XMIN,XMAX
!
        INTEGER I
!------------------------------------------------------------------------------
        IF((N1.LT.1).OR.(N2.GT.N).OR.(N2.LT.N1))THEN
          WRITE(*,101)'ERROR: limits out of range in FINDMML'
          WRITE(*,101)'=> Returned values: XMIN = XMAX = 0'
          XMIN=0.
          XMAX=0.
          RETURN
        END IF
!
        XMIN=X(N1)
        XMAX=XMIN
        DO I=N1+1,N2
          IF(X(I).LT.XMIN) XMIN=X(I)
          IF(X(I).GT.XMAX) XMAX=X(I)
        END DO
101     FORMAT(A)
        END
