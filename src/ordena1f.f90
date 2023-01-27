!------------------------------------------------------------------------------
! Version 15-June-1998                                         File: ordena1f.f
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
! SUBROUTINE ORDENA1F(N,A)
!
! Input: N,A
! Output: A
!
! Sorts the real array A(N) into ascending numerical order. Note that the input
! array is returned rearranged. We follow the Heapsort method, as described by
! D. Knuth, The Art of Computer Programming (pag.146, 5.2.3.).
!
! INTEGER N -> input number of data in A
! REAL    A(N) -> data matrix to be sorted 
!
!omment
!------------------------------------------------------------------------------
        SUBROUTINE ORDENA1F(N,A)
        IMPLICIT NONE
!
        INTEGER N
        REAL A(N)
! local variables
        INTEGER L,R,I,J
        REAL RR
!------------------------------------------------------------------------------
        IF(N.EQ.1)RETURN                                              !evidente
!
! H1: initialize
        L=N/2+1
        R=N
!
! H2: decrease L or R
10      IF(L.GT.1)THEN
          L=L-1
          RR=A(L)
        ELSE
          RR=A(R)
          A(R)=A(1)
          R=R-1
          IF(R.EQ.1)THEN
            A(1)=RR
            RETURN
          END IF
        END IF
!
! H3: prepare for sift-up
        J=L
!
! H4: advance downward
20      I=J
        J=2*J
        IF(J.GT.R)THEN
          GOTO 30
        END IF
!
! H5: find "larger" son
        IF(J+1.LE.R)THEN
          IF(A(J).LT.A(J+1))THEN
            J=J+1
          END IF
        END IF
!
! H6: larger than A(J)?
        IF(RR.GE.A(J))THEN
          GOTO 30
        END IF
!
! H7: move it up
        A(I)=A(J)
        GOTO 20
!
! H8: store RR
30      A(I)=RR
        GOTO 10                                              !return to step H2
!
        END
