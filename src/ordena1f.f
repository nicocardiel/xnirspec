C------------------------------------------------------------------------------
C Version 15-June-1998                                         File: ordena1f.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This routine is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C SUBROUTINE ORDENA1F(N,A)
C
C Input: N,A
C Output: A
C
C Sorts the real array A(N) into ascending numerical order. Note that the input
C array is returned rearranged. We follow the Heapsort method, as described by
C D. Knuth, The Art of Computer Programming (pag.146, 5.2.3.).
C
C INTEGER N -> input number of data in A
C REAL    A(N) -> data matrix to be sorted 
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE ORDENA1F(N,A)
        IMPLICIT NONE
C
        INTEGER N
        REAL A(N)
C local variables
        INTEGER L,R,I,J
        REAL RR
C------------------------------------------------------------------------------
        IF(N.EQ.1)RETURN                                              !evidente
C
C H1: initialize
        L=N/2+1
        R=N
C
C H2: decrease L or R
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
C
C H3: prepare for sift-up
        J=L
C
C H4: advance downward
20      I=J
        J=2*J
        IF(J.GT.R)THEN
          GOTO 30
        END IF
C
C H5: find "larger" son
        IF(J+1.LE.R)THEN
          IF(A(J).LT.A(J+1))THEN
            J=J+1
          END IF
        END IF
C
C H6: larger than A(J)?
        IF(RR.GE.A(J))THEN
          GOTO 30
        END IF
C
C H7: move it up
        A(I)=A(J)
        GOTO 20
C
C H8: store RR
30      A(I)=RR
        GOTO 10                                              !return to step H2
C
        END
