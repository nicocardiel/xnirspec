!------------------------------------------------------------------------------
! Version 20-october-2000                                    File: ordena3f3i.f
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
! SUBROUTINE ORDENA3F3I(N,A,B,C,D,E,F)
!
! Input: N,A,B,C,D,E,F
! Output: A,B,C,D,E,F
!
! Sorts the real array A(N) into ascending numerical order. The additional
! arrays B(N), C(N), D(N), E(N) and F(N) are simultaneously changed in 
! parallel with the array A(N). Note that the five input arrays are returned 
! rearranged. We follow the Heapsort method, as described by D. Knuth, The Art 
! of Computer Programming (pag.146, 5.2.3.).
!
! INTEGER N -> input number of data in A
! REAL    A(N) -> data matrix to be sorted 
! REAL    B(N) -> data matrix to be sorted in parallel with matrix A
! REAL    C(N) -> data matrix to be sorted in parallel with matrix A
! INTEGER D(N) -> data matrix to be sorted in parallel with matrix A
! INTEGER E(N) -> data matrix to be sorted in parallel with matrix A
! INTEGER F(N) -> data matrix to be sorted in parallel with matrix A
!
!omment
!------------------------------------------------------------------------------
        SUBROUTINE ORDENA3F3I(N,A,B,C,D,E,F)
        IMPLICIT NONE
!
        INTEGER N
        REAL A(N)
        REAL B(N)
        REAL C(N)
        INTEGER D(N),E(N),F(N)
! local variables
        INTEGER L,R,I,J
        REAL RRB,RRC
        INTEGER RRD,RRE,RRF
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
          RRB=B(L)
          RRC=C(L)
          RRD=D(L)
          RRE=E(L)
          RRF=F(L)
        ELSE
          RR=A(R)
          RRB=B(R)
          RRC=C(R)
          RRD=D(R)
          RRE=E(R)
          RRF=F(R)
          A(R)=A(1)
          B(R)=B(1)
          C(R)=C(1)
          D(R)=D(1)
          E(R)=E(1)
          F(R)=F(1)
          R=R-1
          IF(R.EQ.1)THEN
            A(1)=RR
            B(1)=RRB
            C(1)=RRC
            D(1)=RRD
            E(1)=RRE
            F(1)=RRF
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
        B(I)=B(J)
        C(I)=C(J)
        D(I)=D(J)
        E(I)=E(J)
        F(I)=F(J)
        GOTO 20
!
! H8: store RR
30      A(I)=RR
        B(I)=RRB
        C(I)=RRC
        D(I)=RRD
        E(I)=RRE
        F(I)=RRF
        GOTO 10                                              !return to step H2
!
        END
