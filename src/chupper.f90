        SUBROUTINE CHUPPER(CADENA)
        IMPLICIT NONE
        CHARACTER*(*) CADENA
!
        INTEGER I,N
!------------------------------------------------------------------------------
        DO I=1,LEN(CADENA)
          N=ICHAR(CADENA(I:I))
          IF((N.GE.97).AND.(N.LE.122)) CADENA(I:I)=CHAR(N-32)
        END DO
        END
