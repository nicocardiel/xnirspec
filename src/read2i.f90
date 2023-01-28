        SUBROUTINE READ2I(CQUESTION,CDEF,N1,N2)
        IMPLICIT NONE
        CHARACTER*(*) CQUESTION
        INTEGER N1,N2
        CHARACTER*(*) CDEF
!
        INTEGER I,L1,L2
        INTEGER TRUEBEG,TRUELEN
        INTEGER NERR
        CHARACTER*1 C
        CHARACTER*255 CADENA
        LOGICAL LECHO
        COMMON/BLKLECHO/LECHO
!------------------------------------------------------------------------------
        NERR=0
10      L2=TRUELEN(CQUESTION)
        IF(L2.NE.0) WRITE(*,100) CQUESTION(1:L2)
        IF(CDEF.NE.'@')THEN
          L1=TRUEBEG(CDEF)
          IF(L1.NE.0)THEN
            L2=TRUELEN(CDEF)
            WRITE(*,100) ' ['
            WRITE(*,100) CDEF(L1:L2)
            WRITE(*,100) '] ? '
          END IF
        ELSE
          WRITE(*,100) '? '
        END IF
        READ(*,'(A)',ERR=20) CADENA
        IF(TRUELEN(CADENA).EQ.0)THEN
          IF(CDEF.EQ.'@')THEN
            GOTO 10
          END IF
          CADENA=CDEF
        END IF
        DO I=1,TRUELEN(CADENA)
          C=CADENA(I:I)
          IF((INDEX('abcdefghijklmnopqrstuvwxyz',C).NE.0).OR. &
             (INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ./',C).NE.0))THEN
            GOTO 20
          END IF
        END DO
        READ(CADENA,*,ERR=20) N1,N2
        IF(LECHO) WRITE(*,*) N1, N2
        RETURN
20      WRITE(*,101) 'ERROR: invalid character(s) found in numbers. Try again.'
        NERR=NERR+1
        IF(NERR.GT.10)THEN
          INCLUDE 'deallocate_arrays.inc'
          STOP 'FATAL ERROR: too many errors in READ2I.'
        END IF
        GOTO 10
100     FORMAT(A,$)
101     FORMAT(A)
        END
