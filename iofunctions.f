C Version October 18, 2000
C******************************************************************************
C******************************************************************************
C Funciones de entrada/salida por teclado             (c) ncl@astrax.fis.ucm.es
C******************************************************************************
C******************************************************************************
C
        CHARACTER*(*) FUNCTION READC(CQUESTION,CDEF,CVAL)
        IMPLICIT NONE
        CHARACTER*(*) CQUESTION
        CHARACTER*(*) CDEF
        CHARACTER*(*) CVAL
C
        INTEGER I,L1,L2
        INTEGER TRUEBEG,TRUELEN
        INTEGER NERR
        CHARACTER*255 CADENA
        LOGICAL LECHO
        COMMON/BLKLECHO/LECHO
C------------------------------------------------------------------------------
        L1=0 !evita un WARNING de compilación
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
        READ(*,101,ERR=20) CADENA
        IF(CVAL.EQ.'@')THEN
          IF(TRUELEN(CADENA).EQ.0)THEN
            IF(CDEF.EQ.'@')THEN
              GOTO 10
            END IF
            CADENA=CDEF(L1:L2)
          END IF
        ELSE
          IF(TRUELEN(CADENA).EQ.0)THEN
            IF(CDEF.EQ.'@')THEN
              GOTO 10
            END IF
            CADENA=CDEF(L1:L2)
          ELSE
            DO I=1,TRUELEN(CADENA)
              IF(INDEX(CVAL,CADENA(I:I)).EQ.0)THEN
                WRITE(*,101) 'ERROR: invalid character(s). Try again.'
                NERR=NERR+1
                IF(NERR.GT.10) 
     +           STOP 'FATAL ERROR: too many errors in READC.'
                GOTO 10
              END IF
            END DO
          END IF
        END IF
        READC=CADENA
        IF(LECHO) WRITE(*,101) CADENA(1:TRUELEN(CADENA))
        RETURN
20      WRITE(*,101) 'ERROR: invalid character(s). Try again.'
        NERR=NERR+1
        IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors in READC.'
        GOTO 10
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
        INTEGER FUNCTION READI(CQUESTION,CDEF)
        IMPLICIT NONE
        CHARACTER*(*) CQUESTION
        CHARACTER*(*) CDEF
C
        INTEGER I,L1,L2
        INTEGER N
        INTEGER NERR
        INTEGER TRUEBEG,TRUELEN
        CHARACTER*1 C
        CHARACTER*255 CADENA
        LOGICAL LECHO
        COMMON/BLKLECHO/LECHO
C------------------------------------------------------------------------------
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
        READ(*,101,ERR=20) CADENA
        IF(TRUELEN(CADENA).EQ.0)THEN
          IF(CDEF.EQ.'@')THEN
            GOTO 10
          END IF
          CADENA=CDEF
        END IF
        DO I=1,TRUELEN(CADENA)
          C=CADENA(I:I)
          IF((INDEX('abcdefghijklmnopqrstuvwxyz',C).NE.0).OR.
     +     (INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ./',C).NE.0))THEN
            GOTO 20
          END IF
        END DO
        READ(CADENA,*,ERR=20) N
        READI=N
        IF(LECHO) WRITE(*,*) N
        RETURN
20      WRITE(*,101) 'ERROR: invalid character(s) found in '//
     +   'number. Try again.'
        NERR=NERR+1
        IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors in READI.'
        GOTO 10
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
        SUBROUTINE READ2I(CQUESTION,CDEF,N1,N2)
        IMPLICIT NONE
        CHARACTER*(*) CQUESTION
        INTEGER N1,N2
        CHARACTER*(*) CDEF
C
        INTEGER I,L1,L2
        INTEGER TRUEBEG,TRUELEN
        INTEGER NERR
        CHARACTER*1 C
        CHARACTER*255 CADENA
        LOGICAL LECHO
        COMMON/BLKLECHO/LECHO
C------------------------------------------------------------------------------
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
          IF((INDEX('abcdefghijklmnopqrstuvwxyz',C).NE.0).OR.
     +     (INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ./',C).NE.0))THEN
            GOTO 20
          END IF
        END DO
        READ(CADENA,*,ERR=20) N1,N2
        IF(LECHO) WRITE(*,*) N1, N2
        RETURN
20      WRITE(*,101) 'ERROR: invalid character(s) found in '//
     +   'numbers. Try again.'
        NERR=NERR+1
        IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors in READ2I.'
        GOTO 10
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
        INTEGER FUNCTION READILIM(CQUESTION,CDEF,N1,N2)
        IMPLICIT NONE
        CHARACTER*(*) CQUESTION
        CHARACTER*(*) CDEF
        INTEGER N1,N2
C
        INTEGER I,L1,L2
        INTEGER N
        INTEGER NERR
        INTEGER TRUEBEG,TRUELEN
        CHARACTER*1 C
        CHARACTER*255 CDUMMY
        CHARACTER*255 CADENA
        LOGICAL LECHO
        COMMON/BLKLECHO/LECHO
C------------------------------------------------------------------------------
        IF(N2.LT.N1) STOP 'ERROR: N2.LT.N1 in function: READILIM'
        NERR=0
10      L2=TRUELEN(CQUESTION)
        IF(L2.NE.0) WRITE(*,100) CQUESTION(1:L2)
        WRITE(CDUMMY,'(A1,I10,A5,I10,A1)') '(',N1,',...,',N2,')'
        CALL RMBLANK(CDUMMY,CDUMMY,L2)
        WRITE(*,100) ' '//CDUMMY(1:L2)
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
        READ(*,101,ERR=20) CADENA
        IF(TRUELEN(CADENA).EQ.0)THEN
          IF(CDEF.EQ.'@')THEN
            GOTO 10
          END IF
          CADENA=CDEF
        END IF
        DO I=1,TRUELEN(CADENA)
          C=CADENA(I:I)
          IF((INDEX('abcdefghijklmnopqrstuvwxyz',C).NE.0).OR.
     +     (INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ./',C).NE.0))THEN
            GOTO 20
          END IF
        END DO
        READ(CADENA,*,ERR=20) N
        READILIM=N
C
        IF((N.LT.N1).OR.(N.GT.N2))THEN
          WRITE(*,100) 'ERROR: invalid number. Valid range is ['
          WRITE(CDUMMY,*) N1
          CALL RMBLANK(CDUMMY,CDUMMY,L2)
          WRITE(*,100) CDUMMY(1:L2)//','
          WRITE(CDUMMY,*) N2
          CALL RMBLANK(CDUMMY,CDUMMY,L2)
          WRITE(*,101) CDUMMY(1:L2)//']. Try again.'
          NERR=NERR+1
          GOTO 22
        END IF
        IF(LECHO) WRITE(*,*) N
        RETURN
C------------------------------------------------------------------------------
20      WRITE(*,101) 'ERROR: invalid character(s) found in '//
     +   'number. Try again.'
        NERR=NERR+1
22      IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors in READILIM.'
        GOTO 10
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
        REAL FUNCTION READF(CQUESTION,CDEF)
        IMPLICIT NONE
        CHARACTER*(*) CQUESTION
        CHARACTER*(*) CDEF
C
        INTEGER I,L1,L2
        INTEGER NERR
        REAL F
        INTEGER TRUEBEG,TRUELEN
        CHARACTER*1 C
        CHARACTER*255 CADENA
        LOGICAL LECHO
        COMMON/BLKLECHO/LECHO
C------------------------------------------------------------------------------
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
        READ(*,101,ERR=20) CADENA
        IF(TRUELEN(CADENA).EQ.0)THEN
          IF(CDEF.EQ.'@')THEN
            GOTO 10
          END IF
          CADENA=CDEF
        END IF
        DO I=1,TRUELEN(CADENA)
          C=CADENA(I:I)
          IF((INDEX('abcfghijklmnoprstuvwxyz',C).NE.0).OR.
     +     (INDEX('ABCFGHIJKLMNOPRSTUVWXYZ/',C).NE.0))THEN
            GOTO 20
          END IF
        END DO
        READ(CADENA,*,ERR=20) F
        READF=F
        IF(LECHO) WRITE(*,*) F
        RETURN
20      WRITE(*,101) 'ERROR: invalid character(s) found in '//
     +   'number. Try again.'
        NERR=NERR+1
        IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors in READF.'
        GOTO 10
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
        SUBROUTINE RMBLANK(C1,C2,L)
        IMPLICIT NONE
        INTEGER L
        CHARACTER*(*) C1,C2
C
        INTEGER I,K,L0
C------------------------------------------------------------------------------
        K=0
        L0=LEN(C1)
        DO I=1,L0
          IF(C1(I:I).NE.CHAR(32))THEN
            K=K+1
            C2(K:K)=C1(I:I)
          END IF
        END DO
        L=K
        L0=LEN(C2)
        IF(L.LT.L0)THEN
          DO I=L+1,L0
            C2(I:I)=' '
          END DO
        END IF
        END
C
C******************************************************************************
C
        INTEGER FUNCTION TRUELEN(CADENA)
        IMPLICIT NONE
        CHARACTER*(*) CADENA
C
        INTEGER I,L
C------------------------------------------------------------------------------
        L=LEN(CADENA)
C
        DO I=L,1,-1
          IF(ICHAR(CADENA(I:I)).GT.32)THEN
            TRUELEN=I
            RETURN
          END IF
        END DO
        TRUELEN=0
        END
C
C******************************************************************************
C
        INTEGER FUNCTION TRUEBEG(CADENA)
        IMPLICIT NONE
        CHARACTER*(*) CADENA
C
        INTEGER I,L
C------------------------------------------------------------------------------
        L=LEN(CADENA)
C
        DO I=1,L
          IF(ICHAR(CADENA(I:I)).GT.32)THEN
            TRUEBEG=I
            RETURN
          END IF
        END DO
        TRUEBEG=0
        END
C
C******************************************************************************
C
        SUBROUTINE CHUPPER(CADENA)
        IMPLICIT NONE
        CHARACTER*(*) CADENA
C
        INTEGER I,N
C------------------------------------------------------------------------------
        DO I=1,LEN(CADENA)
          N=ICHAR(CADENA(I:I))
          IF((N.GE.97).AND.(N.LE.122)) CADENA(I:I)=CHAR(N-32)
        END DO
        END
