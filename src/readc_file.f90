!------------------------------------------------------------------------------
! CHARACTER*255 FUNCTION READC_FILE(CQUESTION,CDEF,CVAL)
!
! Return a file name entered by the user through the keyboard. The interface
! and behaviour are the same as in READC, except that, if the program has
! been compiled with readline support (see configure option --with-readline)
! and both stdin and stdout are terminals, the line is read with readline:
! TAB completes the file name up to the longest unambiguous prefix, a second
! TAB lists the possible completions, and previous entries can be recalled
! with the arrow keys. Otherwise (e.g. input redirected from a file) the line
! is read with a Fortran READ, as in READC.
!
! Wildcards ("*", "?") are not handled here: as with READC, the calling
! routine is responsible for listing the matching files.
!
! Parameters
! ----------
! CQUESTION : CHARACTER*(*)
!     Question displayed before the default value.
! CDEF : CHARACTER*(*)
!     Default value ('@' if there is no default).
! CVAL : CHARACTER*(*)
!     String with the valid characters ('@' if all characters are valid).
!
! Returns
! -------
! READC_FILE : CHARACTER*255
!     String entered by the user (or the default value).
!------------------------------------------------------------------------------
        CHARACTER*255 FUNCTION READC_FILE(CQUESTION,CDEF,CVAL)
        USE ISO_C_BINDING, ONLY: C_CHAR,C_INT,C_NULL_CHAR
        USE ISO_FORTRAN_ENV, ONLY: OUTPUT_UNIT
        IMPLICIT NONE
        CHARACTER*(*) CQUESTION
        CHARACTER*(*) CDEF
        CHARACTER*(*) CVAL
!
! C routine in freadline.c
        INTERFACE
          INTEGER(C_INT) FUNCTION FREADLINE(PROMPT,BUFFER,BUFLEN) &
           BIND(C,NAME='freadline')
          IMPORT C_CHAR,C_INT
          CHARACTER(KIND=C_CHAR) PROMPT(*)
          CHARACTER(KIND=C_CHAR) BUFFER(*)
          INTEGER(C_INT),VALUE :: BUFLEN
          END FUNCTION FREADLINE
        END INTERFACE
!
        INTEGER I,LQ,LD1,LD2,LP
        INTEGER TRUEBEG,TRUELEN
        INTEGER NERR
        INTEGER ISTATUS
        CHARACTER*255 CADENA
        CHARACTER*1024 CPROMPT
        LOGICAL LECHO
        COMMON/BLKLECHO/LECHO
!------------------------------------------------------------------------------
        NERR=0
! build the prompt (same text as in READC); LP is its length, including
! the trailing blank. LD1:LD2 is the non-blank part of the default value
! (LD1=0 if there is no default value).
10      LQ=TRUELEN(CQUESTION)
        IF(LQ.GT.0)THEN
          CPROMPT=CQUESTION(1:LQ)
        ELSE
          CPROMPT=' '
        END IF
        LP=LQ
        LD1=0
        LD2=0
        IF(CDEF.NE.'@')THEN
          LD1=TRUEBEG(CDEF)
          IF(LD1.NE.0)THEN
            LD2=TRUELEN(CDEF)
            CPROMPT(LP+1:)=' ['//CDEF(LD1:LD2)//'] ? '
            LP=LP+LD2-LD1+7
          END IF
        ELSE
          CPROMPT(LP+1:)='? '
          LP=LP+2
        END IF
! flush pending Fortran output (e.g. text written with FORMAT(A,$)) before
! readline writes the prompt through the C standard I/O library
        FLUSH(OUTPUT_UNIT)
        ISTATUS=FREADLINE(CPROMPT(1:LP)//C_NULL_CHAR,CADENA,LEN(CADENA))
        IF(ISTATUS.EQ.0)THEN                  !readline not used: as in READC
          IF(LP.GT.0) WRITE(*,100,ADVANCE='NO') CPROMPT(1:LP)
          READ(*,100,ERR=20) CADENA
        ELSEIF(ISTATUS.LT.0)THEN                       !EOF (Ctrl-D) in a tty
          WRITE(*,*)
          GOTO 20
        END IF
        IF(TRUELEN(CADENA).EQ.0)THEN
          IF(CDEF.EQ.'@')THEN
            GOTO 10
          END IF
          IF(LD1.GT.0)THEN
            CADENA=CDEF(LD1:LD2)
          END IF
        ELSEIF(CVAL.NE.'@')THEN
          DO I=1,TRUELEN(CADENA)
            IF(INDEX(CVAL,CADENA(I:I)).EQ.0)THEN
              WRITE(*,100) 'ERROR: invalid character(s). Try again.'
              NERR=NERR+1
              IF(NERR.GT.10)THEN
                INCLUDE 'deallocate_arrays.inc'
                STOP 'FATAL ERROR: too many errors in READC_FILE.'
              END IF
              GOTO 10
            END IF
          END DO
        END IF
        READC_FILE=CADENA
        IF(LECHO) WRITE(*,100) CADENA(1:TRUELEN(CADENA))
        RETURN
20      WRITE(*,100) 'ERROR: invalid entry. Try again.'
        NERR=NERR+1
        IF(NERR.GT.10)THEN
          INCLUDE 'deallocate_arrays.inc'
          STOP 'FATAL ERROR: too many errors in READC_FILE.'
        END IF
        GOTO 10
100     FORMAT(A)
        END
