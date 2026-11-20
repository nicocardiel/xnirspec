!------------------------------------------------------------------------------
! Driver for test_readc_file.sh: calls READC_FILE three times and prints each
! returned value as "RESULT:<value>|". The script feeds stdin from a file, so
! READC_FILE uses the plain Fortran READ (no readline).
!------------------------------------------------------------------------------
        PROGRAM READC_FILE_DRIVER
        IMPLICIT NONE
!
        CHARACTER*255 READC_FILE
        INTEGER TRUELEN
        CHARACTER*255 C
        LOGICAL LECHO
        COMMON/BLKLECHO/LECHO
!------------------------------------------------------------------------------
        LECHO=.FALSE.
! default value (empty answer)
        C=READC_FILE('First file name','default.fits','@')
        WRITE(*,'(A)') 'RESULT:'//C(1:TRUELEN(C))//'|'
! explicit answer, no default
        C=READC_FILE('Second file name','@','@')
        WRITE(*,'(A)') 'RESULT:'//C(1:TRUELEN(C))//'|'
! restricted set of valid characters: the first answer is rejected
        C=READC_FILE('Third answer','@','abc')
        WRITE(*,'(A)') 'RESULT:'//C(1:TRUELEN(C))//'|'
        END PROGRAM READC_FILE_DRIVER
