!------------------------------------------------------------------------------
! MODULE TEST_UTILS
!
! Minimal helpers for the unit tests run with "make check". Each CHECK_*
! routine prints one PASS/FAIL line and counts the failures; TEST_SUMMARY
! prints the totals and terminates the program with a non-zero exit status
! (ERROR STOP) if any check has failed, which is how the Automake test
! harness detects a failing test. DELETE_FILE removes temporary files.
!------------------------------------------------------------------------------
        MODULE TEST_UTILS
        IMPLICIT NONE
        PRIVATE
        PUBLIC :: CHECK_REAL,CHECK_INT,CHECK_TRUE,CHECK_STR,TEST_SUMMARY
        PUBLIC :: DELETE_FILE
!
        INTEGER, SAVE :: NCHECKS=0
        INTEGER, SAVE :: NFAILS=0
!
        CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE CHECK_REAL(NAME,GOT,EXPECTED,RTOL,ATOL)
!
! Check that ABS(GOT-EXPECTED) <= ATOL + RTOL*ABS(EXPECTED) (the same
! criterion as numpy.isclose).
!
! Parameters
! ----------
! NAME : CHARACTER*(*)
!     Description of the check.
! GOT, EXPECTED : REAL
!     Computed and expected values.
! RTOL, ATOL : REAL
!     Relative and absolute tolerances.
!------------------------------------------------------------------------------
        SUBROUTINE CHECK_REAL(NAME,GOT,EXPECTED,RTOL,ATOL)
        CHARACTER*(*), INTENT(IN) :: NAME
        REAL, INTENT(IN) :: GOT,EXPECTED,RTOL,ATOL
!
        NCHECKS=NCHECKS+1
        IF(ABS(GOT-EXPECTED).LE.ATOL+RTOL*ABS(EXPECTED))THEN
          WRITE(*,'(A)') 'PASS: '//NAME
        ELSE
          NFAILS=NFAILS+1
          WRITE(*,'(A)') 'FAIL: '//NAME
          WRITE(*,'(A,ES16.8,A,ES16.8)') '      got: ',GOT, &
           '   expected: ',EXPECTED
        END IF
        END SUBROUTINE CHECK_REAL
!------------------------------------------------------------------------------
! SUBROUTINE CHECK_INT(NAME,GOT,EXPECTED)
!
! Check that two integers are equal.
!
! Parameters
! ----------
! NAME : CHARACTER*(*)
!     Description of the check.
! GOT, EXPECTED : INTEGER
!     Computed and expected values.
!------------------------------------------------------------------------------
        SUBROUTINE CHECK_INT(NAME,GOT,EXPECTED)
        CHARACTER*(*), INTENT(IN) :: NAME
        INTEGER, INTENT(IN) :: GOT,EXPECTED
!
        NCHECKS=NCHECKS+1
        IF(GOT.EQ.EXPECTED)THEN
          WRITE(*,'(A)') 'PASS: '//NAME
        ELSE
          NFAILS=NFAILS+1
          WRITE(*,'(A)') 'FAIL: '//NAME
          WRITE(*,'(A,I0,A,I0)') '      got: ',GOT,'   expected: ',EXPECTED
        END IF
        END SUBROUTINE CHECK_INT
!------------------------------------------------------------------------------
! SUBROUTINE CHECK_TRUE(NAME,CONDITION)
!
! Check that a logical condition is true.
!
! Parameters
! ----------
! NAME : CHARACTER*(*)
!     Description of the check.
! CONDITION : LOGICAL
!     Condition that must be .TRUE.
!------------------------------------------------------------------------------
        SUBROUTINE CHECK_TRUE(NAME,CONDITION)
        CHARACTER*(*), INTENT(IN) :: NAME
        LOGICAL, INTENT(IN) :: CONDITION
!
        NCHECKS=NCHECKS+1
        IF(CONDITION)THEN
          WRITE(*,'(A)') 'PASS: '//NAME
        ELSE
          NFAILS=NFAILS+1
          WRITE(*,'(A)') 'FAIL: '//NAME
        END IF
        END SUBROUTINE CHECK_TRUE
!------------------------------------------------------------------------------
! SUBROUTINE CHECK_STR(NAME,GOT,EXPECTED)
!
! Check that two strings are equal (trailing blanks are ignored, as in the
! Fortran comparison of character strings).
!
! Parameters
! ----------
! NAME : CHARACTER*(*)
!     Description of the check.
! GOT, EXPECTED : CHARACTER*(*)
!     Computed and expected strings.
!------------------------------------------------------------------------------
        SUBROUTINE CHECK_STR(NAME,GOT,EXPECTED)
        CHARACTER*(*), INTENT(IN) :: NAME
        CHARACTER*(*), INTENT(IN) :: GOT,EXPECTED
!
        NCHECKS=NCHECKS+1
        IF(GOT.EQ.EXPECTED)THEN
          WRITE(*,'(A)') 'PASS: '//NAME
        ELSE
          NFAILS=NFAILS+1
          WRITE(*,'(A)') 'FAIL: '//NAME
          WRITE(*,'(A)') '      got: "'//TRIM(GOT)//'"   expected: "'// &
           TRIM(EXPECTED)//'"'
        END IF
        END SUBROUTINE CHECK_STR
!------------------------------------------------------------------------------
! SUBROUTINE DELETE_FILE(FILENAME)
!
! Delete a file if it exists (nothing is done otherwise).
!
! Parameters
! ----------
! FILENAME : CHARACTER*(*)
!     Name of the file to be deleted.
!------------------------------------------------------------------------------
        SUBROUTINE DELETE_FILE(FILENAME)
        CHARACTER*(*), INTENT(IN) :: FILENAME
!
        INTEGER IUNIT
        LOGICAL LEXIST
!
        INQUIRE(FILE=FILENAME,EXIST=LEXIST)
        IF(LEXIST)THEN
          OPEN(NEWUNIT=IUNIT,FILE=FILENAME,STATUS='OLD')
          CLOSE(IUNIT,STATUS='DELETE')
        END IF
        END SUBROUTINE DELETE_FILE
!------------------------------------------------------------------------------
! SUBROUTINE TEST_SUMMARY()
!
! Print the number of checks and failures. Terminate with ERROR STOP (non-zero
! exit status) if any check has failed.
!------------------------------------------------------------------------------
        SUBROUTINE TEST_SUMMARY()
        WRITE(*,'(A,I0,A,I0,A)') 'Summary: ',NCHECKS,' checks, ', &
         NFAILS,' failures'
        IF(NFAILS.GT.0) ERROR STOP 1
        END SUBROUTINE TEST_SUMMARY
!
        END MODULE TEST_UTILS
