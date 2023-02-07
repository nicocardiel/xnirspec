MODULE Dynamic_Array_COEFFBL
    REAL, DIMENSION(:, :), ALLOCATABLE :: COEFFBL
END MODULE Dynamic_Array_COEFFBL

!------------------------------------------------------------------------------

SUBROUTINE Initialize_Dynamic_Array_COEFFBL
USE Dynamic_Array_COEFFBL

IMPLICIT NONE

! Include dimensions
INCLUDE 'largest.inc'

! Declare local variables
INTEGER :: AllocateStatus


! Allocate storage for array COEFFBL
ALLOCATE (COEFFBL(20, NXYMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array COEFFBL ***"

END SUBROUTINE Initialize_Dynamic_Array_COEFFBL

!------------------------------------------------------------------------------

SUBROUTINE Deallocate_Array_COEFFBL
USE Dynamic_Array_COEFFBL

IMPLICIT NONE

! Declare local variables
INTEGER :: DeAllocateStatus

! Deallocate storage for array COEFFBL
DEALLOCATE(COEFFBL, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array COEFFBL ***"

END SUBROUTINE Deallocate_Array_COEFFBL
