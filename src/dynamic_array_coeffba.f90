MODULE Dynamic_Array_COEFFBA
    REAL, DIMENSION(:, :), ALLOCATABLE :: COEFFBA
END MODULE Dynamic_Array_COEFFBA

!------------------------------------------------------------------------------

SUBROUTINE Initialize_Dynamic_Array_COEFFBA
USE Dynamic_Array_COEFFBA

IMPLICIT NONE

! Include dimensions
INCLUDE 'largest.inc'

! Declare local variables
INTEGER :: AllocateStatus


! Allocate storage for array COEFFBA
ALLOCATE (COEFFBA(20, NXYMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array COEFFBA ***"

END SUBROUTINE Initialize_Dynamic_Array_COEFFBA

!------------------------------------------------------------------------------

SUBROUTINE Deallocate_Array_COEFFBA
USE Dynamic_Array_COEFFBA

IMPLICIT NONE

! Declare local variables
INTEGER :: DeAllocateStatus

! Deallocate storage for array COEFFBA
DEALLOCATE(COEFFBA, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array COEFFBA ***"

END SUBROUTINE Deallocate_Array_COEFFBA
