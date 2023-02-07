MODULE Dynamic_Array_FPIXUSED
    REAL, DIMENSION(:, :), ALLOCATABLE :: FPIXUSED
END MODULE Dynamic_Array_FPIXUSED

!------------------------------------------------------------------------------

SUBROUTINE Initialize_Dynamic_Array_FPIXUSED
USE Dynamic_Array_FPIXUSED

IMPLICIT NONE

! Include dimensions
INCLUDE 'dimensions.inc'

! Declare local variables
INTEGER :: AllocateStatus


! Allocate storage for array FPIXUSED
ALLOCATE (FPIXUSED(NXMAXB9, NYMAXB9), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array FPIXUSED ***"

END SUBROUTINE Initialize_Dynamic_Array_FPIXUSED

!------------------------------------------------------------------------------

SUBROUTINE Deallocate_Array_FPIXUSED
USE Dynamic_Array_FPIXUSED

IMPLICIT NONE

! Declare local variables
INTEGER :: DeAllocateStatus

! Deallocate storage for array FPIXUSED
DEALLOCATE(FPIXUSED, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array FPIXUSED ***"

END SUBROUTINE Deallocate_Array_FPIXUSED
