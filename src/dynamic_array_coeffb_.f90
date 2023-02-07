MODULE Dynamic_Array_COEFFB_
    REAL, DIMENSION(:, :), ALLOCATABLE :: COEFFB_
END MODULE Dynamic_Array_COEFFB_

!------------------------------------------------------------------------------

SUBROUTINE Initialize_Dynamic_Array_COEFFB_
USE Dynamic_Array_COEFFB_

IMPLICIT NONE

! Include dimensions
INCLUDE 'largest.inc'

! Declare local variables
INTEGER :: AllocateStatus


! Allocate storage for array COEFFB_
ALLOCATE (COEFFB_(20, NXYMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array COEFFB_ ***"

END SUBROUTINE Initialize_Dynamic_Array_COEFFB_

!------------------------------------------------------------------------------

SUBROUTINE Deallocate_Array_COEFFB_
USE Dynamic_Array_COEFFB_

IMPLICIT NONE

! Declare local variables
INTEGER :: DeAllocateStatus

! Deallocate storage for array COEFFB_
DEALLOCATE(COEFFB_, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array COEFFB_ ***"

END SUBROUTINE Deallocate_Array_COEFFB_
