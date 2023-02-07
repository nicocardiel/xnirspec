MODULE Dynamic_Array_COEFFBA_
    REAL, DIMENSION(:, :), ALLOCATABLE :: COEFFBA_
END MODULE Dynamic_Array_COEFFBA_

!------------------------------------------------------------------------------

SUBROUTINE Initialize_Dynamic_Array_COEFFBA_
USE Dynamic_Array_COEFFBA_

IMPLICIT NONE

! Include dimensions
INCLUDE 'largest.inc'

! Declare local variables
INTEGER :: AllocateStatus


! Allocate storage for array COEFFBA_
ALLOCATE (COEFFBA_(20, NXYMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array COEFFBA_ ***"

END SUBROUTINE Initialize_Dynamic_Array_COEFFBA_

!------------------------------------------------------------------------------

SUBROUTINE Deallocate_Array_COEFFBA_
USE Dynamic_Array_COEFFBA_

IMPLICIT NONE

! Declare local variables
INTEGER :: DeAllocateStatus

! Deallocate storage for array COEFFBA_
DEALLOCATE(COEFFBA_, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array COEFFBA_ ***"

END SUBROUTINE Deallocate_Array_COEFFBA_
