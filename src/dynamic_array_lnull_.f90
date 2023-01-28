MODULE Dynamic_Array_LNULL_
    LOGICAL, DIMENSION(:, :), ALLOCATABLE :: LNULL_
END MODULE Dynamic_Array_LNULL_

!------------------------------------------------------------------------------

SUBROUTINE Initialize_Dynamic_Array_LNULL_
USE Dynamic_Array_LNULL_

IMPLICIT NONE

! Include dimensions
INCLUDE 'largest.inc'

! Declare local variables
INTEGER :: AllocateStatus


! Allocate storage for array LNULL_
ALLOCATE (LNULL_(NXYMAX, NXYMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array LNULL_ ***"

END SUBROUTINE Initialize_Dynamic_Array_LNULL_

!------------------------------------------------------------------------------

SUBROUTINE Deallocate_Array_LNULL_
USE Dynamic_Array_LNULL_

IMPLICIT NONE

! Declare local variables
INTEGER :: DeAllocateStatus

! Deallocate storage for array LNULL_
DEALLOCATE(LNULL_, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array LNULL_ ***"

END SUBROUTINE Deallocate_Array_LNULL_
