MODULE Dynamic_Array_EPIXEL
    REAL, DIMENSION(:), ALLOCATABLE :: EPIXEL
END MODULE Dynamic_Array_EPIXEL

!------------------------------------------------------------------------------

SUBROUTINE Initialize_Dynamic_Array_EPIXEL
USE Dynamic_Array_EPIXEL

IMPLICIT NONE

! Include dimensions
INCLUDE 'largest.inc'

! Declare local variables
INTEGER :: AllocateStatus


! Allocate storage for array EPIXEL
ALLOCATE (EPIXEL(NXYMAX*NXYMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array EPIXEL ***"

END SUBROUTINE Initialize_Dynamic_Array_EPIXEL

!------------------------------------------------------------------------------

SUBROUTINE Deallocate_Array_EPIXEL
USE Dynamic_Array_EPIXEL

IMPLICIT NONE

! Declare local variables
INTEGER :: DeAllocateStatus

! Deallocate storage for array EPIXEL
DEALLOCATE(EPIXEL, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array EPIXEL ***"

END SUBROUTINE Deallocate_Array_EPIXEL
