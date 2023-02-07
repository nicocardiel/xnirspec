MODULE Dynamic_Array_IPIXEL
    INTEGER, DIMENSION(:), ALLOCATABLE :: IPIXEL
END MODULE Dynamic_Array_IPIXEL

!------------------------------------------------------------------------------

SUBROUTINE Initialize_Dynamic_Array_IPIXEL
USE Dynamic_Array_IPIXEL

IMPLICIT NONE

! Include dimensions
INCLUDE 'largest.inc'

! Declare local variables
INTEGER :: AllocateStatus


! Allocate storage for array IPIXEL
ALLOCATE (IPIXEL(NXYMAX*NXYMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array IPIXEL ***"

END SUBROUTINE Initialize_Dynamic_Array_IPIXEL

!------------------------------------------------------------------------------

SUBROUTINE Deallocate_Array_IPIXEL
USE Dynamic_Array_IPIXEL

IMPLICIT NONE

! Declare local variables
INTEGER :: DeAllocateStatus

! Deallocate storage for array IPIXEL
DEALLOCATE(IPIXEL, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array IPIXEL ***"

END SUBROUTINE Deallocate_Array_IPIXEL
