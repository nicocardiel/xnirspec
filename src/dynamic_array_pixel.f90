MODULE Dynamic_Array_PIXEL
    REAL, DIMENSION(:), ALLOCATABLE :: PIXEL
END MODULE Dynamic_Array_PIXEL

!------------------------------------------------------------------------------

SUBROUTINE Initialize_Dynamic_Array_PIXEL
USE Dynamic_Array_PIXEL

IMPLICIT NONE

! Include dimensions
INCLUDE 'largest.inc'

! Declare local variables
INTEGER :: AllocateStatus


! Allocate storage for array PIXEL
ALLOCATE (PIXEL(NXYMAX*NXYMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array PIXEL ***"

END SUBROUTINE Initialize_Dynamic_Array_PIXEL

!------------------------------------------------------------------------------

SUBROUTINE Deallocate_Array_PIXEL
USE Dynamic_Array_PIXEL

IMPLICIT NONE

! Declare local variables
INTEGER :: DeAllocateStatus

! Deallocate storage for array PIXEL
DEALLOCATE(PIXEL, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array PIXEL ***"

END SUBROUTINE Deallocate_Array_PIXEL
