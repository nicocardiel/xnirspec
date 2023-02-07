MODULE Dynamic_Array_JPIXEL
    INTEGER, DIMENSION(:), ALLOCATABLE :: JPIXEL
END MODULE Dynamic_Array_JPIXEL

!------------------------------------------------------------------------------

SUBROUTINE Initialize_Dynamic_Array_JPIXEL
USE Dynamic_Array_JPIXEL

IMPLICIT NONE

! Include dimensions
INCLUDE 'largest.inc'

! Declare local variables
INTEGER :: AllocateStatus


! Allocate storage for array JPIXEL
ALLOCATE (JPIXEL(NXYMAX*NXYMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array JPIXEL ***"

END SUBROUTINE Initialize_Dynamic_Array_JPIXEL

!------------------------------------------------------------------------------

SUBROUTINE Deallocate_Array_JPIXEL
USE Dynamic_Array_JPIXEL

IMPLICIT NONE

! Declare local variables
INTEGER :: DeAllocateStatus

! Deallocate storage for array JPIXEL
DEALLOCATE(JPIXEL, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array JPIXEL ***"

END SUBROUTINE Deallocate_Array_JPIXEL
