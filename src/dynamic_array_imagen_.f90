MODULE Dynamic_Array_IMAGEN_
    REAL, DIMENSION(:, :), ALLOCATABLE :: IMAGEN_
END MODULE Dynamic_Array_IMAGEN_

!------------------------------------------------------------------------------

SUBROUTINE Initialize_Dynamic_Array_IMAGEN_
USE Dynamic_Array_IMAGEN_

IMPLICIT NONE

! Include dimensions
INCLUDE 'dimensions.inc'

! Declare local variables
INTEGER :: AllocateStatus


! Allocate storage for array IMAGEN_
ALLOCATE (IMAGEN_(NXMAX, NYMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array IMAGEN_ ***"

END SUBROUTINE Initialize_Dynamic_Array_IMAGEN_

!------------------------------------------------------------------------------

SUBROUTINE Deallocate_Array_IMAGEN_
USE Dynamic_Array_IMAGEN_

IMPLICIT NONE

! Declare local variables
INTEGER :: DeAllocateStatus

! Deallocate storage for array IMAGEN
DEALLOCATE(IMAGEN_, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array IMAGEN_ ***"

END SUBROUTINE Deallocate_Array_IMAGEN_
