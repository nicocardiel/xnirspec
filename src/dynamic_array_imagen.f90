MODULE Dynamic_Array_IMAGEN
    REAL, DIMENSION(:, :, :), ALLOCATABLE :: IMAGEN
END MODULE Dynamic_Array_IMAGEN

!------------------------------------------------------------------------------

SUBROUTINE Initialize_Dynamic_Array_IMAGEN
USE Dynamic_Array_IMAGEN

IMPLICIT NONE

! Include dimensions
INCLUDE 'dimensions.inc'

! Declare local variables
INTEGER :: AllocateStatus


! Allocate storage for array IMAGEN
ALLOCATE (IMAGEN(NXMAX, NYMAX, NMAXBUFF), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array IMAGEN ***"

END SUBROUTINE Initialize_Dynamic_Array_IMAGEN

!------------------------------------------------------------------------------

SUBROUTINE Deallocate_Array_IMAGEN
USE Dynamic_Array_IMAGEN

IMPLICIT NONE

! Declare local variables
INTEGER :: DeAllocateStatus

! Deallocate storage for array IMAGEN
DEALLOCATE(IMAGEN, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array IMAGEN ***"

END SUBROUTINE Deallocate_Array_IMAGEN
