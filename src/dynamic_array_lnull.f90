MODULE Dynamic_Array_LNULL
    LOGICAL, DIMENSION(:, :, :), ALLOCATABLE :: LNULL
END MODULE Dynamic_Array_LNULL

!------------------------------------------------------------------------------

SUBROUTINE Initialize_Dynamic_Array_LNULL
USE Dynamic_Array_LNULL

IMPLICIT NONE

! Include dimensions
INCLUDE 'dimensions.inc'

! Declare local variables
INTEGER :: AllocateStatus


! Allocate storage for array LNULL
ALLOCATE (LNULL(NXMAX, NYMAX, NMAXBUFF), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array LNULL ***"

END SUBROUTINE Initialize_Dynamic_Array_LNULL

!------------------------------------------------------------------------------

SUBROUTINE Deallocate_Array_LNULL
USE Dynamic_Array_LNULL

IMPLICIT NONE

! Declare local variables
INTEGER :: DeAllocateStatus

! Deallocate storage for array LNULL
DEALLOCATE(LNULL, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array LNULL ***"

END SUBROUTINE Deallocate_Array_LNULL
