MODULE Dynamic_Arrays_XYPB
    REAL, DIMENSION(:, :), ALLOCATABLE :: XPB
    REAL, DIMENSION(:, :), ALLOCATABLE :: YPB
END MODULE Dynamic_Arrays_XYPB

!------------------------------------------------------------------------------

SUBROUTINE Initialize_Dynamic_Arrays_XYPB
USE Dynamic_Arrays_XYPB

IMPLICIT NONE

! Include dimensions
INCLUDE 'dimensions.inc'
INCLUDE 'largest.inc'

! Declare local variables
INTEGER :: AllocateStatus

! Allocate storage for array XPB
ALLOCATE (XPB(NXYMAX*NOVERSAMPMAX,NMAXBUFF), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array XPB ***"

! Allocate storage for array YPB
ALLOCATE (YPB(NXYMAX*NOVERSAMPMAX,NMAXBUFF), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array YPB ***"

END SUBROUTINE Initialize_Dynamic_Arrays_XYPB

!------------------------------------------------------------------------------

SUBROUTINE Deallocate_Arrays_XYPB
USE Dynamic_Arrays_XYPB

IMPLICIT NONE

! Declare local variables
INTEGER :: DeAllocateStatus

! Deallocate storage for array XPB
DEALLOCATE(XPB, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array XPB ***"

! Deallocate storage for array YPB
DEALLOCATE(YPB, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array YPB ***"

END SUBROUTINE Deallocate_Arrays_XYPB
