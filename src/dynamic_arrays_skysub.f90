MODULE Dynamic_Arrays_SKYSUB
    REAL, DIMENSION(:), ALLOCATABLE :: XCUT
    REAL, DIMENSION(:), ALLOCATABLE :: XCUT2
    REAL, DIMENSION(:), ALLOCATABLE :: XCUT2_ERR
    REAL, DIMENSION(:), ALLOCATABLE :: YCUT
    REAL, DIMENSION(:), ALLOCATABLE :: YCUT_ERR
    REAL, DIMENSION(:), ALLOCATABLE :: XP
END MODULE Dynamic_Arrays_SKYSUB

!------------------------------------------------------------------------------

SUBROUTINE Initialize_Dynamic_Arrays_SKYSUB
USE Dynamic_Arrays_SKYSUB

IMPLICIT NONE

! Include dimensions
INCLUDE 'largest.inc'

! Declare local variables
INTEGER :: AllocateStatus

! Allocate storage for array XCUT
ALLOCATE (XCUT(NXYMAX*NOVERSAMPMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array XCUT ***"

! Allocate storage for array XCUT2
ALLOCATE (XCUT2(NXYMAX*NOVERSAMPMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array XCUT2 ***"

! Allocate storage for array XCUT2_ERR
ALLOCATE (XCUT2_ERR(NXYMAX*NOVERSAMPMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array XCUT2_ERR ***"

! Allocate storage for array YCUT
ALLOCATE (YCUT(NXYMAX*NOVERSAMPMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array YCUT ***"

! Allocate storage for array YCUT_ERR
ALLOCATE (YCUT_ERR(NXYMAX*NOVERSAMPMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array YCUT_ERR ***"

! Allocate storage for array XP
ALLOCATE (XP(NXYMAX*NOVERSAMPMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array XP ***"

END SUBROUTINE Initialize_Dynamic_Arrays_SKYSUB

!------------------------------------------------------------------------------

SUBROUTINE Deallocate_Arrays_SKYSUB
USE Dynamic_Arrays_SKYSUB

IMPLICIT NONE

! Declare local variables
INTEGER :: DeAllocateStatus

! Deallocate storage for array XCUT
DEALLOCATE(XCUT, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array XCUT ***"

! Deallocate storage for array XCUT2
DEALLOCATE(XCUT2, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array XCUT2 ***"

! Deallocate storage for array XCUT2_ERR
DEALLOCATE(XCUT2_ERR, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array XCUT2_ERR ***"

! Deallocate storage for array YCUT
DEALLOCATE(YCUT, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array YCUT ***"

! Deallocate storage for array YCUT_ERR
DEALLOCATE(YCUT_ERR, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array YCUT_ERR ***"

! Deallocate storage for array XP
DEALLOCATE(XP, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array XP ***"

END SUBROUTINE Deallocate_Arrays_SKYSUB
