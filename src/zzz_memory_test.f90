MODULE Dynamic_Arrays
    REAL, DIMENSION(:, :, :), ALLOCATABLE :: IMAGEN
    REAL, DIMENSION(:, :), ALLOCATABLE :: IMAGEN_
END MODULE Dynamic_Arrays

!------------------------------------------------------------------------------

SUBROUTINE Initialize_Dynamic_Arrays
USE Dynamic_Arrays

IMPLICIT NONE

! Include dimensions
INCLUDE 'dimensions.inc'

! Declare local variables
INTEGER :: AllocateStatus

! Allocate storage for array IMAGEN
ALLOCATE (IMAGEN(NXMAX, NYMAX, NMAXBUFF), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array IMAGEN ***"

! Allocate storage for array IMAGEN_
ALLOCATE (IMAGEN_(NXMAX, NYMAX), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory defining the array IMAGEN_ ***"

END SUBROUTINE Initialize_Dynamic_Arrays

!------------------------------------------------------------------------------

SUBROUTINE Deallocate_Arrays
USE Dynamic_Arrays

IMPLICIT NONE

! Declare local variables
INTEGER :: DeAllocateStatus

! Deallocate storage for array IMAGEN
DEALLOCATE(IMAGEN, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array IMAGEN ***"

! Deallocate storage for array IMAGEN_
DEALLOCATE(IMAGEN_, STAT = DeAllocateStatus)
IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating the array IMAGEN_ ***"

END SUBROUTINE Deallocate_Arrays

!------------------------------------------------------------------------------

PROGRAM zzz_MemoryTest
USE Dynamic_Arrays
IMPLICIT NONE

INTERFACE
  SUBROUTINE Initialize_Dynamic_Arrays
  END SUBROUTINE Initialize_Dynamic_Arrays
END INTERFACE

INTERFACE
  SUBROUTINE Deallocate_Arrays
  END SUBROUTINE Deallocate_Arrays
END INTERFACE

integer :: i, j

CALL Initialize_Dynamic_Arrays

write(*,*) 'Details of array IMAGEN:'
write(*,*) 'rank:', rank(IMAGEN)
write(*,*) 'shape:', shape(IMAGEN)
write(*,*) 'total size: ', size(IMAGEN)
write(*,*) 'size along axis 1: ', size(IMAGEN, 1)
write(*,*) 'size along axis 2: ', size(IMAGEN, 2)
write(*,*) 'size along axis 3: ', size(IMAGEN, 3)

write(*,*) 'Details of array IMAGEN_:'
write(*,*) 'rank:', rank(IMAGEN_)
write(*,*) 'shape:', shape(IMAGEN_)
write(*,*) 'total size: ', size(IMAGEN_)
write(*,*) 'size along axis 1: ', size(IMAGEN, 1)
write(*,*) 'size along axis 3: ', size(IMAGEN, 2)

CALL Modify_Arrays

do i=1,2
  do j=1, 2
    write(*,*) 'imagen(',j,',',i,'): ', imagen(j, i, 1)
  end do
end do

CALL Deallocate_Arrays

END PROGRAM zzz_MemoryTest

!------------------------------------------------------------------------------

SUBROUTINE Modify_Arrays
USE Dynamic_Arrays
IMPLICIT NONE

IMAGEN(1, 1, 1) = 10
IMAGEN_(1, 1) = 20

END SUBROUTINE modify_arrays
