        SUBROUTINE TRY_REDUCEME(INFILE,ISTATUS,NAXIS_)
        USE Dynamic_Array_IMAGEN_
        IMPLICIT NONE
        INCLUDE 'interface_imagen_.inc'
! subroutine arguments
        CHARACTER*(*) INFILE
        INTEGER ISTATUS
        INTEGER NAXIS_(0:3)
!
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
!
        INTEGER I,J
        INTEGER NCHAR
        INTEGER NSCAN,NCHAN
        REAL STWV,DISP
        REAL AIRMASS,TIMEXPOS
!delete REAL IMAGEN_(NXYMAX,NXYMAX)
        CHARACTER*12 IDENTIFICATION
        CHARACTER*255 OBJECT,FITSFILE,COMMENT
!
!delete COMMON/BLKIMAGEN1_/IMAGEN_
        COMMON/BLKREDUCEME/STWV,DISP,AIRMASS,TIMEXPOS
!------------------------------------------------------------------------------
        ISTATUS=0 !de momento todo bien
!
        OPEN(88,FILE=INFILE,STATUS='OLD',FORM='UNFORMATTED')
        READ(88) IDENTIFICATION
        IF(IDENTIFICATION.NE.'abcdefghijkl')THEN
          CLOSE(88)
          ISTATUS=1
          RETURN
        END IF
        READ(88) NSCAN,NCHAN
        IF(NSCAN.GT.NYMAX)THEN
          WRITE(*,101) 'ERROR: NSCAN.GT.NYMAX'
          ISTATUS=1
        END IF
        IF(NCHAN.GT.NXMAX)THEN
          WRITE(*,101) 'ERROR: NCHAN.GT.NXMAX'
          ISTATUS=1
        END IF
        READ(88) STWV,DISP
        READ(88) AIRMASS
        READ(88) TIMEXPOS
        READ(88) NCHAR
        IF(NCHAR.GT.0) READ(88) OBJECT(1:NCHAR)
        READ(88) NCHAR
        IF(NCHAR.GT.0) READ(88) FITSFILE(1:NCHAR)
        READ(88) NCHAR
        IF(NCHAR.GT.0) READ(88) COMMENT(1:NCHAR)
! read data frame
        IF(ISTATUS.EQ.0)THEN
          DO I=1,NSCAN
            READ(88) (IMAGEN_(J,I),J=1,NCHAN)
          END DO
        ELSE
          DO I=1,NYMAX
            DO J=1,NXMAX
              IMAGEN_(J,I)=0.0
            END DO
          END DO
        END IF
        CLOSE(88)
        NAXIS_(0)=2
        NAXIS_(1)=NCHAN
        NAXIS_(2)=NSCAN
        NAXIS_(3)=1
!
101     FORMAT(A)
        END
