! Salva una imagen en formato REDUCEME
        SUBROUTINE SESCRREDU(OUTFILE,NCBUFF,LERR)
        USE Dynamic_Array_IMAGEN
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
! subroutine arguments
        CHARACTER*(*) OUTFILE
        INTEGER NCBUFF
        LOGICAL LERR
!
        INCLUDE 'dimensions.inc'
!
!!!        INTEGER TRUELEN
        REAL READF
!
        INTEGER I,J
        INTEGER NCHAR
        INTEGER NAXIS(2,NMAXBUFF)
        REAL STWV,DISP
        REAL AIRMASS,TIMEXPOS
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        CHARACTER*12 IDENTIFICATION
        CHARACTER*50 CDUMMY
!!!        CHARACTER*255 OBJECT
!!!        CHARACTER*255 FITSFILE
!!!        CHARACTER*255 COMMENT
!
        COMMON/BLKNAXIS/NAXIS
!delete COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKREDUCEME/STWV,DISP,AIRMASS,TIMEXPOS
!------------------------------------------------------------------------------
        WRITE(CDUMMY,*) STWV
        STWV=READF('New STWV for output image',CDUMMY)
        WRITE(CDUMMY,*) DISP
        DISP=READF('New DISP for output image',CDUMMY)
! open file
        OPEN(10,FILE=OUTFILE,STATUS='NEW',FORM='UNFORMATTED')
! write header information
        IDENTIFICATION='abcdefghijkl'
        WRITE(10) IDENTIFICATION
        WRITE(10) NAXIS(2,NCBUFF),NAXIS(1,NCBUFF)
        WRITE(10) STWV,DISP
        WRITE(10) AIRMASS
        WRITE(10) TIMEXPOS
!!!        NCHAR=TRUELEN(OBJECT)
        IF(LERR)THEN
          NCHAR=8
          WRITE(10) NCHAR
          WRITE(10) ' @ERROR@'
        ELSE
          NCHAR=0
          WRITE(10) NCHAR
        END IF
!!!        IF(NCHAR.GT.0) WRITE(10) OBJECT(1:NCHAR)
!!!        NCHAR=TRUELEN(FITSFILE)
        NCHAR=0
        WRITE(10) NCHAR
!!!        IF(NCHAR.GT.0) WRITE(10) FITSFILE(1:NCHAR)
!!!        NCHAR=TRUELEN(COMMENT)
        NCHAR=0
        WRITE(10) NCHAR
!!!        IF(NCHAR.GT.0) WRITE(10) COMMENT(1:NCHAR)
! write data frame
        DO I=1,NAXIS(2,NCBUFF)
          WRITE(10) (IMAGEN(J,I,NCBUFF),J=1,NAXIS(1,NCBUFF))
        END DO
        CLOSE(10)
! end of subroutine
        END
