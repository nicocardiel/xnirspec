C Salva una imagen en formato REDUCEME
        SUBROUTINE SESCRREDU(OUTFILE,NCBUFF,LERR)
        IMPLICIT NONE
        CHARACTER*(*) OUTFILE
        INTEGER NCBUFF
        LOGICAL LERR
C
        INCLUDE 'dimensions.inc'
C
ccc        INTEGER TRUELEN
        REAL READF
C
        INTEGER I,J
        INTEGER NCHAR
        INTEGER NAXIS(2,NMAXBUFF)
        REAL STWV,DISP
        REAL AIRMASS,TIMEXPOS
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        CHARACTER*12 IDENTIFICATION
        CHARACTER*50 CDUMMY
ccc        CHARACTER*255 OBJECT
ccc        CHARACTER*255 FITSFILE
ccc        CHARACTER*255 COMMENT
C
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKREDUCEME/STWV,DISP,AIRMASS,TIMEXPOS
C------------------------------------------------------------------------------
        WRITE(CDUMMY,*) STWV
        STWV=READF('New STWV for output image',CDUMMY)
        WRITE(CDUMMY,*) DISP
        DISP=READF('New DISP for output image',CDUMMY)
C open file
        OPEN(10,FILE=OUTFILE,STATUS='NEW',FORM='UNFORMATTED')
C write header information
        IDENTIFICATION='abcdefghijkl'
        WRITE(10) IDENTIFICATION
        WRITE(10) NAXIS(2,NCBUFF),NAXIS(1,NCBUFF)
        WRITE(10) STWV,DISP
        WRITE(10) AIRMASS
        WRITE(10) TIMEXPOS
ccc        NCHAR=TRUELEN(OBJECT)
        IF(LERR)THEN
          NCHAR=8
          WRITE(10) NCHAR
          WRITE(10) ' @ERROR@'
        ELSE
          NCHAR=0
          WRITE(10) NCHAR
        END IF
ccc        IF(NCHAR.GT.0) WRITE(10) OBJECT(1:NCHAR)
ccc        NCHAR=TRUELEN(FITSFILE)
        NCHAR=0
        WRITE(10) NCHAR
ccc        IF(NCHAR.GT.0) WRITE(10) FITSFILE(1:NCHAR)
ccc        NCHAR=TRUELEN(COMMENT)
        NCHAR=0
        WRITE(10) NCHAR
ccc        IF(NCHAR.GT.0) WRITE(10) COMMENT(1:NCHAR)
C write data frame
        DO I=1,NAXIS(2,NCBUFF)
          WRITE(10) (IMAGEN(J,I,NCBUFF),J=1,NAXIS(1,NCBUFF))
        END DO
        CLOSE(10)
C end of subroutine
        END
