C Devuelve la dimension maxima que debe tener cada frame en el mosaico con 
C las imagenes del box9, asi como las dimensiones individuales de cada frame
        SUBROUTINE MAXDIMB9(NCBUFF,INFILE,NSIZE)
        IMPLICIT NONE
C
        INTEGER NCBUFF
        CHARACTER*(*) INFILE
        INTEGER NSIZE
C
        INCLUDE 'dimensions.inc'
C
C
        INTEGER TRUELEN
C
        INTEGER NAXISFRAME(2,9,NMAXBUFF)
        INTEGER I,NFRAMES
        INTEGER IUNIT,IREADWRITE,BLOCKSIZE,ISTATUS
        INTEGER NAXIS_(0:2),NFOUND
        CHARACTER*50 COMMENT
        CHARACTER*80 FILEFITS
        LOGICAL LOGFILE
C
        COMMON/BLKNAXISFRAME/NAXISFRAME
C------------------------------------------------------------------------------
        IUNIT=80
        IREADWRITE=0
        ISTATUS=0
C
        NSIZE=0
        OPEN(11,FILE=INFILE,STATUS='OLD',FORM='FORMATTED')
        READ(11,*) NFRAMES
        DO I=1,NFRAMES
          READ(11,101) FILEFITS
          INQUIRE(FILE=FILEFITS,EXIST=LOGFILE)
          IF(.NOT.LOGFILE)THEN
            WRITE(*,100) 'FATAL ERROR: the file "'
            WRITE(*,100) FILEFITS(1:TRUELEN(FILEFITS))
            WRITE(*,101) '" does not exist.'
            STOP
          END IF
          CALL FTOPEN(IUNIT,FILEFITS,IREADWRITE,BLOCKSIZE,ISTATUS)
          CALL FTGKYJ(IUNIT,'NAXIS',NAXIS_(0),COMMENT,ISTATUS)
          CALL FTGKNJ(IUNIT,'NAXIS',1,2,NAXIS_(1),NFOUND,ISTATUS)
          CALL FTCLOS(IUNIT,ISTATUS)
          IF(ISTATUS.GT.0)THEN
            CALL PRINTERROR(ISTATUS)
          END IF
          IF(NAXIS_(1).GT.NSIZE) NSIZE=NAXIS_(1)
          IF(NAXIS_(2).GT.NSIZE) NSIZE=NAXIS_(2)
          NAXISFRAME(1,I,NCBUFF)=NAXIS_(1)
          NAXISFRAME(2,I,NCBUFF)=NAXIS_(2)
        END DO
        CLOSE(11)
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
