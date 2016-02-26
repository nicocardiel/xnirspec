C
C******************************************************************************
C Subrutina para escribir una image FITS
        SUBROUTINE SESCRFITS(FITSFILE,NCBUFF)
        IMPLICIT NONE
        CHARACTER*(*) FITSFILE
        INTEGER NCBUFF
C
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
C
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
C
        INTEGER TRUELEN
        INTEGER READI
        CHARACTER*255 READC
C
        INTEGER I,J,K,L
        INTEGER BITPIX,NAXIS(2,NMAXBUFF)
        INTEGER NAXIS_(2,NMAXBUFF)
        INTEGER ISTATUS,IUNIT,IUNIT_
        INTEGER BLOCKSIZE,NULLVAL
        INTEGER NAXIS1,NAXIS2
        INTEGER N1,N2,NSKIP,NKEYS,NSPACE
        INTEGER DI(NBOXMAX),DJ(NBOXMAX)
        INTEGER READWRITE
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL IMAGEN_(NXYMAX,NXYMAX)
        CHARACTER*1 CBOX9
        CHARACTER*2 CNIMAGE
        CHARACTER*80 RECORD
        CHARACTER*255 INFILE,CLISTHEAD,INFILE_(NBOXMAX)
        LOGICAL LOGFILE
        LOGICAL SIMPLE,EXTEND
C
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKIMAGEN1_/IMAGEN_
        COMMON/BLKNAXIS/NAXIS
C------------------------------------------------------------------------------
C Note: the pattern of the frames in box-9 is the following:
C       6 9 4
C       3 1 7
C       8 5 2
C offsets of each frame in the 768x768 composite mosaic
        DATA (DI(K),DJ(K),K=1,NBOXMAX) /
     +   256,256,  !1
     +   000,512,  !2
     +   256,000,  !3
     +   512,512,  !4
     +   000,256,  !5
     +   512,000,  !6
     +   256,512,  !7
     +   000,000,  !8
     +   512,256/  !9
C inicializamos variables
        ISTATUS=0               !controla posibles errores durante la ejecucion
        NULLVAL=-999
        BLOCKSIZE=1           !normalmente usaremos un valor de 1 (=2880 bytes)
        NSKIP=0 !evita WARNING de compilacion
C------------------------------------------------------------------------------
C pedimos nombre de la imagen
        L=TRUELEN(FITSFILE)
        INQUIRE(FILE=FITSFILE(1:L),EXIST=LOGFILE)
        IF(LOGFILE)THEN
          WRITE(*,101) '***FATAL ERROR***'
          WRITE(*,101) '=> The file '
          WRITE(*,100) FITSFILE(1:L)
          WRITE(*,101) ' already exists.'
          STOP
        END IF
        IF((NAXIS(1,NCBUFF).EQ.768).AND.(NAXIS(2,NCBUFF).EQ.768))THEN
          WRITE(*,*)
          CBOX9(1:1)=READC('Save 9 256x256 images (y/n)','y',
     +     'yn')
          IF(CBOX9.EQ.'y')THEN
            LOGFILE=.FALSE.
            DO WHILE(.NOT.LOGFILE)
              CLISTHEAD=READC('File list to copy FITS headers','@','@')
              INQUIRE(FILE=CLISTHEAD,EXIST=LOGFILE)
              IF(.NOT.LOGFILE)THEN
                WRITE(*,101) 'ERROR: this file does not exist. '//
     +           'Try again.'
                WRITE(*,100) 'Pres <CR> to continue...'
                READ(*,*)
              END IF
            END DO
            OPEN(35,FILE=CLISTHEAD,STATUS='OLD',FORM='FORMATTED')
            READ(35,*) !skip first line with number of files
            DO K=1,NBOXMAX
              READ(35,101) INFILE_(K)
              WRITE(*,100) 'Header will be copied from file: '
              WRITE(*,101) INFILE_(K)(1:TRUELEN(INFILE_(K)))
            END DO
            CLOSE(35)
          END IF
        ELSE
          CBOX9='n'
        END IF
        IF(CBOX9.EQ.'y')THEN
          N1=1
          N2=9
          WRITE(*,*)
          WRITE(*,101)'NOTE: the first keywords of the header templates'
          WRITE(*,101) 'will be skipped since they have already been'
          WRITE(*,101) 'written by this program. In particular:'
          WRITE(*,101) 'SIMPLE  =                    T'
          WRITE(*,101) 'BITPIX  =                  -32'
          WRITE(*,101) 'NAXIS   =                    2'
          WRITE(*,101) 'NAXIS1  =                  256'
          WRITE(*,101) 'NAXIS2  =                  256'
          WRITE(*,101) ' '
          WRITE(*,101) 'So, you must know the actual number of'
          WRITE(*,101) 'keywords including the NAXIS2 description.'
          NSKIP=READI('No. of initial keywords to skip','5')
        ELSE
          N1=1
          N2=1
        END IF
C localizamos un numero de unidad de fichero no utilizada
ccc        CALL FTGIOU(IUNIT,ISTATUS)
        IUNIT=82 !ojo, IUNIT=99 entra en conflicto con el Postcript y es
                 !precisamente este numero el que toma esta funcion
        IUNIT_=83
        DO K=N1,N2
C creamos un nuevo fichero FITS vacio
          IF(CBOX9.EQ.'y')THEN
            WRITE(CNIMAGE,'(I2.2)') K
            IF(FITSFILE(L-4:L).EQ.'.fits')THEN
              INFILE=FITSFILE(1:L-5)//CNIMAGE//'.fits'
            ELSE
              INFILE=FITSFILE(1:L)//CNIMAGE//'.fits'
            END IF
            CALL FTINIT(IUNIT,INFILE,BLOCKSIZE,ISTATUS)
            NAXIS_(1,NCBUFF)=256
            NAXIS_(2,NCBUFF)=256
          ELSE
            DI(1)=0
            DJ(1)=0
            NAXIS_(1,NCBUFF)=NAXIS(1,NCBUFF)
            NAXIS_(2,NCBUFF)=NAXIS(2,NCBUFF)
            CALL FTINIT(IUNIT,FITSFILE,BLOCKSIZE,ISTATUS)
          END IF
C inicializamos los parametros
          SIMPLE=.TRUE.
          BITPIX=-32
          EXTEND=.FALSE.
C escribimos los keywords indispensables
          CALL FTPHPR(IUNIT,SIMPLE,BITPIX,2,NAXIS_(1,NCBUFF),0,1,
     +     EXTEND,ISTATUS)
C si estamos en Box-9, copiamos los header de las imagenes de referencia
          IF(CBOX9.EQ.'y')THEN
            READWRITE=0 !solo lectura
            CALL FTOPEN(IUNIT_,INFILE_(K),READWRITE,BLOCKSIZE,ISTATUS)
            CALL FTGHSP(IUNIT_,NKEYS,NSPACE,ISTATUS)
            DO I=1,NSKIP
              CALL FTGREC(IUNIT_,I,RECORD,ISTATUS)
            END DO
            DO I=NSKIP+1,NKEYS
              CALL FTGREC(IUNIT_,I,RECORD,ISTATUS)
              CALL FTPREC(IUNIT,RECORD,ISTATUS) 
            END DO
            CALL FTCLOS(IUNIT_,ISTATUS)
          END IF
C salvamos la imagen
          DO I=1,NAXIS_(2,NCBUFF)
            DO J=1,NAXIS_(1,NCBUFF)
              IMAGEN_(J,I)=IMAGEN(J+DJ(K),I+DI(K),NCBUFF)
            END DO
          END DO
          NAXIS1=NAXIS_(1,NCBUFF)
          NAXIS2=NAXIS_(2,NCBUFF)
          IF(BITPIX.EQ.-32)THEN
            CALL FTP2DE(IUNIT,1,NXYMAX,NAXIS1,NAXIS2,IMAGEN_,ISTATUS)
          ELSE
            WRITE(*,101) '***FATAL ERROR***'
            WRITE(*,100) '=> subroutine ESCFITS: BITPIX ='
            WRITE(*,*) BITPIX
            CALL FTCLOS(IUNIT,ISTATUS)
            STOP
          END IF
C cerramos el fichero
          CALL FTCLOS(IUNIT,ISTATUS)
C liberamos el numero de unidad del fichero utilizado
ccc        CALL FTFIOU(IUNIT,ISTATUS)
C chequeamos si se ha producido algun error
          IF(ISTATUS.GT.0)THEN
            CALL PRINTERROR(ISTATUS)
          END IF
        END DO
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
