!
!******************************************************************************
! Subrutina para escribir una image FITS
!
! Write buffer NCBUFF in the FITS file FITSFILE (BITPIX=-32), optionally
! copying the header of the FITS file HDRFILE ('none': no template). A
! 768x768 image can be saved as the 9 individual 256x256 frames of a box-9
! mosaic. If the buffer has wavelength calibration (LWAVECAL), CTYPE1,
! CUNIT1, CRPIX1, CRVAL1 and CDELT1 are written (see WRITEWAVECAL; not for
! the box-9 frames); otherwise, a wavelength calibration copied from the
! header template is removed (see DELWAVECAL).
        SUBROUTINE SESCRFITS(FITSFILE,NCBUFF,HDRFILE)
        USE Dynamic_Array_IMAGEN
        USE Dynamic_Array_IMAGEN_
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
        INCLUDE 'interface_imagen_.inc'
! subroutine arguments
        CHARACTER*(*) FITSFILE
        INTEGER NCBUFF
        CHARACTER*(*) HDRFILE
!
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
!
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
!
        INTEGER TRUELEN
        INTEGER READI
        CHARACTER*255 READC
        CHARACTER*255 READC_FILE
!
        INTEGER I,J,K,L
        INTEGER BITPIX,NAXIS(2,NMAXBUFF)
        INTEGER NAXIS_(2,NMAXBUFF)
        INTEGER ISTATUS,IUNIT,IUNIT_
        INTEGER BLOCKSIZE,NULLVAL
        INTEGER NAXIS1,NAXIS2
        INTEGER N1,N2,NSKIP,NKEYS,NSPACE
        INTEGER DI(NBOXMAX),DJ(NBOXMAX)
        INTEGER DI0,DJ0           !offsets of the region of IMAGEN to be saved
        INTEGER READWRITE
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
!delete REAL IMAGEN_(NXYMAX,NXYMAX)
        CHARACTER*1 CBOX9
        CHARACTER*2 CNIMAGE
        CHARACTER*255 INFILE,CLISTHEAD,INFILE_(NBOXMAX)
        CHARACTER*255 C255
        LOGICAL LOGFILE
        LOGICAL SIMPLE,EXTEND
        LOGICAL LDELETED
        REAL CRPIX1(NMAXBUFF),CRVAL1(NMAXBUFF),CDELT1(NMAXBUFF)
        LOGICAL LWAVECAL(NMAXBUFF)
        CHARACTER*20 CTYPE1(NMAXBUFF),CUNIT1(NMAXBUFF)
!
!delete COMMON/BLKIMAGEN1/IMAGEN
!delete COMMON/BLKIMAGEN1_/IMAGEN_
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKWAVECAL1/CRPIX1,CRVAL1,CDELT1         !wavelength calibration
        COMMON/BLKWAVECAL2/LWAVECAL
        COMMON/BLKWAVECAL3/CTYPE1,CUNIT1
!------------------------------------------------------------------------------
! Note: the pattern of the frames in box-9 is the following:
!       6 9 4
!       3 1 7
!       8 5 2
! offsets of each frame in the 768x768 composite mosaic
        DATA (DI(K),DJ(K),K=1,NBOXMAX) / &
         256,256, & !1
         000,512, & !2
         256,000, & !3
         512,512, & !4
         000,256, & !5
         512,000, & !6
         256,512, & !7
         000,000, & !8
         512,256/  !9
! inicializamos variables
        ISTATUS=0               !controla posibles errores durante la ejecucion
        NULLVAL=-999
        BLOCKSIZE=1           !normalmente usaremos un valor de 1 (=2880 bytes)
        NSKIP=0 !evita WARNING de compilacion
!------------------------------------------------------------------------------
! pedimos nombre de la imagen
        L=TRUELEN(FITSFILE)
        INQUIRE(FILE=FITSFILE(1:L),EXIST=LOGFILE)
        IF(LOGFILE)THEN
          WRITE(*,101) '***FATAL ERROR***'
          WRITE(*,101) '=> The file '
          WRITE(*,100) FITSFILE(1:L)
          WRITE(*,101) ' already exists.'
          INCLUDE 'deallocate_arrays.inc'
          STOP
        END IF
        IF((NAXIS(1,NCBUFF).EQ.768).AND.(NAXIS(2,NCBUFF).EQ.768))THEN
          WRITE(*,*)
          C255=READC('Save 9 256x256 images (y/n)','y','yn')
          CBOX9=C255(1:1)
          IF(CBOX9.EQ.'y')THEN
            LOGFILE=.FALSE.
            DO WHILE(.NOT.LOGFILE)
              CLISTHEAD=READC_FILE('File list to copy FITS headers','@','@')
              INQUIRE(FILE=CLISTHEAD,EXIST=LOGFILE)
              IF(.NOT.LOGFILE)THEN
                WRITE(*,101) 'ERROR: this file does not exist. Try again.'
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
! localizamos un numero de unidad de fichero no utilizada
!!!        CALL FTGIOU(IUNIT,ISTATUS)
        IUNIT=82 !ojo, IUNIT=99 entra en conflicto con el Postcript y es
                 !precisamente este numero el que toma esta funcion
        IUNIT_=83
        DO K=N1,N2
! creamos un nuevo fichero FITS vacio
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
            DI0=DI(K)
            DJ0=DJ(K)
          ELSE
! do not modify DI(1), DJ(1): they are initialized with DATA (and are
! therefore SAVEd), so changing them would corrupt later box-9 saves
            DI0=0
            DJ0=0
            NAXIS_(1,NCBUFF)=NAXIS(1,NCBUFF)
            NAXIS_(2,NCBUFF)=NAXIS(2,NCBUFF)
            CALL FTINIT(IUNIT,FITSFILE,BLOCKSIZE,ISTATUS)
          END IF
! inicializamos los parametros
          SIMPLE=.TRUE.
          BITPIX=-32
          EXTEND=.FALSE.
! escribimos los keywords indispensables
          CALL FTPHPR(IUNIT,SIMPLE,BITPIX,2,NAXIS_(1,NCBUFF),0,1,EXTEND,ISTATUS)
! si estamos en Box-9, copiamos los header de las imagenes de referencia
          IF(CBOX9.EQ.'y')THEN
            READWRITE=0 !solo lectura
            CALL FTOPEN(IUNIT_,INFILE_(K),READWRITE,BLOCKSIZE,ISTATUS)
            CALL FTGHSP(IUNIT_,NKEYS,NSPACE,ISTATUS)
            CALL SESCRFITS_COPYHEADER(IUNIT_,IUNIT,NSKIP,NKEYS,ISTATUS)
            CALL FTCLOS(IUNIT_,ISTATUS)
          ELSE
! si HDRFILE no es 'none' salvamos el header de una imagen externa
            IF(HDRFILE.NE.'none')THEN
              READWRITE=0 !solo lectura
              CALL FTOPEN(IUNIT_,HDRFILE,READWRITE,BLOCKSIZE,ISTATUS)
              CALL FTGHSP(IUNIT_,NKEYS,NSPACE,ISTATUS)
              WRITE(*,100) 'KEYS='
              WRITE(*,*) NKEYS
              WRITE(*,101) 'NOTE: the first keywords of the header'
              WRITE(*,101) 'template file will be skipped since they'
              WRITE(*,101) 'have already been written by this program.'
              WRITE(*,101) 'In particular:'
              WRITE(*,101) 'SIMPLE  =                    T'
              WRITE(*,101) 'BITPIX  =                  -32'
              WRITE(*,101) 'NAXIS   =                    2'
              WRITE(*,101) 'NAXIS1  =                 ????'
              WRITE(*,101) 'NAXIS2  =                 ????'
              WRITE(*,101) ' '
              WRITE(*,101) 'So, you must know the actual number of'
              WRITE(*,101) 'keywords including the NAXIS2 description.'
              NSKIP=READI('No. of initial keywords to be skipped','5')
              CALL SESCRFITS_COPYHEADER(IUNIT_,IUNIT,NSKIP,NKEYS,ISTATUS)
              CALL FTCLOS(IUNIT_,ISTATUS)
            END IF
          END IF
! wavelength calibration of the buffer (not for the individual frames of a
! box-9 mosaic); keywords copied from a header template are updated. If the
! buffer is not calibrated, a wavelength calibration copied from the
! template is removed (it does not correspond to the saved image)
          IF(CBOX9.NE.'y')THEN
            IF(LWAVECAL(NCBUFF))THEN
              CALL WRITEWAVECAL(IUNIT,CRPIX1(NCBUFF),CRVAL1(NCBUFF), &
               CDELT1(NCBUFF),CTYPE1(NCBUFF),CUNIT1(NCBUFF),ISTATUS)
            ELSE
              CALL DELWAVECAL(IUNIT,LDELETED,ISTATUS)
              IF(LDELETED)THEN
                WRITE(*,101) 'WARNING: the wavelength calibration '// &
                 'copied from the header template has been removed'
                WRITE(*,101) '         (the image in the buffer is '// &
                 'not wavelength calibrated)'
              END IF
            END IF
          END IF
! salvamos la imagen
          DO I=1,NAXIS_(2,NCBUFF)
            DO J=1,NAXIS_(1,NCBUFF)
              IMAGEN_(J,I)=IMAGEN(J+DJ0,I+DI0,NCBUFF)
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
            INCLUDE 'deallocate_arrays.inc'
            STOP
          END IF
! cerramos el fichero
          CALL FTCLOS(IUNIT,ISTATUS)
! liberamos el numero de unidad del fichero utilizado
!!!        CALL FTFIOU(IUNIT,ISTATUS)
! chequeamos si se ha producido algun error
          IF(ISTATUS.GT.0)THEN
            CALL PRINTERROR(ISTATUS)
          END IF
        END DO
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
!
!******************************************************************************
! SUBROUTINE SESCRFITS_COPYHEADER(IUNIT_,IUNIT,NSKIP,NKEYS,ISTATUS)
!
! Copy the header records NSKIP+1,...,NKEYS of the template file opened in
! unit IUNIT_ to the header of the file being written in unit IUNIT. Records
! identical to any of those already written in IUNIT (e.g. the two COMMENT
! lines with the reference to the FITS standard that CFITSIO writes after
! NAXIS2) are not copied, to avoid duplicating them.
!
! Parameters
! ----------
! IUNIT_ : INTEGER
!     Unit of the template file (opened read-only).
! IUNIT : INTEGER
!     Unit of the output file.
! NSKIP : INTEGER
!     Number of initial records of the template that are skipped.
! NKEYS : INTEGER
!     Number of records in the header of the template.
! ISTATUS : INTEGER
!     CFITSIO status (input and output).
!------------------------------------------------------------------------------
        SUBROUTINE SESCRFITS_COPYHEADER(IUNIT_,IUNIT,NSKIP,NKEYS,ISTATUS)
        IMPLICIT NONE
        INTEGER IUNIT_,IUNIT,NSKIP,NKEYS,ISTATUS
!
        INTEGER NOWNMAX
        PARAMETER (NOWNMAX=100)
!
        INTEGER I,K,NOWN,NSPACE
        CHARACTER*80 RECORD
        CHARACTER*80 OWNREC(NOWNMAX)
        LOGICAL LDUPLICATED
!------------------------------------------------------------------------------
! records already written in the output header
        CALL FTGHSP(IUNIT,NOWN,NSPACE,ISTATUS)
        NOWN=MIN(NOWN,NOWNMAX)
        DO K=1,NOWN
          CALL FTGREC(IUNIT,K,OWNREC(K),ISTATUS)
        END DO
! copy the remaining records of the template, skipping duplicates
        DO I=NSKIP+1,NKEYS
          CALL FTGREC(IUNIT_,I,RECORD,ISTATUS)
          LDUPLICATED=.FALSE.
          DO K=1,NOWN
            IF(RECORD.EQ.OWNREC(K)) LDUPLICATED=.TRUE.
          END DO
          IF(.NOT.LDUPLICATED) CALL FTPREC(IUNIT,RECORD,ISTATUS)
        END DO
        END
