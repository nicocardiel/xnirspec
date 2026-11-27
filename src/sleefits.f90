!
!******************************************************************************
! Subrutina para leer imagenes FITS.
! FITSFILE: nombre del fichero a leer
! LSHOW: if .TRUE. mostramos las keywords
! IROTATE: rotation angle (only 0, +90, +180, +270, -90, -180, -270)
! LBOX9: if .TRUE. la subrutina devuelve la imagen en IMAGEN_ y no en la
!        variable IMAGEN(J,I,NEWBUFF) -> los valores de IROTATE y de NEWBUFF no 
!        tienen ninguna importancia en ese caso
! IBUFF: if <> 0 means that the photometry keywords of the HST are searched
!        and, if found, stored in the corresponding array (buffer IBUFF).
!
! Tile-compressed images are supported (image parameters are obtained
! with FTGIDT, FTGIDM and FTGISZ).
!
! Files with more than one HDU: the user is asked for the extension number,
! counted from 0 (primary HDU) as in astropy. The default is the first image
! HDU with 1 <= NAXIS <= 3; the answer -1 prints a list of the HDUs (see
! FITSINFO) and asks again, and HDUs that are not such images are rejected.
!
! Wavelength calibration of FITS files (CRPIX1, CRVAL1, CDELT1, with
! CTYPE1 and CUNIT1) is read with FITSWAVECAL, and it is only accepted when
! CTYPE1 is a linear wavelength axis and CDELT1 (or CD1_1) is not zero;
! REDUCEME files are calibrated when DISP is not zero. When the
! image is rotated: with IROTATE=+-180 it is transformed to the reversed X
! axis (CRPIX1 -> NAXIS1-CRPIX1+1, CDELT1 -> -CDELT1); with IROTATE=+-90 or
! +-270 the original NAXIS1 becomes the Y axis, so the calibration is not
! applicable and LWAVECAL is set to .FALSE. The calibration of files in
! REDUCEME format (STWV, DISP) is not modified.

!------------------------------------------------------------------------------
        SUBROUTINE SLEEFITS(FITSFILE,LSHOW,IROTATE,NEWBUFF,LBOX9,IBUFF)
        USE Dynamic_Array_IMAGEN
        USE Dynamic_Array_IMAGEN_
        USE Dynamic_Array_LNULL
        USE Dynamic_Array_LNULL_
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
        INCLUDE 'interface_imagen_.inc'
        INCLUDE 'interface_lnull.inc'
        INCLUDE 'interface_lnull_.inc'
! subroutine arguments
        CHARACTER*(*) FITSFILE
        LOGICAL LSHOW
        INTEGER IROTATE
        INTEGER NEWBUFF
        LOGICAL LBOX9
        INTEGER IBUFF
! parametros
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
! funciones auxiliares
        INTEGER READILIM
        INTEGER TRUELEN
! variables globales (COMMONs)
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
!delete REAL IMAGEN_(NXYMAX,NXYMAX)
        REAL CRPIX1(NMAXBUFF),CRVAL1(NMAXBUFF),CDELT1(NMAXBUFF)
        LOGICAL LWAVECAL(NMAXBUFF)
        CHARACTER*20 CTYPE1(NMAXBUFF),CUNIT1(NMAXBUFF)
        REAL STWV,DISP
        REAL AIRMASS,TIMEXPOS
!delete LOGICAL LNULL(NXMAX,NYMAX,NMAXBUFF)
        LOGICAL ANYNULL
        INTEGER NAXIS(2,NMAXBUFF)
! variables locales
        INTEGER NEW_HDU,HDUTYPE
        INTEGER NHDU,IEXT,IEXTDEF,NDIMHDU
        CHARACTER*10 CDEF
        INTEGER JROW(NXYMAX)
        INTEGER I,J,K,L
        INTEGER FIRSTPIX
        INTEGER BITPIX
        INTEGER ISTATUS,IREADWRITE,IUNIT
        INTEGER BLOCKSIZE,NULLVAL
        INTEGER NKEYS,NSPACE
        INTEGER NAXES
        INTEGER NAXIS_(3)                                  !OJO: el limite es 3
        INTEGER NAXIS3FIXED1,NAXIS3FIXED2
        INTEGER NANYNULLS
        REAL FROW(NXYMAX)
        CHARACTER*20 EXTNAME
        CHARACTER*50 COMMENT
        CHARACTER*80 CLINEA
        REAL PHOTFLAM(NMAXBUFF),PHOTZPT(NMAXBUFF),EXPTIME(NMAXBUFF)
        LOGICAL EXTEND
        LOGICAL ZIMAGE
!delete LOGICAL LNULL_(NXYMAX,NXYMAX)
        LOGICAL LOGFILE,LANYNULL
        LOGICAL LROW(NXYMAX)
        LOGICAL L_PHOTFLAM(NMAXBUFF),L_PHOTZPT(NMAXBUFF)
        LOGICAL L_EXPTIME(NMAXBUFF)
        LOGICAL LERROR
        LOGICAL LREDUCEME                  !.TRUE. if file has REDUCEME format
!
!delete COMMON/BLKIMAGEN1/IMAGEN             !imagen FITS leida en formato REAL
!delete COMMON/BLKIMAGEN1_/IMAGEN_              !es global para ahorrar memoria
        COMMON/BLKWAVECAL1/CRPIX1,CRVAL1,CDELT1         !wavelength calibration
        COMMON/BLKWAVECAL2/LWAVECAL                     !wavelength calibration
        COMMON/BLKWAVECAL3/CTYPE1,CUNIT1                !wavelength calibration
        COMMON/BLKREDUCEME/STWV,DISP,AIRMASS,TIMEXPOS     !Reduceme header info
!delete COMMON/BLKLNULL/LNULL,ANYNULL   !mascara que indica si existen NaN, etc
        COMMON/BLKANYNULL/ANYNULL                   !indica si existen NaN, etc
        COMMON/BLKNAXIS/NAXIS                                      !dimensiones
        COMMON/BLK_HSTFLUX1/L_PHOTFLAM,L_PHOTZPT,L_EXPTIME !keywords with info
        COMMON/BLK_HSTFLUX2/PHOTFLAM,PHOTZPT,EXPTIME !concerning HST flux cal.
!------------------------------------------------------------------------------
        CALL Initialize_Dynamic_Array_LNULL_
!------------------------------------------------------------------------------
! inicializamos variables
        ISTATUS=0               !controla posibles errores durante la ejecucion
        IREADWRITE=0                      !la imagen se abrira en modo READONLY
        NULLVAL=-999
        LANYNULL=.FALSE.
        NANYNULLS=0
        LREDUCEME=.FALSE.
!------------------------------------------------------------------------------
! chequeamos si existe la imagen
        INQUIRE(FILE=FITSFILE,EXIST=LOGFILE)
        IF(LOGFILE)THEN
        ELSE
          WRITE(*,100) 'FATAL ERROR: the file "'
          WRITE(*,100) FITSFILE(1:TRUELEN(FITSFILE))
          WRITE(*,101) '" does not exist.'
          INCLUDE 'deallocate_arrays.inc'
          CALL Deallocate_Array_LNULL_
          STOP
        END IF
! localizamos un numero de unidad de fichero no utilizada
!!!        CALL FTGIOU(IUNIT,ISTATUS)
        IUNIT=80 !ojo, IUNIT=99 entra en conflicto con el Postcript y es
                 !precisamente este numero el que toma esta funcion
! abrimos el fichero
        CALL FTOPEN(IUNIT,FITSFILE,IREADWRITE,BLOCKSIZE,ISTATUS)
! si la imagen no es FITS, miramos a ver si es de formato REDUCEME
        IF(ISTATUS.NE.0)THEN
          CALL FTCLOS(IUNIT,ISTATUS)
          CALL TRY_REDUCEME(FITSFILE,ISTATUS,NAXIS_)
          IF(ISTATUS.EQ.0)THEN
            DO I=1,NAXIS_(2)
              DO J=1,NAXIS_(1)
                LNULL_(J,I)=.FALSE.
              END DO
            END DO
! REDUCEME files have no CTYPE1: they are considered wavelength calibrated
! when DISP is not zero (CTYPE1='WAVE', units unknown)
            LWAVECAL(NEWBUFF)=(DISP.NE.0.0)
            CRPIX1(NEWBUFF)=1.0
            CRVAL1(NEWBUFF)=STWV
            CDELT1(NEWBUFF)=DISP
            CTYPE1(NEWBUFF)='WAVE'
            CUNIT1(NEWBUFF)=' '
            LREDUCEME=.TRUE.
            WRITE(*,100) 'WARNING: reading file with REDUCEME format'
            GOTO 777
          ELSE
            WRITE(*,100) 'ERROR: REDUCEME format not found'
            CALL Deallocate_Array_LNULL_
            RETURN
          END IF
        END IF
! miramos si la imagen tiene extensiones y, en caso afirmativo, preguntamos
! cual queremos leer
        CALL FTGKYL(IUNIT,'EXTEND',EXTEND,COMMENT,ISTATUS)
        IF(ISTATUS.EQ.202)THEN
          EXTEND=.FALSE.
          ISTATUS=0
        END IF
! The user numbers the HDUs from 0 (primary HDU), as in astropy, whereas
! CFITSIO (FTMAHD) numbers them from 1: NEW_HDU is always the CFITSIO
! number (NEW_HDU=IEXT+1). EXTEND=T does not guarantee that the file has
! extensions, so the question is only asked when there is more than one HDU.
        NHDU=1
        IF(EXTEND) CALL FTTHDU(IUNIT,NHDU,ISTATUS)
        IF(NHDU.GT.1)THEN
! default: first image HDU that can be read (1 <= NAXIS <= 3)
          IEXTDEF=-1
          DO I=1,NHDU
            CALL FTMAHD(IUNIT,I,HDUTYPE,ISTATUS)
            IF(HDUTYPE.EQ.0)THEN
              CALL FTGIDM(IUNIT,NDIMHDU,ISTATUS)
              IF((NDIMHDU.GE.1).AND.(NDIMHDU.LE.3))THEN
                IEXTDEF=I-1
                EXIT
              END IF
            END IF
          END DO
          IF((IEXTDEF.EQ.-1).OR.(ISTATUS.GT.0))THEN
            IF(ISTATUS.GT.0) CALL PRINTERROR(ISTATUS)
            CALL FITSINFO(IUNIT,FITSFILE,NHDU)
            WRITE(*,101) '***FATAL ERROR***'
            WRITE(*,101) '=> no image HDU with 1 <= NAXIS <= 3 found'
            ISTATUS=0
            CALL FTCLOS(IUNIT,ISTATUS)
            INCLUDE 'deallocate_arrays.inc'
            CALL Deallocate_Array_LNULL_
            STOP
          END IF
          WRITE(*,101) '***WARNING***'
          WRITE(*,101) '=> this file contains extensions'
          WRITE(CDEF,'(I10)') IEXTDEF
10        IEXT=READILIM('Extension number to be read (0=primary, -1=list)', &
           CDEF,-1,NHDU-1)
          IF(IEXT.EQ.-1)THEN
            CALL FITSINFO(IUNIT,FITSFILE,NHDU)
            GOTO 10
          END IF
          NEW_HDU=IEXT+1
          CALL FTMAHD(IUNIT,NEW_HDU,HDUTYPE,ISTATUS)
          NDIMHDU=0
          IF(HDUTYPE.EQ.0) CALL FTGIDM(IUNIT,NDIMHDU,ISTATUS)
          IF((HDUTYPE.NE.0).OR.(NDIMHDU.LT.1).OR.(NDIMHDU.GT.3))THEN
            WRITE(*,101) 'ERROR: this HDU is not an image with 1 <= NAXIS <= 3'// &
             ' (use -1 to list the extensions)'
            GOTO 10
          END IF
        ELSE
          NEW_HDU=1
        END IF
! determinamos el numero de keywords en la cabecera y se muestran
        CALL FTGHSP(IUNIT,NKEYS,NSPACE,ISTATUS)
        IF(LSHOW)THEN
          DO I=1,NKEYS
            CALL FTGREC(IUNIT,I,CLINEA,ISTATUS)
            L=TRUELEN(CLINEA)
            WRITE(*,101) CLINEA(1:L)
          END DO
          IF(ISTATUS.EQ.0)THEN                                !todo ha ido bien
            WRITE(*,101) 'END'
            WRITE(*,*)
          END IF
        END IF
! leemos BITPIX. FTGIDT, FTGIDM and FTGISZ are used instead of reading the
! keywords BITPIX, NAXIS and NAXISn: for tile-compressed images these
! keywords describe the binary table that stores the compressed data, and
! CFITSIO returns the values of the uncompressed image (ZBITPIX, ZNAXIS...)
        CALL FTGIDT(IUNIT,BITPIX,ISTATUS)
        IF(.NOT.LSHOW)THEN
          WRITE(*,*)
          IF(NEW_HDU.EQ.1)THEN
            WRITE(*,101) 'CFITSIO> Primary HDU'
          ELSE
            WRITE(*,100) 'CFITSIO> Extension number: '
            WRITE(*,*) NEW_HDU-1
            CALL FTGKYS(IUNIT,'EXTNAME',EXTNAME,COMMENT,ISTATUS)
            IF(ISTATUS.EQ.202)THEN              !EXTNAME is optional
              ISTATUS=0
              EXTNAME=' '
! default name given by astropy to tile-compressed images (ZIMAGE=T)
              CALL FTGKYL(IUNIT,'ZIMAGE',ZIMAGE,COMMENT,ISTATUS)
              IF(ISTATUS.EQ.202)THEN
                ISTATUS=0
                ZIMAGE=.FALSE.
              END IF
              IF(ZIMAGE) EXTNAME='COMPRESSED_IMAGE'
            END IF
            WRITE(*,100) 'CFITSIO> EXTNAME: '
            WRITE(*,101) EXTNAME
          ENDIF
          WRITE(*,100) 'CFITSIO> BITPIX : '
          WRITE(*,*) BITPIX
        END IF
! leemos NAXIS, NAXIS1, NAXIS2 y NAXIS3
        CALL FTGIDM(IUNIT,NAXES,ISTATUS)
        IF(.NOT.LSHOW)THEN
          WRITE(*,100) 'CFITSIO> NAXIS  : '
          WRITE(*,*) NAXES 
        END IF
        LERROR=.FALSE.
        IF(NAXES.EQ.1)THEN
          CALL FTGISZ(IUNIT,1,NAXIS_,ISTATUS)
          NAXIS_(2)=1
          NAXIS_(3)=1
          NAXIS3FIXED1=1
          NAXIS3FIXED2=1
        ELSEIF(NAXES.EQ.2)THEN
          CALL FTGISZ(IUNIT,2,NAXIS_,ISTATUS)
          NAXIS_(3)=1
          NAXIS3FIXED1=1
          NAXIS3FIXED2=1
        ELSEIF(NAXES.EQ.3)THEN
          CALL FTGISZ(IUNIT,3,NAXIS_,ISTATUS)
          IF(NAXIS_(3).EQ.1)THEN
            WRITE(*,101) 'WARNING: NAXIS=3 with NAXIS3=1. Reading anyway!'
            NAXIS3FIXED1=1
            NAXIS3FIXED2=1
          END IF
        ELSE
          LERROR=.TRUE.
        END IF
!
        IF(LERROR)THEN
          WRITE(*,101) '***FATAL ERROR***'
          WRITE(*,100) '=> Invalid NAXIS='
          WRITE(*,*) NAXES
          CALL FTCLOS(IUNIT,ISTATUS)
          INCLUDE 'deallocate_arrays.inc'
          CALL Deallocate_Array_LNULL_
          STOP
        END IF
!
        IF(.NOT.LSHOW)THEN
          WRITE(*,100) 'CFITSIO> NAXIS1 : '
          WRITE(*,*) NAXIS_(1)
          WRITE(*,100) 'CFITSIO> NAXIS2 : '
          WRITE(*,*) NAXIS_(2)
          IF(NAXES.EQ.3)THEN
            WRITE(*,100) 'CFITSIO> NAXIS3 : '
            WRITE(*,*) NAXIS_(3)
          END IF
        END IF
        IF(NAXIS_(3).GT.1)THEN
          WRITE(*,101) 'Select option to collapse data along NAXIS3:'
          WRITE(*,101) '        0 : collapse whole 3D data cube'
          WRITE(*,101) '       -1 : collapse slice range'
          WRITE(*,101) ' <number> : read single slice <number>'
          K=READILIM('Option','0',-1,NAXIS_(3))
          IF(K.EQ.0)THEN
            NAXIS3FIXED1=1
            NAXIS3FIXED2=NAXIS_(3)
          ELSEIF(K.EQ.-1)THEN
            NAXIS3FIXED1=READILIM('Fist slice number','@',1,NAXIS_(3))
            NAXIS3FIXED2=READILIM('Last slice number','@',NAXIS3FIXED1,NAXIS_(3))
          ELSE
            NAXIS3FIXED1=K
            NAXIS3FIXED2=K
          END IF
        END IF
!
        IF((IROTATE.EQ.0).OR.(IROTATE.EQ.180).OR.(IROTATE.EQ.-180))THEN
          IF(NAXIS_(1).GT.NXMAX)THEN
            WRITE(*,100) 'NAXIS(1), NXMAX: '
            WRITE(*,*) NAXIS_(1),NXMAX
            WRITE(*,101) '* FATAL ERROR in subroutine LEEFITS:'
            WRITE(*,101) 'NAXIS(1) > NXMAX'
            CALL FTCLOS(IUNIT,ISTATUS)
            INCLUDE 'deallocate_arrays.inc'
            CALL Deallocate_Array_LNULL_
            STOP
          END IF
        ELSE
          IF(NAXIS_(1).GT.NYMAX)THEN
            WRITE(*,100) 'NAXIS(1), NYMAX: '
            WRITE(*,*) NAXIS_(1),NYMAX
            WRITE(*,101) '* FATAL ERROR in subroutine LEEFITS:'
            WRITE(*,101) 'NAXIS(1) > NYMAX'
            CALL FTCLOS(IUNIT,ISTATUS)
            INCLUDE 'deallocate_arrays.inc'
            CALL Deallocate_Array_LNULL_
            STOP
          END IF
        END IF
        IF((IROTATE.EQ.0).OR.(IROTATE.EQ.180).OR.(IROTATE.EQ.-180))THEN
          IF(NAXIS_(2).GT.NYMAX)THEN
            WRITE(*,100) 'NAXIS(2), NYMAX: '
            WRITE(*,*) NAXIS_(2),NYMAX
            WRITE(*,101) '* FATAL ERROR in subroutine LEEFITS:'
            WRITE(*,101) 'NAXIS(2) > NYMAX'
            CALL FTCLOS(IUNIT,ISTATUS)
            INCLUDE 'deallocate_arrays.inc'
            CALL Deallocate_Array_LNULL_
            STOP
          END IF
        ELSE
          IF(NAXIS_(2).GT.NXMAX)THEN
            WRITE(*,100) 'NAXIS(2), NXMAX: '
            WRITE(*,*) NAXIS_(2),NXMAX
            WRITE(*,101) '* FATAL ERROR in subroutine LEEFITS:'
            WRITE(*,101) 'NAXIS(2) > NXMAX'
            CALL FTCLOS(IUNIT,ISTATUS)
            INCLUDE 'deallocate_arrays.inc'
            CALL Deallocate_Array_LNULL_
            STOP
          END IF
        END IF
! comprobamos si es una imagen HST con informacion sobre calibracion en
! flujo; de existir, esta informacion puede ser utilizada para calcular
! las magnitudes con SExtractor
        IF(IBUFF.NE.0)THEN
          CALL FTGKYE(IUNIT,'PHOTFLAM',PHOTFLAM(IBUFF),COMMENT,ISTATUS)
          IF(ISTATUS.GT.0)THEN
            L_PHOTFLAM(IBUFF)=.FALSE.
            PHOTFLAM(IBUFF)=0.0
            ISTATUS=0
          ELSE
            L_PHOTFLAM(IBUFF)=.TRUE.
          END IF
          CALL FTGKYE(IUNIT,'PHOTZPT',PHOTZPT(IBUFF),COMMENT,ISTATUS)
          IF(ISTATUS.GT.0)THEN
            L_PHOTZPT(IBUFF)=.FALSE.
            PHOTZPT(IBUFF)=0.0
            ISTATUS=0
          ELSE
            L_PHOTZPT(IBUFF)=.TRUE.
          END IF
          CALL FTGKYE(IUNIT,'EXPTIME',EXPTIME(IBUFF),COMMENT,ISTATUS)
          IF(ISTATUS.GT.0)THEN
            L_EXPTIME(IBUFF)=.FALSE.
            EXPTIME(IBUFF)=0.0
            ISTATUS=0
          ELSE
            L_EXPTIME(IBUFF)=.TRUE.
          END IF
        END IF
! wavelength calibration (CTYPE1 must be a linear wavelength axis; see
! FITSWAVECAL); with LBOX9 the image is not stored in a buffer, and NEWBUFF
! may not be a valid buffer number
        IF((NEWBUFF.GE.1).AND.(NEWBUFF.LE.NMAXBUFF))THEN
          CALL FITSWAVECAL(IUNIT,LWAVECAL(NEWBUFF),CRPIX1(NEWBUFF), &
           CRVAL1(NEWBUFF),CDELT1(NEWBUFF),CTYPE1(NEWBUFF),CUNIT1(NEWBUFF))
          IF(LWAVECAL(NEWBUFF))THEN
            WRITE(*,100) 'CFITSIO> CTYPE1 : '
            WRITE(*,101) TRIM(CTYPE1(NEWBUFF))
            IF(LEN_TRIM(CUNIT1(NEWBUFF)).GT.0)THEN
              WRITE(*,100) 'CFITSIO> CUNIT1 : '
              WRITE(*,101) TRIM(CUNIT1(NEWBUFF))
            END IF
            WRITE(*,100) 'CFITSIO> CRPIX1 : '
            WRITE(*,*) CRPIX1(NEWBUFF)
            WRITE(*,100) 'CFITSIO> CRVAL1 : '
            WRITE(*,*) CRVAL1(NEWBUFF)
            WRITE(*,100) 'CFITSIO> CDELT1 : '
            WRITE(*,*) CDELT1(NEWBUFF)
          END IF
        END IF
! leemos la imagen
        DO I=1,NAXIS_(2)
          DO J=1,NAXIS_(1)
            IMAGEN_(J,I)=0.0
            LNULL_(J,I)=.FALSE.
          END DO
        END DO
        IF((BITPIX.EQ.8).OR.(BITPIX.EQ.16))THEN
          DO K=NAXIS3FIXED1,NAXIS3FIXED2
            DO I=1,NAXIS_(2)
              FIRSTPIX=(K-1)*(NAXIS_(1)*NAXIS_(2))+(I-1)*NAXIS_(1)+1
              CALL FTGPFJ(IUNIT,1,FIRSTPIX,NAXIS_(1),JROW(1),LROW(1),ANYNULL,ISTATUS)
              IF(ANYNULL)THEN
                DO J=1,NAXIS_(1)
                  LNULL_(J,I)=LNULL_(J,I).OR.LROW(J)
                  IF(LNULL_(J,I))THEN
                    NANYNULLS=NANYNULLS+1
                    IF(NANYNULLS.LT.10)THEN
                      print*,j,i,k,jrow(j)
                    END IF
                    JROW(J)=0
                  END IF
                END DO
                LANYNULL=.TRUE.
              END IF
              DO J=1,NAXIS_(1)
                IMAGEN_(J,I)=IMAGEN_(J,I)+REAL(JROW(J))
              END DO
            END DO
          END DO
        ELSEIF((BITPIX.EQ.32).OR.(BITPIX.EQ.-32).OR.(BITPIX.EQ.-64))THEN
          DO K=NAXIS3FIXED1,NAXIS3FIXED2
            DO I=1,NAXIS_(2)
              FIRSTPIX=(K-1)*(NAXIS_(1)*NAXIS_(2))+(I-1)*NAXIS_(1)+1
              CALL FTGPFE(IUNIT,1,FIRSTPIX,NAXIS_(1),FROW(1),LROW(1),ANYNULL,ISTATUS)
              IF(ANYNULL)THEN
                DO J=1,NAXIS_(1)
                  LNULL_(J,I)=LNULL_(J,I).OR.LROW(J)
                  IF(LNULL_(J,I))THEN
                    NANYNULLS=NANYNULLS+1
                    IF(NANYNULLS.LT.10)THEN
                      print*,j,i,k,frow(j)
                    END IF
                    FROW(J)=0.0
                  END IF
                END DO
                LANYNULL=.TRUE.
              END IF
              DO J=1,NAXIS_(1)
                IMAGEN_(J,I)=IMAGEN_(J,I)+FROW(J)
              END DO
            END DO
          END DO
        ELSE
          WRITE(*,100) 'FATAL ERROR in subroutine LEEFITS: BITPIX ='
          WRITE(*,*) BITPIX
          CALL FTCLOS(IUNIT,ISTATUS)
          INCLUDE 'deallocate_arrays.inc'
          CALL Deallocate_Array_LNULL_
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
        ANYNULL=LANYNULL                  !basta que haya ocurrido una sola vez
        IF(LBOX9)THEN
          CALL Deallocate_Array_LNULL_
          RETURN
        END IF
!------------------------------------------------------------------------------
777     CONTINUE
! tratamos la orientacion de la imagen
        IF(IROTATE.EQ.0)THEN
          DO I=1,NAXIS_(2)
            DO J=1,NAXIS_(1)
              IMAGEN(J,I,NEWBUFF)=IMAGEN_(J,I)
              LNULL(J,I,NEWBUFF)=LNULL_(J,I)
            END DO
          END DO
          NAXIS(1,NEWBUFF)=NAXIS_(1)
          NAXIS(2,NEWBUFF)=NAXIS_(2)
        ELSEIF((IROTATE.EQ.90).OR.(IROTATE.EQ.-270))THEN
          DO I=1,NAXIS_(2)
            DO J=1,NAXIS_(1)
              IMAGEN(NAXIS_(2)-I+1,J,NEWBUFF)=IMAGEN_(J,I)
              LNULL(NAXIS_(2)-I+1,J,NEWBUFF)=LNULL_(J,I)
            END DO
          END DO
          NAXIS(1,NEWBUFF)=NAXIS_(2)
          NAXIS(2,NEWBUFF)=NAXIS_(1)
          CALL SLEEFITS_NOWAVECAL(LREDUCEME,LWAVECAL(NEWBUFF))
        ELSEIF((IROTATE.EQ.180).OR.(IROTATE.EQ.-180))THEN
          DO I=1,NAXIS_(2)
            DO J=1,NAXIS_(1)
              IMAGEN(J,I,NEWBUFF)=IMAGEN_(NAXIS_(1)-J+1,NAXIS_(2)-I+1)
              LNULL(J,I,NEWBUFF)=LNULL_(NAXIS_(1)-J+1,NAXIS_(2)-I+1)
            END DO
          END DO
          NAXIS(1,NEWBUFF)=NAXIS_(1)
          NAXIS(2,NEWBUFF)=NAXIS_(2)
! the X axis is reversed: pixel X' corresponds to the original pixel
! NAXIS1-X'+1, so the same wavelength is obtained with the new CRPIX1 and
! CDELT1 (CRVAL1 is unchanged)
          IF(LWAVECAL(NEWBUFF).AND.(.NOT.LREDUCEME))THEN
            CRPIX1(NEWBUFF)=REAL(NAXIS_(1))-CRPIX1(NEWBUFF)+1.0
            CDELT1(NEWBUFF)=-CDELT1(NEWBUFF)
          END IF
        ELSEIF((IROTATE.EQ.270).OR.(IROTATE.EQ.-90))THEN
          DO I=1,NAXIS_(2)
            DO J=1,NAXIS_(1)
              IMAGEN(I,NAXIS_(1)-J+1,NEWBUFF)=IMAGEN_(J,I)
              LNULL(I,NAXIS_(1)-J+1,NEWBUFF)=LNULL_(J,I)
            END DO
          END DO
          NAXIS(1,NEWBUFF)=NAXIS_(2)
          NAXIS(2,NEWBUFF)=NAXIS_(1)
          CALL SLEEFITS_NOWAVECAL(LREDUCEME,LWAVECAL(NEWBUFF))
        ELSE
          WRITE(*,101) '***FATAL ERROR***'
          WRITE(*,100) '=> IROTATE='
          WRITE(*,*) IROTATE
          WRITE(*,101) '=> Invalid IROTATE in subroutine SLEEFITS.'
          INCLUDE 'deallocate_arrays.inc'
          CALL Deallocate_Array_LNULL_
          STOP
        END IF
        IF(NANYNULLS.GT.0)THEN
          WRITE(*,100) 'Total number of ANYNULL pixels: '
          WRITE(*,*) NANYNULLS
          IF(NAXES.EQ.3)THEN
            WRITE(*,101) 'WARNING: null values found!'
            WRITE(*,101) 'WARNING: spaxel values set to zero for all the slices'
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END IF
        CALL Deallocate_Array_LNULL_
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
!
!******************************************************************************
! SUBROUTINE SLEEFITS_NOWAVECAL(LREDUCEME,LWAVECAL)
!
! Invalidate the wavelength calibration of a FITS image read with
! IROTATE=+-90 or +-270: the original NAXIS1 (to which CRPIX1, CRVAL1 and
! CDELT1 refer) becomes the Y axis. Files in REDUCEME format are not
! modified.
!
! Parameters
! ----------
! LREDUCEME : LOGICAL
!     .TRUE. if the file has REDUCEME format (nothing is done).
! LWAVECAL : LOGICAL
!     Wavelength calibration flag of the buffer; set to .FALSE. (with a
!     warning) if it was .TRUE. and the file is not in REDUCEME format.
!------------------------------------------------------------------------------
        SUBROUTINE SLEEFITS_NOWAVECAL(LREDUCEME,LWAVECAL)
        IMPLICIT NONE
        LOGICAL LREDUCEME
        LOGICAL LWAVECAL
!
        IF(LWAVECAL.AND.(.NOT.LREDUCEME))THEN
          LWAVECAL=.FALSE.
          WRITE(*,'(A)') 'WARNING: wavelength calibration (CRPIX1, CRVAL1, '// &
           'CDELT1) ignored: NAXIS1 is not the X axis after rotation'
        END IF
        END
!
!******************************************************************************
!
        SUBROUTINE PRINTERROR(ISTATUS)
! Print out the FITSIO error messages to the user
        INTEGER ISTATUS
        CHARACTER ERRTEXT*30,ERRMESSAGE*80
! Check if status is OK (no error); if so, simply return
        IF(ISTATUS.LE.0) RETURN
! Get the text string which describes the error
        CALL FTGERR(ISTATUS,ERRTEXT)
        WRITE(*,'(A,$)') 'FITSIO Error Status = '
        WRITE(*,*) ISTATUS
        WRITE(*,'(A)') ERRTEXT
! Read and print out all the error messages on the FITSIO stack
        CALL FTGMSG(ERRMESSAGE)
        DO WHILE(ERRMESSAGE.NE.' ')
          WRITE(*,'(A)') ERRMESSAGE
          CALL FTGMSG(ERRMESSAGE)
        END DO
        END
