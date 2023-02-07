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
        INTEGER READI
        INTEGER TRUELEN
! variables globales (COMMONs)
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
!delete REAL IMAGEN_(NXYMAX,NXYMAX)
        REAL CRPIX1(NMAXBUFF),CRVAL1(NMAXBUFF),CDELT1(NMAXBUFF)
        LOGICAL LWAVECAL(NMAXBUFF)
        REAL STWV,DISP
        REAL AIRMASS,TIMEXPOS
!delete LOGICAL LNULL(NXMAX,NYMAX,NMAXBUFF)
        LOGICAL ANYNULL
        INTEGER NAXIS(2,NMAXBUFF)
! variables locales
        INTEGER NEW_HDU,HDUTYPE
        INTEGER JROW(NXYMAX)
        INTEGER I,J,L
        INTEGER FIRSTPIX
        INTEGER BITPIX
        INTEGER ISTATUS,IREADWRITE,IUNIT
        INTEGER BLOCKSIZE,NULLVAL
        INTEGER NKEYS,NSPACE,NFOUND
        INTEGER NAXES
        INTEGER NAXIS_(3)                                  !OJO: el limite es 2
        INTEGER NANYNULLS
        REAL FROW(NXYMAX)
        CHARACTER*50 COMMENT
        CHARACTER*80 CLINEA
        REAL PHOTFLAM(NMAXBUFF),PHOTZPT(NMAXBUFF),EXPTIME(NMAXBUFF)
        LOGICAL EXTEND
!delete LOGICAL LNULL_(NXYMAX,NXYMAX)
        LOGICAL LOGFILE,LANYNULL
        LOGICAL LROW(NXYMAX)
        LOGICAL L_PHOTFLAM(NMAXBUFF),L_PHOTZPT(NMAXBUFF)
        LOGICAL L_EXPTIME(NMAXBUFF)
        LOGICAL LERROR
!
!delete COMMON/BLKIMAGEN1/IMAGEN             !imagen FITS leida en formato REAL
!delete COMMON/BLKIMAGEN1_/IMAGEN_              !es global para ahorrar memoria
        COMMON/BLKWAVECAL1/CRPIX1,CRVAL1,CDELT1         !wavelength calibration
        COMMON/BLKWAVECAL2/LWAVECAL                     !wavelength calibration
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
            LWAVECAL(NEWBUFF)=.TRUE.
            CRPIX1(NEWBUFF)=1.0
            CRVAL1(NEWBUFF)=STWV
            CDELT1(NEWBUFF)=DISP
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
        IF(EXTEND)THEN
          WRITE(*,101) '***WARNING***'
          WRITE(*,101) '=> this file contains extensions'
          NEW_HDU=READI('Extension number to be read (1=primary)','1')
          CALL FTMAHD(IUNIT,NEW_HDU,HDUTYPE,ISTATUS)
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
! leemos BITPIX
        CALL FTGKYJ(IUNIT,'BITPIX',BITPIX,COMMENT,ISTATUS)
        IF(.NOT.LSHOW)THEN
          WRITE(*,*)
          WRITE(*,100) 'CFITSIO> BITPIX: '
          WRITE(*,*) BITPIX
        END IF
! leemos NAXIS, NAXIS1 y NAXIS2
        CALL FTGKYJ(IUNIT,'NAXIS',NAXES,COMMENT,ISTATUS)
        IF(.NOT.LSHOW)THEN
          WRITE(*,100) 'CFITSIO> NAXIS : '
          WRITE(*,*) NAXES 
        END IF
        LERROR=.FALSE.
        IF(NAXES.EQ.1)THEN
          CALL FTGKYJ(IUNIT,'NAXIS1',NAXIS_(1),COMMENT,ISTATUS)
          NAXIS_(2)=1
        ELSEIF(NAXES.EQ.2)THEN
          CALL FTGKYJ(IUNIT,'NAXIS1',NAXIS_(1),COMMENT,ISTATUS)
          CALL FTGKYJ(IUNIT,'NAXIS2',NAXIS_(2),COMMENT,ISTATUS)
        ELSEIF(NAXES.EQ.3)THEN
          CALL FTGKYJ(IUNIT,'NAXIS1',NAXIS_(1),COMMENT,ISTATUS)
          CALL FTGKYJ(IUNIT,'NAXIS2',NAXIS_(2),COMMENT,ISTATUS)
          CALL FTGKYJ(IUNIT,'NAXIS3',NAXIS_(3),COMMENT,ISTATUS)
          IF(NAXIS_(3).EQ.1)THEN
            WRITE(*,101) 'WARNING: NAXIS=3 with NAXIS3=1. Reading anyway!'
          ELSE
            LERROR=.TRUE.
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
          WRITE(*,100) 'CFITSIO> NAXIS1: '
          WRITE(*,*) NAXIS_(1)
          WRITE(*,100) 'CFITSIO> NAXIS2: '
          WRITE(*,*) NAXIS_(2)
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
        IF((IROTATE.EQ.0).AND.(IROTATE.EQ.180).OR.(IROTATE.EQ.-180))THEN
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
! leemos calibracion en longitud de onda
        LWAVECAL(NEWBUFF)=.TRUE.           !salvo que se demuestre lo contrario
        CALL FTGKYE(IUNIT,'CRPIX1',CRPIX1(NEWBUFF),COMMENT,ISTATUS)
        IF(ISTATUS.GT.0)THEN
          CRPIX1(NEWBUFF)=1.0
          ISTATUS=0
        END IF
        CALL FTGKYE(IUNIT,'CRVAL1',CRVAL1(NEWBUFF),COMMENT,ISTATUS)
        IF(ISTATUS.GT.0)THEN
          LWAVECAL(NEWBUFF)=.FALSE.
          CRVAL1(NEWBUFF)=0.0
          ISTATUS=0
        END IF
        CALL FTGKYE(IUNIT,'CDELT1',CDELT1(NEWBUFF),COMMENT,ISTATUS)
        IF(ISTATUS.GT.0)THEN
          LWAVECAL(NEWBUFF)=.FALSE.
          CDELT1(NEWBUFF)=0.0
          ISTATUS=0
        END IF
        IF(LWAVECAL(NEWBUFF))THEN
          WRITE(*,100) 'CFITSIO> CRPIX1: '
          WRITE(*,*) CRPIX1(NEWBUFF)
          WRITE(*,100) 'CFITSIO> CRVAL1: '
          WRITE(*,*) CRVAL1(NEWBUFF)
          WRITE(*,100) 'CFITSIO> CDELT1: '
          WRITE(*,*) CDELT1(NEWBUFF)
        END IF
! leemos la imagen
        IF((BITPIX.EQ.8).OR.(BITPIX.EQ.16))THEN
          DO I=1,NAXIS_(2)
            FIRSTPIX=(I-1)*NAXIS_(1)+1
            CALL FTGPFJ(IUNIT,1,FIRSTPIX,NAXIS_(1),JROW(1),LROW(1),ANYNULL,ISTATUS)
            IF(ANYNULL)THEN
              DO J=1,NAXIS_(1)
                LNULL_(J,I)=LROW(J)
                IF(LNULL_(J,I))THEN
                  NANYNULLS=NANYNULLS+1
                  IF(NANYNULLS.LT.10)THEN
                    print*,j,i,jrow(j)
                  END IF
                  JROW(J)=0
                END IF
              END DO
              LANYNULL=.TRUE.
            END IF
            DO J=1,NAXIS_(1)
              IMAGEN_(J,I)=REAL(JROW(J))
            END DO
          END DO
!        ELSEIF(BITPIX.EQ.32)THEN
!          DO I=1,NAXIS_(2)
!            FIRSTPIX=(I-1)*NAXIS_(1)+1
!            CALL FTGPFJ(IUNIT,1,FIRSTPIX,NAXIS_(1),JROW(1),LROW(1),
!     +       ANYNULL,ISTATUS)
!            DO J=1,NAXIS_(1)
!              IMAGEN_(J,I)=REAL(JROW(J))
!            END DO
!            IF(ANYNULL)THEN
!              DO J=1,NAXIS_(1)
!                LNULL_(J,I)=LROW(J)
!              END DO
!              LANYNULL=.TRUE.
!            END IF
!          END DO
        ELSEIF((BITPIX.EQ.32).OR.(BITPIX.EQ.-32).OR.(BITPIX.EQ.-64))THEN
          DO I=1,NAXIS_(2)
            FIRSTPIX=(I-1)*NAXIS_(1)+1
            CALL FTGPFE(IUNIT,1,FIRSTPIX,NAXIS_(1),FROW(1),LROW(1),ANYNULL,ISTATUS)
            IF(ANYNULL)THEN
              DO J=1,NAXIS_(1)
                LNULL_(J,I)=LROW(J)
                IF(LNULL_(J,I))THEN
                  NANYNULLS=NANYNULLS+1
                  IF(NANYNULLS.LT.10)THEN
                    print*,j,i,frow(j)
                  END IF
                  FROW(J)=0.0
                END IF
              END DO
              LANYNULL=.TRUE.
            END IF
            DO J=1,NAXIS_(1)
              IMAGEN_(J,I)=FROW(J)
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
        ELSEIF((IROTATE.EQ.180).OR.(IROTATE.EQ.-180))THEN
          DO I=1,NAXIS_(2)
            DO J=1,NAXIS_(1)
              IMAGEN(J,I,NEWBUFF)=IMAGEN_(NAXIS_(1)-J+1,NAXIS_(2)-I+1)
              LNULL(J,I,NEWBUFF)=LNULL_(J,I)
            END DO
          END DO
          NAXIS(1,NEWBUFF)=NAXIS_(1)
          NAXIS(2,NEWBUFF)=NAXIS_(2)
        ELSEIF((IROTATE.EQ.270).OR.(IROTATE.EQ.-90))THEN
          DO I=1,NAXIS_(2)
            DO J=1,NAXIS_(1)
              IMAGEN(I,NAXIS_(1)-J+1,NEWBUFF)=IMAGEN_(J,I)
              LNULL(I,NAXIS_(1)-J+1,NEWBUFF)=LNULL_(J,I)
            END DO
          END DO
          NAXIS(1,NEWBUFF)=NAXIS_(2)
          NAXIS(2,NEWBUFF)=NAXIS_(1)
          CRPIX1(NEWBUFF)=1.0
          CRVAL1(NEWBUFF)=STWV
          CDELT1(NEWBUFF)=DISP
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
        END IF
        CALL Deallocate_Array_LNULL_
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
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
