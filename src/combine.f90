! Version 6-November-2000
!
! Nota: en la combinacion comprobamos que el fichero de errores no contiene
! sen~al negativa (la sen~al negativa en la imagen de errores se emplea, por 
! tanto, como mascara para indicar que pixels pueden sumarse; esto es esencial
! cuando sumamos imagenes que contienen huecos y queremos que el numero de
! cuentas se conserve adecuadamente).
        SUBROUTINE COMBINE(NCBUFF)
        USE Dynamic_Array_IMAGEN
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
        INTEGER NCBUFF
!
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
!
        INCLUDE 'dimensions.inc'
!
        INTEGER READILIM
        INTEGER TRUEBEG,TRUELEN
        REAL READF
        REAL FMEAN0
        REAL FMEDIAN1
        REAL FTSTUDENTI
        CHARACTER*255 READC
        LOGICAL INSIDE
!
        INTEGER I,J,K,II,JJ,KK,II1,II2,JJ1,JJ2,IJ
        INTEGER NK,NN
        INTEGER L,IDUM
        INTEGER L1,L2
        INTEGER DI(9),DJ(9)
        INTEGER NEWBUFF1,NEWBUFF2
        INTEGER NFRAMES(NMAXBUFF)
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER NAXISFRAME(2,9,NMAXBUFF)
        INTEGER NSIZEFB9(NMAXBUFF)
        INTEGER NF,NFRAMES_
        INTEGER NEXTRA_XMINB9,NEXTRA_XMAXB9
        INTEGER NEXTRA_YMINB9,NEXTRA_YMAXB9
        INTEGER NREMOV_XMINB9,NREMOV_XMAXB9
        INTEGER NREMOV_YMINB9,NREMOV_YMAXB9
        INTEGER J0(NBOXMAX),I0(NBOXMAX)
        INTEGER NXMAXB9_,NYMAXB9_
        INTEGER NXMAXB9_EFF,NYMAXB9_EFF
        INTEGER NPIXREGION
        INTEGER K_LAST(9*4),II_LAST(9*4),JJ_LAST(9*4)
        INTEGER NPIXELSVECINOS
        INTEGER NEXTINFO
        INTEGER NX1,NX2,NY1,NY2
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL FPIXUSED(NXMAXB9,NYMAXB9)
        REAL FOFFSETX(NBOXMAX),FOFFSETY(NBOXMAX)
        REAL FOFFSETX_,FOFFSETY_
        REAL EFOFFSETX_,EFOFFSETY_
        REAL FJ0(NBOXMAX),FI0(NBOXMAX)
        REAL FPIXELTHRESHOLDCR
        REAL FRACCION_PIXEL
        REAL PIXEL(9*4),EPIXEL(9*4),PIXELF(9*4)
        REAL PIXELW(9*4),SUMPIXELW
        REAL PIXELVECINO(24)
        REAL FMEAN_PIX,FSIGMA_PIX
        REAL TSHOT
        REAL FMEAN,FSIGMA,FMEDIAN,FMIN,FMAX
        REAL BG,FG
        REAL FERRORLIMIT
        CHARACTER*1 CREJECT,CWEIGHT,CWMODE,CMASK
        CHARACTER*50 CDUMMY
        CHARACTER*80 COFFSETFILE
        LOGICAL LOGFILE
        LOGICAL LZERO
        LOGICAL LMORE
        LOGICAL LANYCR
!
!delete COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKNFRAMES/NFRAMES
        COMMON/BLKXYLIMPLOT/NX1,NX2,NY1,NY2
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKNAXISFRAME/NAXISFRAME
        COMMON/BLKNSIZEFB9/NSIZEFB9
        COMMON/BLKESTADISTICA/FMEAN,FSIGMA,FMEDIAN,FMIN,FMAX
!------------------------------------------------------------------------------
! Note: the pattern of the frames in box-9 is the following:
!       6 9 4
!       3 1 7
!       8 5 2
! offsets of each frame in the 768x768 composite mosaic
        DATA (DI(K),DJ(K),K=1,NBOXMAX) / &
!!!     +   256,256,  !1
!!!     +   000,512,  !2
!!!     +   256,000,  !3
!!!     +   512,512,  !4
!!!     +   000,256,  !5
!!!     +   512,000,  !6
!!!     +   256,512,  !7
!!!     +   000,000,  !8
!!!     +   512,256/  !9
         1,1, & !1
         0,2, & !2
         1,0, & !3
         2,2, & !4
         0,1, & !5
         2,0, & !6
         1,2, & !7
         0,0, & !8
         2,1/  !9
!
        II2=0   !evita un WARNING de compilacion
        JJ2=0   !evita un WARNING de compilacion
        FRACCION_PIXEL=0.0 !evita un WARNING de compilacion
        FERRORLIMIT=0.0 !evita un WARNING de compilacion
!------------------------------------------------------------------------------
        CREJECT='n'
        CWEIGHT='n'
        FPIXELTHRESHOLDCR=1.0
!
        CWMODE='X'
        CMASK='X'
!------------------------------------------------------------------------------
! pedimos offsets
        LOGFILE=.FALSE.
        DO WHILE(.NOT.LOGFILE)
          COFFSETFILE(1:80)=READC('Name of the file with offsets (none=exit)','offsets.dat','@')
          L1=TRUEBEG(COFFSETFILE)
          L2=TRUELEN(COFFSETFILE)
          IF(COFFSETFILE(L1:L2).EQ.'none') RETURN !permitimos escapar
          INQUIRE(FILE=COFFSETFILE(L1:L2),EXIST=LOGFILE)
          IF(.NOT.LOGFILE)THEN
            WRITE(*,101) 'ERROR: this file does not exist. Try again.'
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END DO
        WRITE(*,101) 'Reading offsets...'
        OPEN(20,FILE=COFFSETFILE,STATUS='OLD',FORM='FORMATTED')
        NF=0
5        READ(20,*,END=6) IDUM,FOFFSETX_,EFOFFSETX_,FOFFSETY_,EFOFFSETY_
        NF=NF+1
        IF(NF.LE.NBOXMAX)THEN
          FOFFSETX(NF)=FOFFSETX_
          FOFFSETY(NF)=FOFFSETY_
          WRITE(*,'(I1,2(2X,F8.3))') NF,FOFFSETX(NF),FOFFSETY(NF)
        ELSE
          CLOSE(20)
          WRITE(*,101) 'ERROR: file with offsets contains too many lines (>9)'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
        GOTO 5
6        CLOSE(20)
        NFRAMES_=NF
        WRITE(*,101) '...OK! File read and closed.'
        WRITE(*,100) '> No. of frames/image: '
        WRITE(*,*) NFRAMES_
        WRITE(*,100) '> NSIZEFB9...........: '
        WRITE(*,*) NSIZEFB9(NCBUFF)
! proteccion #1
        IF(NFRAMES_.NE.NFRAMES(NCBUFF))THEN
          WRITE(*,*)
          WRITE(*,100) '> NFRAMES(NCBUFF),NFRAMES_: '
          WRITE(*,*) NFRAMES(NCBUFF),NFRAMES_
          WRITE(*,100) 'ERROR: number of frames does not correspond'
          WRITE(*,101) ' with expected value.'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
! proteccion #2
        DO NF=1,NFRAMES_
          IF((NAXISFRAME(1,NF,NCBUFF).EQ.0).OR.(NAXISFRAME(2,NF,NCBUFF).EQ.0))THEN
            WRITE(*,*)
            WRITE(*,100) '> NAXISFRAME(1,NF,NCBUFF): '
            WRITE(*,*) NAXISFRAME(1,NF,NCBUFF)
            WRITE(*,100) '> NAXISFRAME(2,NF,NCBUFF): '
            WRITE(*,*) NAXISFRAME(2,NF,NCBUFF)
            WRITE(*,101) 'ERROR: invalid dimensions of individual frame'
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
            RETURN
          END IF
        END DO
! estimamos partes enteras y fraccionarias de los offsets, y determinamos la
! dimension maxima que vamos a necesitar para la imagen combinada
        NEXTRA_XMINB9=0
        NEXTRA_XMAXB9=0
        NEXTRA_YMINB9=0
        NEXTRA_YMAXB9=0
        DO NF=1,NFRAMES_
          J0(NF)=INT(FOFFSETX(NF))
          FJ0(NF)=FOFFSETX(NF)-REAL(J0(NF))
          IF(FOFFSETX(NF).GT.0.0)THEN
            IF(FJ0(NF).EQ.0.0)THEN
              NEXTRA_XMINB9=AMIN0(NEXTRA_XMINB9,-J0(NF))
            ELSE
              NEXTRA_XMINB9=AMIN0(NEXTRA_XMINB9,-J0(NF)-1)
            END IF
          ELSEIF(FOFFSETX(NF).LE.0.0)THEN
            IF(FJ0(NF).EQ.0.0)THEN
              NEXTRA_XMAXB9=AMAX0(NEXTRA_XMAXB9,-J0(NF))
            ELSE
              NEXTRA_XMAXB9=AMAX0(NEXTRA_XMAXB9,-J0(NF)+1)
            END IF
          END IF
          I0(NF)=INT(FOFFSETY(NF))
          FI0(NF)=FOFFSETY(NF)-REAL(I0(NF))
          IF(FOFFSETY(NF).GT.0.0)THEN
            IF(FI0(NF).EQ.0.0)THEN
              NEXTRA_YMINB9=AMIN0(NEXTRA_YMINB9,-I0(NF))
            ELSE
              NEXTRA_YMINB9=AMIN0(NEXTRA_YMINB9,-I0(NF)-1)
            END IF
          ELSEIF(FOFFSETY(NF).LE.0.0)THEN
            IF(FI0(NF).EQ.0.0)THEN
              NEXTRA_YMAXB9=AMAX0(NEXTRA_YMAXB9,-I0(NF))
            ELSE
              NEXTRA_YMAXB9=AMAX0(NEXTRA_YMAXB9,-I0(NF)+1)
            END IF
          END IF
        END DO
        WRITE(*,*)
        WRITE(*,100) '>>> Initial extra limits of combined frame: '
        WRITE(*,*) NEXTRA_XMINB9,NEXTRA_XMAXB9,NEXTRA_YMINB9,NEXTRA_YMAXB9
        NXMAXB9_=NSIZEFB9(NCBUFF)-NEXTRA_XMINB9+NEXTRA_XMAXB9
        NYMAXB9_=NSIZEFB9(NCBUFF)-NEXTRA_YMINB9+NEXTRA_YMAXB9
! mostramos taman~o
        WRITE(*,100) '>>> New dimensions of combined frame: '
        WRITE(CDUMMY,'(I10,A1,I10)') NXMAXB9_,'x',NYMAXB9_
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101) CDUMMY(1:L)
        NPIXREGION=NXMAXB9_*NYMAXB9_
! proteccion
        IF(NXMAXB9_.GT.NXMAXB9)THEN
          WRITE(*,101) 'ERROR: X dimension > NXMAXB9'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
        IF(NYMAXB9_.GT.NYMAXB9)THEN
          WRITE(*,101) 'ERROR: Y dimension > NYMAXB9'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
!------------------------------------------------------------------------------
! definimos los buffers de trabajo
        NEWBUFF1=NCBUFF+1
        IF(NEWBUFF1.GT.NMAXBUFF/2) NEWBUFF1=1
        WRITE(CDUMMY,*) NEWBUFF1
        NEWBUFF1=READILIM('First  auxiliary buffer # to store manipulated data..',CDUMMY,1,NMAXBUFF/2)
!
        NEWBUFF2=NEWBUFF1+1
        IF(NEWBUFF2.GT.NMAXBUFF/2) NEWBUFF2=1
        WRITE(CDUMMY,*) NEWBUFF2
        NEWBUFF2=READILIM('Second auxiliary buffer # to store manipulated data..',CDUMMY,1,NMAXBUFF/2)
!------------------------------------------------------------------------------
! decidimos si eliminamos pixels extran~os
        WRITE(*,*)
        CREJECT(1:1)=READC('Reject hot/cool pixels (y/n)',CREJECT,'yn')
        IF(CREJECT.EQ.'y')THEN
! numero de pixels umbral para deteccion de pixels calientes
          WRITE(CDUMMY,*) FPIXELTHRESHOLDCR
          FPIXELTHRESHOLDCR=READF('Threshold in no. of pixels > Times_Sigma for hot/cool pixels',CDUMMY)
        END IF
!------------------------------------------------------------------------------
! decidimos si pesamos la imagen suma con los errores
        CWEIGHT(1:1)=READC('Perform weighted sum (y/n)',CWEIGHT,'yn')
        IF(CWEIGHT.EQ.'y')THEN
          CWMODE(1:1)=READC('Weigth with [e]rrors or [s]ignal/noise','e','es')
          WRITE(*,101) '* Note: data can be masked if error is below a given threshold error'
          CMASK(1:1)=READC('Do you want to mask data (y/n)','y','yn')
          IF(CMASK.EQ.'y')THEN
            FERRORLIMIT=READF('Threshold error to mask data','0.0')
          END IF
        END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! it is time to coadd the frames into a single image
        NAXIS(1,NEWBUFF1)=NXMAXB9_
        NAXIS(2,NEWBUFF1)=NYMAXB9_
        NFRAMES(NEWBUFF1)=1
        NAXIS(1,NEWBUFF1+NMAXBUFF/2)=NXMAXB9_
        NAXIS(2,NEWBUFF1+NMAXBUFF/2)=NYMAXB9_
        NFRAMES(NEWBUFF1+NMAXBUFF/2)=1
        DO I=1,NAXIS(2,NEWBUFF1)
          DO J=1,NAXIS(1,NEWBUFF1)
            IMAGEN(J,I,NEWBUFF1)=0. !inicializamos a cero el nuevo frame
            FPIXUSED(J,I)=0. !aqui almacenamos el numero de pixels usados en
                             !el calculo de cada nuevo pixel en la imagen final
            IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=0. !imagen de errores
          END DO
        END DO
!------------------------------------------------------------------------------
! recorremos cada pixel en la nueva imagen y determinamos cuantos pixels de 
! los frames iniciales se pueden usar en cada caso; aqui hay que tener cuidado
! porque, al contrario de lo que ocurre al combinar un box-9 formado por 
! imagenes de 256x256 pixels, en este caso tenemos que combinar varias imagenes
! que inicialmente pueden tener taman~os distintos!
        LANYCR=.FALSE.
        NEXTINFO=0
        DO I=1,NAXIS(2,NEWBUFF1)
          DO J=1,NAXIS(1,NEWBUFF1)
            KK=0 !numero total de pixels (o fracciones de pixels) a utilizar
            DO K=1,NFRAMES_
              JJ1=J+NEXTRA_XMINB9+J0(K) !ojo: NEXTRA_XMINB9 es <=0
              IF(FJ0(K).EQ.0.0)THEN
                JJ2=0
              ELSEIF(FJ0(K).GT.0.0)THEN
                JJ2=JJ1+1
              ELSEIF(FJ0(K).LT.0.0)THEN
                JJ2=JJ1-1
              END IF
              II1=I+NEXTRA_YMINB9+I0(K) !ojo: NEXTRA_YMINB9 es <=0
              IF(FI0(K).EQ.0.0)THEN
                II2=0
              ELSEIF(FI0(K).GT.0.0)THEN
                II2=II1+1
              ELSEIF(FI0(K).LT.0.0)THEN
                II2=II1-1
              END IF
              DO IJ=1,4 !recorremos los 4 posibles pixels
                IF(IJ.EQ.1)THEN
                  JJ=JJ1
                  II=II1
                  FRACCION_PIXEL=(1.-ABS(FJ0(K)))*(1.-ABS(FI0(K)))
                ELSEIF(IJ.EQ.2)THEN
                  JJ=JJ1
                  II=II2
                  FRACCION_PIXEL=(1.-ABS(FJ0(K)))*ABS(FI0(K))
                ELSEIF(IJ.EQ.3)THEN
                  JJ=JJ2
                  II=II1
                  FRACCION_PIXEL=ABS(FJ0(K))*(1.-ABS(FI0(K)))
                ELSEIF(IJ.EQ.4)THEN
                  JJ=JJ2
                  II=II2
                  FRACCION_PIXEL=ABS(FJ0(K))*ABS(FI0(K))
                END IF
                IF(INSIDE(JJ,1,NAXISFRAME(1,K,NCBUFF),II,1,NAXISFRAME(2,K,NCBUFF)))THEN
                  IF(IMAGEN(JJ+DJ(K)*NSIZEFB9(NCBUFF),II+DI(K)*NSIZEFB9(NCBUFF),NCBUFF+NMAXBUFF/2).GE.0.0)THEN
                    KK=KK+1
                    PIXEL(KK)=IMAGEN(JJ+DJ(K)*NSIZEFB9(NCBUFF),II+DI(K)*NSIZEFB9(NCBUFF),NCBUFF)
                    EPIXEL(KK)=IMAGEN(JJ+DJ(K)*NSIZEFB9(NCBUFF),II+DI(K)*NSIZEFB9(NCBUFF),NCBUFF+NMAXBUFF/2)
                    PIXELF(KK)=FRACCION_PIXEL
                    K_LAST(KK)=K !necesario para eliminar hot/cool pixels
                    JJ_LAST(KK)=JJ !idem
                    II_LAST(KK)=II !idem
                  END IF
                END IF
              END DO
            END DO
!..............................................................................
!..............................................................................
! En caso necesario, eliminamos hot/cool pixels.
            IF(CREJECT.EQ.'y')THEN
!----->
              IF(KK.LE.4)THEN !si tenemos menos de cinco pixels
                DO NK=1,KK
                  NPIXELSVECINOS=0
                  DO II=II_LAST(NK)-2,II_LAST(NK)+2
                    IF((II.GE.1).AND.(II.LE.NAXISFRAME(2,K_LAST(NK),NCBUFF)))THEN
                      DO JJ=JJ_LAST(NK)-2,JJ_LAST(NK)+2
                        IF((JJ.GE.1).AND.(JJ.LE.NAXISFRAME(1,K_LAST(NK),NCBUFF)))THEN
                          IF((II.NE.II_LAST(NK)).OR.(JJ.NE.JJ_LAST(NK)))THEN
                            NPIXELSVECINOS=NPIXELSVECINOS+1
                            PIXELVECINO(NPIXELSVECINOS)= &
                             IMAGEN(JJ+DJ(K_LAST(NK))*NSIZEFB9(NCBUFF), &
                                    II+DI(K_LAST(NK))*NSIZEFB9(NCBUFF),NCBUFF)
                          END IF
                        END IF
                      END DO
                    END IF
                  END DO
                  CALL ORDENA1F(NPIXELSVECINOS,PIXELVECINO)
                  FMEAN_PIX=FMEAN0(NPIXELSVECINOS,PIXELVECINO,FSIGMA_PIX)
                  TSHOT=FTSTUDENTI(NPIXELSVECINOS-2, & !grados de libertad
                   FPIXELTHRESHOLDCR/REAL(NPIXREGION)) !probabilidad
                  IF(ABS(PIXEL(NK)-FMEAN_PIX).GT.TSHOT*FSIGMA_PIX)THEN
                    PIXEL(NK)=FMEDIAN1(NPIXELSVECINOS,PIXELVECINO)
                    EPIXEL(NK)=FSIGMA_PIX
                    CALL PGSCI(4)
                    CALL PGPOINT(1, &
                     REAL(JJ_LAST(NK)+DJ(K_LAST(NK))*NSIZEFB9(NCBUFF)), &
                     REAL(II_LAST(NK)+DI(K_LAST(NK))*NSIZEFB9(NCBUFF)),21)
                    CALL PGSCI(1)
                    LANYCR=.TRUE.
                  END IF
                END DO
!----->
              ELSE !si tenemos mas de 4 pixels
! ordenamos los pixels por sen~al creciente
                CALL ORDENA3F3I(KK,PIXEL,EPIXEL,PIXELF,K_LAST,II_LAST,JJ_LAST)
! calculamos la media y la dispersion eliminando el pixel mas brillante y
! el mas debil
                FMEAN_PIX=FMEAN0(KK-2,PIXEL(2),FSIGMA_PIX)
                TSHOT=FTSTUDENTI(KK-4, &!grados de libertad
                 FPIXELTHRESHOLDCR/REAL(NPIXREGION)) !/ probabilidad
! eliminamos el pixel mas brillante si excede en TSHOT veces sigma
                IF(PIXEL(KK)-FMEAN_PIX.GT.TSHOT*FSIGMA_PIX)THEN
                  CALL PGSCI(5)
                  CALL PGPOINT(1, &
                   REAL(JJ_LAST(KK)+DJ(K_LAST(KK))*NSIZEFB9(NCBUFF)), &
                   REAL(II_LAST(KK)+DI(K_LAST(KK))*NSIZEFB9(NCBUFF)),21)
                  CALL PGSCI(1)
                  KK=KK-1
                  LANYCR=.TRUE.
                END IF
! eliminamos el pixel mas debil si excede en -TSHOT veces sigma
                IF(FMEAN_PIX-PIXEL(1).GT.TSHOT*FSIGMA_PIX)THEN
                  CALL PGSCI(6)
                  CALL PGPOINT(1, &
                   REAL(JJ_LAST(1)+DJ(K_LAST(1))*NSIZEFB9(NCBUFF)), &
                   REAL(II_LAST(1)+DI(K_LAST(1))*NSIZEFB9(NCBUFF)),21)
                  CALL PGSCI(1)
                  DO NN=1,KK-1
                    PIXEL(NN)=PIXEL(NN+1)
                    EPIXEL(NN)=EPIXEL(NN+1)
                    PIXELF(NN)=PIXELF(NN+1)
                  END DO
                  KK=KK-1
                  LANYCR=.TRUE.
                END IF
              END IF
            END IF
!..............................................................................
!..............................................................................
! normalizamos cada pixel por el numero total de frames que hemos utilizado
! para calcular la sen~al total en dicho pixel
            IF(KK.EQ.0)THEN !lo dejamos por claridad
!!!              IMAGEN(J,I,NEWBUFF1)=0.
!!!              IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=0.
!!!              FPIXUSED(J,I)=0.
            ELSE
              IF(CWEIGHT.EQ.'y')THEN                   !suma pesada con errores
                LZERO=.FALSE.
                IF(CMASK.EQ.'n')THEN !chequeamos si hay algun error que sea 0
                  DO K=1,KK
                    IF(EPIXEL(K).EQ.0.0)THEN
                      LZERO=.TRUE.
                    END IF
                  END DO
                ELSEIF(CMASK.EQ.'y')THEN
                ELSE
                  WRITE(*,101) 'FATAL ERROR: unexpected CMASK value'
                  CALL Deallocate_Array_IMAGEN
                  CALL Deallocate_Array_IMAGEN_
                  STOP
                END IF
                IF(LZERO)THEN    !si hay algun error nulo, hacemos pesos de los
                                 !datos con error nulo igual a 1, y los pesos 
                                 !de los demas datos igual a 0
                  DO K=1,KK
                    IF(EPIXEL(K).EQ.0.0)THEN
                      PIXELW(K)=1.0
                    ELSE
                      PIXELW(K)=0.0
                    END IF
                  END DO
                ELSE !si no hay errores nulos, tomamos los pesos normales
                  IF(CWMODE.EQ.'e')THEN
                    DO K=1,KK !calculamos pesos iniciales
                      IF(CMASK.EQ.'y')THEN
                        IF(EPIXEL(K).LE.FERRORLIMIT)THEN
                          PIXELW(K)=0.
                        ELSE
                          PIXELW(K)=1./(EPIXEL(K)*EPIXEL(K))
                        END IF
                      ELSE
                        PIXELW(K)=1./(EPIXEL(K)*EPIXEL(K))
                      END IF
                    END DO
                  ELSEIF(CWMODE.EQ.'s')THEN
                    DO K=1,KK !calculamos pesos iniciales
                      IF(CMASK.EQ.'y')THEN
                        IF(EPIXEL(K).LE.FERRORLIMIT)THEN
                          PIXELW(K)=0.
                        ELSE
                          PIXELW(K)=PIXEL(K)*PIXEL(K)/(EPIXEL(K)*EPIXEL(K))
                        END IF
                      ELSE
                        PIXELW(K)=PIXEL(K)*PIXEL(K)/(EPIXEL(K)*EPIXEL(K))
                      END IF
                    END DO
                  ELSE
                    WRITE(*,101) 'FATAL ERROR: invalid CWMODE'
                    CALL Deallocate_Array_IMAGEN
                    CALL Deallocate_Array_IMAGEN_
                    STOP
                  END IF
                END IF
              ELSE                                  !suma no pesada con errores
                DO K=1,KK !no pesamos <==> hacemos todos pesos iguales
                  PIXELW(K)=1.0
                END DO
              END IF
              DO K=1,KK !an~adimos al peso la contribucion de la fraccion de
                        !pixel utilizada
                PIXELW(K)=PIXELW(K)*PIXELF(K)
              END DO
              SUMPIXELW=0.0 !normalizamos los pesos al numero total de pixels
                            !empleados para calcular el nuevo pixel (J,I)
              DO K=1,KK
                SUMPIXELW=SUMPIXELW+PIXELW(K)
                FPIXUSED(J,I)=FPIXUSED(J,I)+PIXELF(K)
              END DO
              IF(SUMPIXELW.EQ.0.0)THEN
                DO K=1,KK
                  PIXELW(K)=PIXELF(K)
                  SUMPIXELW=SUMPIXELW+PIXELW(K)
                END DO
              END IF
              SUMPIXELW=SUMPIXELW/FPIXUSED(J,I)
              DO K=1,KK
                PIXELW(K)=PIXELW(K)/SUMPIXELW
              END DO
              DO K=1,KK
!!!                IMAGEN(J,I,NEWBUFF1)=IMAGEN(J,I,NEWBUFF1)+
!!!     +           PIXEL(K)*PIXELW(K)*PIXELF(K)
                IMAGEN(J,I,NEWBUFF1)=IMAGEN(J,I,NEWBUFF1)+PIXEL(K)*PIXELW(K)
!!!                IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=
!!!     +           IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)+
!!!     +           EPIXEL(K)*EPIXEL(K)*PIXELW(K)*PIXELW(K)*
!!!     +           PIXELF(K)*PIXELF(K)
                IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)+EPIXEL(K)*EPIXEL(K)*PIXELW(K)*PIXELW(K)
!!!                FPIXUSED(J,I)=FPIXUSED(J,I)+PIXELF(K)
              END DO
              IMAGEN(J,I,NEWBUFF1)=IMAGEN(J,I,NEWBUFF1)/FPIXUSED(J,I)
              IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=SQRT(IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2))/FPIXUSED(J,I)
            END IF
          END DO
          CALL SHOWPERC(1,NAXIS(2,NEWBUFF1),1,I,NEXTINFO)
        END DO
!..............................................................................
! salvamos en otro buffer el numero de pixels utilizados para el calculo de 
! cada pixel (ademas, en aquellos pixels en los que no hemos sumado nada, 
! introducimos un error negativo)
        NAXIS(1,NEWBUFF2)=NXMAXB9_
        NAXIS(2,NEWBUFF2)=NYMAXB9_
        NFRAMES(NEWBUFF2)=1
        DO I=1,NAXIS(2,NEWBUFF2)
          DO J=1,NAXIS(1,NEWBUFF2)
            IMAGEN(J,I,NEWBUFF2)=FPIXUSED(J,I)
            IF(FPIXUSED(J,I).LE.0.0)THEN
              IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=-1.
            END IF
          END DO
        END DO

!------------------------------------------------------------------------------
! refinamos las dimensiones iniciales y calculamos unas dimensiones mas
! ajustadas al taman~o real de la imagen combinada (notar que el limite en la
! busqueda lo establezco a -0.5, cuando sabemos que los errores se han
! establecido en -1.0 para indicar pixels no utilizados; de esta forma
! garantizamos que no cometemos ningun error de redondeo al comparar numeros
! reales)
        NREMOV_XMINB9=0
        J=0
        LMORE=.TRUE.
        DO WHILE((J.LT.NXMAXB9_).AND.(LMORE))
          J=J+1
          DO I=1,NYMAXB9_
            IF(IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2).GT.-0.5) LMORE=.FALSE.
          END DO
          IF(LMORE) NREMOV_XMINB9=NREMOV_XMINB9+1
        END DO
!
        NREMOV_XMAXB9=0
        J=NXMAXB9_+1
        LMORE=.TRUE.
        DO WHILE((J.GT.1).AND.(LMORE))
          J=J-1
          DO I=1,NYMAXB9_
            IF(IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2).GT.-0.5) LMORE=.FALSE.
          END DO
          IF(LMORE) NREMOV_XMAXB9=NREMOV_XMAXB9+1
        END DO
!
        NREMOV_YMINB9=0
        I=0
        LMORE=.TRUE.
        DO WHILE((I.LT.NYMAXB9_).AND.(LMORE))
          I=I+1
          DO J=1,NXMAXB9_
            IF(IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2).GT.-0.5) LMORE=.FALSE.
          END DO
          IF(LMORE) NREMOV_YMINB9=NREMOV_YMINB9+1
        END DO
!
        NREMOV_YMAXB9=0
        I=NYMAXB9_+1
        LMORE=.TRUE.
        DO WHILE((I.GT.1).AND.(LMORE))
          I=I-1
          DO J=1,NXMAXB9_
            IF(IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2).GT.-0.5) LMORE=.FALSE.
          END DO
          IF(LMORE) NREMOV_YMAXB9=NREMOV_YMAXB9+1
        END DO
!
        WRITE(*,100) '> No. pixels to remove in Xmin, Xmax, Ymin, Ymax: '
        WRITE(*,*) NREMOV_XMINB9,NREMOV_XMAXB9,NREMOV_YMINB9,NREMOV_YMAXB9
        NXMAXB9_EFF=NXMAXB9_-NREMOV_XMINB9-NREMOV_XMAXB9
        NYMAXB9_EFF=NYMAXB9_-NREMOV_YMINB9-NREMOV_YMAXB9
        WRITE(*,100) '>>> Previous dimensions of  combined frame: '
        WRITE(CDUMMY,'(I10,A1,I10)') NXMAXB9_,'x',NYMAXB9_
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101) CDUMMY(1:L)
        WRITE(*,100) '>>> New dimensions of final combined frame: '
        WRITE(CDUMMY,'(I10,A1,I10)') NXMAXB9_EFF,'x',NYMAXB9_EFF
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101) CDUMMY(1:L)
!..............................................................................
! movemos imagenes y recortamos al taman~o final
        IF((NREMOV_XMINB9.GT.0).OR.(NREMOV_YMINB9.GT.0))THEN
          DO I=1,NYMAXB9_EFF
            DO J=1,NXMAXB9_EFF
              IMAGEN(J,I,NEWBUFF1)=IMAGEN(J+NREMOV_XMINB9,I+NREMOV_YMINB9,NEWBUFF1)
              IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=IMAGEN(J+NREMOV_XMINB9,I+NREMOV_YMINB9,NEWBUFF1+NMAXBUFF/2)
              IMAGEN(J,I,NEWBUFF2)=IMAGEN(J+NREMOV_XMINB9,I+NREMOV_YMINB9,NEWBUFF2)
            END DO
          END DO
        END IF
        NAXIS(1,NEWBUFF1)=NXMAXB9_EFF
        NAXIS(2,NEWBUFF1)=NYMAXB9_EFF
        NAXIS(1,NEWBUFF1+NMAXBUFF/2)=NXMAXB9_EFF
        NAXIS(2,NEWBUFF1+NMAXBUFF/2)=NYMAXB9_EFF
        NAXIS(1,NEWBUFF2)=NXMAXB9_EFF
        NAXIS(2,NEWBUFF2)=NYMAXB9_EFF
!------------------------------------------------------------------------------
        NX1=1
        NX2=NAXIS(1,NEWBUFF1)
        NY1=1
        NY2=NAXIS(2,NEWBUFF1)
        CALL STATISTIC(NEWBUFF1,NX1,NX2,NY1,NY2,.FALSE.,.TRUE.,.TRUE.,0.0,.FALSE.)
        IF(FSIGMA.GT.0.0)THEN
          BG=FMEAN-5.*FSIGMA
          FG=FMEAN+5.*FSIGMA
        ELSE
          BG=FMEAN-1.0
          FG=FMEAN+1.0
        END IF
        CALL HISTOGRAM(NEWBUFF1)
        CALL SUBLOOK(.FALSE.,NEWBUFF1,.FALSE.)
!
        WRITE(*,100) 'Press <CR> to return to initial buffer...'
        READ(*,*)
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
