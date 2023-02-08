! Version 6-November-2000
!
! Escribe en ficheros las imagenes individuales del box-9 desplazadas en
! funcion de los offsets introducidos a traves de un fichero. De esta forma
! las imagenes individuales pueden sumarse posteriormente.
        SUBROUTINE SHIFTB9(NCBUFF)
        USE Dynamic_Array_IMAGEN
        USE Dynamic_Array_FPIXUSED
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
        INCLUDE 'interface_fpixused.inc'
! subroutine argument
        INTEGER NCBUFF
!
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
!
        INCLUDE 'dimensions.inc'
!
        INTEGER READILIM
        INTEGER TRUEBEG,TRUELEN
        CHARACTER*255 READC
        LOGICAL INSIDE
!
        INTEGER I,J,K,II,JJ,KK,II1,II2,JJ1,JJ2,IJ
        INTEGER L1,L2,KFRAME
        INTEGER L,IDUM
        INTEGER DI(9),DJ(9)
        INTEGER NEWBUFF1
        INTEGER NFRAMES(NMAXBUFF)
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER NAXISFRAME(2,9,NMAXBUFF)
        INTEGER NSIZEFB9(NMAXBUFF)
        INTEGER NF,NFRAMES_
        INTEGER NEXTRA_XMINB9,NEXTRA_XMAXB9
        INTEGER NEXTRA_YMINB9,NEXTRA_YMAXB9
        INTEGER J0(NBOXMAX),I0(NBOXMAX)
        INTEGER NXMAXB9_,NYMAXB9_
        INTEGER NPIXREGION
        INTEGER NX1,NX2,NY1,NY2
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
!delete REAL FPIXUSED(NXMAXB9,NYMAXB9)
        REAL FOFFSETX(NBOXMAX),FOFFSETY(NBOXMAX)
        REAL FOFFSETX_,FOFFSETY_
        REAL EFOFFSETX_,EFOFFSETY_
        REAL FJ0(NBOXMAX),FI0(NBOXMAX)
        REAL FRACCION_PIXEL
        REAL PIXEL(9*4),EPIXEL(9*4),PIXELF(9*4)
        REAL FMEAN,FSIGMA,FMEDIAN,FMIN,FMAX
        REAL BG,FG
        CHARACTER*255 CDUMMY
        CHARACTER*255 OUTFILE,ERRFILE,OUTFILEBASE
        CHARACTER*255 COFFSETFILE
        LOGICAL LOGFILE
!
!delete COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKNFRAMES/NFRAMES
        COMMON/BLKXYLIMPLOT/NX1,NX2,NY1,NY2
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKNAXISFRAME/NAXISFRAME
        COMMON/BLKNSIZEFB9/NSIZEFB9
        COMMON/BLKESTADISTICA/FMEAN,FSIGMA,FMEDIAN,FMIN,FMAX
!------------------------------------------------------------------------------
        CALL Initialize_Dynamic_Array_FPIXUSED
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
        II2=0 !evita WARNING de compilacion
        JJ2=0 !evita WARNING de compilacion
        FRACCION_PIXEL=0.0 !evita WARNING de compilacion
!------------------------------------------------------------------------------
! pedimos offsets
        LOGFILE=.FALSE.
        DO WHILE(.NOT.LOGFILE)
          COFFSETFILE=READC('Name of the file with offsets (none=exit)','offsets.dat','@')
          L1=TRUEBEG(COFFSETFILE)
          L2=TRUELEN(COFFSETFILE)
          IF(COFFSETFILE(L1:L2).EQ.'none')THEN
            CALL Deallocate_Array_FPIXUSED
            RETURN !permitimos escapar
          END IF
          INQUIRE(FILE=COFFSETFILE(L1:L2),EXIST=LOGFILE)
          IF(.NOT.LOGFILE)THEN
            WRITE(*,101) 'ERROR: this file does not exist. Try again.'
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END DO
        WRITE(*,101) 'Reading offsets...'
        OPEN(20,FILE=COFFSETFILE(L1:L2),STATUS='OLD',FORM='FORMATTED')
        NF=0
5       READ(20,*,END=6) IDUM,FOFFSETX_,EFOFFSETX_,FOFFSETY_,EFOFFSETY_
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
          CALL Deallocate_Array_FPIXUSED
          RETURN
        END IF
        GOTO 5
6       CLOSE(20)
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
          CALL Deallocate_Array_FPIXUSED
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
            CALL Deallocate_Array_FPIXUSED
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
              NEXTRA_XMINB9=MIN0(NEXTRA_XMINB9,-J0(NF))
            ELSE
              NEXTRA_XMINB9=MIN0(NEXTRA_XMINB9,-J0(NF)-1)
            END IF
          ELSEIF(FOFFSETX(NF).LE.0.0)THEN
            IF(FJ0(NF).EQ.0.0)THEN
              NEXTRA_XMAXB9=MAX0(NEXTRA_XMAXB9,-J0(NF))
            ELSE
              NEXTRA_XMAXB9=MAX0(NEXTRA_XMAXB9,-J0(NF)+1)
            END IF
          END IF
          I0(NF)=INT(FOFFSETY(NF))
          FI0(NF)=FOFFSETY(NF)-REAL(I0(NF))
          IF(FOFFSETY(NF).GT.0.0)THEN
            IF(FI0(NF).EQ.0.0)THEN
              NEXTRA_YMINB9=MIN0(NEXTRA_YMINB9,-I0(NF))
            ELSE
              NEXTRA_YMINB9=MIN0(NEXTRA_YMINB9,-I0(NF)-1)
            END IF
          ELSEIF(FOFFSETY(NF).LE.0.0)THEN
            IF(FI0(NF).EQ.0.0)THEN
              NEXTRA_YMAXB9=MAX0(NEXTRA_YMAXB9,-I0(NF))
            ELSE
              NEXTRA_YMAXB9=MAX0(NEXTRA_YMAXB9,-I0(NF)+1)
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
          CALL Deallocate_Array_FPIXUSED
          RETURN
        END IF
        IF(NYMAXB9_.GT.NYMAXB9)THEN
          WRITE(*,101) 'ERROR: Y dimension > NYMAXB9'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          CALL Deallocate_Array_FPIXUSED
          RETURN
        END IF
!------------------------------------------------------------------------------
! definimos el buffer de trabajo
        NEWBUFF1=NCBUFF+1
        IF(NEWBUFF1.GT.NMAXBUFF/2) NEWBUFF1=1
        WRITE(CDUMMY,*) NEWBUFF1
        NEWBUFF1=READILIM('Auxiliary buffer # to store manipulated data',CDUMMY,1,NMAXBUFF/2)
!------------------------------------------------------------------------------
! indicamos el nombre base de los ficheros de salida
        LOGFILE=.TRUE.
        DO WHILE(LOGFILE)
          OUTFILEBASE=READC('Base name for output files','@','@')
          LOGFILE=.FALSE.
          DO K=1,NFRAMES_
            IF(.NOT.LOGFILE)THEN
              L1=TRUEBEG(OUTFILEBASE)
              L2=TRUELEN(OUTFILEBASE)
              WRITE(CDUMMY,'(I1)') K
              OUTFILE=OUTFILEBASE(L1:L2)//'_'//CDUMMY(1:1)//'.fits'
              INQUIRE(FILE=OUTFILE,EXIST=LOGFILE)
              IF(LOGFILE)THEN
                WRITE(*,100) 'ERROR: '
                WRITE(*,101) OUTFILE(TRUEBEG(OUTFILE):TRUELEN(OUTFILE))
              END IF
            END IF
          END DO
          IF(LOGFILE)THEN
            WRITE(*,101) 'The previous files already exist. Try again.'
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END DO
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
        NAXIS(1,NEWBUFF1)=NXMAXB9_
        NAXIS(2,NEWBUFF1)=NYMAXB9_
        NFRAMES(NEWBUFF1)=1
        NAXIS(1,NEWBUFF1+NMAXBUFF/2)=NXMAXB9_
        NAXIS(2,NEWBUFF1+NMAXBUFF/2)=NYMAXB9_
        NFRAMES(NEWBUFF1+NMAXBUFF/2)=1
!------------------------------------------------------------------------------
        DO KFRAME=1,NFRAMES_
          DO I=1,NAXIS(2,NEWBUFF1)
            DO J=1,NAXIS(1,NEWBUFF1)
              IMAGEN(J,I,NEWBUFF1)=0. !inicializamos a cero el nuevo frame
              FPIXUSED(J,I)=0. !aqui almacenamos el numero de pixels usados en
                               !el calculo de cada nuevo pixel en la imagen final
              IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=0. !imagen de errores
            END DO
          END DO
!..............................................................................
! recorremos cada pixel en la nueva imagen y determinamos cuantos pixels del
! frame considerado se pueden usar
          DO I=1,NAXIS(2,NEWBUFF1)
            DO J=1,NAXIS(1,NEWBUFF1)
              KK=0 !numero total de pixels (o fracciones de pixels) a utilizar
              JJ1=J+NEXTRA_XMINB9+J0(KFRAME) !ojo: NEXTRA_XMINB9 es <=0
              IF(FJ0(KFRAME).EQ.0.0)THEN
                JJ2=0
              ELSEIF(FJ0(KFRAME).GT.0.0)THEN
                JJ2=JJ1+1
              ELSEIF(FJ0(KFRAME).LT.0.0)THEN
                JJ2=JJ1-1
              END IF
              II1=I+NEXTRA_YMINB9+I0(KFRAME) !ojo: NEXTRA_YMINB9 es <=0
              IF(FI0(KFRAME).EQ.0.0)THEN
                II2=0
              ELSEIF(FI0(KFRAME).GT.0.0)THEN
                II2=II1+1
              ELSEIF(FI0(KFRAME).LT.0.0)THEN
                II2=II1-1
              END IF
              DO IJ=1,4 !recorremos los 4 posibles pixels
                IF(IJ.EQ.1)THEN
                  JJ=JJ1
                  II=II1
                  FRACCION_PIXEL=(1.-ABS(FJ0(KFRAME)))*(1.-ABS(FI0(KFRAME)))
                ELSEIF(IJ.EQ.2)THEN
                  JJ=JJ1
                  II=II2
                  FRACCION_PIXEL=(1.-ABS(FJ0(KFRAME)))*ABS(FI0(KFRAME))
                ELSEIF(IJ.EQ.3)THEN
                  JJ=JJ2
                  II=II1
                  FRACCION_PIXEL=ABS(FJ0(KFRAME))*(1.-ABS(FI0(KFRAME)))
                ELSEIF(IJ.EQ.4)THEN
                  JJ=JJ2
                  II=II2
                  FRACCION_PIXEL=ABS(FJ0(KFRAME))*ABS(FI0(KFRAME))
                END IF
                IF(INSIDE(JJ,1,NAXISFRAME(1,KFRAME,NCBUFF),II,1,NAXISFRAME(2,KFRAME,NCBUFF)))THEN
                  IF(IMAGEN(JJ+DJ(KFRAME)*NSIZEFB9(NCBUFF),II+DI(KFRAME)*NSIZEFB9(NCBUFF),NCBUFF+NMAXBUFF/2).GE.0.0)THEN
                    KK=KK+1
                    PIXEL(KK)=IMAGEN(JJ+DJ(KFRAME)*NSIZEFB9(NCBUFF),II+DI(KFRAME)*NSIZEFB9(NCBUFF),NCBUFF)
                    EPIXEL(KK)=IMAGEN(JJ+DJ(KFRAME)*NSIZEFB9(NCBUFF),II+DI(KFRAME)*NSIZEFB9(NCBUFF),NCBUFF+NMAXBUFF/2)
                    PIXELF(KK)=FRACCION_PIXEL
                  END IF
                END IF
              END DO
!..............................................................................
! normalizamos cada pixel por el numero total de frames que hemos utilizado
! para calcular la sen~al total en dicho pixel
              IF(KK.EQ.0)THEN !lo dejamos por claridad
!!!                IMAGEN(J,I,NEWBUFF1)=0.
!!!                IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=0.
!!!                FPIXUSED(J,I)=0.
              ELSE
                DO K=1,KK
                  IMAGEN(J,I,NEWBUFF1)=IMAGEN(J,I,NEWBUFF1)+PIXEL(K)*PIXELF(K)
                  IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)+EPIXEL(K)*EPIXEL(K)*PIXELF(K)*PIXELF(K)
                  FPIXUSED(J,I)=FPIXUSED(J,I)+PIXELF(K)
                END DO
                IMAGEN(J,I,NEWBUFF1)=IMAGEN(J,I,NEWBUFF1)/FPIXUSED(J,I)
                IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=SQRT(IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2))/FPIXUSED(J,I)
              END IF
            END DO
          END DO
! en aquellos pixels en los que no hemos sumado nada, introducimos un error 
! negativo
          DO I=1,NYMAXB9_
            DO J=1,NXMAXB9_
              IF(FPIXUSED(J,I).LE.0.0)THEN
                IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=-1.
              END IF
            END DO
          END DO
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
!..............................................................................
          L1=TRUEBEG(OUTFILEBASE)
          L2=TRUELEN(OUTFILEBASE)
          WRITE(CDUMMY,'(I1)') KFRAME
          OUTFILE=OUTFILEBASE(L1:L2)//'_'//CDUMMY(1:1)//'.fits'
          L1=TRUEBEG(OUTFILE)
          L2=TRUELEN(OUTFILE)
          WRITE(*,100) 'Creating file: '
          WRITE(*,100) OUTFILE(L1:L2)
          WRITE(*,100) '...'
          CALL SESCRFITS(OUTFILE(L1:L2),NEWBUFF1,'none')
          CALL GUESSEF(OUTFILE(L1:L2),ERRFILE)
          L1=TRUEBEG(ERRFILE)
          L2=TRUELEN(ERRFILE)
          CALL SESCRFITS(ERRFILE(L1:L2),NEWBUFF1+NMAXBUFF/2,'none')
          WRITE(*,101) 'OK!'
        END DO
!
        WRITE(*,100) 'Press <CR> to return to initial buffer...'
        READ(*,*)
!------------------------------------------------------------------------------
        CALL Deallocate_Array_FPIXUSED
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
