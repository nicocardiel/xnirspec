        SUBROUTINE MIDEOFFSET(NFRAMES_,NEWBUFF1,NEWBUFF2,NEWBUFF3,NEWBUFF4,FOFFSETX,FOFFSETY)
        USE Dynamic_Array_IMAGEN
        USE Dynamic_Array_IMAGEN_
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
        INCLUDE 'interface_imagen_.inc'
! subroutine arguments
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
        INTEGER NOBJFITMAX
        PARAMETER (NOBJFITMAX=20)
        INTEGER NSIMULMAX
        PARAMETER (NSIMULMAX=1000)
        REAL PI
        PARAMETER (PI=3.141592654)
!
        INTEGER NFRAMES_
        INTEGER NEWBUFF1,NEWBUFF2,NEWBUFF3,NEWBUFF4
        REAL FOFFSETX(NBOXMAX),FOFFSETY(NBOXMAX)
        REAL RANDOMNUMBER
!
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
!
        INTEGER READI,READILIM
        INTEGER TRUEBEG,TRUELEN
        REAL READF
        REAL FMEAN0,FMEAN0W
        CHARACTER*255 READC
        LOGICAL LPOST
        LOGICAL INSIDE
!
        INTEGER I,J,K,NO
        INTEGER L1,L2
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER NSIZEFB9(NMAXBUFF)
        INTEGER NAXISFRAME(2,NBOXMAX,NMAXBUFF)
        INTEGER IX1,IX2,IY1,IY2
        INTEGER DII,DJJ
        INTEGER IDOLD,IDNEW,NPOST
        INTEGER DI(NBOXMAX),DJ(NBOXMAX)
        INTEGER NX1,NX2,NY1,NY2
        INTEGER NOBJFIT,NOBJEFF
        INTEGER NBOX_CENTROID,NS,NSIMUL,IMODEG2D
        INTEGER NPLOT,NPLOT2
        INTEGER NSEED
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
!delete REAL IMAGEN_(NXMAX,NYMAX)
        REAL BG,FG
        REAL FMEAN,FSIGMA,FMEDIAN,FMIN,FMAX
        REAL XC,YC,TR(6)
        REAL XC_,YC_
        REAL FACTOR,FACTOR1,FACTOR2,FACTOR3
        REAL X0(NOBJFITMAX,NBOXMAX),Y0(NOBJFITMAX,NBOXMAX)
        REAL EX0(NOBJFITMAX,NBOXMAX),EY0(NOBJFITMAX,NBOXMAX)
        REAL SIGMAX(NOBJFITMAX,NBOXMAX),SIGMAY(NOBJFITMAX,NBOXMAX)
        REAL BETA(NOBJFITMAX,NBOXMAX),AMP(NOBJFITMAX,NBOXMAX)
        REAL CTE(NOBJFITMAX,NBOXMAX)
        REAL X0_,Y0_,SIGMAX_,SIGMAY_,BETA_,AMP_,CTE_
        REAL EX0_,EY0_,ESIGMAX_,ESIGMAY_,EBETA_,EAMP_,ECTE_
        REAL EEX0_,EEY0_,EESIGMAX_,EESIGMAY_,EEBETA_,EEAMP_,EECTE_
        REAL X0__,Y0__,SIGMAX__,SIGMAY__,AMPX__,AMPY__,CTEX__,CTEY__
        REAL AREA1,AREA2,AREA3,EAREA1,EAREA2,EAREA3
        REAL FOFFSETX_(NOBJFITMAX,9),FOFFSETY_(NOBJFITMAX,9)
        REAL EFOFFSETX_(NOBJFITMAX,9),EFOFFSETY_(NOBJFITMAX,9)
        REAL FOFFSETX_M(9),FOFFSETY_M(9) !mean values
        REAL FOFFSETX_S(9),FOFFSETY_S(9) !standard deviation
        REAL FOFFSETX_WM(9),FOFFSETY_WM(9) !mean values
        REAL FOFFSETX_WS(9),FOFFSETY_WS(9) !standard deviation
        REAL YRMSTOL
        REAL XV1,XV2,YV1,YV2,XW1,XW2,YW1,YW2
        REAL OLD_CH
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XMINOFF,XMAXOFF,YMINOFF,YMAXOFF,DYOFF
        REAL XDATA(NOBJFITMAX),EXDATA(NOBJFITMAX)
        REAL YDATA(NOBJFITMAX),EYDATA(NOBJFITMAX)
        REAL UDATA(NOBJFITMAX),EUDATA(NOBJFITMAX)
        REAL VDATA(NOBJFITMAX),EVDATA(NOBJFITMAX)
        REAL XDATA_(NOBJFITMAX),YDATA_(NOBJFITMAX)
        REAL UDATA_(NOBJFITMAX),VDATA_(NOBJFITMAX)
        REAL ERRORX,ERRORY,ERRORU,ERRORV,R1,R2
        REAL A11,A12,A21,A22,A31,A32
        REAL A11_(NSIMULMAX),A12_(NSIMULMAX),A21_(NSIMULMAX)
        REAL A22_(NSIMULMAX),A31_(NSIMULMAX),A32_(NSIMULMAX)
        REAL A11_M,A12_M,A21_M,A22_M,A31_M,A32_M
        REAL A11_S,A12_S,A21_S,A22_S,A31_S,A32_S
        REAL CHISQR,CHISQR_(NSIMULMAX),CHISQR_M,CHISQR_S
        CHARACTER*1 CPOST,CH,CCONT,CSAVE,CMODE,CFITOK
        CHARACTER*1 CCHANGE,CWEIGHT
        CHARACTER*5 CAXISX(9),CAXISY(9)
        CHARACTER*50 CBASEPOST,CDUMMY
        CHARACTER*255 OUTFILE
        CHARACTER*255 C255
        LOGICAL LOOP
        LOGICAL LOGFILE
        LOGICAL LOK
        LOGICAL LOUTSIDE
        LOGICAL LFITOK(NOBJFITMAX,NBOXMAX)
!
!delete COMMON/BLKIMAGEN1/IMAGEN
!delete COMMON/BLKIMAGEN1_/IMAGEN_
        COMMON/BLKNAXIS/NAXIS                                      !dimensiones
        COMMON/BLKNAXISFRAME/NAXISFRAME
        COMMON/BLKNSIZEFB9/NSIZEFB9
        COMMON/BLKBGFG/BG,FG
        COMMON/BLKB9POST1/CPOST,CBASEPOST
        COMMON/BLKXYLIMPLOT/NX1,NX2,NY1,NY2
        COMMON/BLKB9POST2/NPOST,IDOLD
        COMMON/BLKESTADISTICA/FMEAN,FSIGMA,FMEDIAN,FMIN,FMAX
        COMMON/BLKDEFAULTS8/NBOX_CENTROID
        COMMON/BLKDEFAULTS9/YRMSTOL
        COMMON/BLKDEFAULTS10/NSIMUL
        COMMON/BLKDEFAULTS11/IMODEG2D
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
        DATA (CAXISX(K),CAXISY(K),K=1,NBOXMAX) / &
         'BCTS ','BCTS ', & !1
         'BCTSN','BCTS ', & !2
         'BCTS ','BCTSN', & !3
         'BCTS ','BCTS ', & !4
         'BCTSN','BCTS ', & !5
         'BCTS ','BCTSN', & !6
         'BCTS ','BCTS ', & !7
         'BCTSN','BCTSN', & !8
         'BCTS ','BCTS '/ !9
!------------------------------------------------------------------------------
        TR(1)=0.
        TR(2)=1.
        TR(3)=0.
        TR(4)=0.
        TR(5)=0.
        TR(6)=1.
        NSEED=-1
!------------------------------------------------------------------------------
! inicializamos buffer con posibles ajustes aceptados
        DO I=1,NYMAX
          DO J=1,NXMAX
            IMAGEN(J,I,NEWBUFF4)=0.0
          END DO
        END DO
        NAXIS(1,NEWBUFF4)=0
        NAXIS(2,NEWBUFF4)=0
!------------------------------------------------------------------------------
        LOOP=.TRUE.
        NOBJFIT=0
        DO WHILE(LOOP)
!..............................................................................
! mostramos parametros para los ajustes
          WRITE(*,*)
          WRITE(*,100) 'Box size (pixels)..........: '
          WRITE(*,*) NBOX_CENTROID
          WRITE(*,100) 'YRMSTOL for DOWNHILL.......: '
          WRITE(*,*) YRMSTOL
          WRITE(*,100) 'NSIMUL.....................: '
          WRITE(*,*) NSIMUL
          WRITE(*,100) '2D Gaussian fit mode.......: '
          WRITE(*,*) IMODEG2D
!..............................................................................
          LOK=.FALSE.
          DO WHILE(.NOT.LOK)
            CALL STATISTICB9(NFRAMES_,NEWBUFF1,FMEAN,FSIGMA)
            IF(FSIGMA.GT.0.0)THEN
              BG=FMEAN-5.*FSIGMA
              FG=FMEAN+5.*FSIGMA
            ELSE
              BG=FMEAN-1.0
              FG=FMEAN+1.0
            END IF
            CALL HISTOGRAM(NEWBUFF1)
! hacemos un zoom alrededor del frame central
            NX1=DJ(1)*NSIZEFB9(NEWBUFF1)+1-NBOX_CENTROID
            NX2=DJ(1)*NSIZEFB9(NEWBUFF1)+NSIZEFB9(NEWBUFF1)-1+NBOX_CENTROID
            NY1=DI(1)*NSIZEFB9(NEWBUFF1)+1-NBOX_CENTROID
            NY2=DI(1)*NSIZEFB9(NEWBUFF1)+NSIZEFB9(NEWBUFF1)-1+NBOX_CENTROID
            CALL SUBLOOK(.FALSE.,NEWBUFF1,.FALSE.)
            CALL DRAWSEPB9(0,0,NSIZEFB9(NEWBUFF1))
            DO K=1,NFRAMES_
              IF(NOBJFIT.GT.0)THEN
                DO NO=1,NOBJFIT
                  CALL PGSCI(MOD(NO-1,5)+2)
                  IF(LFITOK(NO,K)) CALL PGPOINT(1,X0(NO,K),Y0(NO,K),2)
                END DO
              END IF
            END DO
            CALL PGSCI(1)
            IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
              CALL SUBLOOK(.FALSE.,NEWBUFF1,.TRUE.)
              CALL DRAWSEPB9(0,0,NSIZEFB9(NEWBUFF1))
              DO K=1,NFRAMES_
                IF(NOBJFIT.GT.0)THEN
                  DO NO=1,NOBJFIT
                    CALL PGSCI(MOD(NO-1,5)+2)
                    IF(LFITOK(NO,K)) CALL PGPOINT(1,X0(NO,K),Y0(NO,K),2)
                  END DO
                END IF
              END DO
              CALL PGSCI(1)
              CALL PGCLOS(IDNEW)
              CALL PGSLCT(IDOLD)
            END IF
! seleccionamos punto inicial para ajustar un objeto del frame #1
            LOUTSIDE=.TRUE.
            DO WHILE(LOUTSIDE)
              WRITE(*,*)
              WRITE(*,101) ' <left button>: select object'
              WRITE(*,101) '<right button>: change fitting options/exit'
              WRITE(*,100) 'Select object from frame #1...'
              CALL PGSCI(7)
              CALL RPGBAND(7,0,0.,0.,XC,YC,CH)
              CALL PGSCI(1)
              WRITE(*,101) '  ...OK!'
              IF(CH.EQ.'X')THEN
                LOUTSIDE=.FALSE.
              ELSE
                LOUTSIDE=.NOT.INSIDE( &
                 NINT(XC), &
                 1+DJ(1)*NSIZEFB9(NEWBUFF1)+NBOX_CENTROID/2+1, &
                 NAXISFRAME(1,1,NEWBUFF1)+ &
                  DJ(1)*NSIZEFB9(NEWBUFF1)-NBOX_CENTROID/2-1, &
                 NINT(YC), &
                 1+DI(1)*NSIZEFB9(NEWBUFF1)+NBOX_CENTROID/2+1, &
                 NAXISFRAME(2,1,NEWBUFF1)+ &
                  DI(1)*NSIZEFB9(NEWBUFF1)-NBOX_CENTROID/2-1)
              END IF
              IF(LOUTSIDE)THEN
                WRITE(*,100) 'ERROR: cursor out of frame #1. Try again.'
              END IF
            END DO
            IF(CH.EQ.'X')THEN
              C255=READC('[c]hange fitting parameters or e[x]it (c/x)','c','cx')
              CCHANGE(1:1)=C255(1:1)
              IF(CCHANGE.EQ.'x')THEN
                RETURN
              END IF
              WRITE(CDUMMY,*) NBOX_CENTROID
              NBOX_CENTROID=READILIM('Box size in pixels (odd)',CDUMMY,5,MIN(NAXIS(1,NEWBUFF1),NAXIS(2,NEWBUFF1)))
              IF(MOD(NBOX_CENTROID,2).EQ.0)THEN
                NBOX_CENTROID=NBOX_CENTROID+1
                IF(NBOX_CENTROID.GT.MIN(NAXIS(1,NEWBUFF1),NAXIS(2,NEWBUFF1))) NBOX_CENTROID=NBOX_CENTROID-2
              END IF
              WRITE(CDUMMY,*) YRMSTOL
              YRMSTOL=READF('YRMSTOL for DOWNHILL',CDUMMY)
              WRITE(CDUMMY,*) NSIMUL
              NSIMUL=READI('Number of simulations',CDUMMY)
              WRITE(CDUMMY,*) IMODEG2D
              IMODEG2D=READI('2D Gaussian fit mode (0,1)',CDUMMY)
            ELSE
! sumamos todas las regiones disponibles para hacer un ajuste inicial que
! nos sirva para luego hacer los ajustes individuales de forma mas sencilla
! (nota: aqui hay que tener cuidado porque en los bordes del frame#1, hay
! zonas que no estan presentes en todas las imagenes del box-9)
              NAXIS(1,NEWBUFF2)=NBOX_CENTROID
              NAXIS(2,NEWBUFF2)=NBOX_CENTROID
              DO I=1,NAXIS(2,NEWBUFF2)
                DO J=1,NAXIS(1,NEWBUFF2)
                  IMAGEN(J,I,NEWBUFF2)=0.
                  IMAGEN_(J,I)=0. !almacenamos aqui numero de pixels sumados
                END DO
              END DO
              DO K=1,NFRAMES_
                IX1=NINT(XC+FOFFSETX(K))-NBOX_CENTROID/2+DJ(K)*NSIZEFB9(NEWBUFF1)-DJ(1)*NSIZEFB9(NEWBUFF1)
                IX2=IX1+NBOX_CENTROID-1
                IY1=NINT(YC+FOFFSETY(K))-NBOX_CENTROID/2+DI(K)*NSIZEFB9(NEWBUFF1)-DI(1)*NSIZEFB9(NEWBUFF1)
                IY2=IY1+NBOX_CENTROID-1
                CALL PGSCI(4)
                CALL PGMOVE(REAL(IX1),REAL(IY1))
                CALL PGDRAW(REAL(IX2),REAL(IY1))
                CALL PGDRAW(REAL(IX2),REAL(IY2))
                CALL PGDRAW(REAL(IX1),REAL(IY2))
                CALL PGDRAW(REAL(IX1),REAL(IY1))
                CALL PGSCI(1)
                DO I=IY1,IY2
                  DO J=IX1,IX2
                    IF(INSIDE(J,1+DJ(K)*NSIZEFB9(NEWBUFF1), &
                     NAXISFRAME(1,K,NEWBUFF1)+DJ(K)*NSIZEFB9(NEWBUFF1), &
                     I,1+DI(K)*NSIZEFB9(NEWBUFF1), &
                     NAXISFRAME(2,K,NEWBUFF1)+ &
                     DI(K)*NSIZEFB9(NEWBUFF1)))THEN
                      IMAGEN(J-IX1+1,I-IY1+1,NEWBUFF2)=IMAGEN(J-IX1+1,I-IY1+1,NEWBUFF2)+IMAGEN(J,I,NEWBUFF1)
                      IMAGEN_(J-IX1+1,I-IY1+1)=IMAGEN_(J-IX1+1,I-IY1+1)+1.0
                    END IF
                  END DO
                END DO
              END DO
              DO I=1,NBOX_CENTROID
                DO J=1,NBOX_CENTROID
                  IF(IMAGEN_(J,I).GT.0.0)THEN
                    IMAGEN(J,I,NEWBUFF2)=IMAGEN(J,I,NEWBUFF2)/IMAGEN_(J,I)
                  END IF
                END DO
              END DO
              CALL STATISTIC(NEWBUFF2,1,NBOX_CENTROID,1,NBOX_CENTROID,.FALSE.,.FALSE.,.TRUE.,0.0,.FALSE.)
              IF(FSIGMA.GT.0.0)THEN
                BG=FMEAN-5.*FSIGMA
                FG=FMEAN+5.*FSIGMA
              ELSE
                BG=FMEAN-1.0
                FG=FMEAN+1.0
              END IF
!
              CALL PGQVP(0,XV1,XV2,YV1,YV2)
              CALL PGQWIN(XW1,XW2,YW1,YW2)
              CALL PGQCH(OLD_CH)
              CALL PGSCH(0.7)
              CALL RPGERASW(0.00,0.40,0.39,0.80,0)
              CALL PGSVP(0.01,0.19,0.40,0.79)
              CALL PGWNAD(REAL(1)-0.6,REAL(NBOX_CENTROID)+0.6,REAL(1)-0.6,REAL(NBOX_CENTROID)+0.6)
              CALL PGIMAG(IMAGEN(1,1,NEWBUFF2),NXMAX,NYMAX,1,NBOX_CENTROID,1,NBOX_CENTROID,FG,BG,TR)
              CALL PGBOX('BCTS',0.0,0,'BCTS',0.0,0)
              CALL PGSCI(7)
              CALL PGMTXT('T',0.5,0.5,0.5,'averaged frame')
              CALL PGSCI(1)
!
              CALL GAUSS2DFIT(NEWBUFF2,1,NBOX_CENTROID,1,NBOX_CENTROID, &
               .FALSE., &
               X0_,Y0_,SIGMAX_,SIGMAY_,BETA_,AMP_,CTE_, &
               EX0_,EY0_,ESIGMAX_,ESIGMAY_,EBETA_,EAMP_,ECTE_, &
               EEX0_,EEY0_,EESIGMAX_,EESIGMAY_,EEBETA_,EEAMP_,EECTE_, &
               X0__,Y0__,SIGMAX__,SIGMAY__,AMPX__,AMPY__,CTEX__,CTEY__, &
               AREA1,AREA2,AREA3,EAREA1,EAREA2,EAREA3, &
               YRMSTOL,0,IMODEG2D,.FALSE.)
!
              DO I=1,NBOX_CENTROID
                DO J=1,NBOX_CENTROID
                  FACTOR1=(REAL(J)-X0_)*(REAL(J)-X0_)/(2.*SIGMAX_*SIGMAX_)
                  FACTOR2=(REAL(I)-Y0_)*(REAL(I)-Y0_)/(2.*SIGMAY_*SIGMAY_)
                  FACTOR3=BETA_*(REAL(J)-X0_)*(REAL(I)-Y0_)/(SIGMAX_*SIGMAY_)
                  FACTOR=FACTOR1+FACTOR2+FACTOR3
                  IF(FACTOR.GT.60.)THEN
                    IMAGEN_(J,I)=CTE_
                  ELSE
                    IMAGEN_(J,I)=CTE_+AMP_*EXP(-FACTOR)
                  END IF
                END DO
              END DO
              CALL PGSVP(0.21,0.39,0.40,0.79)
              CALL PGWNAD(REAL(1)-0.6,REAL(NBOX_CENTROID)+0.6,REAL(1)-0.6,REAL(NBOX_CENTROID)+0.6)
              CALL PGIMAG(IMAGEN_,NXYMAX,NXYMAX,1,NBOX_CENTROID,1,NBOX_CENTROID,FG,BG,TR)
              CALL PGBOX('BCTS',0.0,0,'BCTS',0.0,0)
              CALL PGSCI(7)
              CALL PGMTXT('T',0.5,0.5,0.5,'fit to averaged frame')
              CALL PGSCI(1)
              CALL PGSVP(XV1,XV2,YV1,YV2)
              CALL PGSWIN(XW1,XW2,YW1,YW2)
              CALL PGSCH(OLD_CH)
              WRITE(*,100) 'Press <left/right> button to <confirm/repeat>...'
              CALL RPGBAND(0,0,0.,0.,XC_,YC_,CH)
              LOK=(CH.NE.'X')
            END IF
          END DO
!..............................................................................
! posible nuevo objeto
          NOBJFIT=NOBJFIT+1
!..............................................................................
          CALL STATISTICB9(NFRAMES_,NEWBUFF1,FMEAN,FSIGMA)
          IF(FSIGMA.GT.0.0)THEN
            BG=FMEAN-5.*FSIGMA
            FG=FMEAN+5.*FSIGMA
          ELSE
            BG=FMEAN-1.0
            FG=FMEAN+1.0
          END IF
          CALL HISTOGRAM(NEWBUFF1)
! deshacemos un zoom y mostramos de nuevo toda la imagen
          NX1=1
          NX2=3*NSIZEFB9(NEWBUFF1)
          NY1=1
          NY2=3*NSIZEFB9(NEWBUFF1)
          CALL SUBLOOK(.FALSE.,NEWBUFF1,.FALSE.)
          CALL DRAWSEPB9(0,0,NSIZEFB9(NEWBUFF1))
          DO K=1,NFRAMES_
            IF(NOBJFIT-1.GT.0)THEN
              DO NO=1,NOBJFIT-1
                CALL PGSCI(MOD(NO-1,5)+2)
                IF(LFITOK(NO,K)) CALL PGPOINT(1,X0(NO,K),Y0(NO,K),2)
              END DO
            END IF
          END DO
          CALL PGSCI(1)
          DO K=1,NFRAMES_
            IX1=NINT(XC+FOFFSETX(K))-NBOX_CENTROID/2+DJ(K)*NSIZEFB9(NEWBUFF1)-DJ(1)*NSIZEFB9(NEWBUFF1)
            IX2=IX1+NBOX_CENTROID-1
            IY1=NINT(YC+FOFFSETY(K))-NBOX_CENTROID/2+DI(K)*NSIZEFB9(NEWBUFF1)-DI(1)*NSIZEFB9(NEWBUFF1)
            IY2=IY1+NBOX_CENTROID-1
            IF( &
             INSIDE(IX1,1+DJ(K)*NSIZEFB9(NEWBUFF1), &
                    NAXISFRAME(1,K,NEWBUFF1)+DJ(K)*NSIZEFB9(NEWBUFF1), &
                    IY1,1+DI(K)*NSIZEFB9(NEWBUFF1), &
                    NAXISFRAME(2,K,NEWBUFF1)+DI(K)*NSIZEFB9(NEWBUFF1)) &
             .AND. &
             INSIDE(IX2,1+DJ(K)*NSIZEFB9(NEWBUFF1), &
                    NAXISFRAME(1,K,NEWBUFF1)+DJ(K)*NSIZEFB9(NEWBUFF1), &
                    IY1,1+DI(K)*NSIZEFB9(NEWBUFF1), &
                    NAXISFRAME(2,K,NEWBUFF1)+DI(K)*NSIZEFB9(NEWBUFF1)) &
             .AND. &
             INSIDE(IX2,1+DJ(K)*NSIZEFB9(NEWBUFF1), &
                    NAXISFRAME(1,K,NEWBUFF1)+DJ(K)*NSIZEFB9(NEWBUFF1), &
                    IY2,1+DI(K)*NSIZEFB9(NEWBUFF1), &
                    NAXISFRAME(2,K,NEWBUFF1)+DI(K)*NSIZEFB9(NEWBUFF1)) &
             .AND. &
             INSIDE(IX1,1+DJ(K)*NSIZEFB9(NEWBUFF1), &
                    NAXISFRAME(1,K,NEWBUFF1)+DJ(K)*NSIZEFB9(NEWBUFF1), &
                    IY2,1+DI(K)*NSIZEFB9(NEWBUFF1), &
                    NAXISFRAME(2,K,NEWBUFF1)+DI(K)*NSIZEFB9(NEWBUFF1)) &
             )THEN
              CALL PGSCI(4)
              LFITOK(NOBJFIT,K)=.TRUE.
            ELSE
              CALL PGSCI(2)
              LFITOK(NOBJFIT,K)=.FALSE.
            END IF
            CALL PGMOVE(REAL(IX1),REAL(IY1))
            CALL PGDRAW(REAL(IX2),REAL(IY1))
            CALL PGDRAW(REAL(IX2),REAL(IY2))
            CALL PGDRAW(REAL(IX1),REAL(IY2))
            CALL PGDRAW(REAL(IX1),REAL(IY1))
            CALL PGSCI(1)
          END DO
          IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
            CALL SUBLOOK(.FALSE.,NEWBUFF1,.TRUE.)
            CALL DRAWSEPB9(0,0,NSIZEFB9(NEWBUFF1))
            DO K=1,NFRAMES_
              IF(NOBJFIT-1.GT.0)THEN
                DO NO=1,NOBJFIT-1
                  CALL PGSCI(MOD(NO-1,5)+2)
                  IF(LFITOK(NO,K)) CALL PGPOINT(1,X0(NO,K),Y0(NO,K),2)
                END DO
              END IF
            END DO
            CALL PGSCI(1)
            DO K=1,NFRAMES_
              IX1=NINT(XC+FOFFSETX(K))-NBOX_CENTROID/2+DJ(K)*NSIZEFB9(NEWBUFF1)-DJ(1)*NSIZEFB9(NEWBUFF1)
              IX2=IX1+NBOX_CENTROID-1
              IY1=NINT(YC+FOFFSETY(K))-NBOX_CENTROID/2+DI(K)*NSIZEFB9(NEWBUFF1)-DI(1)*NSIZEFB9(NEWBUFF1)
              IY2=IY1+NBOX_CENTROID-1
              IF(LFITOK(NOBJFIT,K))THEN
                CALL PGSCI(4)
              ELSE
                CALL PGSCI(2)
              END IF
              CALL PGMOVE(REAL(IX1),REAL(IY1))
              CALL PGDRAW(REAL(IX2),REAL(IY1))
              CALL PGDRAW(REAL(IX2),REAL(IY2))
              CALL PGDRAW(REAL(IX1),REAL(IY2))
              CALL PGDRAW(REAL(IX1),REAL(IY1))
              CALL PGSCI(1)
            END DO
            CALL PGCLOS(IDNEW)
            CALL PGSLCT(IDOLD)
          END IF
!..............................................................................
! ajustamos centroides
          DO K=1,NFRAMES_
            IF(LFITOK(NOBJFIT,K))THEN
              CALL PGSCI(2)
              CALL PGPOINT(1, &
               XC+FOFFSETX(K)+DJ(K)*NSIZEFB9(NEWBUFF1)-DJ(1)*NSIZEFB9(NEWBUFF1), &
               YC+FOFFSETY(K)+DI(K)*NSIZEFB9(NEWBUFF1)-DI(1)*NSIZEFB9(NEWBUFF1),2)
              CALL PGSCI(1)
              IF(K.EQ.9)THEN
                CMODE='D'
              ELSE
                CMODE='d'
              END IF
              X0(NOBJFIT,K)=X0_-REAL(NBOX_CENTROID/2)+XC+FOFFSETX(K)+DJ(K)*NSIZEFB9(NEWBUFF1)-DJ(1)*NSIZEFB9(NEWBUFF1)
              Y0(NOBJFIT,K)=Y0_-REAL(NBOX_CENTROID/2)+YC+FOFFSETY(K)+DI(K)*NSIZEFB9(NEWBUFF1)-DI(1)*NSIZEFB9(NEWBUFF1)
              SIGMAX(NOBJFIT,K)=SIGMAX_
              SIGMAY(NOBJFIT,K)=SIGMAY_
              BETA(NOBJFIT,K)=BETA_
              AMP(NOBJFIT,K)=AMP_
              CTE(NOBJFIT,K)=CTE_
              CALL CENTROID(NEWBUFF1,CMODE, &
               XC+FOFFSETX(K)+ &
               DJ(K)*NSIZEFB9(NEWBUFF1)-DJ(1)*NSIZEFB9(NEWBUFF1), &
               YC+FOFFSETY(K)+ &
               DI(K)*NSIZEFB9(NEWBUFF1)-DI(1)*NSIZEFB9(NEWBUFF1), &
               .TRUE., &
               X0(NOBJFIT,K),Y0(NOBJFIT,K), &
               SIGMAX(NOBJFIT,K),SIGMAY(NOBJFIT,K), &
               BETA(NOBJFIT,K),AMP(NOBJFIT,K),CTE(NOBJFIT,K), &
               EX0(NOBJFIT,K),EY0(NOBJFIT,K), &
               ESIGMAX_,ESIGMAY_,EBETA_,EAMP_,ECTE_,.FALSE.)
            END IF
          END DO
!..............................................................................
! mostramos ajustes
          NAXIS(1,NEWBUFF2)=3*NBOX_CENTROID
          NAXIS(2,NEWBUFF2)=3*NBOX_CENTROID
          DO I=1,NAXIS(2,NEWBUFF2)
            DO J=1,NAXIS(1,NEWBUFF2)
              IMAGEN(J,I,NEWBUFF2)=0.
            END DO
          END DO
          DO K=1,NFRAMES_
            IF(LFITOK(NOBJFIT,K))THEN
              IX1=NINT(XC+FOFFSETX(K))-NBOX_CENTROID/2+DJ(K)*NSIZEFB9(NEWBUFF1)-DJ(1)*NSIZEFB9(NEWBUFF1)
              IX2=IX1+NBOX_CENTROID-1
              IY1=NINT(YC+FOFFSETY(K))-NBOX_CENTROID/2+DI(K)*NSIZEFB9(NEWBUFF1)-DI(1)*NSIZEFB9(NEWBUFF1)
              IY2=IY1+NBOX_CENTROID-1
              DII=DI(K)*NBOX_CENTROID
              DJJ=DJ(K)*NBOX_CENTROID
              DO I=IY1,IY2
                DO J=IX1,IX2
                  IF(INSIDE(J,1,NAXIS(1,NEWBUFF1),I,1,NAXIS(2,NEWBUFF1)))THEN
                    IMAGEN(J-IX1+1+DJJ,I-IY1+1+DII,NEWBUFF2)=IMAGEN(J,I,NEWBUFF1)
                  END IF
                END DO
              END DO
            END IF
          END DO
          CALL STATISTIC(NEWBUFF2,1,3*NBOX_CENTROID,1,3*NBOX_CENTROID,.FALSE.,.FALSE.,.TRUE.,0.0,.FALSE.)
          IF(FSIGMA.GT.0.0)THEN
            BG=FMEAN-5.*FSIGMA
            FG=FMEAN+5.*FSIGMA
          ELSE
            BG=FMEAN-1.0
            FG=FMEAN+1.0
          END IF
          XMIN=1.0-0.6
          XMAX=REAL(NBOX_CENTROID)*3+0.6
          YMIN=1.0-0.6
          YMAX=REAL(NBOX_CENTROID)*3+0.6
          CALL PGQVP(0,XV1,XV2,YV1,YV2)
          CALL PGQWIN(XW1,XW2,YW1,YW2)
          CALL PGQCH(OLD_CH)
          CALL PGSCH(0.7)
          CALL RPGERASW(0.00,0.40,0.39,0.80,0)
          CALL PGSVP(0.01,0.19,0.40,0.79)
          CALL PGWNAD(XMIN,XMAX,YMIN,YMAX)
          CALL PGIMAG(IMAGEN(1,1,NEWBUFF2),NXMAX,NYMAX,1,NAXIS(1,NEWBUFF2),1,NAXIS(2,NEWBUFF2),FG,BG,TR)
          CALL PGBOX('BCTS',0.0,0,'BCTS',0.0,0)
          CALL PGSCI(7)
          CALL PGMTXT('T',0.5,0.5,0.5,'original frame')
          CALL PGSCI(1)
          NAXIS(1,NEWBUFF3)=3*NBOX_CENTROID
          NAXIS(2,NEWBUFF3)=3*NBOX_CENTROID
          DO I=1,NAXIS(2,NEWBUFF3)
            DO J=1,NAXIS(1,NEWBUFF3)
              IMAGEN(J,I,NEWBUFF3)=0.
            END DO
          END DO
          DO K=1,NFRAMES_
            IF(LFITOK(NOBJFIT,K))THEN
              IX1=NINT(XC+FOFFSETX(K))-NBOX_CENTROID/2+DJ(K)*NSIZEFB9(NEWBUFF1)-DJ(1)*NSIZEFB9(NEWBUFF1)
              IX2=IX1+NBOX_CENTROID-1
              IY1=NINT(YC+FOFFSETY(K))-NBOX_CENTROID/2+DI(K)*NSIZEFB9(NEWBUFF1)-DI(1)*NSIZEFB9(NEWBUFF1)
              IY2=IY1+NBOX_CENTROID-1
              DII=DI(K)*NBOX_CENTROID
              DJJ=DJ(K)*NBOX_CENTROID
              DO I=IY1,IY2
                DO J=IX1,IX2
                  FACTOR1=(REAL(J)-X0(NOBJFIT,K))*(REAL(J)-X0(NOBJFIT,K))/(2.*SIGMAX(NOBJFIT,K)*SIGMAX(NOBJFIT,K))
                  FACTOR2=(REAL(I)-Y0(NOBJFIT,K))*(REAL(I)-Y0(NOBJFIT,K))/(2.*SIGMAY(NOBJFIT,K)*SIGMAY(NOBJFIT,K))
                  FACTOR3=BETA(NOBJFIT,K)*(REAL(J)-X0(NOBJFIT,K))*(REAL(I)-Y0(NOBJFIT,K))/(SIGMAX(NOBJFIT,K)*SIGMAY(NOBJFIT,K))
                  FACTOR=FACTOR1+FACTOR2+FACTOR3
                  IF(FACTOR.GT.60.)THEN
                    IMAGEN(J-IX1+1+DJJ,I-IY1+1+DII,NEWBUFF3)=CTE(NOBJFIT,K)
                  ELSE
                    IMAGEN(J-IX1+1+DJJ,I-IY1+1+DII,NEWBUFF3)=CTE(NOBJFIT,K)+AMP(NOBJFIT,K)*EXP(-FACTOR)
                  END IF
                END DO
              END DO
            END IF
          END DO
          CALL PGSVP(0.21,0.39,0.40,0.79)
          CALL PGWNAD(XMIN,XMAX,YMIN,YMAX)
          CALL PGIMAG(IMAGEN(1,1,NEWBUFF3),NXMAX,NYMAX,1,NAXIS(1,NEWBUFF3),1,NAXIS(2,NEWBUFF3),FG,BG,TR)
          CALL PGBOX('BCTS',0.0,0,'BCTS',0.0,0)
          CALL PGSCI(7)
          CALL PGMTXT('T',0.5,0.5,0.5,'fit to original frame')
          CALL PGSCI(1)
          CALL PGSVP(XV1,XV2,YV1,YV2)
          CALL PGSWIN(XW1,XW2,YW1,YW2)
          CALL PGSCH(OLD_CH)
!..............................................................................
! calculamos offsets medios y dispersion
          DO K=1,NFRAMES_
            IF(LFITOK(NOBJFIT,K))THEN
              FOFFSETX_(NOBJFIT,K)=X0(NOBJFIT,K)-X0(NOBJFIT,1)-DJ(K)*NSIZEFB9(NEWBUFF1)+DJ(1)*NSIZEFB9(NEWBUFF1)-FOFFSETX(K)
              EFOFFSETX_(NOBJFIT,K)=EX0(NOBJFIT,K)
              FOFFSETY_(NOBJFIT,K)=Y0(NOBJFIT,K)-Y0(NOBJFIT,1)-DI(K)*NSIZEFB9(NEWBUFF1)+DI(1)*NSIZEFB9(NEWBUFF1)-FOFFSETY(K)
              EFOFFSETY_(NOBJFIT,K)=EY0(NOBJFIT,K)
            END IF
          END DO
          DO K=1,NFRAMES_
            NOBJEFF=0
            DO NO=1,NOBJFIT
              IF(LFITOK(NO,K))THEN
                NOBJEFF=NOBJEFF+1
                XDATA(NOBJEFF)=FOFFSETX_(NO,K)
                EXDATA(NOBJEFF)=EFOFFSETX_(NO,K)
                YDATA(NOBJEFF)=FOFFSETY_(NO,K)
                EYDATA(NOBJEFF)=EFOFFSETY_(NO,K)
              END IF
            END DO
            IF(NOBJEFF.GT.0)THEN
              FOFFSETX_M(K)=FMEAN0(NOBJEFF,XDATA,FOFFSETX_S(K))
              FOFFSETX_WM(K)=FMEAN0W(NOBJEFF,XDATA,EXDATA,FOFFSETX_WS(K))
              FOFFSETY_M(K)=FMEAN0(NOBJEFF,YDATA,FOFFSETY_S(K))
              FOFFSETY_WM(K)=FMEAN0W(NOBJEFF,YDATA,EYDATA,FOFFSETY_WS(K))
              WRITE(*,100) 'Offsets > '
              WRITE(*,*) K,FOFFSETX_M(K)+FOFFSETX(K),FOFFSETX_S(K),FOFFSETY_M(K)+FOFFSETY(K),FOFFSETY_S(K)
              WRITE(*,100) 'OffsetsW> '
              WRITE(*,*) K,FOFFSETX_WM(K)+FOFFSETX(K),FOFFSETX_WS(K),FOFFSETY_WM(K)+FOFFSETY(K),FOFFSETY_WS(K)
            ELSE
              FOFFSETX_M(K)=0.
              FOFFSETX_S(K)=0.
              FOFFSETX_WM(K)=0.
              FOFFSETX_WS(K)=0.
              WRITE(*,101) 'Offsets > unknown'
              WRITE(*,101) 'OffsetsW> unknown'
            END IF
          END DO
! dibujamos offsets medidos hasta ahora
          CALL PGQVP(0,XV1,XV2,YV1,YV2)
          CALL PGQWIN(XW1,XW2,YW1,YW2)
          CALL PGQCH(OLD_CH)
          IF(CPOST.EQ.'y')THEN
            NPLOT2=2
          ELSE
            NPLOT2=1
          END IF
          DO NPLOT=1,NPLOT2
            IF(NPLOT.EQ.1)THEN !en X11
              CALL PGSCH(0.7)
              CALL RPGERASW(0.40,1.00,0.00,0.80,0)
              YMINOFF=FOFFSETX_(1,1)
              YMAXOFF=YMINOFF
              DO K=1,NFRAMES_
                DO NO=1,NOBJFIT
                  IF(LFITOK(NO,K))THEN
                    IF(YMINOFF.GT.FOFFSETX_(NO,K)-EFOFFSETX_(NO,K)) YMINOFF=FOFFSETX_(NO,K)-EFOFFSETX_(NO,K)
                    IF(YMINOFF.GT.FOFFSETY_(NO,K)-EFOFFSETY_(NO,K)) YMINOFF=FOFFSETY_(NO,K)-EFOFFSETY_(NO,K)
                    IF(YMAXOFF.LT.FOFFSETX_(NO,K)+EFOFFSETX_(NO,K)) YMAXOFF=FOFFSETX_(NO,K)+EFOFFSETX_(NO,K)
                    IF(YMAXOFF.LT.FOFFSETY_(NO,K)+EFOFFSETY_(NO,K)) YMAXOFF=FOFFSETY_(NO,K)+EFOFFSETY_(NO,K)
                  END IF
                END DO
                IF(YMINOFF.GT.FOFFSETX_M(K)-FOFFSETX_S(K)) YMINOFF=FOFFSETX_M(K)-FOFFSETX_S(K)
                IF(YMINOFF.GT.FOFFSETX_WM(K)-FOFFSETX_WS(K)) YMINOFF=FOFFSETX_WM(K)-FOFFSETX_WS(K)
                IF(YMAXOFF.LT.FOFFSETX_M(K)+FOFFSETX_S(K)) YMAXOFF=FOFFSETX_M(K)+FOFFSETX_S(K)
                IF(YMAXOFF.LT.FOFFSETX_WM(K)+FOFFSETX_WS(K)) YMAXOFF=FOFFSETX_WM(K)+FOFFSETX_WS(K)
              END DO
              DYOFF=YMAXOFF-YMINOFF
              YMINOFF=YMINOFF-DYOFF/20.
              YMAXOFF=YMAXOFF+DYOFF/20.
            ELSE !en PostScript
              IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN !abrimos fichero
              END IF !OJO: no borrar aunque el IF-ENDIF este vacio!
            END IF
            DO K=1,NBOXMAX !dibujamos las 9 cajas aunque solo tengamos NFRAMES_
              XMIN=REAL(DJ(K))/3.
              XMAX=XMIN+1./3.
              XMIN=0.48+XMIN*(0.98-0.48)
              XMAX=0.48+XMAX*(0.98-0.48)
              YMIN=REAL(DI(K))/3.
              YMAX=YMIN+1./3.
              YMIN=0.08+YMIN*(0.78-0.08)
              YMAX=0.08+YMAX*(0.78-0.08)
              CALL PGSVP(XMIN,XMAX,YMIN,YMAX)
              XMAXOFF=REAL(NOBJFIT+1)-0.01
              IF(NOBJFIT.EQ.1)THEN
                XMINOFF=0.01
                CALL PGSWIN(XMINOFF,XMAXOFF,YMINOFF,YMAXOFF)
              ELSE
                XMINOFF=-2.99
                CALL PGSWIN(XMINOFF,XMAXOFF,YMINOFF,YMAXOFF)
                CALL PGSLS(2)
                CALL PGSCI(4)
                CALL PGMOVE(0.,YMINOFF)
                CALL PGDRAW(0.,YMAXOFF)
                CALL PGSCI(1)
                CALL PGSLS(1)
              END IF
              CALL PGSLS(4)
              CALL PGMOVE(XMINOFF,0.)
              CALL PGDRAW(XMAXOFF,0.)
              CALL PGSLS(1)
              CALL PGBOX(CAXISX(K),0.0,0,CAXISY(K),0.0,0)
              IF(INDEX(CAXISX(K),'N').NE.0)THEN
                CALL PGSCI(1)
                CALL PGMTXT('B',2.5,0.5,0.5,'object #')
              END IF
              IF(INDEX(CAXISY(K),'N').NE.0)THEN
                CALL PGSCI(3)
                CALL PGMTXT('L',2.5,0.49,1.00,'\\gDx')
                CALL PGSCI(2)
                CALL PGMTXT('L',2.5,0.51,0.00,'\\gDy')
              END IF
              CALL PGSCI(1)
              IF(K.LE.NFRAMES_)THEN
                DO NO=1,NOBJFIT
                  IF(LFITOK(NO,K))THEN
                    CALL PGSCI(3)
                    CALL PGPOINT(1,REAL(NO)-0.2,FOFFSETX_(NO,K),1)
                    CALL PGERRY(1,REAL(NO)-0.2,FOFFSETX_(NO,K)-EFOFFSETX_(NO,K),FOFFSETX_(NO,K)+EFOFFSETX_(NO,K),1.0)
                    CALL PGSCI(2)
                    CALL PGPOINT(1,REAL(NO)+0.2,FOFFSETY_(NO,K),1)
                    CALL PGERRY(1,REAL(NO)+0.2,FOFFSETY_(NO,K)-EFOFFSETY_(NO,K),FOFFSETY_(NO,K)+EFOFFSETY_(NO,K),1.0)
                    CALL PGSCI(1)
                  END IF
                END DO
                IF(NOBJFIT.GT.1)THEN
                  CALL PGSCI(3)
                  CALL PGPOINT(1,-2.0-0.2,FOFFSETX_M(K),1)
                  CALL PGERRY(1,-2.0-0.2,FOFFSETX_M(K)-FOFFSETX_S(K),FOFFSETX_M(K)+FOFFSETX_S(K),1.0)
                  CALL PGPOINT(1,-1.0-0.2,FOFFSETX_WM(K),1)
                  CALL PGERRY(1,-1.0-0.2,FOFFSETX_WM(K)-FOFFSETX_WS(K),FOFFSETX_WM(K)+FOFFSETX_WS(K),1.0)
                  CALL PGSCI(2)
                  CALL PGPOINT(1,-2.0+0.2,FOFFSETY_M(K),1)
                  CALL PGERRY(1,-2.0+0.2,FOFFSETY_M(K)-FOFFSETY_S(K),FOFFSETY_M(K)+FOFFSETY_S(K),1.0)
                  CALL PGPOINT(1,-1.0+0.2,FOFFSETY_WM(K),1)
                  CALL PGERRY(1,-1.0+0.2,FOFFSETY_WM(K)-FOFFSETY_WS(K),FOFFSETY_WM(K)+FOFFSETY_WS(K),1.0)
                  CALL PGSCI(1)
                END IF
              END IF
            END DO
            IF(NPLOT.EQ.2)THEN
              CALL PGCLOS(IDNEW)
              CALL PGSLCT(IDOLD)
            END IF
          END DO
          CALL PGSVP(XV1,XV2,YV1,YV2)
          CALL PGSWIN(XW1,XW2,YW1,YW2)
          CALL PGSCH(OLD_CH)
!..............................................................................
! si el numero de objetos ajustados es mayor o igual a 3, podemos calcular
! las correspondientes transformaciones afines (de esta manera podemos estimar
! la importancia de una posible distorsion/rotacion). Tomamos como referencia
! los objetos del frame #1 (notar que el error en los puntos del frame #1 se
! los pasamos a la subrutina de calculo de la transformacion afin como si
! fueran errores de la imagen K-esima, con K>1).
          IF(NOBJFIT.GE.1)THEN
            DO K=2,NFRAMES_ !para los demas frames que no son el #1
              WRITE(*,*)
              WRITE(*,100) '> Frame #'
              WRITE(*,*) K
              NOBJEFF=0 !calculamos el numero real de objetos en comun con #1
              DO NO=1,NOBJFIT
                IF(LFITOK(NO,K))THEN
                  NOBJEFF=NOBJEFF+1
                  XDATA(NOBJEFF)=X0(NO,1)
                  EXDATA(NOBJEFF)=EX0(NO,1)
                  YDATA(NOBJEFF)=Y0(NO,1)
                  EYDATA(NOBJEFF)=EY0(NO,1)
                  UDATA(NOBJEFF)=X0(NO,K)-DJ(K)*NSIZEFB9(NEWBUFF1)+DJ(1)*NSIZEFB9(NEWBUFF1)-FOFFSETX(K)
                  EUDATA(NOBJEFF)=EX0(NO,K)
                  VDATA(NOBJEFF)=Y0(NO,K)-DI(K)*NSIZEFB9(NEWBUFF1)+DI(1)*NSIZEFB9(NEWBUFF1)-FOFFSETY(K)
                  EVDATA(NOBJEFF)=EY0(NO,K)
                END IF
              END DO
              IF(NOBJEFF.GE.1)THEN !en la imagen K-esima hay al menos 1 punto
                CALL AFINE_TRA(NOBJEFF,XDATA,YDATA,UDATA,VDATA,A11,A21,A31,A12,A22,A32,CHISQR)
                IF(NSIMUL.GT.0)THEN
                  DO NS=1,NSIMUL
                    DO NO=1,NOBJEFF
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORX=1.41421356*EXDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      XDATA_(NO)=XDATA(NO)+ERRORX
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORY=1.41421356*EYDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      YDATA_(NO)=YDATA(NO)+ERRORY
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORU=1.41421356*EUDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      UDATA_(NO)=UDATA(NO)+ERRORU
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORV=1.41421356*EVDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      VDATA_(NO)=VDATA(NO)+ERRORV
                    END DO
                    CALL AFINE_TRA(NOBJEFF,XDATA_,YDATA_,UDATA_,VDATA_, &
                     A11_(NS),A21_(NS),A31_(NS),A12_(NS),A22_(NS),A32_(NS),CHISQR_(NS))
                  END DO
                  A31_M=FMEAN0(NSIMUL,A31_,A31_S)
                  A32_M=FMEAN0(NSIMUL,A32_,A32_S)
                  CHISQR_M=FMEAN0(NSIMUL,CHISQR_,CHISQR_S)
                END IF
                WRITE(*,100) 'Afine transformation (translation)'
                WRITE(*,100) '.........:'
                WRITE(*,*) CHISQR,CHISQR_S
                WRITE(*,100) 'a31,a32: '
                WRITE(*,*) A31,A32,A31_S,A32_S
!
                CALL AFINE_ESC(NOBJEFF,XDATA,YDATA,UDATA,VDATA,A11,A21,A31,A12,A22,A32,CHISQR)
                IF(NSIMUL.GT.0)THEN
                  DO NS=1,NSIMUL
                    DO NO=1,NOBJEFF
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORX=1.41421356*EXDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      XDATA_(NO)=XDATA(NO)+ERRORX
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORY=1.41421356*EYDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      YDATA_(NO)=YDATA(NO)+ERRORY
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORU=1.41421356*EUDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      UDATA_(NO)=UDATA(NO)+ERRORU
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORV=1.41421356*EVDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      VDATA_(NO)=VDATA(NO)+ERRORV
                    END DO
                    CALL AFINE_ESC(NOBJEFF,XDATA_,YDATA_,UDATA_,VDATA_, &
                     A11_(NS),A21_(NS),A31_(NS),A12_(NS),A22_(NS),A32_(NS),CHISQR_(NS))
                  END DO
                  A11_M=FMEAN0(NSIMUL,A11_,A11_S)
                  A22_M=FMEAN0(NSIMUL,A22_,A22_S)
                  CHISQR_M=FMEAN0(NSIMUL,CHISQR_,CHISQR_S)
                END IF
                WRITE(*,100) 'Afine transformation (scale)......'
                WRITE(*,100) '.........:'
                WRITE(*,*) CHISQR,CHISQR_S
                WRITE(*,100) 'a11,a22: '
                WRITE(*,*) A11,A22,A11_S,A22_S
!
                CALL AFINE_SHX(NOBJEFF,XDATA,YDATA,UDATA,VDATA,A11,A21,A31,A12,A22,A32,CHISQR)
                IF(NSIMUL.GT.0)THEN
                  DO NS=1,NSIMUL
                    DO NO=1,NOBJEFF
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORX=1.41421356*EXDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      XDATA_(NO)=XDATA(NO)+ERRORX
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORY=1.41421356*EYDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      YDATA_(NO)=YDATA(NO)+ERRORY
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORU=1.41421356*EUDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      UDATA_(NO)=UDATA(NO)+ERRORU
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORV=1.41421356*EVDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      VDATA_(NO)=VDATA(NO)+ERRORV
                    END DO
                    CALL AFINE_SHX(NOBJEFF,XDATA_,YDATA_,UDATA_,VDATA_, &
                     A11_(NS),A21_(NS),A31_(NS),A12_(NS),A22_(NS),A32_(NS),CHISQR_(NS))
                  END DO
                  A12_M=FMEAN0(NSIMUL,A12_,A12_S)
                  CHISQR_M=FMEAN0(NSIMUL,CHISQR_,CHISQR_S)
                END IF
                WRITE(*,100) 'Afine transformation (shear X)....'
                WRITE(*,100) '.........:'
                WRITE(*,*) CHISQR,CHISQR_S
                WRITE(*,100) 'a12: '
                WRITE(*,*) A12,A12_S
!
                CALL AFINE_SHY(NOBJEFF,XDATA,YDATA,UDATA,VDATA,A11,A21,A31,A12,A22,A32,CHISQR)
                IF(NSIMUL.GT.0)THEN
                  DO NS=1,NSIMUL
                    DO NO=1,NOBJEFF
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORX=1.41421356*EXDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      XDATA_(NO)=XDATA(NO)+ERRORX
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORY=1.41421356*EYDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      YDATA_(NO)=YDATA(NO)+ERRORY
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORU=1.41421356*EUDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      UDATA_(NO)=UDATA(NO)+ERRORU
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORV=1.41421356*EVDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      VDATA_(NO)=VDATA(NO)+ERRORV
                    END DO
                    CALL AFINE_SHY(NOBJEFF,XDATA_,YDATA_,UDATA_,VDATA_, &
                     A11_(NS),A21_(NS),A31_(NS),A12_(NS),A22_(NS),A32_(NS),CHISQR_(NS))
                  END DO
                  A21_M=FMEAN0(NSIMUL,A21_,A21_S)
                  CHISQR_M=FMEAN0(NSIMUL,CHISQR_,CHISQR_S)
                END IF
                WRITE(*,100) 'Afine transformation (shear Y)....'
                WRITE(*,100) '.........:'
                WRITE(*,*) CHISQR,CHISQR_S
                WRITE(*,100) 'a21: '
                WRITE(*,*) A21,A21_S
              END IF
!
              IF(NOBJEFF.GE.2)THEN !en la imagen K-esima hay al menos 2 puntos
                CALL AFINE_ROT(NOBJEFF,XDATA,YDATA,UDATA,VDATA,A11,A21,A31,A12,A22,A32,CHISQR)
                IF(NSIMUL.GT.0)THEN
                  DO NS=1,NSIMUL
                    DO NO=1,NOBJEFF
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORX=1.41421356*EXDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      XDATA_(NO)=XDATA(NO)+ERRORX
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORY=1.41421356*EYDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      YDATA_(NO)=YDATA(NO)+ERRORY
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORU=1.41421356*EUDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      UDATA_(NO)=UDATA(NO)+ERRORU
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORV=1.41421356*EVDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      VDATA_(NO)=VDATA(NO)+ERRORV
                    END DO
                    CALL AFINE_ROT(NOBJEFF,XDATA_,YDATA_,UDATA_,VDATA_, &
                     A11_(NS),A21_(NS),A31_(NS),A12_(NS),A22_(NS),A32_(NS),CHISQR_(NS))
                  END DO
                  A11_M=FMEAN0(NSIMUL,A11_,A11_S)
                  A12_M=FMEAN0(NSIMUL,A12_,A12_S)
                  CHISQR_M=FMEAN0(NSIMUL,CHISQR_,CHISQR_S)
                END IF
                WRITE(*,100) 'Afine transformation (rotation)...'
                WRITE(*,100) '.........:'
                WRITE(*,*) CHISQR,CHISQR_S
                WRITE(*,100) 'a11,a12: '
                WRITE(*,*) A11,A12,A11_S,A12_S
              END IF
!
              IF(NOBJEFF.GE.3)THEN !en la imagen K-esima hay al menos 3 puntos
                CALL AFINE_GEN(NOBJEFF,XDATA,YDATA,UDATA,VDATA,A11,A21,A31,A12,A22,A32,CHISQR)
                IF(NSIMUL.GT.0)THEN
                  DO NS=1,NSIMUL
                    DO NO=1,NOBJEFF
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORX=1.41421356*EXDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      XDATA_(NO)=XDATA(NO)+ERRORX
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORY=1.41421356*EYDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      YDATA_(NO)=YDATA(NO)+ERRORY
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORU=1.41421356*EUDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      UDATA_(NO)=UDATA(NO)+ERRORU
                      R1=RANDOMNUMBER(NSEED)
                      R2=RANDOMNUMBER(NSEED)
                      ERRORV=1.41421356*EVDATA(NO)*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
                      VDATA_(NO)=VDATA(NO)+ERRORV
                    END DO
                    CALL AFINE_GEN(NOBJEFF,XDATA_,YDATA_,UDATA_,VDATA_, &
                     A11_(NS),A21_(NS),A31_(NS),A12_(NS),A22_(NS),A32_(NS),CHISQR_(NS))
                  END DO
                  A11_M=FMEAN0(NSIMUL,A11_,A11_S)
                  A21_M=FMEAN0(NSIMUL,A21_,A21_S)
                  A31_M=FMEAN0(NSIMUL,A31_,A31_S)
                  A12_M=FMEAN0(NSIMUL,A12_,A12_S)
                  A22_M=FMEAN0(NSIMUL,A22_,A22_S)
                  A32_M=FMEAN0(NSIMUL,A32_,A32_S)
                  CHISQR_M=FMEAN0(NSIMUL,CHISQR_,CHISQR_S)
                END IF
                WRITE(*,100) 'Afine transformation (general)....'
                WRITE(*,100) '.........:'
                WRITE(*,*) CHISQR,CHISQR_S
                WRITE(*,100) 'a11,a12: '
                WRITE(*,*) A11,A12,A11_S,A12_S
                WRITE(*,100) 'a21,a22: '
                WRITE(*,*) A21,A22,A21_S,A22_S
                WRITE(*,100) 'a31,a32: '
                WRITE(*,*) A31,A32,A31_S,A32_S
              END IF
            END DO
          END IF
!..............................................................................
          WRITE(*,*)
          WRITE(*,100) 'Is the last fit OK (use mouse...)'
          CALL RPGBAND(0,0,0.,0.,XC_,YC_,CFITOK)
          IF(CFITOK.EQ.'X')THEN
            CFITOK='n'
          ELSE
            CFITOK='y'
          END IF
          IF(CFITOK.EQ.'y')THEN
            IF(9*NBOX_CENTROID.LE.NXMAX/3)THEN
              DO I=1,3*NBOX_CENTROID
                DO J=1,3*NBOX_CENTROID
                  IMAGEN(J,I+NAXIS(2,NEWBUFF4),NEWBUFF4)=IMAGEN(J,I,NEWBUFF2)
                END DO
              END DO
              DO I=1,3*NBOX_CENTROID
                DO J=1,3*NBOX_CENTROID
                  IMAGEN(J+3*NBOX_CENTROID,I+NAXIS(2,NEWBUFF4),NEWBUFF4)=IMAGEN(J,I,NEWBUFF3)
                END DO
              END DO
              DO I=1,3*NBOX_CENTROID
                DO J=1,3*NBOX_CENTROID
                  IMAGEN(J+6*NBOX_CENTROID,I+NAXIS(2,NEWBUFF4),NEWBUFF4)=IMAGEN(J,I,NEWBUFF2)-IMAGEN(J,I,NEWBUFF3)
                END DO
              END DO
              NAXIS(1,NEWBUFF4)=MAX0(NAXIS(1,NEWBUFF4),9*NBOX_CENTROID)
              NAXIS(2,NEWBUFF4)=NAXIS(2,NEWBUFF4)+3*NBOX_CENTROID
            END IF
          ELSE
            NOBJFIT=NOBJFIT-1
            IF(NOBJFIT.GT.0)THEN
              WRITE(*,*)
              WRITE(*,101) 'Recomputing offsets after removing last fit...'
              DO K=1,NFRAMES_
                NOBJEFF=0
                DO NO=1,NOBJFIT
                  IF(LFITOK(NO,K))THEN
                    NOBJEFF=NOBJEFF+1
                    XDATA(NOBJEFF)=FOFFSETX_(NO,K)
                    EXDATA(NOBJEFF)=EFOFFSETX_(NO,K)
                    YDATA(NOBJEFF)=FOFFSETY_(NO,K)
                    EYDATA(NOBJEFF)=EFOFFSETY_(NO,K)
                  END IF
                END DO
                IF(NOBJEFF.GT.0)THEN
                  FOFFSETX_M(K)=FMEAN0(NOBJEFF,XDATA,FOFFSETX_S(K))
                  FOFFSETX_WM(K)=FMEAN0W(NOBJEFF,XDATA,EXDATA,FOFFSETX_WS(K))
                  FOFFSETY_M(K)=FMEAN0(NOBJEFF,YDATA,FOFFSETY_S(K))
                  FOFFSETY_WM(K)=FMEAN0W(NOBJEFF,YDATA,EYDATA,FOFFSETY_WS(K))
                  WRITE(*,100) 'Offsets > '
                  WRITE(*,*) K,FOFFSETX_M(K)+FOFFSETX(K),FOFFSETX_S(K),FOFFSETY_M(K)+FOFFSETY(K),FOFFSETY_S(K)
                  WRITE(*,100) 'OffsetsW> '
                  WRITE(*,*) K,FOFFSETX_WM(K)+FOFFSETX(K),FOFFSETX_WS(K),FOFFSETY_WM(K)+FOFFSETY(K),FOFFSETY_WS(K)
                ELSE
                  FOFFSETX_M(K)=0.
                  FOFFSETX_S(K)=0.
                  FOFFSETX_WM(K)=0.
                  FOFFSETX_WS(K)=0.
                  WRITE(*,101) 'Offsets > unknown'
                  WRITE(*,101) 'OffsetsW> unknown'
                END IF
              END DO
            END IF
          END IF
!..............................................................................
          WRITE(*,*)
          WRITE(*,100) 'Add more offsets measurement (use mouse...)'
          CALL RPGBAND(0,0,0.,0.,XC_,YC_,CCONT)
          IF(CCONT.EQ.'X')THEN
            CCONT='n'
          ELSE
            CCONT='y'
          END IF
          IF(CCONT.EQ.'y')THEN
            LOOP=.TRUE.
          ELSE
            LOOP=.FALSE.
          END IF
!..............................................................................
        END DO
!------------------------------------------------------------------------------
! si no tenemos ningun ajuste, regresamos
        IF(NOBJFIT.EQ.0)THEN
          RETURN
        END IF
!------------------------------------------------------------------------------
! damos la posibilidad de salvar los offsets calculados
        DO K=1,NFRAMES_
          WRITE(*,100) 'Final Offsets > '
          WRITE(*,*) K,FOFFSETX_M(K)+FOFFSETX(K),FOFFSETX_S(K),FOFFSETY_M(K)+FOFFSETY(K),FOFFSETY_S(K)
          WRITE(*,100) 'Final OffsetsW> '
          WRITE(*,*) K,FOFFSETX_WM(K)+FOFFSETX(K),FOFFSETX_WS(K),FOFFSETY_WM(K)+FOFFSETY(K),FOFFSETY_WS(K)
        END DO
        C255=READC('Save these measured offsets into a file (y/n)','n','yn')
        CSAVE=C255(1:1)
        IF(CSAVE.EQ.'y')THEN
          C255=READC('Save [w]eighted or [u]nweighted fits (w/u)','w','wu')
          CWEIGHT=C255(1:1)
          LOGFILE=.TRUE.
          DO WHILE(LOGFILE)
            OUTFILE=READC('Output file name','@','@')
            L1=TRUEBEG(OUTFILE)
            L2=TRUELEN(OUTFILE)
            INQUIRE(FILE=OUTFILE(L1:L2),EXIST=LOGFILE)
            IF(LOGFILE)THEN
              WRITE(*,100) 'ERROR: this file already exist.'
              WRITE(*,101) ' Try again.'
              WRITE(*,100) 'Press <CR> to continue...'
              READ(*,*)
            END IF
          END DO
          OPEN(30,FILE=OUTFILE(L1:L2),STATUS='NEW',FORM='FORMATTED')
          IF(CWEIGHT.EQ.'w')THEN
            DO K=1,NFRAMES_
              WRITE(30,*) K,FOFFSETX_WM(K)+FOFFSETX(K),FOFFSETX_WS(K),FOFFSETY_WM(K)+FOFFSETY(K),FOFFSETY_WS(K)
            END DO
          ELSE
            DO K=1,NFRAMES_
              WRITE(30,*) K,FOFFSETX_M(K)+FOFFSETX(K),FOFFSETX_S(K),FOFFSETY_M(K)+FOFFSETY(K),FOFFSETY_S(K)
            END DO
          END IF
          CLOSE(30)
        END IF
!------------------------------------------------------------------------------
! redibujamos botones del histograma
        CALL RPGERASW(0.00,0.40,0.32,0.38,0)
        CALL RPGERASW(0.00,0.40,0.00,0.32,0)
        CALL BUTTON(161,'zoom',0)
        CALL BUTTON(161,'zoom',-5)
        CALL BUTTON(162,'min[,]max',0)
        CALL BUTTON(162,'min[,]max',-5)
        CALL BUTTON(163,'z1[/]z2',0)
        CALL BUTTON(163,'z1[/]z2',-5)
        CALL BUTTON(164,'BG[:]FG',0)
        CALL BUTTON(164,'BG[:]FG',-5)
        CALL HISTOGRAM(NEWBUFF2)
        CPOST='n'
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
!
!******************************************************************************
! Calcula la transformacion afin general (translation+rotation+scale+shear) 
! usando NOBJEFF puntos. La solucion se obtiene resolviendo dos sistemas de 
! ecuaciones con tres incognitas (para NOBJEFF>3 el sistema esta 
! sobredeterminado).
        SUBROUTINE AFINE_GEN(N,X,Y,U,V,A11,A21,A31,A12,A22,A32,CHISQR)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL Y(N)
        REAL U(N)
        REAL V(N)
        REAL A11,A21,A31
        REAL A12,A22,A32
        REAL CHISQR
!
        REAL DETER3P3
!
        INTEGER I
        REAL SUMU,SUMV,SUMUU,SUMVV,SUMUV
        REAL SUMXU,SUMXV,SUMX
        REAL SUMYU,SUMYV,SUMY
        REAL DETERMINANTE
!------------------------------------------------------------------------------
        A11=1.
        A21=0.
        A31=0.
        A12=0.
        A22=1.
        A32=0.
        CHISQR=0.
! precaucion
        IF(N.LT.3)THEN
          WRITE(*,100) 'N = '
          WRITE(*,*) N
          WRITE(*,101) 'ERROR in AFINE_GEN: number of points < 3'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
!------------------------------------------------------------------------------
        SUMU=0.
        SUMV=0.
        SUMUU=0.
        SUMVV=0.
        SUMUV=0.
        SUMXU=0.
        SUMXV=0.
        SUMX=0.
        SUMYU=0.
        SUMYV=0.
        SUMY=0.
        DO I=1,N
          SUMU=SUMU+U(I)
          SUMV=SUMV+V(I)
          SUMUU=SUMUU+U(I)*U(I)
          SUMVV=SUMVV+V(I)*V(I)
          SUMUV=SUMUV+U(I)*V(I)
          SUMXU=SUMXU+U(I)*X(I)
          SUMXV=SUMXV+V(I)*X(I)
          SUMX=SUMX+X(I)
          SUMYU=SUMYU+U(I)*Y(I)
          SUMYV=SUMYV+V(I)*Y(I)
          SUMY=SUMY+Y(I)
        END DO
! resolvemos el problema aplicando la regla de Cramer
        DETERMINANTE=DETER3P3(SUMUU,SUMUV,SUMU, SUMUV,SUMVV,SUMV, SUMU,SUMV,REAL(N))
        IF(DETERMINANTE.GE.1.E-6)THEN
          A11=DETER3P3(SUMXU,SUMUV,SUMU, SUMXV,SUMVV,SUMV, SUMX,SUMV,REAL(N))/DETERMINANTE
          A21=DETER3P3(SUMUU,SUMXU,SUMU, SUMUV,SUMXV,SUMV, SUMU,SUMX,REAL(N))/DETERMINANTE
          A31=DETER3P3(SUMUU,SUMUV,SUMXU, SUMUV,SUMVV,SUMXV, SUMU,SUMV,SUMX)/DETERMINANTE
          A12=DETER3P3(SUMYU,SUMUV,SUMU, SUMYV,SUMVV,SUMV, SUMY,SUMV,REAL(N))/DETERMINANTE
          A22=DETER3P3(SUMUU,SUMYU,SUMU, SUMUV,SUMYV,SUMV, SUMU,SUMY,REAL(N))/DETERMINANTE
          A32=DETER3P3(SUMUU,SUMUV,SUMYU, SUMUV,SUMVV,SUMYV, SUMU,SUMV,SUMY)/DETERMINANTE
        ELSE
          WRITE(*,100) 'WARNING in AFINE_GEN: Determinant ='
          WRITE(*,*) DETERMINANTE
        END IF
!
        DO I=1,N
          CHISQR=CHISQR+ &
           (A11*U(I)+A21*V(I)+A31-X(I))*(A11*U(I)+A21*V(I)+A31-X(I))+ &
           (A12*U(I)+A22*V(I)+A32-Y(I))*(A12*U(I)+A22*V(I)+A32-Y(I))
        END DO
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
!
!******************************************************************************
! Calcula la transformacion afin traslacion pura.
        SUBROUTINE AFINE_TRA(N,X,Y,U,V,A11,A21,A31,A12,A22,A32,CHISQR)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL Y(N)
        REAL U(N)
        REAL V(N)
        REAL A11,A21,A31
        REAL A12,A22,A32
        REAL CHISQR
!
        INTEGER I
        REAL SUMU,SUMV
        REAL SUMX,SUMY
!------------------------------------------------------------------------------
        A11=1.
        A21=0.
        A31=0.
        A12=0.
        A22=1.
        A32=0.
        CHISQR=0.
! precaucion
        IF(N.LT.1)THEN
          WRITE(*,100) 'N = '
          WRITE(*,*) N
          WRITE(*,101) 'ERROR in AFINE_TRA: number of points < 1'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
!------------------------------------------------------------------------------
        SUMU=0.
        SUMV=0.
        SUMX=0.
        SUMY=0.
        DO I=1,N
          SUMU=SUMU+U(I)
          SUMV=SUMV+V(I)
          SUMX=SUMX+X(I)
          SUMY=SUMY+Y(I)
        END DO
        A31=SUMX/REAL(N)-SUMU/REAL(N)
        A32=SUMY/REAL(N)-SUMV/REAL(N)
!
        DO I=1,N
          CHISQR=CHISQR+ &
           (A11*U(I)+A21*V(I)+A31-X(I))*(A11*U(I)+A21*V(I)+A31-X(I))+ &
           (A12*U(I)+A22*V(I)+A32-Y(I))*(A12*U(I)+A22*V(I)+A32-Y(I))
        END DO
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
!
!******************************************************************************
! Calcula la transformacion afin rotacion pura.
        SUBROUTINE AFINE_ROT(N,X,Y,U,V,A11,A21,A31,A12,A22,A32,CHISQR)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL Y(N)
        REAL U(N)
        REAL V(N)
        REAL A11,A21,A31
        REAL A12,A22,A32
        REAL CHISQR
!
        INTEGER I
        REAL SUMUU,SUMVV
        REAL SUMXV,SUMXU,SUMYU,SUMYV
        REAL COSTHETA,SENTHETA
        REAL LAMBDA1,LAMBDA2
        REAL A,B,C,DELTA
        REAL A11_,A12_,A21_,A22_,CHISQR_
!------------------------------------------------------------------------------
        A11=1.
        A21=0.
        A31=0.
        A12=0.
        A22=1.
        A32=0.
        CHISQR=0.
! precaucion
        IF(N.LT.2)THEN
          WRITE(*,100) 'N = '
          WRITE(*,*) N
          WRITE(*,101) 'ERROR in AFINE_ROT: number of points < 2'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
!------------------------------------------------------------------------------
        SUMUU=0.
        SUMVV=0.
        SUMXU=0.
        SUMXV=0.
        SUMYU=0.
        SUMYV=0.
        DO I=1,N
          SUMUU=SUMUU+U(I)*U(I)
          SUMVV=SUMVV+V(I)*V(I)
          SUMXU=SUMXU+X(I)*U(I)
          SUMXV=SUMXV+X(I)*V(I)
          SUMYU=SUMYU+Y(I)*U(I)
          SUMYV=SUMYV+Y(I)*V(I)
        END DO
!
        A=1.0
        B=2.0*(SUMUU+SUMVV)
        C=(SUMUU+SUMVV)*(SUMUU+SUMVV)-(SUMXU+SUMYV)*(SUMXU+SUMYV)-(SUMYU-SUMXV)*(SUMYU-SUMXV)
        DELTA=B*B-4.0*A*C
        IF(DELTA.LT.0.0)THEN
          WRITE(*,101) 'ERROR in AFINE_ROT: no solution'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        ELSE
          LAMBDA1=(-B+SQRT(DELTA))/(2.0*A)
          LAMBDA2=(-B-SQRT(DELTA))/(2.0*A)
        END IF
!
        COSTHETA=(SUMXU+SUMYV)/(SUMUU+SUMVV+LAMBDA1)
        SENTHETA=(SUMYU-SUMXV)/(SUMUU+SUMVV+LAMBDA1)
        A11=COSTHETA
        A12=SENTHETA
        A21=-SENTHETA
        A22=COSTHETA
        DO I=1,N
          CHISQR=CHISQR+ &
           (A11*U(I)+A21*V(I)+A31-X(I))*(A11*U(I)+A21*V(I)+A31-X(I))+ &
           (A12*U(I)+A22*V(I)+A32-Y(I))*(A12*U(I)+A22*V(I)+A32-Y(I))
        END DO
!
        COSTHETA=(SUMXU+SUMYV)/(SUMUU+SUMVV+LAMBDA2)
        SENTHETA=(SUMYU-SUMXV)/(SUMUU+SUMVV+LAMBDA2)
        A11_=COSTHETA
        A12_=SENTHETA
        A21_=-SENTHETA
        A22_=COSTHETA
        CHISQR_=0.
        DO I=1,N
          CHISQR_=CHISQR_+ &
           (A11_*U(I)+A21_*V(I)+A31-X(I))*(A11_*U(I)+A21_*V(I)+A31-X(I))+ &
           (A12_*U(I)+A22_*V(I)+A32-Y(I))*(A12_*U(I)+A22_*V(I)+A32-Y(I))
        END DO
!
        IF(CHISQR_.LT.CHISQR)THEN
          A11=A11_
          A12=A12_
          A21=A21_
          A22=A22_
          CHISQR=CHISQR_
        END IF
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
!
!******************************************************************************
! Calcula la transformacion afin de escala pura
        SUBROUTINE AFINE_ESC(N,X,Y,U,V,A11,A21,A31,A12,A22,A32,CHISQR)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL Y(N)
        REAL U(N)
        REAL V(N)
        REAL A11,A21,A31
        REAL A12,A22,A32
        REAL CHISQR
!
        INTEGER I
        REAL SUMUU,SUMVV
        REAL SUMXU,SUMYV
!------------------------------------------------------------------------------
        A11=1.
        A21=0.
        A31=0.
        A12=0.
        A22=1.
        A32=0.
        CHISQR=0.
! precaucion
        IF(N.LT.1)THEN
          WRITE(*,100) 'N = '
          WRITE(*,*) N
          WRITE(*,101) 'ERROR in AFINE_ESC: number of points < 1'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
!------------------------------------------------------------------------------
        SUMUU=0.
        SUMVV=0.
        SUMXU=0.
        SUMYV=0.
        DO I=1,N
          SUMUU=SUMUU+U(I)*U(I)
          SUMVV=SUMVV+V(I)*V(I)
          SUMXU=SUMXU+X(I)*U(I)
          SUMYV=SUMYV+Y(I)*V(I)
        END DO
        A11=SUMXU/SUMUU
        A22=SUMYV/SUMVV
!
        DO I=1,N
          CHISQR=CHISQR+ &
           (A11*U(I)+A21*V(I)+A31-X(I))*(A11*U(I)+A21*V(I)+A31-X(I))+ &
           (A12*U(I)+A22*V(I)+A32-Y(I))*(A12*U(I)+A22*V(I)+A32-Y(I))
        END DO
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
!
!******************************************************************************
! Calcula la transformacion afin de shear en X
        SUBROUTINE AFINE_SHX(N,X,Y,U,V,A11,A21,A31,A12,A22,A32,CHISQR)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL Y(N)
        REAL U(N)
        REAL V(N)
        REAL A11,A21,A31
        REAL A12,A22,A32
        REAL CHISQR
!
        INTEGER I
        REAL SUMUU
        REAL SUMYU,SUMUV
!------------------------------------------------------------------------------
        A11=1.
        A21=0.
        A31=0.
        A12=0.
        A22=1.
        A32=0.
        CHISQR=0.
! precaucion
        IF(N.LT.1)THEN
          WRITE(*,100) 'N = '
          WRITE(*,*) N
          WRITE(*,101) 'ERROR in AFINE_SHX: number of points < 1'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
!------------------------------------------------------------------------------
        SUMUU=0.
        SUMYU=0.
        SUMUV=0.
        DO I=1,N
          SUMUU=SUMUU+U(I)*U(I)
          SUMYU=SUMYU+Y(I)*U(I)
          SUMUV=SUMUV+U(I)*V(I)
        END DO
        A12=(SUMYU-SUMUV)/SUMUU
!
        DO I=1,N
          CHISQR=CHISQR+ &
           (A11*U(I)+A21*V(I)+A31-X(I))*(A11*U(I)+A21*V(I)+A31-X(I))+ &
           (A12*U(I)+A22*V(I)+A32-Y(I))*(A12*U(I)+A22*V(I)+A32-Y(I))
        END DO
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
!
!******************************************************************************
! Calcula la transformacion afin de shear en Y
        SUBROUTINE AFINE_SHY(N,X,Y,U,V,A11,A21,A31,A12,A22,A32,CHISQR)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL Y(N)
        REAL U(N)
        REAL V(N)
        REAL A11,A21,A31
        REAL A12,A22,A32
        REAL CHISQR
!
        INTEGER I
        REAL SUMVV
        REAL SUMXV,SUMUV
!------------------------------------------------------------------------------
        A11=1.
        A21=0.
        A31=0.
        A12=0.
        A22=1.
        A32=0.
        CHISQR=0.
! precaucion
        IF(N.LT.1)THEN
          WRITE(*,100) 'N = '
          WRITE(*,*) N
          WRITE(*,101) 'ERROR in AFINE_SHY: number of points < 1'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
!------------------------------------------------------------------------------
        SUMVV=0.
        SUMXV=0.
        SUMUV=0.
        DO I=1,N
          SUMVV=SUMVV+V(I)*V(I)
          SUMXV=SUMXV+X(I)*V(I)
          SUMUV=SUMUV+U(I)*V(I)
        END DO
        A21=(SUMXV-SUMUV)/SUMVV
!
        DO I=1,N
          CHISQR=CHISQR+ &
           (A11*U(I)+A21*V(I)+A31-X(I))*(A11*U(I)+A21*V(I)+A31-X(I))+ &
           (A12*U(I)+A22*V(I)+A32-Y(I))*(A12*U(I)+A22*V(I)+A32-Y(I))
        END DO
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
!
!******************************************************************************
!
        REAL FUNCTION DETER3P3(A11,A12,A13,A21,A22,A23,A31,A32,A33)
        IMPLICIT NONE
        REAL A11,A12,A13
        REAL A21,A22,A23
        REAL A31,A32,A33
!
        DETER3P3=A11*A22*A33+A21*A32*A13+A12*A23*A31-A13*A22*A31-A23*A32*A11-A21*A12*A33
!
        END
