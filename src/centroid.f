C Calcula el centroide en una region seleccionada con el raton
C Si CH='d', la posicion no se lee del raton sino que se toma de los
C parametros de entrada XC,YC. El resultado del ajuste y los errores
C correspondientes son parametros de salida de la subrutina.
C Si LGUESS=.TRUE., los parametros X0_F,Y0_F,SIGMAX_F,SIGMAY_F,BETA_F,AMP_F,
C y CTE_F son un ajuste preliminar de los datos.
C Si LSHOW=.TRUE. muestra resultados intermedios.
        SUBROUTINE CENTROID(NBUFF,CH,XC,YC,LGUESS,
     +   X0_F,Y0_F,SIGMAX_F,SIGMAY_F,BETA_F,AMP_F,CTE_F,
     +   EX0_F,EY0_F,ESIGMAX_F,ESIGMAY_F,EBETA_F,EAMP_F,ECTE_F,LSHOW)
        IMPLICIT NONE
        INTEGER NBUFF
        CHARACTER*1 CH
        REAL XC,YC
        LOGICAL LGUESS
        REAL X0_F,Y0_F,SIGMAX_F,SIGMAY_F,BETA_F,AMP_F,CTE_F
        REAL EX0_F,EY0_F,ESIGMAX_F,ESIGMAY_F,EBETA_F,EAMP_F,ECTE_F
        LOGICAL LSHOW
C parameters
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
C functions
        INTEGER READI
        INTEGER READILIM
        REAL READF
C
        INTEGER I,J
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER NBOX_CENTROID
        INTEGER IX1,IX2,IY1,IY2
        INTEGER NCOLOR3
        INTEGER NSIMUL
        INTEGER NPFITX,NPFITY
        INTEGER IMODEG2D
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL IMAGEN_(NXMAX,NYMAX)
        REAL BG,FG,BG_,FG_
        REAL XV1,XV2,YV1,YV2,XW1,XW2,YW1,YW2
        REAL OLD_CH
        REAL TR(6)
        REAL YRMSTOL
        REAL X0,Y0,EX0,EY0,EEX0,EEY0
        REAL SIGMAX,SIGMAY,ESIGMAX,ESIGMAY,EESIGMAX,EESIGMAY
        REAL BETA,EBETA,EEBETA
        REAL AMP,EAMP,EEAMP
        REAL CTE,ECTE,EECTE
        REAL FACTOR,FACTOR1,FACTOR2,FACTOR3
        REAL X0_,Y0_,SIGMAX_,SIGMAY_,AMPX_,AMPY_,CTEX_,CTEY_
        REAL YMIN,YMAX,YMIN_,YMAX_,DY
        REAL XP(NXYMAX),YP(NXYMAX)
        REAL XP_(NXYMAX),YP_(NXYMAX)
        REAL FMEAN,FSIGMA
        REAL Y__,X0__,Y0__,SIGMAX__,SIGMAY__
        REAL AREA1,AREA2,AREA3,EAREA1,EAREA2,EAREA3
        REAL XC_,YC_
        REAL XOFFSET_ORIGEN,YOFFSET_ORIGEN
        REAL M00, M10, M01
        CHARACTER*1 CH_
        CHARACTER*50 CDUMMY
        LOGICAL LEXIT
C
        COMMON/BLKIMAGEN1/IMAGEN             !imagen FITS leida en formato REAL
        COMMON/BLKIMAGEN1_/IMAGEN_
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKDEFAULTS8/NBOX_CENTROID
        COMMON/BLKBGFG/BG,FG
        COMMON/BLKDEFAULTS9/YRMSTOL
        COMMON/BLKDEFAULTS10/NSIMUL
        COMMON/BLKDEFAULTS11/IMODEG2D
        COMMON/BLKDEFAULTS12/XOFFSET_ORIGEN,YOFFSET_ORIGEN
        COMMON/BLKFUNK2DG/Y__,CTE,AMP,X0__,Y0__,SIGMAX__,SIGMAY__,BETA
C------------------------------------------------------------------------------
        TR(1)=0.
        TR(2)=1.
        TR(3)=0.
        TR(4)=0.
        TR(5)=0.
        TR(6)=1.
C
        IF(CH.EQ.'C')THEN
          WRITE(CDUMMY,*) NBOX_CENTROID
          NBOX_CENTROID=READILIM('Box size in pixels (odd)',CDUMMY,
     +     5,MIN(NAXIS(1,NBUFF),NAXIS(2,NBUFF)))
          IF(MOD(NBOX_CENTROID,2).EQ.0)THEN
            NBOX_CENTROID=NBOX_CENTROID+1
            IF(NBOX_CENTROID.GT.MIN(NAXIS(1,NBUFF),NAXIS(2,NBUFF)))
     +       NBOX_CENTROID=NBOX_CENTROID-2
          END IF
          WRITE(CDUMMY,*) XOFFSET_ORIGEN
          XOFFSET_ORIGEN=READF('Xoffset for centroids',CDUMMY)
          WRITE(CDUMMY,*) YOFFSET_ORIGEN
          YOFFSET_ORIGEN=READF('Yoffset for centroids',CDUMMY)
          WRITE(CDUMMY,*) YRMSTOL
          YRMSTOL=READF('YRMSTOL for DOWNHILL',CDUMMY)
          WRITE(CDUMMY,*) NSIMUL
          NSIMUL=READI('Number of simulations',CDUMMY)
          WRITE(CDUMMY,*) IMODEG2D
          IMODEG2D=READI('2D Gaussian fit mode (0,1)',CDUMMY)
        END IF
C
        LEXIT=.FALSE.
        DO WHILE(.NOT.LEXIT)
          IF((CH.EQ.'d').OR.(CH.EQ.'D'))THEN
            XC_=XC
            YC_=YC
            CH_='A'
            LEXIT=.TRUE.
          ELSE
            WRITE(*,*)
            WRITE(*,101) 'Click around region to be fitted...'
            CALL RPGBAND(7,0,0.,0.,XC_,YC_,CH_)
          END IF
          IF(CH_.NE.'X')THEN
            IX1=NINT(XC_)-NBOX_CENTROID/2
            IF(IX1.LT.1) IX1=1
            IX2=IX1+NBOX_CENTROID-1
            IF(IX2.GT.NAXIS(1,NBUFF))THEN
              IX2=NAXIS(1,NBUFF)
              IX1=IX2-NBOX_CENTROID+1
            END IF
            IY1=NINT(YC_)-NBOX_CENTROID/2
            IF(IY1.LT.1) IY1=1
            IY2=IY1+NBOX_CENTROID-1
            IF(IY2.GT.NAXIS(2,NBUFF))THEN
              IY2=NAXIS(1,NBUFF)
              IY1=IY2-NBOX_CENTROID+1
            END IF
            CALL PGSCI(3)
            CALL PGMOVE(REAL(IX1),REAL(IY1))
            CALL PGDRAW(REAL(IX2),REAL(IY1))
            CALL PGDRAW(REAL(IX2),REAL(IY2))
            CALL PGDRAW(REAL(IX1),REAL(IY2))
            CALL PGDRAW(REAL(IX1),REAL(IY1))
            CALL PGSCI(1)
            DO I=IY1,IY2
              DO J=IX1,IX2
                IMAGEN_(J,I)=IMAGEN(J,I,NBUFF)
              END DO
            END DO
C calculo sencillo del centro de masas (usando momentos)
C (ver https://en.wikipedia.org/wiki/Image_moment)
            M00 = 0.0
            M10 = 0.0
            M01 = 0.0
            DO I=IY1,IY2
              DO J=IX1,IX2
                M00 = M00 + IMAGEN_(J,I)
                M10 = M10 + REAL(J)*IMAGEN_(J,I)
                M01 = M01 + REAL(I)*IMAGEN_(J,I)
              END DO
            END DO
            WRITE(*,100) 'X_cm (center of mass): '
            WRITE(*,*) M10/M00
            WRITE(*,100) 'Y_cm (center of mass): '
            WRITE(*,*) M01/M00
C comenzamos buffering
            CALL PGBBUF
C almacenamos region de dibujo actual
            CALL PGQVP(0,XV1,XV2,YV1,YV2)
            CALL PGQWIN(XW1,XW2,YW1,YW2)
            CALL PGQCH(OLD_CH)
            CALL PGSCH(0.7)
C borramos el plot previo
            CALL RPGERASW(0.00,0.40,0.00,0.80,0)
C definimos nueva region de dibujo para el dibujo 3D y lo hacemos
            CALL PGSVP(0.01,0.19,0.39,0.79)
            CALL PLOT3DBARS(IX1,IX2,IY1,IY2,0,.TRUE.,BG,FG)
C definimos nueva region de dibujo para el dibujo 2D y lo hacemos
            CALL PGSVP(0.0275,0.1425,0.21,0.38)
            CALL PGWNAD(REAL(IX1)-0.6,REAL(IX2)+0.6,
     +                  REAL(IY1)-0.6,REAL(IY2)+0.6)
            CALL PGBOX('BC',0.0,0,'BC',0.0,0)
            CALL PGIMAG(IMAGEN(1,1,NBUFF),NXMAX,NYMAX,
     +       IX1,IX2,IY1,IY2,FG,BG,TR)
            CALL PGSCI(7)
            CALL PGMTXT('T',-1.5,0.5,0.5,'original data')
            CALL PGSCI(1)
C si tenemos un ajuste inicial, lo pasamos a GAUSS2DFIT
            IF(LGUESS)THEN
              X0=X0_F
              Y0=Y0_F
              SIGMAX=SIGMAX_F
              SIGMAY=SIGMAY_F
              BETA=BETA_F
              AMP=AMP_F
              CTE=CTE_F
            END IF
C realizamos ajuste
            CALL GAUSS2DFIT(NBUFF,IX1,IX2,IY1,IY2,LGUESS,
     +       X0,Y0,SIGMAX,SIGMAY,BETA,AMP,CTE,
     +       EX0,EY0,ESIGMAX,ESIGMAY,EBETA,EAMP,ECTE,
     +       EEX0,EEY0,EESIGMAX,EESIGMAY,EEBETA,EEAMP,EECTE,
     +       X0_,Y0_,SIGMAX_,SIGMAY_,AMPX_,AMPY_,CTEX_,CTEY_,
     +       AREA1,AREA2,AREA3,EAREA1,EAREA2,EAREA3,
     +       YRMSTOL,NSIMUL,IMODEG2D,LSHOW)
C calculamos nueva region ajustada
            DO I=IY1,IY2
              DO J=IX1,IX2
                FACTOR1=(REAL(J)-X0)*(REAL(J)-X0)/(2.*SIGMAX*SIGMAX)
                FACTOR2=(REAL(I)-Y0)*(REAL(I)-Y0)/(2.*SIGMAY*SIGMAY)
                FACTOR3=BETA*(REAL(J)-X0)*(REAL(I)-Y0)/(SIGMAX*SIGMAY)
                FACTOR=FACTOR1+FACTOR2+FACTOR3
                IF(FACTOR.GT.60.)THEN
                  IMAGEN_(J,I)=CTE
                ELSE
                  IMAGEN_(J,I)=CTE+AMP*EXP(-FACTOR)
                END IF
              END DO
            END DO
C dibujamos resultados del ajuste 1D: calculamos limites
            NPFITX=IX2-IX1+1
            DO J=IX1,IX2
              XP(J-IX1+1)=REAL(J)
              YP(J-IX1+1)=0.
              DO I=IY1,IY2
                YP(J-IX1+1)=YP(J-IX1+1)+IMAGEN(J,I,NBUFF)
              END DO
              YP(J-IX1+1)=YP(J-IX1+1)/REAL(NPFITX)
            END DO
            CALL FINDMML(NPFITX,1,NPFITX,YP,YMIN,YMAX)
            NPFITY=IY2-IY1+1
            DO I=IY1,IY2
              XP_(I-IY1+1)=REAL(I)
              YP_(I-IY1+1)=0.
              DO J=IX1,IX2
                YP_(I-IY1+1)=YP_(I-IY1+1)+IMAGEN(J,I,NBUFF)
              END DO
              YP_(I-IY1+1)=YP_(I-IY1+1)/REAL(NPFITY)
            END DO
            CALL FINDMML(NPFITY,1,NPFITY,YP_,YMIN_,YMAX_)
            IF(YMIN_.LT.YMIN) YMIN=YMIN_
            IF(YMAX_.GT.YMAX) YMAX=YMAX_
            DY=YMAX-YMIN
            YMIN=YMIN-DY/20.
            YMAX=YMAX+DY/20.
c corte promedio en X
            CALL PGSVP(0.04,0.20,0.04,0.19)
            CALL PGSWIN(REAL(IX1)-0.6,REAL(IX2)+0.6,YMIN,YMAX)
            CALL PGBOX('BCTSN',0.0,0,'BCTSN',0.0,0)
            CALL PGSCI(5)
            CALL PGMTXT('L',+2.5,0.5,0.5,'mean X cut')
            CALL PGSCI(3)
            CALL PGMTXT('T',-1.5,0.05,0.0,'1D fit')
            CALL PGSCI(2)
            CALL PGMTXT('T',-1.5,0.95,1.0,'2D fit')
            CALL PGSCI(1)
            CALL PGBIN(NPFITX,XP,YP,.TRUE.)
            DO I=1,NXYMAX
              XP(I)=REAL(IX1)+REAL(I-1)/REAL(NXYMAX-1)*REAL(IX2-IX1)
              FACTOR=(XP(I)-X0_)*(XP(I)-X0_)/(2.*SIGMAX_*SIGMAX_)
              IF(FACTOR.GT.60.)THEN
                YP(I)=CTEX_
              ELSE
                YP(I)=CTEX_+AMPX_*EXP(-FACTOR)
              END IF
            END DO
            CALL PGSCI(3)
            CALL PGLINE(NXYMAX,XP,YP)
            DO J=IX1,IX2
              XP(J-IX1+1)=REAL(J)
              YP(J-IX1+1)=0.0
              DO I=IY1,IY2
                YP(J-IX1+1)=YP(J-IX1+1)+IMAGEN_(J,I)
              END DO
              YP(J-IX1+1)=YP(J-IX1+1)/REAL(NPFITY)
            END DO
            CALL PGSCI(2)
            CALL PGBIN(NPFITX,XP,YP,.TRUE.)
            CALL PGSCI(1)
c corte promedio en Y
            CALL PGSVP(0.20,0.36,0.04,0.19)
            CALL PGSWIN(REAL(IY1)-0.6,REAL(IY2)+0.6,YMIN,YMAX)
            CALL PGBOX('BCTSN',0.0,0,'BCTSM',0.0,0)
            CALL PGSCI(5)
            CALL PGMTXT('R',+2.5,0.5,0.5,'mean Y cut')
            CALL PGSCI(3)
            CALL PGMTXT('T',-1.5,0.05,0.0,'1D fit')
            CALL PGSCI(2)
            CALL PGMTXT('T',-1.5,0.95,1.0,'2D fit')
            CALL PGSCI(1)
            CALL PGBIN(NPFITY,XP_,YP_,.TRUE.)
            DO I=1,NXYMAX
              XP_(I)=REAL(IY1)+REAL(I-1)/REAL(NXYMAX-1)*REAL(IY2-IY1)
              FACTOR=(XP_(I)-Y0_)*(XP_(I)-Y0_)/(2.*SIGMAY_*SIGMAY_)
              IF(FACTOR.GT.60.)THEN
                YP_(I)=CTEY_
              ELSE
                YP_(I)=CTEY_+AMPY_*EXP(-FACTOR)
              END IF
            END DO
            CALL PGSCI(3)
            CALL PGLINE(NXYMAX,XP_,YP_)
            DO I=IY1,IY2
              XP_(I-IY1+1)=REAL(I)
              YP_(I-IY1+1)=0.0
              DO J=IX1,IX2
                YP_(I-IY1+1)=YP_(I-IY1+1)+IMAGEN_(J,I)
              END DO
              YP_(I-IY1+1)=YP_(I-IY1+1)/REAL(NPFITX)
            END DO
            CALL PGSCI(2)
            CALL PGBIN(NPFITY,XP_,YP_,.TRUE.)
            CALL PGSCI(1)
C mostramos resultados del ajuste
            IF(LSHOW)THEN
              WRITE(*,100) 'NPFIT: '
              WRITE(*,*) (IX2-IX1+1)*(IY2-IY1+1)
              WRITE(*,100) 'X0,EX0,EEX0............: '
              WRITE(*,*) X0,EX0,EEX0,X0_
              WRITE(*,100) 'Y0,EY0,EEY0............: '
              WRITE(*,*) Y0,EY0,EEY0,Y0_
              WRITE(*,100) 'X0-XOFFSET,Y0-YOFFSET..: '
              WRITE(*,*) X0-XOFFSET_ORIGEN,Y0-YOFFSET_ORIGEN
              WRITE(*,100) 'XOFFSET-X0,YOFFSET-Y0..: '
              WRITE(*,*) XOFFSET_ORIGEN-X0,YOFFSET_ORIGEN-Y0
              WRITE(*,100) 'SIGMAX,ESIGMAX,EESIGMAX: '
              WRITE(*,*) SIGMAX,ESIGMAX,EESIGMAX,SIGMAX_
              WRITE(*,100) 'SIGMAY,ESIGMAY,EESIGMAY: '
              WRITE(*,*) SIGMAY,ESIGMAY,EESIGMAY,SIGMAY_
              WRITE(*,100) 'BETA,EBETA,EEBETA......: '
              WRITE(*,*) BETA,EBETA,EEBETA
              WRITE(*,100) 'AMP,EAMP,EEAMP.........: '
              WRITE(*,*) AMP,EAMP,EEAMP,
     +         (AMPX_*2.5*SIGMAY_+AMPY_*2.5*SIGMAX_)/2.
              WRITE(*,100) 'CTE,ECTE,EECTE.........: '
              WRITE(*,*) CTE,ECTE,EECTE,(CTEX_+CTEY_)/2.
              IF(EAREA1.GT.0.0)THEN
                WRITE(*,100) 'AREA1,EAREA1,ratio.....: '
                WRITE(*,*) AREA1,EAREA1,AREA1/EAREA1
              ELSE
                WRITE(*,100) 'AREA1,EAREA1...........: '
                WRITE(*,*) AREA1,EAREA1
              END IF
              IF(EAREA2.GT.0.0)THEN
                WRITE(*,100) 'AREA2,EAREA2,ratio.....: '
                WRITE(*,*) AREA2,EAREA2,AREA2/EAREA2
              ELSE
                WRITE(*,100) 'AREA2,EAREA2...........: '
                WRITE(*,*) AREA2,EAREA2
              END IF
              IF(EAREA3.GT.0.0)THEN
                WRITE(*,100) 'AREA3,EAREA3,ratio.....: '
                WRITE(*,*) AREA3,EAREA3,AREA3/EAREA3
              ELSE
                WRITE(*,100) 'AREA3,EAREA3...........: '
                WRITE(*,*) AREA3,EAREA3
              END IF
            END IF
C almacenamos valores ajustados como parametros de salida de la subrutina
            X0_F=X0
            Y0_F=Y0
            SIGMAX_F=SIGMAX
            SIGMAY_F=SIGMAY
            BETA_F=BETA
            AMP_F=AMP
            CTE_F=CTE
            EX0_F=EX0
            EY0_F=EY0
            ESIGMAX_F=ESIGMAX
            ESIGMAY_F=ESIGMAY
            EBETA_F=EBETA
            EAMP_F=EAMP
            ECTE_F=ECTE
C definimos nueva region de dibujo para el dibujo 3D y lo hacemos
            CALL PGSVP(0.21,0.39,0.39,0.79)
            CALL PLOT3DBARS(IX1,IX2,IY1,IY2,0,.TRUE.,BG,FG)
C definimos nueva region de dibujo para el dibujo 2D y lo hacemos
            CALL PGSVP(0.1425,0.2575,0.21,0.38)
            CALL PGWNAD(REAL(IX1)-0.6,REAL(IX2)+0.6,
     +                  REAL(IY1)-0.6,REAL(IY2)+0.6)
            CALL PGBOX('BC',0.0,0,'BC',0.0,0)
            CALL PGIMAG(IMAGEN_,NXMAX,NYMAX,
     +       IX1,IX2,IY1,IY2,FG,BG,TR)
            CALL PGSCI(7)
            CALL PGMTXT('T',-1.5,0.5,0.5,'2D fit')
            CALL PGSCI(1)
C calculamos residuos
            DO I=IY1,IY2
              DO J=IX1,IX2
                IMAGEN_(J,I)=IMAGEN(J,I,NBUFF)-IMAGEN_(J,I)
              END DO
            END DO
C calculamos media, r.m.s., y dibujamos residuos
            FMEAN=0.
            DO I=IY1,IY2
              DO J=IX1,IX2
                FMEAN=FMEAN+IMAGEN_(J,I)
              END DO
            END DO
            FMEAN=FMEAN/REAL(NPFITX*NPFITY)
            FSIGMA=0.
            DO I=IY1,IY2
              DO J=IX1,IX2
                FSIGMA=FSIGMA+(FMEAN-IMAGEN_(J,I))*(FMEAN-IMAGEN_(J,I))
              END DO
            END DO
            FSIGMA=SQRT(FSIGMA/REAL(NPFITX*NPFITY-1))
            BG_=FMEAN-5*FSIGMA
            FG_=FMEAN+5*FSIGMA
            IF(LSHOW)THEN
              WRITE(*,100) 'Residuals mean,sigma...: '
              WRITE(*,*) FMEAN,FSIGMA
              WRITE(*,100) 'Residuals BG,FG........: '
              WRITE(*,*) BG_,FG_
            END IF
            CALL PGSVP(0.2575,0.3725,0.21,0.38)
            CALL PGWNAD(REAL(IX1)-0.6,REAL(IX2)+0.6,
     +                  REAL(IY1)-0.6,REAL(IY2)+0.6)
            CALL PGBOX('BC',0.0,0,'BC',0.0,0)
            CALL PGIMAG(IMAGEN_,NXMAX,NYMAX,
     +       IX1,IX2,IY1,IY2,FG_,BG_,TR)
            CALL PGSCI(7)
            CALL PGMTXT('T',-1.5,0.5,0.5,'residuals')
            CALL PGSCI(1)
C recuperamos region de dibujo inicial
            CALL PGSVP(XV1,XV2,YV1,YV2)
            CALL PGSWIN(XW1,XW2,YW1,YW2)
C terminamos buffering
            CALL PGSCH(OLD_CH)
            CALL PGEBUF
          ELSE
            LEXIT=.TRUE.
          END IF
        END DO
C------------------------------------------------------------------------------
C recuperamos histograma
        IF(CH.NE.'d')THEN
          NCOLOR3=-4-1
          CALL RPGERASW(0.00,0.40,0.32,0.38,0)
          CALL RPGERASW(0.00,0.40,0.00,0.32,0)
          CALL BUTTON(161,'zoom',0)
          CALL BUTTON(161,'zoom',NCOLOR3)
          CALL BUTTON(162,'min[,]max',0)
          CALL BUTTON(162,'min[,]max',NCOLOR3)
          CALL BUTTON(163,'z1[/]z2',0)
          CALL BUTTON(163,'z1[/]z2',NCOLOR3)
          CALL BUTTON(164,'BG[:]FG',0)
          CALL BUTTON(164,'BG[:]FG',NCOLOR3)
          CALL HISTOGRAM(NBUFF)
        END IF
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C Subrutina del Numerical Recipes
      SUBROUTINE QSIMP2DG(A,B,S)
      PARAMETER (EPS=1.E-6, JMAX=20)
      OST=-1.E30
      OS= -1.E30
      DO 11 J=1,JMAX
        CALL TRAPZD2DG(A,B,ST,J)
        S=(4.*ST-OST)/3.
        IF (ABS(S-OS).LT.EPS*ABS(OS)) RETURN
        OS=S
        OST=ST
11    CONTINUE
ccc      PAUSE 'Too many steps.'
      END
C******************************************************************************
C Subrutina del Numerical Recipes
      SUBROUTINE TRAPZD2DG(A,B,S,N)
      it=1 !introducido por NCL para evitar un warning
      IF (N.EQ.1) THEN
        S=0.5*(B-A)*(FUNK2DG(A)+FUNK2DG(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNK2DG(X)
          X=X+DEL
11      CONTINUE
        S=0.5*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END
C
C******************************************************************************
C Funcion a integrar (aqui Y es constante, y la unica variable es X)
        REAL FUNCTION FUNK2DG(X)
        IMPLICIT NONE
        REAL X
C
        REAL Y,CTE,AMP,X0,Y0,SIGMAX,SIGMAY,BETA
        COMMON/BLKFUNK2DG/Y,CTE,AMP,X0,Y0,SIGMAX,SIGMAY,BETA
C
        REAL FACTOR1,FACTOR2,FACTOR3,FACTOR
C------------------------------------------------------------------------------
        FACTOR1=(X-X0)*(X-X0)/(2.*SIGMAX*SIGMAX)
        FACTOR2=(Y-Y0)*(Y-Y0)/(2.*SIGMAY*SIGMAY)
        FACTOR3=BETA*(X-X0)*(Y-Y0)/(SIGMAX*SIGMAY)
        FACTOR=FACTOR1+FACTOR2+FACTOR3
        IF(FACTOR.GT.60.)THEN
          FUNK2DG=CTE
        ELSE
          FUNK2DG=CTE+AMP*EXP(-FACTOR)
        END IF
        END
