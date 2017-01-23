C
C******************************************************************************
C Dibuja un histograma de la señal en el rectangulo definido por NX1,NX2,
C NY1,NY2, con limites en X entre BG y FG
        SUBROUTINE HISTOGRAM(NBUFF)
        IMPLICIT NONE
        INTEGER NBUFF
C
        INCLUDE 'dimensions.inc'
        INTEGER NBINMAX
        PARAMETER(NBINMAX=100)                 !numero de bins en el histograma
        INTEGER NDEGPOLYMAX
        PARAMETER (NDEGPOLYMAX=3) !grado maximo del polinomio ajustado al hist.
C
        INTEGER TRUEBEG
        INTEGER TRUELEN
C
        INTEGER NX1,NX2,NY1,NY2
        INTEGER I,J,K
        INTEGER NPIX(NBINMAX)
        INTEGER ISUM
        INTEGER NBIN
        INTEGER L1,L2
        INTEGER KMODA,KMAX
        INTEGER NDEGPOLY
        INTEGER ITER,ITERMAX
        REAL BG,FG
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF),DBIN,DBIN2
        REAL XPLOT(NBINMAX)
        REAL FPIX(NBINMAX)
        REAL YGAUSS(NBINMAX),AMP
        REAL YPOLY(NBINMAX),YPOLYMAX
!       REAL DYPOLY(NBINMAX),DDYPOLY(NBINMAX)
        REAL XV1,XV2,YV1,YV2
        REAL XW1,XW2,YW1,YW2
        REAL XMIN,XMAX,DX,YMIN,YMAX,DY
        REAL FMEAN,FSIGMA,FMEDIAN,FMIN,FMAX
        REAL OLD_CH
        REAL FDUM,FMODA
        REAL COEF(NDEGPOLYMAX+1),CHISQR
        REAL XX0,XX1,YFUN,YDER
        CHARACTER*50 CDUMMY
        LOGICAL LOOP
C
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKBGFG/BG,FG
        COMMON/BLKXYLIMPLOT/NX1,NX2,NY1,NY2
        COMMON/BLKESTADISTICA/FMEAN,FSIGMA,FMEDIAN,FMIN,FMAX
        COMMON/BLKHLIMITS/XMIN,XMAX,YMIN,YMAX             !local con ZOOMHISTOG
C------------------------------------------------------------------------------
ccc     IF(INT(ABS(FG-BG)).LT.NBINMAX)THEN
ccc       NBIN=INT(ABS(FG-BG))
ccc       IF(NBIN.EQ.0) NBIN=1
ccc     ELSE
          NBIN=NBINMAX
ccc     END IF
C calculamos el histograma
        DO K=1,NBIN
          NPIX(K)=0
        END DO
        DBIN=(FG-BG)/REAL(NBIN-1)                          !anchura de cada bin
        DO K=1,NBIN
          XPLOT(K)=BG+REAL(K-1)*DBIN
        END DO
C
        DO I=NY1,NY2
          DO J=NX1,NX2
            K=NINT((IMAGEN(J,I,NBUFF)-BG)/DBIN)+1
            IF((K.GE.1).AND.(K.LE.NBIN)) NPIX(K)=NPIX(K)+1
          END DO
        END DO
C calculamos los limites para dibujar el histograma
        XMIN=BG
        XMAX=FG
        DX=XMAX-XMIN
        XMIN=XMIN-DX/20.
        XMAX=XMAX+DX/20.
C
        YMIN=-1.
        YMAX=YMIN
        DO K=1,NBIN
          IF(NPIX(K).EQ.0)THEN
            FPIX(K)=-1
          ELSE
            FPIX(K)=ALOG10(REAL(NPIX(K)))
          END IF
          IF(FPIX(K).GT.YMAX) YMAX=FPIX(K)
        END DO
        IF(YMAX.EQ.0.0) YMAX=1.0
        DY=YMAX-YMIN
        IF(DY.EQ.0.0)THEN
          YMAX=1.1
        ELSE
          YMAX=YMAX+DY/20.
        END IF
C------------------------------------------------------------------------------
C almacenamos region de dibujo actual
        CALL PGQVP(0,XV1,XV2,YV1,YV2)
        CALL PGQWIN(XW1,XW2,YW1,YW2)
        CALL PGQCH(OLD_CH)
C------------------------------------------------------------------------------
        CALL PGBBUF
C borramos histograma previo
        CALL RPGERASW(0.00,0.40,0.00,0.317,0)
C definimos nueva region de dibujo para el histograma
        CALL PGSVP(0.05,0.37,0.08,0.30)
        CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
        CALL PGSCH(0.8)
        CALL PGBOX('BCTSN',0.0,0,'BCTSN',0.0,0)
        CALL PGLABEL('signal','Log[No.pixels]',' ')
        CALL PGSCI(7)
        DO K=1,NBIN
ccc       write(69,*) k,xplot(k),fpix(k)
          CALL PGMOVE(XPLOT(K)-0.5*DBIN,-1.0)
          CALL PGDRAW(XPLOT(K)-0.5*DBIN,FPIX(K))
          CALL PGDRAW(XPLOT(K)+0.5*DBIN,FPIX(K))
          CALL PGDRAW(XPLOT(K)+0.5*DBIN,-1.0)
        END DO
C dibujamos gaussiana
        IF(FSIGMA.GT.0.0)THEN
          ISUM=0
          DO K=1,NBIN
            ISUM=ISUM+NPIX(K)
          END DO
          AMP=REAL(ISUM)*DBIN/(FSIGMA*SQRT(2.*3.141593))
          DBIN2=(FG-BG)/REAL(NBIN-1)
          DO K=1,NBIN
            XPLOT(K)=BG+REAL(K-1)*DBIN2
          END DO
          DO K=1,NBIN
            FDUM=(XPLOT(K)-FMEAN)*(XPLOT(K)-FMEAN)/(2.*FSIGMA*FSIGMA)
            IF(FDUM.GT.70)THEN !evitamos underflow (IEEE)
              YGAUSS(K)=0.
            ELSE
              YGAUSS(K)=AMP*EXP(-FDUM)
            END IF
            IF(YGAUSS(K).LT.0.1)THEN
              YGAUSS(K)=-1.
            ELSE
              YGAUSS(K)=ALOG10(YGAUSS(K))
            END IF
          END DO
          CALL PGSCI(2)
          CALL PGLINE(NBIN,XPLOT,YGAUSS)
        END IF
        CALL PGSCI(1)
C mostramos BG y FG
        CALL PGSCI(3)
        WRITE(CDUMMY,*) BG
        L1=TRUEBEG(CDUMMY)
        L2=TRUELEN(CDUMMY)
        CALL PGMTXT('B',3.0,0.0,0.0,CDUMMY(L1:L2))
        WRITE(CDUMMY,*) FG
        L1=TRUEBEG(CDUMMY)
        L2=TRUELEN(CDUMMY)
        CALL PGMTXT('B',3.0,1.0,1.0,CDUMMY(L1:L2))
        CALL PGSCH(OLD_CH)
        CALL PGSCI(1)
        CALL PGEBUF
C calculamos la moda
        KMODA=1
        FMODA=FPIX(1)
        DO K=2,NBIN
          IF(FPIX(K).GT.FMODA)THEN
            FMODA=FPIX(K)
            KMODA=K
          END IF
        END DO
        WRITE(*,*)
        WRITE(*,100) '>>> Mean......: '
        WRITE(*,*) FMEAN
        WRITE(*,100) '>>> Sigma.....: '
        WRITE(*,*) FSIGMA
        WRITE(*,100) '>>> Mode......: '
        WRITE(*,*) XPLOT(KMODA)
C ajustamos un polinomio de grado alto a los datos del histograma en 
C escala logaritmica
        CALL POLFIT(XPLOT,FPIX,FPIX,NBIN,NDEGPOLYMAX+1,0,COEF,CHISQR)
        KMAX=1 !(lo pongo aqui para evitar un WARNING pero hace falta)
        YPOLYMAX=0.0 !no hace falta; evita un WARNING de compilacion
        DO K=1,NBIN
          YPOLY(K)=COEF(NDEGPOLYMAX+1)
          DO NDEGPOLY=NDEGPOLYMAX,1,-1
            YPOLY(K)=COEF(NDEGPOLY)+YPOLY(K)*XPLOT(K)
          END DO
!         DYPOLY(K)=REAL(NDEGPOLYMAX)*COEF(NDEGPOLYMAX+1)
!         DO NDEGPOLY=NDEGPOLYMAX,2,-1
!           DYPOLY(K)=REAL(NDEGPOLY-1)*COEF(NDEGPOLY)+DYPOLY(K)*XPLOT(K)
!         END DO
!         DDYPOLY(K)=REAL(NDEGPOLYMAX-1)*REAL(NDEGPOLYMAX)*
!    +     COEF(NDEGPOLYMAX+1)
!         DO NDEGPOLY=NDEGPOLYMAX,3,-1
!           DDYPOLY(K)=REAL(NDEGPOLY-2)*REAL(NDEGPOLY-1)*COEF(NDEGPOLY)+
!    +       DDYPOLY(K)*XPLOT(K)
!         END DO
          IF(K.EQ.1)THEN
            YPOLYMAX=YPOLY(1)
          ELSE
            IF(YPOLY(K).GT.YPOLYMAX)THEN
              KMAX=K
              YPOLYMAX=YPOLY(K)
            END IF
          END IF
        END DO
C usamos Newton-Raphson para calcular el máximo con precisión: ojo, calculamos
C el punto en el que la derivada primera se hace cero, por lo que en el metodo
C de Newton-Raphson aparecen las derivadas primera y segunda, y no la funcion
C y su derivada. De todas formas, para evitar confusion, en el codigo usamos
C YFUN e YDER.
        XX0=XPLOT(KMAX)
        LOOP=.TRUE.
        ITERMAX=500
        ITER=0
        DO WHILE(LOOP)
          YFUN=REAL(NDEGPOLYMAX)*COEF(NDEGPOLYMAX+1)
          DO NDEGPOLY=NDEGPOLYMAX,2,-1
            YFUN=REAL(NDEGPOLY-1)*COEF(NDEGPOLY)+YFUN*XX0
          END DO
          YDER=REAL(NDEGPOLYMAX-1)*REAL(NDEGPOLYMAX)*
     +     COEF(NDEGPOLYMAX+1)
          DO NDEGPOLY=NDEGPOLYMAX,3,-1
            YDER=REAL(NDEGPOLY-2)*REAL(NDEGPOLY-1)*COEF(NDEGPOLY)+
     +       YDER*XX0
          END DO
          XX1=XX0-YFUN/YDER
          IF(ABS(XX1-XX0).LE.1.0E-4)THEN
            LOOP=.FALSE.
          ELSE
            ITER=ITER+1
            XX0=XX1
          END IF
          IF(ITER.GE.ITERMAX)THEN
            LOOP=.FALSE.
          END IF
        END DO
        CALL PGSCI(5)
        CALL PGLINE(NBIN,XPLOT,YPOLY)
        WRITE(*,100) '>>> Max.polyn.: '
        WRITE(*,*) XX1
        CALL PGSCI(1)
C------------------------------------------------------------------------------
C recuperamos region de dibujo inicial
        CALL PGSVP(XV1,XV2,YV1,YV2)
        CALL PGSWIN(XW1,XW2,YW1,YW2)
C------------------------------------------------------------------------------
100     FORMAT(A,$)
        END
