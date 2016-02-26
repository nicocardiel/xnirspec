C Version 26-Oct-2000
C manipula imagenes (hasta un box-9)
C NFRAMES es el numero real de imagenes en el box-9
        SUBROUTINE OPERBOX9(NCBUFF)
        IMPLICIT NONE
        INTEGER NCBUFF
C
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
C
        INCLUDE 'dimensions.inc'
C
        INTEGER NBINMAX
        PARAMETER (NBINMAX=100)
        INTEGER MAXNCOEFF
        PARAMETER (MAXNCOEFF=16)
        INTEGER NOPERMAX_STACK
        PARAMETER (NOPERMAX_STACK=100)
C
        INTEGER READI,READILIM
        REAL FMEAN0,FMEDIAN1,FMEDIAN1E
        REAL READF
        REAL FTSTUDENTI,FNORTIPI
        CHARACTER*255 READC
        LOGICAL INSIDE
        LOGICAL LPOST
C
        INTEGER I,J,II,JJ,IJ,K,KK,KK_,NC,NN,NK
        INTEGER II1,II2,JJ1,JJ2
        INTEGER M,M_,L,L_
        INTEGER K_LAST(9*4),II_LAST(9*4),JJ_LAST(9*4),NPIXELSVECINOS
        INTEGER I1(4),I2(4),J1(4),J2(4)
        INTEGER NQUAD,KCOEFF
        INTEGER NM1,NM2
        INTEGER DI(NBOXMAX),DJ(NBOXMAX)
        INTEGER NFRAMES(NMAXBUFF)
        INTEGER NF,NFRAMES_
        INTEGER J0(NBOXMAX),I0(NBOXMAX)
        INTEGER FLATBUFF,NEWBUFF1,NEWBUFF2,NEWBUFF3,NEWBUFF4
        INTEGER NX1,NX2,NY1,NY2
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER NAXISFRAME(2,9,NMAXBUFF)
        INTEGER NSIZEFB9(NMAXBUFF)
        INTEGER NPIX(NBINMAX,9)
        INTEGER NP(9),SUMNP
        INTEGER NXMAXB9_,NYMAXB9_
        INTEGER NEXTRA_XMINB9,NEXTRA_XMAXB9
        INTEGER NEXTRA_YMINB9,NEXTRA_YMAXB9
        INTEGER NIMSTAT,NSUBQUA
        INTEGER NPIXREGION
        INTEGER NEXTINFO
        INTEGER GX1,GY1,GX2,GY2
        INTEGER GXFILL,GYFILL
        INTEGER NIMA,NQUA
        INTEGER ITER,NITER
        INTEGER NCOEFF
        INTEGER NOPER,NOPER_STACK,TYPEOPER_STACK(NOPERMAX_STACK)
        INTEGER NIMA_STACK(NOPERMAX_STACK)
        INTEGER NQUA_STACK(NOPERMAX_STACK)
        INTEGER GX_STACK(NOPERMAX_STACK)
        INTEGER GY_STACK(NOPERMAX_STACK)
        INTEGER IDOLD,IDNEW,NPOST
        INTEGER JUST
        INTEGER IDUM
        REAL FOFFSETX(NBOXMAX),FOFFSETY(NBOXMAX)
        REAL FOFFSETX_,FOFFSETY_
        REAL EFOFFSETX_,EFOFFSETY_
        REAL FRACCION_PIXEL
        REAL FJ0(NBOXMAX),FI0(NBOXMAX)
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL IMAGEN_(NXMAX,NYMAX)
        REAL PIXEL(9*4),EPIXEL(9*4),PIXELF(9*4)
        REAL PIXELVECINO(24)
        REAL FPIXUSED(NXMAXB9,NYMAXB9)
        REAL FMEAN,FSIGMA,FMEDIAN,FMIN,FMAX
        REAL BG,FG,TR(6)
        REAL BG_,FG_
        REAL XPLOT(NBINMAX),DBIN
        REAL XMIN,XMAX,DX
        REAL FFMEAN(9),FFSIGMA(9)
        REAL FFSIGMA1
        REAL FFSIGMANOR(9)
        REAL TIMESSIGMA(9)
        REAL FMEAN_PIX,FSIGMA_PIX
        REAL HMIN,HMAX
        REAL XP(72),YP(72) !72=36x2 (donde 36=numero de quadrantes)
        REAL FPIXELTHRESHOLD
        REAL FPIXELTHRESHOLDCR
        REAL FTAILS
        REAL TSHOT
        REAL YRMSTOL1,YRMSTOL2,YRMSTOLFILL
        REAL QMIN(36),QMAX(36)
        REAL FMEDIAN_STACK(36,NOPERMAX_STACK)
        REAL FCOEFF_STACK(MAXNCOEFF,36,NOPERMAX_STACK)
        REAL FFACTOR,X(128),Y(128),XSOL(MAXNCOEFF,36)
        REAL QMIN_STACK(36,NOPERMAX_STACK),QMAX_STACK(36,NOPERMAX_STACK)
        REAL XW1,XW2,YW1,YW2
        CHARACTER*1 CINTERACTIVE,CCONT,CREJECT,CERRSKY
        CHARACTER*1 CMEDIANA1,CMEDIANA2,CREFINE1,CREFINE2
        CHARACTER*1 CMEDIANAFILL,CREFINEFILL
        CHARACTER*1 CSCALE1,CSCALE2
        CHARACTER*1 CSCALE1OPER,CSCALE2OPER
        CHARACTER*1 CNOR1,CNOR2
        CHARACTER*1 CMORE,CPAUSE
        CHARACTER*1 CFILL
        CHARACTER*1 CNOR_STACK(NOPERMAX_STACK),CSTACK
        CHARACTER*1 CSCALEOPER_STACK(NOPERMAX_STACK)
        CHARACTER*1 CPOST
        CHARACTER*50 CDUMMY
        CHARACTER*50 CBASEPOST
        CHARACTER*80 COFFSETFILE
        LOGICAL LOGFILE
        LOGICAL MASKBOX9(NXMAXB9,NYMAXB9)
        LOGICAL MASKBOX9_(NXMAXB9,NYMAXB9)
        LOGICAL SUBMASKBOX9(256,256,9)
        LOGICAL LANY,LANYCR
        LOGICAL LASK
        LOGICAL LECHO
C
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKIMAGEN1_/IMAGEN_
        COMMON/BLKNAXIS/NAXIS                                      !dimensiones
        COMMON/BLKNAXISFRAME/NAXISFRAME
        COMMON/BLKNSIZEFB9/NSIZEFB9
        COMMON/BLKNFRAMES/NFRAMES
        COMMON/BLKBGFG/BG,FG
        COMMON/BLKESTADISTICA/FMEAN,FSIGMA,FMEDIAN,FMIN,FMAX
        COMMON/BLKXYLIMPLOT/NX1,NX2,NY1,NY2
        COMMON/BLKMASKBOX9/MASKBOX9
        COMMON/BLKFPIXUSED/FPIXUSED
        COMMON/BLKSUBMASKBOX9/SUBMASKBOX9
        COMMON/BLKJUST/JUST
        COMMON/BLKLECHO/LECHO
        COMMON/BLKB9POST1/CPOST,CBASEPOST
        COMMON/BLKB9POST2/NPOST,IDOLD
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
C Note: the old numeration of the frames was the following:
C        2  5  8
C        7  1  3
C        4  9  6
ccc     +   256,256,  !1
ccc     +   512,000,  !2
ccc     +   256,512,  !3
ccc     +   000,000,  !4
ccc     +   512,256,  !5
ccc     +   000,512,  !6
ccc     +   256,000,  !7
ccc     +   512,512,  !8
ccc     +   000,256/  !9
C
C limites en pixels de cada cuadrante
        DATA (I1(NQUAD),I2(NQUAD),J1(NQUAD),J2(NQUAD),NQUAD=1,4) /
     +   001,128,001,128,
     +   129,256,001,128,
     +   001,128,129,256,
     +   129,256,129,256/
C
        TR(1)=0.
        TR(2)=1.
        TR(3)=0.
        TR(4)=0.
        TR(5)=0.
        TR(6)=1.
C
        II2=0 !evita un WARNING de compilacion
        JJ2=0 !evita un WARNING de compilacion
        FRACCION_PIXEL=0.0 !evita un WARNING de compilacion
C
        FLATBUFF=0
C
        FPIXELTHRESHOLD=1.0
        FPIXELTHRESHOLDCR=10.0
        FTAILS=0.5
        CREJECT='n'
        CSCALE1='n'
        CSCALE1OPER='s'
        CNOR1='n'
        CSCALE2='n'
        CSCALE2OPER='s'
        CNOR2='n'
        CFILL='n'
        GX1=1
        GY1=1
        GX2=1
        GY2=1
        GXFILL=1
        GYFILL=1
        CMEDIANA1='y'
        CMEDIANA2='y'
        CMEDIANAFILL='n'
        CERRSKY='n'
        YRMSTOL1=1.E-3
        YRMSTOL2=1.E-3
        YRMSTOLFILL=1.E-3
        CREFINE1='y'
        CREFINE2='y'
        CREFINEFILL='n'
        CMORE='n'
        CINTERACTIVE='y'
        NM2=0
        NPOST=0
        CBASEPOST='pgplot'
C------------------------------------------------------------------------------
C La subrutina almacena las operaciones realizadas en un STACK de operaciones.
C Como hay varias operaciones posibles (ajustar la mediana, calculo de
C superficie), cada cual con un numero de parametros diferente, tenemos que
C almacenar esos parametros en funcion del numero de operacion. He declarado un
C stack para "guardar" dichos parametros, pero por ahorrar un par de variables,
C el stack de cada tipo de operacion puede contener huecos. Esto facilita
C tambien un poco la programacion.
        NOPER_STACK=0
        CSTACK='n'
C------------------------------------------------------------------------------
C pedimos offsets
        LOGFILE=.FALSE.
        DO WHILE(.NOT.LOGFILE)
          COFFSETFILE(1:80)=
     +     READC('Name of the file with offsets (none=exit)',
     +     'offsets.dat','@')
          IF(COFFSETFILE.EQ.'none') RETURN !permitimos escapar
          INQUIRE(FILE=COFFSETFILE,EXIST=LOGFILE)
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
          WRITE(*,101) 'ERROR: file with offsets contains too many'//
     +     'lines (>9)'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
        GOTO 5
6        CLOSE(20)
        NFRAMES_=NF
        WRITE(*,101) '...OK! File read and closed.'
C proteccion
        IF(NFRAMES_.NE.NFRAMES(NCBUFF))THEN
          WRITE(*,100) 'NFRAMES(NCBUFF),NFRAMES_: '
          WRITE(*,*) NFRAMES(NCBUFF),NFRAMES_
          WRITE(*,100) 'ERROR: number of frames does not correspond'
          WRITE(*,101) ' with expected value.'
          CCONT(1:1)=
     +     READC('Do you want to continue anyway (y/n)','n','yn')
          IF(CCONT.EQ.'n') RETURN
        END IF
C estimamos partes enteras y fraccionarias de los offsets, y determinamos la
C dimension maxima que vamos a necesitar para la imagen combinada
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
        WRITE(*,*) NEXTRA_XMINB9,NEXTRA_XMAXB9,
     +   NEXTRA_YMINB9,NEXTRA_YMAXB9
        NXMAXB9_=256-NEXTRA_XMINB9+NEXTRA_XMAXB9
        NYMAXB9_=256-NEXTRA_YMINB9+NEXTRA_YMAXB9
C mostramos taman~o
        WRITE(*,100) '>>> New dimensions of combined frame: '
        WRITE(CDUMMY,'(I10,A1,I10)') NXMAXB9_,'x',NYMAXB9_
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101) CDUMMY(1:L)
        NPIXREGION=NXMAXB9_*NYMAXB9_
C proteccion
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
C------------------------------------------------------------------------------
C pedimos el numero de pixels umbral para el cual vamos a calcular el numero 
C de veces sigma (Times_Sigma) tal que en cada region ese sea el numero de 
C pixels que, estadisticamente, exceden el valor Times_Sigma.
        WRITE(CDUMMY,*) FPIXELTHRESHOLD
        FPIXELTHRESHOLD=READF(
     +   'Threshold in no. of pixels > Times_Sigma for mask',CDUMMY)
C factor para aplicar FPIXELTHRESHOLD a la busqueda de las regiones mas
C debiles alrededor de los objetos que claramente hay que "tapar" con la
C mascara
        WRITE(CDUMMY,'(F8.6)') FTAILS
        FTAILS=READF(
     +   'Factor to search for extended mask pixels (0 <= factor <= 1)',
     +    CDUMMY)
C------------------------------------------------------------------------------
C definimos los buffers de trabajo
        WRITE(CDUMMY,*) FLATBUFF
        FLATBUFF=READILIM(
     +   'Enter buffer # with flatfield image, 0=no flatfield..',
     +   CDUMMY,0,NMAXBUFF/2)
        IF(FLATBUFF.GT.0)THEN
          IF((NAXIS(1,FLATBUFF).NE.256).OR.
     +       (NAXIS(2,FLATBUFF).NE.256))THEN
            WRITE(*,100) 'ERROR: image dimensions of flatfield buffer: '
            WRITE(*,*) NAXIS(1,FLATBUFF),NAXIS(2,FLATBUFF)
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
            IF(LECHO) WRITE(*,*)
            RETURN
          END IF
        END IF
c
        NEWBUFF1=NCBUFF+1
        IF(NEWBUFF1.GT.NMAXBUFF/2) NEWBUFF1=1
        WRITE(CDUMMY,*) NEWBUFF1
        NEWBUFF1=READILIM(
     +   'First  auxiliary buffer # to store manipulated data..',
     +   CDUMMY,1,NMAXBUFF/2)
c
        NEWBUFF2=NEWBUFF1+1
        IF(NEWBUFF2.GT.NMAXBUFF/2) NEWBUFF2=1
        WRITE(CDUMMY,*) NEWBUFF2
        NEWBUFF2=READILIM(
     +   'Second auxiliary buffer # to store manipulated data..',
     +   CDUMMY,1,NMAXBUFF/2)
c
        NEWBUFF3=NEWBUFF2+1
        IF(NEWBUFF3.GT.NMAXBUFF/2) NEWBUFF3=1
        WRITE(CDUMMY,*) NEWBUFF3
        NEWBUFF3=READILIM(
     +   'Third auxiliary buffer # to store manipulated data...',
     +   CDUMMY,1,NMAXBUFF/2)
c
        NEWBUFF4=NEWBUFF3+1
        IF(NEWBUFF4.GT.NMAXBUFF/2) NEWBUFF4=1
        WRITE(CDUMMY,*) NEWBUFF4
        NEWBUFF4=READILIM(
     +   'Fourth auxiliary buffer # to store manipulated data..',
     +   CDUMMY,1,NMAXBUFF/2)
C------------------------------------------------------------------------------
        CINTERACTIVE(1:1)=
     +   READC('Run procedure in interactive mode (y/m/n)',
     +   CINTERACTIVE,'ynm')
        IF(CINTERACTIVE.NE.'y')THEN
          CSCALE1(1:1)=
     +     READC('Fit#1 [c]onstant, [s]urface or [n]one (c/s/n)',
     +     CSCALE1,'csn')
          IF(CSCALE1.EQ.'s')THEN
            WRITE(CDUMMY,*) GX1
            GX1=READILIM('Polynomial degree in X direction',CDUMMY,
     +       0,INT(SQRT(REAL(MAXNCOEFF)))-1)
            WRITE(CDUMMY,*) GY1
            GY1=READILIM('Polynomial degree in Y direction',CDUMMY,
     +       0,INT(SQRT(REAL(MAXNCOEFF)))-1)
            CMEDIANA1(1:1)=READC('Refine fit minimizing median (y/n)',
     +       CMEDIANA1,'yn')
            IF(CMEDIANA1.EQ.'y')THEN
              WRITE(CDUMMY,*) YRMSTOL1
              YRMSTOL1=READF('YRMSTOL for DOWNHILL',CDUMMY)
            END IF
            CREFINE1(1:1)=
     +       READC('Refine fit removing deviating points (y/n)',
     +       CREFINE1,'yn')
          END IF
          IF(CSCALE1.NE.'n')THEN
            CSCALE1OPER(1:1)=
     +       READC('[s]ubtract fit or [d]ivide by fit (s/d)',
     +       CSCALE1OPER,'sd')
            IF(CSCALE1OPER.EQ.'d')THEN
              WRITE(*,101) 'WARNING: normalization is recomended to '//
     +         'preserve number of counts'
              CNOR1(1:1)=READC('Normalize fit (y/n)','y','yn')
            ELSE
              CNOR1='n'
            END IF
          END IF
c
          CSCALE2(1:1)=
     +     READC('Fit#2 [c]onstant, [s]urface or [n]one (c/s/n)',
     +     CSCALE2,'csn')
          IF(CSCALE2.EQ.'s')THEN
            WRITE(CDUMMY,*) GX2
            GX2=READILIM('Polynomial degree in X direction',CDUMMY,
     +       0,INT(SQRT(REAL(MAXNCOEFF)))-1)
            WRITE(CDUMMY,*) GY2
            GY2=READILIM('Polynomial degree in Y direction',CDUMMY,
     +       0,INT(SQRT(REAL(MAXNCOEFF)))-1)
            CMEDIANA2(1:1)=READC('Refine fit minimizing median (y/n)',
     +       CMEDIANA2,'yn')
            IF(CMEDIANA2.EQ.'y')THEN
              WRITE(CDUMMY,*) YRMSTOL2
              YRMSTOL2=READF('YRMSTOL for DOWNHILL',CDUMMY)
            END IF
            CREFINE2(1:1)=
     +       READC('Refine fit removing deviating points (y/n)',
     +       CREFINE2,'yn')
          END IF
          IF(CSCALE2.NE.'n')THEN
            CSCALE2OPER(1:1)=
     +       READC('[s]ubtract fit or [d]ivide by fit (s/d)',
     +       CSCALE2OPER,'sd')
            IF(CSCALE2OPER.EQ.'d')THEN
              WRITE(*,101) 'WARNING: normalization is recomended to '//
     +         'preserve number of counts'
              CNOR2(1:1)=READC('Normalize fit (y/n)','y','yn')
            ELSE
              CNOR2='n'
            END IF
          END IF
C determinamos si en el calculo del cielo rellenamos las regines tapadas por
C la mascara
          CFILL(1:1)=READC('Fill masked regions with 2D polynomial '//
     +     'fits (y/n)',CFILL,'yn')
          IF(CFILL.EQ.'y')THEN
            WRITE(CDUMMY,*) GXFILL
            GXFILL=READILIM('Polynomial degree in X direction',CDUMMY,
     +       0,INT(SQRT(REAL(MAXNCOEFF)))-1)
            WRITE(CDUMMY,*) GYFILL
            GYFILL=READILIM('Polynomial degree in Y direction',CDUMMY,
     +       0,INT(SQRT(REAL(MAXNCOEFF)))-1)
            CMEDIANAFILL(1:1)=
     +       READC('Refine fit minimizing median (y/n)',
     +       CMEDIANAFILL,'yn')
            IF(CMEDIANAFILL.EQ.'y')THEN
              WRITE(CDUMMY,*) YRMSTOLFILL
              YRMSTOLFILL=READF('YRMSTOL for DOWNHILL',CDUMMY)
            END IF
            CREFINEFILL(1:1)=
     +       READC('Refine fit removing deviating points'//
     +       ' (y/n)',CREFINEFILL,'yn')
          END IF
C decidimos si calculamos errores en el calculo del cielo
          CERRSKY(1:1)=READC('Compute error in median sky (y/n)',
     +     CERRSKY,'yn')
C decidimos si eliminamos pixels extran~os
          IF(CINTERACTIVE.EQ.'n')THEN
            CREJECT(1:1)=READC('Reject hot/cool pixels (y/n)',
     +       CREJECT,'yn')
            IF(CREJECT.EQ.'y')THEN
C numero de pixels umbral para deteccion de pixels calientes
              WRITE(CDUMMY,*) FPIXELTHRESHOLDCR
              FPIXELTHRESHOLDCR=READF('Threshold in no. of pixels >'//
     +         ' Times_Sigma for hot/cool pixels',CDUMMY)
            END IF
          END IF
        END IF
C------------------------------------------------------------------------------
C inicializamos la mascara (inicialmente usamos todos los pixels)
        NX1=1
        NX2=NXMAXB9_
        NY1=1
        NY2=NYMAXB9_
        DO I=NY1,NY2
          DO J=NX1,NX2
            MASKBOX9(J,I)=.TRUE.
          END DO
        END DO
C inicializamos el valor de FFSIGMANOR, que en el primer paso (antes de ser
C calculada) intenta ser utilizada por el programa para detectar hot/cool pixels
        DO NC=1,NFRAMES_
          FFSIGMANOR(NC)=0.
          TIMESSIGMA(NC)=0.
        END DO
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C aqui comienza la parte de la subrutina que se emplea recursivamente al iterar
10      CONTINUE
C------------------------------------------------------------------------------
        WRITE(*,*)
        CPOST(1:1)=READC('Create PostScript files (y/n)','n','yn')
ccc        CPOST='n'
        IF(CPOST.EQ.'y')THEN
          CBASEPOST(1:50)=READC('Base for postscript file names',
     +     CBASEPOST,'@')
          WRITE(CDUMMY,*) NPOST+1
          NPOST=READI('Number of first pgplot????.ps file',CDUMMY)
          NPOST=NPOST-1
          CALL PGQID(IDOLD)
        END IF
C------------------------------------------------------------------------------
C a partir de la ultima mascara generada, calculamos las submascaras para cada
C frame que constituye el box 9
        CALL GIVESUBMASK(NFRAMES_,NEXTRA_XMINB9,NEXTRA_YMINB9,
     +   I0,J0,FI0,FJ0,NXMAXB9_,NYMAXB9_)
C Duplicamos buffer para no modificar NCBUFF
        NAXIS(1,NEWBUFF1)=768
        NAXIS(2,NEWBUFF1)=768
        NFRAMES(NEWBUFF1)=NFRAMES(NCBUFF)
        DO I=1,NAXIS(2,NEWBUFF1)
          DO J=1,NAXIS(1,NEWBUFF1)
            IMAGEN(J,I,NEWBUFF1)=IMAGEN(J,I,NCBUFF)
          END DO
        END DO
        NAXIS(1,NEWBUFF1+NMAXBUFF/2)=768
        NAXIS(2,NEWBUFF1+NMAXBUFF/2)=768
        NFRAMES(NEWBUFF1+NMAXBUFF/2)=NFRAMES(NCBUFF+NMAXBUFF/2)
        DO I=1,NAXIS(2,NEWBUFF1+NMAXBUFF/2)
          DO J=1,NAXIS(1,NEWBUFF1+NMAXBUFF/2)
            IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=
     +       IMAGEN(J,I,NCBUFF+NMAXBUFF/2)
          END DO
        END DO
C mas cosas a duplicar por si medimos offsets (ahora o luego)
        DO K=1,NFRAMES_
          NAXISFRAME(1,K,NEWBUFF1)=NAXISFRAME(1,K,NCBUFF)
          NAXISFRAME(2,K,NEWBUFF1)=NAXISFRAME(2,K,NCBUFF)
        END DO
        NSIZEFB9(NEWBUFF1)=NSIZEFB9(NCBUFF)
C flatfielding (if required)
        IF(FLATBUFF.GT.0)THEN
          DO I=1,NAXIS(2,NEWBUFF1)
            II=MOD(I,256)
            IF(II.EQ.0) II=256
            DO J=1,NAXIS(1,NEWBUFF1)
              JJ=MOD(J,256)
              IF(JJ.EQ.0) JJ=256
              IF(IMAGEN(JJ,II,FLATBUFF).GT.0.0)THEN !ignore <=0
                IMAGEN(J,I,NEWBUFF1)=IMAGEN(J,I,NEWBUFF1)/
     +           IMAGEN(JJ,II,FLATBUFF)
                IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=SQRT(
     +           IMAGEN(J,I,NEWBUFF1)*IMAGEN(J,I,NEWBUFF1)*
     +           IMAGEN(JJ,II,FLATBUFF+NMAXBUFF/2)*
     +           IMAGEN(JJ,II,FLATBUFF+NMAXBUFF/2)+
     +           IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)*
     +           IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)*
     +           IMAGEN(JJ,II,FLATBUFF)*IMAGEN(JJ,II,FLATBUFF))/
     +           (IMAGEN(JJ,II,FLATBUFF)*IMAGEN(JJ,II,FLATBUFF))
              END IF
            END DO
          END DO
        END IF
C mostramos imagen a utilizar
        NX1=1
        NX2=768
        NY1=1
        NY2=768
        CALL STATISTICB9(NFRAMES_,NEWBUFF1,FMEAN,FSIGMA)
        IF(FSIGMA.GT.0.0)THEN
          BG=FMEAN-5.*FSIGMA
          FG=FMEAN+5.*FSIGMA
        ELSE
          BG=FMEAN-1.0
          FG=FMEAN+1.0
        END IF
        CALL HISTOGRAM(NEWBUFF1)
        CALL SUBLOOK(.FALSE.,NEWBUFF1,.FALSE.)
        CALL DRAWSEPB9(0,0,256)
        IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
          CALL SUBLOOK(.FALSE.,NEWBUFF1,.TRUE.)
          CALL DRAWSEPB9(0,0,256)
          CALL PGCLOS(IDNEW)
          CALL PGSLCT(IDOLD)
        END IF
C------------------------------------------------------------------------------
C Si tenemos operaciones almacenadas en el stack, preguntamos si queremos
C aplicarlas antes de continuar
        WRITE(*,*)
        WRITE(*,100) '>>> Total no. of operations in stack: '
        WRITE(*,*) NOPER_STACK
        IF(NOPER_STACK.GT.0)THEN
          IF(CINTERACTIVE.EQ.'y')THEN
            CSTACK(1:1)=READC('Apply operations in stack (y/n)',
     +       CSTACK,'yn')
          END IF
          IF(CSTACK.EQ.'y')THEN
            CPAUSE(1:1)=READC('Pause between intermediate plots (y/n)',
     +       'n','yn')
            WRITE(*,*)
            WRITE(*,101) 'Applying operations in stack...'
C..............................................................................
            DO NOPER=1,NOPER_STACK
c
              WRITE(*,*)
              WRITE(*,100) '>>> Operation #'
              WRITE(*,*) NOPER
              WRITE(*,100) '>>> NIMA, NQUA: '
              WRITE(*,*) NIMA_STACK(NOPER),NQUA_STACK(NOPER)
              WRITE(*,100) '>>> Type of operation: '
              WRITE(*,*) TYPEOPER_STACK(NOPER)
              WRITE(*,100) '>>> CSCALE, CNOR: '
              WRITE(*,100) CSCALEOPER_STACK(NOPER)
              WRITE(*,100) ', '
              WRITE(*,101) CNOR_STACK(NOPER)
c
              IF((NIMA_STACK(NOPER).NE.0).OR.
     +         (NQUA_STACK(NOPER).NE.0).OR.
     +         (NFRAMES_.NE.NBOXMAX))THEN
                DO K=1,NBOXMAX
                  DO NQUAD=1,4
                    DO I=I1(NQUAD),I2(NQUAD)
                      DO J=J1(NQUAD),J2(NQUAD)
                        IMAGEN_(J+DJ(K),I+DI(K))=0.0
                      END DO
                    END DO
                  END DO
                END DO
              END IF
c
              IF(TYPEOPER_STACK(NOPER).EQ.1)THEN !mediana
                DO K=1,NFRAMES_
                  IF((NIMA_STACK(NOPER).EQ.0).OR.
     +             (NIMA_STACK(NOPER).EQ.K))THEN
                    DO NQUAD=1,4
                      IF((NQUA_STACK(NOPER).EQ.0).OR.
     +                 (NQUA_STACK(NOPER).EQ.NQUAD))THEN
                        DO I=I1(NQUAD),I2(NQUAD)
                          DO J=J1(NQUAD),J2(NQUAD)
                            IMAGEN_(J+DJ(K),I+DI(K))=
     +                       FMEDIAN_STACK((K-1)*4+NQUAD,NOPER)
                          END DO
                        END DO
                      END IF
                    END DO
                  END IF
                END DO
              ELSEIF(TYPEOPER_STACK(NOPER).EQ.2)THEN !superficie polinomica
                DO J=1,128
                  X(J)=(REAL(J)-64.5)/63.5
                END DO
                DO I=1,128
                  Y(I)=(REAL(I)-64.5)/63.5
                END DO
                DO K=1,NFRAMES_
                  IF((NIMA_STACK(NOPER).EQ.0).OR.
     +             (NIMA_STACK(NOPER).EQ.K))THEN
                    DO NQUAD=1,4
                      IF((NQUA_STACK(NOPER).EQ.0).OR.
     +                 (NQUA_STACK(NOPER).EQ.NQUAD))THEN
                        DO M=I1(NQUAD),I2(NQUAD)
                          M_=M-I1(NQUAD)+1
                          DO L=J1(NQUAD),J2(NQUAD)
                            L_=L-J1(NQUAD)+1
                            IMAGEN_(L+DJ(K),M+DI(K))=0.
                            KK=0
                            DO I=0,GX_STACK(NOPER)
                              DO J=0,GY_STACK(NOPER)
                                KK=KK+1
                                FFACTOR=
     +                           FCOEFF_STACK(KK,(K-1)*4+NQUAD,NOPER)
                                IF(I.NE.0) FFACTOR=FFACTOR*(X(L_)**(I))
                                IF(J.NE.0) FFACTOR=FFACTOR*(Y(M_)**(J))
                                IMAGEN_(L+DJ(K),M+DI(K))=
     +                           IMAGEN_(L+DJ(K),M+DI(K))+FFACTOR
                              END DO
                                END DO
                          END DO
                        END DO
                      END IF
                    END DO
                  END IF
                END DO
              ELSE !fatal error
                WRITE(*,100) 'NOPER, TYPEOPER_STACK: '
                WRITE(*,*) NOPER,TYPEOPER_STACK(NOPER)
                WRITE(*,101) 'FATAL ERROR: invalid operation type.'
                STOP
              END IF
c
              IF(CPAUSE.EQ.'y')THEN
                WRITE(*,100) 'Press <CR> to display fits...'
                READ(*,*)
                IF(LECHO) WRITE(*,*)
              ELSE
                WRITE(*,101) 'Displaying fits...'
              END IF
              CALL GIVEME_BGFG(NFRAMES_,NIMA_STACK(NOPER),
     +         NQUA_STACK(NOPER),BG_,FG_)
              BG=BG_
              FG=FG_
              CALL HISTOGRAM(NEWBUFF1)
              CALL PGIMAG(IMAGEN_,NXMAX,NYMAX,NX1,NX2,NY1,NY2,FG,BG,TR)
              CALL DRAWSEPB9(NIMA_STACK(NOPER),NQUA_STACK(NOPER),256)
              DO I=1,4*NFRAMES_
                XP(2*(I-1)+1)=REAL(I)/4.
                XP(2*I)=REAL(I)/4.
                YP(2*(I-1)+1)=QMIN_STACK(I,NOPER)
                YP(2*I)=QMAX_STACK(I,NOPER)
              END DO
              CALL SUBPLOT(8*NFRAMES_,1,8*NFRAMES_,XP,YP,XP,YP,
     +                     .TRUE.,.TRUE.,.FALSE.,.FALSE.,
     +                     'frame #','signal',' ',
     +                     1,-1,1.0)
              DO I=1,4*NFRAMES_
                IF(MOD((I-1)/4,2).EQ.0)THEN
                  CALL PGSCI(3)
                ELSE
                  CALL PGSCI(5)
                END IF
                CALL SUBPLOTBIS(2,(I-1)*2+1,I*2,XP,YP,XP,YP,
     +           .FALSE.,.FALSE.,MOD(I-1,4)+2,17,1.0)
                CALL SUBPLOTBIS(2,(I-1)*2+1,I*2,XP,YP,XP,YP,
     +           .FALSE.,.FALSE.,MOD(I-1,4)+2,101,1.0)
                WRITE(*,'(A,I2,A,$)') '> Stack #',NOPER,' (residuals): '
                WRITE(*,*) I,YP((I-1)*2+1),YP(2*I)
              END DO
c
              IF(CNOR_STACK(NOPER).EQ.'y')THEN
                IF(CSCALEOPER_STACK(NOPER).EQ.'s')THEN
                  CALL NORMSUMA(NFRAMES_,NIMA_STACK(NOPER),
     +             NQUA_STACK(NOPER))
                ELSE
                  CALL NORMDIVI(NFRAMES_,NIMA_STACK(NOPER),
     +             NQUA_STACK(NOPER))
                END IF
              END IF
              IF(CSCALEOPER_STACK(NOPER).EQ.'s')THEN
                DO I=NY1,NY2
                  DO J=NX1,NX2
                    IMAGEN(J,I,NEWBUFF1)=
     +               IMAGEN(J,I,NEWBUFF1)-IMAGEN_(J,I)
                  END DO
                END DO
              ELSE
                DO I=NY1,NY2
                  DO J=NX1,NX2
                    IMAGEN(J,I,NEWBUFF1)=
     +               IMAGEN(J,I,NEWBUFF1)/IMAGEN_(J,I)
                    IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=
     +               IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)/IMAGEN_(J,I)
                  END DO
                END DO
              END IF
            END DO
C..............................................................................
            IF(CINTERACTIVE.EQ.'y')THEN
              WRITE(*,100) 'Press <CR> to plot image after applying'//
     +         ' operations in stack...'
              READ(*,*)
              IF(LECHO) WRITE(*,*)
            ELSE
              WRITE(*,101) 'Plotting image after applying operations'//
     +         ' in stack...'
            END IF
            CALL STATISTICB9(NFRAMES_,NEWBUFF1,FMEAN,FSIGMA)
            IF(FSIGMA.GT.0.0)THEN
              BG=FMEAN-5.*FSIGMA
              FG=FMEAN+5.*FSIGMA
            ELSE
              BG=FMEAN-1.0
              FG=FMEAN+1.0
            END IF
            CALL HISTOGRAM(NEWBUFF1)
            CALL SUBLOOK(.FALSE.,NEWBUFF1,.FALSE.)
            CALL DRAWSEPB9(0,0,256)
            IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
              CALL SUBLOOK(.FALSE.,NEWBUFF1,.TRUE.)
                CALL DRAWSEPB9(0,0,256)
              CALL PGCLOS(IDNEW)
              CALL PGSLCT(IDOLD)
            END IF
C..............................................................................
          ELSE
            NOPER_STACK=0
          END IF
C
        END IF
C------------------------------------------------------------------------------
C superponemos mascara
        IF(CINTERACTIVE.EQ.'y')THEN
          WRITE(*,100) 'Press <CR> to overplot mask...'
          READ(*,*)
          IF(LECHO) WRITE(*,*)
        ELSE
          WRITE(*,101) 'Overplotting mask...'
        END IF
        CALL PGBBUF
        CALL PGSCI(2)
        DO K=1,NFRAMES_
          DO I=1,256
            DO J=1,256
              IF(.NOT.SUBMASKBOX9(J,I,K))THEN
                CALL PGRECT(REAL(J+DJ(K))-0.5,REAL(J+DJ(K))+0.5,
     +                      REAL(I+DI(K))-0.5,REAL(I+DI(K))+0.5)
              END IF
            END DO
          END DO
        END DO
        CALL PGSCI(1)
        CALL PGEBUF
        IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
          CALL SUBLOOK(.FALSE.,NEWBUFF1,.TRUE.)
          CALL PGSCI(2)
          DO K=1,NFRAMES_
            DO I=1,256
              DO J=1,256
                IF(.NOT.SUBMASKBOX9(J,I,K))THEN
                  CALL PGRECT(REAL(J+DJ(K))-0.5,REAL(J+DJ(K))+0.5,
     +                        REAL(I+DI(K))-0.5,REAL(I+DI(K))+0.5)
                END IF
              END DO
            END DO
          END DO
          CALL PGSCI(1)
            CALL DRAWSEPB9(0,0,256)
          CALL PGCLOS(IDNEW)
          CALL PGSLCT(IDOLD)
        END IF
C------------------------------------------------------------------------------
C Estudio estadistico del ruido frente a la sen~al en cada subcuadrante.
        IF(CINTERACTIVE.EQ.'y')THEN
          CCONT(1:1)=READC('Perform a detailed statistical'//
     +     ' analysis of any subquadrant (y/n)','n','yn')
          IF(CCONT.EQ.'y')THEN
            IF(CCONT.EQ.'y')THEN
              NIMSTAT=10
              NSUBQUA=0
              DO WHILE(NIMSTAT.NE.0)
                IF(NIMSTAT.LE.9)THEN
                  WRITE(CDUMMY,*) NIMSTAT
                ELSE
                  CDUMMY='@'
                END IF
                WRITE(*,101) ' 6  9  4'
                WRITE(*,101) ' 3  1  7'
                WRITE(*,101) ' 8  5  2'
                NIMSTAT=READILIM('Image number (0=exit)',CDUMMY,
     +           0,NFRAMES_)
                IF(NIMSTAT.NE.0) CALL DETAILB9(NEWBUFF1,NCBUFF,
     +           NIMSTAT)
              END DO
            END IF
          END IF
        END IF
C------------------------------------------------------------------------------
C Scale quadrants
        IF(CINTERACTIVE.EQ.'y')THEN
          IF(CSCALE1.NE.'n')THEN
            CCONT(1:1)=
     +       READC('Do you want to scale each quadrant (y/n/x)',
     +       'y','ynx')
          ELSE
            CCONT(1:1)=
     +       READC('Do you want to scale each quadrant (y/n/x)',
     +       'n','ynx')
          END IF
          IF(CCONT.EQ.'x')THEN
            CPOST='n'
            RETURN
          END IF
          IF(CCONT.EQ.'n') CSCALE1='n'
        ELSE
          CCONT='y'
          CPAUSE='n'
        END IF
C------------------------------------------------------------------------------
        IF(CCONT.EQ.'y')THEN
          NIMA=0
          NQUA=0
          NITER=1
          LASK=.TRUE.
          DO WHILE(LASK)
            IF(CINTERACTIVE.EQ.'y')THEN
              CSCALE1(1:1)=
     +         READC('Fit#1 [c]onstant, [s]urface or [n]one '//
     +         '(c/s/n)',CSCALE1,'csn')
              IF(CSCALE1.EQ.'s')THEN
                WRITE(CDUMMY,*) GX1
                GX1=READILIM('Polynomial degree in X direction',CDUMMY,
     +           0,INT(SQRT(REAL(MAXNCOEFF)))-1)
                WRITE(CDUMMY,*) GY1
                GY1=READILIM('Polynomial degree in Y direction',CDUMMY,
     +           0,INT(SQRT(REAL(MAXNCOEFF)))-1)
                CMEDIANA1(1:1)=
     +           READC('Refine fit minimizing median (y/n)',
     +           CMEDIANA1,'yn')
                IF(CMEDIANA1.EQ.'y')THEN
                  WRITE(CDUMMY,*) YRMSTOL1
                  YRMSTOL1=READF('YRMSTOL for DOWNHILL',CDUMMY)
                END IF
                CREFINE1(1:1)=
     +           READC('Refine fit removing deviating points'//
     +           ' (y/n)',CREFINE1,'yn')
              END IF
              IF(CSCALE1.NE.'n')THEN
                CSCALE1OPER(1:1)=
     +           READC('[s]ubtract fit or [d]ivide by fit (s/d)',
     +           CSCALE1OPER,'sd')
                IF(CSCALE1OPER.EQ.'d')THEN
                  WRITE(*,101) 'WARNING: normalization is recomended '//
     +             'to preserve number of counts'
                  CNOR1(1:1)=READC('Normalize fit (y/n)','y','yn')
                ELSE
                  CNOR1='n'
                END IF
                WRITE(*,101) ' 6  9  4'
                WRITE(*,101) ' 3  1  7'
                WRITE(*,101) ' 8  5  2'
                WRITE(CDUMMY,*) NIMA
                NIMA=READILIM('Image number (0=all)',CDUMMY,0,NFRAMES_)
                WRITE(*,101) ' 2  4'
                WRITE(*,101) ' 1  3'
                WRITE(CDUMMY,*) NQUA
                NQUA=READILIM('Quadrant (0=all)',CDUMMY,0,4)
                WRITE(CDUMMY,*) NITER
                NITER=READILIM('No. of iterations',CDUMMY,1,10)
                IF(CINTERACTIVE.EQ.'y')THEN
                  CPAUSE(1:1)=
     +             READC('Pause between intermediate plots (y/n)',
     +             'n','yn')
                ELSE
                  CPAUSE='n'
                END IF
                WRITE(*,101) 'Working...'
              END IF
            END IF
C..............................................................................
            IF(CSCALE1.NE.'n')THEN
              DO ITER=1,NITER
                NOPER_STACK=NOPER_STACK+1
                WRITE(*,100) '>>> Stacking operation #'
                WRITE(*,*) NOPER_STACK
                NIMA_STACK(NOPER_STACK)=NIMA
                NQUA_STACK(NOPER_STACK)=NQUA
                CNOR_STACK(NOPER_STACK)=CNOR1
                CSCALEOPER_STACK(NOPER_STACK)=CSCALE1OPER
                IF(CSCALE1.EQ.'c')THEN
                  TYPEOPER_STACK(NOPER_STACK)=1
                  CALL SCALE_CONSTANT(NEWBUFF1,NFRAMES_,
     +             NIMA,NQUA,QMIN,QMAX)
                  DO NQUAD=1,4*NFRAMES_
                    FMEDIAN_STACK(NQUAD,NOPER_STACK)=QMIN(NQUAD)
                  END DO
                ELSE
                  TYPEOPER_STACK(NOPER_STACK)=2
                  CALL SCALE_SURFACE(NEWBUFF1,NFRAMES_,
     +             NIMA,NQUA,GX1,GY1,
     +             FFSIGMANOR(1),TIMESSIGMA(1),
     +             (CMEDIANA1.EQ.'y'),YRMSTOL1,(CREFINE1.EQ.'y'),
     +             QMIN,QMAX,XSOL)
                  GX_STACK(NOPER_STACK)=GX1
                  GY_STACK(NOPER_STACK)=GX1
                  NCOEFF=(GX1+1)*(GY1+1)
                  DO NQUAD=1,4*NFRAMES_
                    DO KCOEFF=1,NCOEFF
                      FCOEFF_STACK(KCOEFF,NQUAD,NOPER_STACK)=
     +                 XSOL(KCOEFF,NQUAD)
                    END DO
                  END DO
                END IF
                DO I=1,4*NFRAMES_
                  QMIN_STACK(I,NOPER_STACK)=QMIN(I)
                  QMAX_STACK(I,NOPER_STACK)=QMAX(I)
                  XP(2*(I-1)+1)=REAL(I)/4.
                  XP(2*I)=REAL(I)/4.
                  YP(2*(I-1)+1)=QMIN(I)
                  YP(2*I)=QMAX(I)
                END DO
                CALL SUBPLOT(8*NFRAMES_,1,8*NFRAMES_,XP,YP,XP,YP,
     +                       .TRUE.,.TRUE.,.FALSE.,.FALSE.,
     +                       'frame #','signal',' ',
     +                       1,-1,1.0)
                DO I=1,4*NFRAMES_
                  IF(MOD((I-1)/4,2).EQ.0)THEN
                    CALL PGSCI(3)
                  ELSE
                    CALL PGSCI(5)
                  END IF
                  CALL SUBPLOTBIS(2,(I-1)*2+1,I*2,XP,YP,XP,YP,
     +             .FALSE.,.FALSE.,MOD(I-1,4)+2,17,1.0)
                  CALL SUBPLOTBIS(2,(I-1)*2+1,I*2,XP,YP,XP,YP,
     +             .FALSE.,.FALSE.,MOD(I-1,4)+2,101,1.0)
                  WRITE(*,'(A,I2,A,$)') '> Stack #',NOPER_STACK,
     +             ' (residuals): '
                  WRITE(*,*) I,YP((I-1)*2+1),YP(2*I)
                END DO
C..............................................................................
C dibujamos el ajuste y su diferencia
                WRITE(*,101) 'Displaying image with refined cuts...'
                CALL GIVEME_BGFG(NFRAMES_,NIMA,NQUA,BG_,FG_)
                BG=BG_
                FG=FG_
                CALL HISTOGRAM(NEWBUFF1)
                CALL SUBLOOK(.TRUE.,NEWBUFF1,.FALSE.)
                CALL DRAWSEPB9(NIMA,NQUA,256)
                IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
                  CALL SUBLOOK(.FALSE.,NEWBUFF1,.TRUE.)
                    CALL DRAWSEPB9(NIMA,NQUA,256)
                  CALL PGCLOS(IDNEW)
                  CALL PGSLCT(IDOLD)
                END IF
                IF(CPAUSE.EQ.'y')THEN
                  WRITE(*,100) 'Press <CR> to display fits...'
                  READ(*,*)
                  IF(LECHO) WRITE(*,*)
                ELSE
                  WRITE(*,101) 'Displaying fits...'
                END IF
                CALL PGIMAG(IMAGEN_,NXMAX,NYMAX,NX1,NX2,NY1,NY2,
     +           FG,BG,TR)
                CALL DRAWSEPB9(0,0,256)
                CALL PGQWIN(XW1,XW2,YW1,YW2)
                IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
                  CALL PGENV(XW1,XW2,YW1,YW2,JUST,-2)
                  CALL PGBOX('IBCTSN',0.0,0,'IBCTSN',0.0,0)
                  CALL PGIMAG(IMAGEN_,NXMAX,NYMAX,NX1,NX2,NY1,NY2,
     +             BG,FG,TR)
                    CALL DRAWSEPB9(0,0,256)
                  CALL PGCLOS(IDNEW)
                  CALL PGSLCT(IDOLD)
                END IF
                IF(CPAUSE.EQ.'y')THEN
                  WRITE(*,100) 'Press <CR> to apply fits...'
                  READ(*,*)
                  IF(LECHO) WRITE(*,*)
                ELSE
                  WRITE(*,101) 'Applying fits...'
                END IF
                IF(CNOR1.EQ.'y')THEN
                  IF(CSCALE1OPER.EQ.'s')THEN
                    CALL NORMSUMA(NFRAMES_,NIMA,NQUA)
                  ELSE
                    CALL NORMDIVI(NFRAMES_,NIMA,NQUA)
                  END IF
                END IF
                IF(CSCALE1OPER.EQ.'s')THEN
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NEWBUFF1)=
     +                 IMAGEN(J,I,NEWBUFF1)-IMAGEN_(J,I)
                    END DO
                  END DO
                ELSE
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NEWBUFF1)=
     +                 IMAGEN(J,I,NEWBUFF1)/IMAGEN_(J,I)
                      IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=
     +                 IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)/IMAGEN_(J,I)
                    END DO
                END DO
                END IF
                CALL BGFG_MASK(NFRAMES_,NEWBUFF1,BG_,FG_)
                BG=BG_
                FG=FG_
                CALL HISTOGRAM(NEWBUFF1)
                CALL SUBLOOK(.TRUE.,NEWBUFF1,.FALSE.)
                CALL DRAWSEPB9(0,0,256)
                IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
                  CALL SUBLOOK(.FALSE.,NEWBUFF1,.TRUE.)
                    CALL DRAWSEPB9(0,0,256)
                  CALL PGCLOS(IDNEW)
                  CALL PGSLCT(IDOLD)
                END IF
              END DO
              IF(CINTERACTIVE.EQ.'y')THEN
                CMORE(1:1)=READC('More fits (y/n)',CMORE,'yn')
                LASK=(CMORE.EQ.'y')
              ELSE
                LASK=.FALSE.
              END IF
            ELSE
              LASK=.FALSE.
            END IF
          END DO
        END IF
C------------------------------------------------------------------------------
C Estudio estadistico del ruido frente a la sen~al en cada subcuadrante.
        IF((CINTERACTIVE.EQ.'y').AND.(CCONT.EQ.'y'))THEN
          CCONT(1:1)=READC('Perform a detailed statistical'//
     +     ' analysis of any subquadrant (y/n)','n','yn')
          IF(CCONT.EQ.'y')THEN
            IF(CCONT.EQ.'y')THEN
              NIMSTAT=10
              NSUBQUA=0
              DO WHILE(NIMSTAT.NE.0)
                IF(NIMSTAT.LE.9)THEN
                  WRITE(CDUMMY,*) NIMSTAT
                ELSE
                  CDUMMY='@'
                END IF
                WRITE(*,101) ' 6  9  4'
                WRITE(*,101) ' 3  1  7'
                WRITE(*,101) ' 8  5  2'
                NIMSTAT=READILIM('Image number (0=exit)',CDUMMY,
     +           0,NFRAMES_)
                IF(NIMSTAT.NE.0) CALL DETAILB9(NEWBUFF1,NCBUFF,
     +           NIMSTAT)
              END DO
            END IF
          END IF
        END IF
C------------------------------------------------------------------------------
C Median sky
        IF(CINTERACTIVE.EQ.'y')THEN
          CCONT(1:1)=READC('Compute median sky for each frame (y/x)',
     +     'y','yx')
          IF(CCONT.EQ.'x')THEN
            CPOST='n'
            RETURN
          END IF
        END IF
C preguntamos si queremos rellenar regiones de la mascara con ajustes a 
C superficies polinomicas suaves
        IF(CINTERACTIVE.EQ.'y')THEN
          CFILL(1:1)=READC('Fill masked regions with 2D polynomial '//
     +     'fits (y/n)',CFILL,'yn')
          IF(CFILL.EQ.'y')THEN
            WRITE(CDUMMY,*) GXFILL
            GXFILL=READILIM('Polynomial degree in X direction',CDUMMY,
     +       0,INT(SQRT(REAL(MAXNCOEFF)))-1)
            WRITE(CDUMMY,*) GYFILL
            GYFILL=READILIM('Polynomial degree in Y direction',CDUMMY,
     +       0,INT(SQRT(REAL(MAXNCOEFF)))-1)
            CMEDIANAFILL(1:1)=
     +       READC('Refine fit minimizing median (y/n)',
     +       CMEDIANAFILL,'yn')
            IF(CMEDIANAFILL.EQ.'y')THEN
              WRITE(CDUMMY,*) YRMSTOLFILL
              YRMSTOLFILL=READF('YRMSTOL for DOWNHILL',CDUMMY)
            END IF
            CREFINEFILL(1:1)=
     +       READC('Refine fit removing deviating points'//
     +       ' (y/n)',CREFINEFILL,'yn')
          END IF
          CERRSKY(1:1)=READC('Compute error in median sky (y/n)',
     +     CERRSKY,'yn')
        END IF
        IF(CFILL.EQ.'y')THEN
          CALL SCALE_SURFACE(NEWBUFF1,NFRAMES_,
     +     0,0,GXFILL,GYFILL,
     +     FFSIGMANOR(1),TIMESSIGMA(1),
     +     (CMEDIANAFILL.EQ.'y'),YRMSTOLFILL,(CREFINEFILL.EQ.'y'),
     +     QMIN,QMAX,XSOL)
          DO I=1,4*NFRAMES_
            XP(2*(I-1)+1)=REAL(I)/4.
            XP(2*I)=REAL(I)/4.
            YP(2*(I-1)+1)=QMIN(I)
            YP(2*I)=QMAX(I)
          END DO
          CALL SUBPLOT(8*NFRAMES_,1,8*NFRAMES_,XP,YP,XP,YP,
     +                 .TRUE.,.TRUE.,.FALSE.,.FALSE.,
     +                 'frame #','signal',' ',
     +                 1,-1,1.0)
          DO I=1,4*NFRAMES_
            IF(MOD((I-1)/4,2).EQ.0)THEN
              CALL PGSCI(3)
            ELSE
              CALL PGSCI(5)
            END IF
            CALL SUBPLOTBIS(2,(I-1)*2+1,I*2,XP,YP,XP,YP,
     +       .FALSE.,.FALSE.,MOD(I-1,4)+2,17,1.0)
            CALL SUBPLOTBIS(2,(I-1)*2+1,I*2,XP,YP,XP,YP,
     +       .FALSE.,.FALSE.,MOD(I-1,4)+2,101,1.0)
            WRITE(*,100)'> residuals: '
            WRITE(*,*) I,YP((I-1)*2+1),YP(2*I)
          END DO
          WRITE(*,101) 'Displaying fitted mask regions...'
          CALL PGIMAG(IMAGEN_,NXMAX,NYMAX,NX1,NX2,NY1,NY2,FG,BG,TR)
          CALL DRAWSEPB9(0,0,256)
          CALL PGQWIN(XW1,XW2,YW1,YW2)
          IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
            CALL PGENV(XW1,XW2,YW1,YW2,JUST,-2)
            CALL PGBOX('IBCTSN',0.0,0,'IBCTSN',0.0,0)
            CALL PGIMAG(IMAGEN_,NXMAX,NYMAX,NX1,NX2,NY1,NY2,
     +       BG,FG,TR)
              CALL DRAWSEPB9(0,0,256)
            CALL PGCLOS(IDNEW)
            CALL PGSLCT(IDOLD)
          END IF
          IF(CINTERACTIVE.EQ.'y')THEN
            WRITE(*,100) 'Press <CR> to proceed with sky subtraction...'
            READ(*,*)
          ELSE
            WRITE(*,101) 'Computing median sky for each frame...'
          END IF
        ELSE
          WRITE(*,101) 'Computing median sky for each frame...'
        END IF
C------------------------------------------------------------------------------
        DO K=1,NFRAMES_
          WRITE(*,'(A,I1,A,$)') 'Working with frame #',K,'...'
          DO I=1,256
            DO J=1,256
              KK_=0
              DO KK=1,NFRAMES_
                IF(KK.NE.K)THEN
                  IF(SUBMASKBOX9(J,I,KK))THEN
                    KK_=KK_+1
                    PIXEL(KK_)=IMAGEN(J+DJ(KK),I+DI(KK),NEWBUFF1)
                    EPIXEL(KK_)=IMAGEN(J+DJ(KK),I+DI(KK),
     +               NEWBUFF1+NMAXBUFF/2)
                  END IF
                END IF
              END DO
              IF(KK_.GT.0)THEN
                IF(KK_.EQ.1)THEN
                  IMAGEN(J+DJ(K),I+DI(K),NEWBUFF4)=PIXEL(1)
                  IF(CERRSKY.EQ.'y')THEN
                    IMAGEN(J+DJ(K),I+DI(K),NEWBUFF4+NMAXBUFF/2)=
     +               EPIXEL(1)
                  ELSE
                    IMAGEN(J+DJ(K),I+DI(K),NEWBUFF4+NMAXBUFF/2)=0.
                  END IF
                ELSE
                  IF(CERRSKY.EQ.'y')THEN
                    IMAGEN(J+DJ(K),I+DI(K),NEWBUFF4)=FMEDIAN1E(KK_,
     +               PIXEL,EPIXEL,IMAGEN(J+DJ(K),I+DI(K),
     +               NEWBUFF4+NMAXBUFF/2))
                  ELSE
                    IMAGEN(J+DJ(K),I+DI(K),NEWBUFF4)=FMEDIAN1(KK_,PIXEL)
                    IMAGEN(J+DJ(K),I+DI(K),NEWBUFF4+NMAXBUFF/2)=0.
                  END IF
                END IF
              ELSE
                WRITE(*,100) '> I,J,K: '
                WRITE(*,*) I,J,K
                WRITE(*,101) '> No. of pixels to compute median sky: 0'
                KK_=0
                DO KK=1,NFRAMES_
                  IF(KK.NE.K)THEN
                    IF(CFILL.EQ.'y')THEN
                      KK_=KK_+1
                      PIXEL(KK_)=IMAGEN_(J+DJ(KK),I+DI(KK))
                      EPIXEL(KK_)=0. !asumo error despreciable en el ajuste
                    ELSE
                      KK_=KK_+1
                      PIXEL(KK_)=IMAGEN(J+DJ(KK),I+DI(KK),NEWBUFF1)
                      EPIXEL(KK_)=IMAGEN(J+DJ(KK),I+DI(KK),
     +                 NEWBUFF1+NMAXBUFF/2)
                    END IF
                  END IF
                END DO
                IF(KK_.EQ.1)THEN
                  IMAGEN(J+DJ(K),I+DI(K),NEWBUFF4)=PIXEL(1)
                  IF(CERRSKY.EQ.'y')THEN
                    IMAGEN(J+DJ(K),I+DI(K),NEWBUFF4+NMAXBUFF/2)=
     +               EPIXEL(1)
                  ELSE
                    IMAGEN(J+DJ(K),I+DI(K),NEWBUFF4+NMAXBUFF/2)=0.
                  END IF
                ELSE
                  IF(CERRSKY.EQ.'y')THEN
                    IMAGEN(J+DJ(K),I+DI(K),NEWBUFF4)=FMEDIAN1E(KK_,
     +               PIXEL,EPIXEL,IMAGEN(J+DJ(K),I+DI(K),
     +               NEWBUFF4+NMAXBUFF/2))
                  ELSE
                    IMAGEN(J+DJ(K),I+DI(K),NEWBUFF4)=FMEDIAN1(KK_,PIXEL)
                    IMAGEN(J+DJ(K),I+DI(K),NEWBUFF4+NMAXBUFF/2)=0.
                  END IF
                END IF
              END IF
            END DO
          END DO
          WRITE(*,101) 'OK!'
        END DO
        IF(NFRAMES_.LT.NBOXMAX)THEN
          DO K=NFRAMES_+1,NBOXMAX
            DO I=1,256
              DO J=1,256
                IMAGEN(J+DJ(K),I+DI(K),NEWBUFF4)=0.
                IMAGEN(J+DJ(K),I+DI(K),NEWBUFF4+NMAXBUFF/2)=0.
              END DO
            END DO
          END DO
        END IF
        NX1=1
        NX2=768
        NY1=1
        NY2=768
        CALL PGIMAG(IMAGEN(1,1,NEWBUFF4),NXMAX,NYMAX,
     +   NX1,NX2,NY1,NY2,FG,BG,TR)
        CALL DRAWSEPB9(0,0,256)
        CALL PGQWIN(XW1,XW2,YW1,YW2)
        IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
          CALL PGENV(XW1,XW2,YW1,YW2,JUST,-2)
          CALL PGBOX('IBCTSN',0.0,0,'IBCTSN',0.0,0)
          CALL PGIMAG(IMAGEN(1,1,NEWBUFF4),NXMAX,NYMAX,
     +     NX1,NX2,NY1,NY2,BG,FG,TR)
            CALL DRAWSEPB9(0,0,256)
          CALL PGCLOS(IDNEW)
          CALL PGSLCT(IDOLD)
        END IF
C------------------------------------------------------------------------------
        IF(CINTERACTIVE.EQ.'y')THEN
          CCONT(1:1)=READC('Do you want to subtract median sky (y/n/x)',
     +     'y','ynx')
          IF(CCONT.EQ.'x')THEN
            DO I=NY1,NY2
              DO J=NX1,NX2
                IMAGEN(J,I,NEWBUFF1)=IMAGEN(J,I,NEWBUFF4)
                IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=
     +           IMAGEN(J,I,NEWBUFF4+NMAXBUFF/2)
              END DO
            END DO
            CPOST='n'
            RETURN
          END IF
        ELSE
          CCONT='y'
        END IF
C
        IF(CCONT.EQ.'y')THEN
          WRITE(*,101) 'Subtracting median sky...'
        ELSE
          WRITE(*,101) 
     +     'WARNING: median sky is not going to be subtracted!'
        END IF
C------------------------------------------------------------------------------
        IF(CCONT.EQ.'y')THEN
          DO I=NY1,NY2
            DO J=NX1,NX2
              IMAGEN(J,I,NEWBUFF1)=
     +         IMAGEN(J,I,NEWBUFF1)-IMAGEN(J,I,NEWBUFF4)
              IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=SQRT(
     +         IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)*
     +         IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)+
     +         IMAGEN(J,I,NEWBUFF4+NMAXBUFF/2)*
     +         IMAGEN(J,I,NEWBUFF4+NMAXBUFF/2))
            END DO
          END DO
        END IF
        CALL STATISTICB9(NFRAMES_,NEWBUFF1,FMEAN,FSIGMA)
        IF(FSIGMA.GT.0.0)THEN
          BG=FMEAN-5.*FSIGMA
          FG=FMEAN+5.*FSIGMA
        ELSE
          BG=FMEAN-1.0
          FG=FMEAN+1.0
        END IF
        CALL HISTOGRAM(NEWBUFF1)
        CALL SUBLOOK(.TRUE.,NEWBUFF1,.FALSE.)
        CALL DRAWSEPB9(0,0,256)
        IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
          CALL SUBLOOK(.FALSE.,NEWBUFF1,.TRUE.)
            CALL DRAWSEPB9(0,0,256)
          CALL PGCLOS(IDNEW)
          CALL PGSLCT(IDOLD)
        END IF
C------------------------------------------------------------------------------
C Estudio estadistico del ruido frente a la sen~al en cada subcuadrante.
        IF(CINTERACTIVE.EQ.'y')THEN
          CCONT(1:1)=READC('Perform a detailed statistical'//
     +     ' analysis of any subquadrant (y/n)','n','yn')
          IF(CCONT.EQ.'y')THEN
            IF(CCONT.EQ.'y')THEN
              NIMSTAT=10
              NSUBQUA=0
              DO WHILE(NIMSTAT.NE.0)
                IF(NIMSTAT.LE.9)THEN
                  WRITE(CDUMMY,*) NIMSTAT
                ELSE
                  CDUMMY='@'
                END IF
                WRITE(*,101) ' 6  9  4'
                WRITE(*,101) ' 3  1  7'
                WRITE(*,101) ' 8  5  2'
                NIMSTAT=READILIM('Image number (0=exit)',CDUMMY,
     +           0,NFRAMES_)
                IF(NIMSTAT.NE.0) CALL DETAILB9(NEWBUFF1,NCBUFF,
     +           NIMSTAT)
              END DO
            END IF
          END IF
        END IF
C------------------------------------------------------------------------------
C Scale quadrants
        IF(CINTERACTIVE.EQ.'y')THEN
          IF(CSCALE2.NE.'n')THEN
            CCONT(1:1)=
     +       READC('Do you want to scale each quadrant (y/n/x)',
     +       'y','ynx')
          ELSE
            CCONT(1:1)=
     +       READC('Do you want to scale each quadrant (y/n/x)',
     +       'n','ynx')
          END IF
          IF(CCONT.EQ.'x')THEN
            CPOST='n'
            RETURN
          END IF
        ELSE
          CCONT='y'
          IF(CCONT.EQ.'n') CSCALE2='n'
        END IF
C------------------------------------------------------------------------------
        IF(CCONT.EQ.'y')THEN
          NIMA=0
          NQUA=0
          NITER=1
          LASK=.TRUE.
          DO WHILE(LASK)
            IF(CINTERACTIVE.EQ.'y')THEN
              CSCALE2(1:1)=
     +         READC('Fit#2 [c]onstant, [s]urface or [n]one '//
     +         '(c/s/n)',CSCALE2,'csn')
              IF(CSCALE2.EQ.'s')THEN
                WRITE(CDUMMY,*) GX2
                GX2=READILIM('Polynomial degree in X direction',CDUMMY,
     +           0,INT(SQRT(REAL(MAXNCOEFF)))-1)
                WRITE(CDUMMY,*) GY2
                GY2=READILIM('Polynomial degree in Y direction',CDUMMY,
     +           0,INT(SQRT(REAL(MAXNCOEFF)))-1)
                CMEDIANA2(1:1)=
     +           READC('Refine fit minimizing median (y/n)',
     +           CMEDIANA2,'yn')
                IF(CMEDIANA2.EQ.'y')THEN
                  WRITE(CDUMMY,*) YRMSTOL2
                  YRMSTOL2=READF('YRMSTOL for DOWNHILL',CDUMMY)
                END IF
                CREFINE2(1:1)=
     +           READC('Refine fit removing deviating points'//
     +           ' (y/n)',CREFINE2,'yn')
              END IF
              IF(CSCALE2.NE.'n')THEN
                CSCALE2OPER(1:1)=
     +           READC('[s]ubtract fit or [d]ivide by fit (s/d)',
     +           CSCALE2OPER,'sd')
                IF(CSCALE2OPER.EQ.'d')THEN
                  WRITE(*,101) 'WARNING: normalization is recomended '//
     +             'to preserve number of counts'
                  CNOR2(1:1)=READC('Normalize fit (y/n)','y','yn')
                ELSE
                  CNOR2='n'
                END IF
                WRITE(*,101) ' 6  9  4'
                WRITE(*,101) ' 3  1  7'
                WRITE(*,101) ' 8  5  2'
                WRITE(CDUMMY,*) NIMA
                NIMA=READILIM('Image number (0=all)',CDUMMY,0,NFRAMES_)
                WRITE(*,101) ' 2  4'
                WRITE(*,101) ' 1  3'
                WRITE(CDUMMY,*) NQUA
                NQUA=READILIM('Quadrant (0=all)',CDUMMY,0,4)
                WRITE(CDUMMY,*) NITER
                NITER=READILIM('No. of iterations',CDUMMY,1,10)
                IF(CINTERACTIVE.EQ.'y')THEN
                  CPAUSE(1:1)=
     +             READC('Pause between intermediate plots (y/n)',
     +             'n','yn')
                ELSE
                  CPAUSE='n'
                END IF
                WRITE(*,101) 'Working...'
              END IF
            END IF
C..............................................................................
            IF(CSCALE2.NE.'n')THEN
              DO ITER=1,NITER
                NOPER_STACK=NOPER_STACK+1
                WRITE(*,100) '>>> Stacking operation #'
                WRITE(*,*) NOPER_STACK
                NIMA_STACK(NOPER_STACK)=NIMA
                NQUA_STACK(NOPER_STACK)=NQUA
                CNOR_STACK(NOPER_STACK)=CNOR2
                CSCALEOPER_STACK(NOPER_STACK)=CSCALE2OPER
                IF(CSCALE2.EQ.'c')THEN
                  TYPEOPER_STACK(NOPER_STACK)=1
                  CALL SCALE_CONSTANT(NEWBUFF1,NFRAMES_,
     +             NIMA,NQUA,QMIN,QMAX)
                  DO NQUAD=1,4*NFRAMES_
                    FMEDIAN_STACK(NQUAD,NOPER_STACK)=QMIN(NQUAD)
                  END DO
                ELSE
                  TYPEOPER_STACK(NOPER_STACK)=2
                  CALL SCALE_SURFACE(NEWBUFF1,NFRAMES_,
     +             NIMA,NQUA,GX2,GY2,
     +             FFSIGMANOR(1),TIMESSIGMA(1),
     +             (CMEDIANA2.EQ.'y'),YRMSTOL2,(CREFINE2.EQ.'y'),
     +             QMIN,QMAX,XSOL)
                  GX_STACK(NOPER_STACK)=GX2
                  GY_STACK(NOPER_STACK)=GY2
                  NCOEFF=(GX2+1)*(GY2+1)
                  DO NQUAD=1,4*NFRAMES_
                    DO KCOEFF=1,NCOEFF
                      FCOEFF_STACK(KCOEFF,NQUAD,NOPER_STACK)=
     +                 XSOL(KCOEFF,NQUAD)
                    END DO
                  END DO
                END IF
                DO I=1,4*NFRAMES_
                  QMIN_STACK(I,NOPER_STACK)=QMIN(I)
                  QMAX_STACK(I,NOPER_STACK)=QMAX(I)
                  XP(2*(I-1)+1)=REAL(I)/4.
                  XP(2*I)=REAL(I)/4.
                  YP(2*(I-1)+1)=QMIN(I)
                  YP(2*I)=QMAX(I)
                END DO
                CALL SUBPLOT(8*NFRAMES_,1,8*NFRAMES_,XP,YP,XP,YP,
     +                       .TRUE.,.TRUE.,.FALSE.,.FALSE.,
     +                       'frame #','signal',' ',
     +                       1,-1,1.0)
                DO I=1,4*NFRAMES_
                  IF(MOD((I-1)/4,2).EQ.0)THEN
                    CALL PGSCI(3)
                  ELSE
                    CALL PGSCI(5)
                  END IF
                  CALL SUBPLOTBIS(2,(I-1)*2+1,I*2,XP,YP,XP,YP,
     +             .FALSE.,.FALSE.,MOD(I-1,4)+2,17,1.0)
                  CALL SUBPLOTBIS(2,(I-1)*2+1,I*2,XP,YP,XP,YP,
     +             .FALSE.,.FALSE.,MOD(I-1,4)+2,101,1.0)
                  WRITE(*,'(A,I2,A,$)') '> Stack #',NOPER_STACK,
     +             ' (residuals): '
                  WRITE(*,*) I,YP((I-1)*2+1),YP(2*I)
                END DO
C..............................................................................
C dibujamos el ajuste y su diferencia
                WRITE(*,101) 'Displaying image with refined cuts...'
                CALL GIVEME_BGFG(NFRAMES_,NIMA,NQUA,BG_,FG_)
                BG=BG_
                FG=FG_
                CALL HISTOGRAM(NEWBUFF1)
                CALL SUBLOOK(.TRUE.,NEWBUFF1,.FALSE.)
                CALL DRAWSEPB9(NIMA,NQUA,256)
                IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
                  CALL SUBLOOK(.FALSE.,NEWBUFF1,.TRUE.)
                    CALL DRAWSEPB9(NIMA,NQUA,256)
                  CALL PGCLOS(IDNEW)
                  CALL PGSLCT(IDOLD)
                END IF
                IF(CPAUSE.EQ.'y')THEN
                  WRITE(*,100) 'Press <CR> to display fits...'
                  READ(*,*)
                  IF(LECHO) WRITE(*,*)
                ELSE
                  WRITE(*,101) 'Displaying fits...'
                END IF
                CALL PGIMAG(IMAGEN_,NXMAX,NYMAX,NX1,NX2,NY1,NY2,
     +           FG,BG,TR)
                CALL DRAWSEPB9(0,0,256)
                CALL PGQWIN(XW1,XW2,YW1,YW2)
                IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
                  CALL PGENV(XW1,XW2,YW1,YW2,JUST,-2)
                  CALL PGBOX('IBCTSN',0.0,0,'IBCTSN',0.0,0)
                  CALL PGIMAG(IMAGEN_,NXMAX,NYMAX,NX1,NX2,NY1,NY2,
     +             BG,FG,TR)
                    CALL DRAWSEPB9(0,0,256)
                  CALL PGCLOS(IDNEW)
                  CALL PGSLCT(IDOLD)
                END IF
                IF(CPAUSE.EQ.'y')THEN
                  WRITE(*,100) 'Press <CR> to apply fits...'
                  READ(*,*)
                  IF(LECHO) WRITE(*,*)
                ELSE
                  WRITE(*,101) 'Applying fits...'
                END IF
                IF(CNOR2.EQ.'y')THEN
                  IF(CSCALE2OPER.EQ.'s')THEN
                    CALL NORMSUMA(NFRAMES_,NIMA,NQUA)
                  ELSE
                    CALL NORMDIVI(NFRAMES_,NIMA,NQUA)
                  END IF
                END IF
                IF(CSCALE2OPER.EQ.'s')THEN
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NEWBUFF1)=
     +                 IMAGEN(J,I,NEWBUFF1)-IMAGEN_(J,I)
                    END DO
                  END DO
                ELSE
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NEWBUFF1)=
     +                 IMAGEN(J,I,NEWBUFF1)/IMAGEN_(J,I)
                      IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=
     +                 IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)/IMAGEN_(J,I)
                    END DO
                  END DO
                END IF
                CALL BGFG_MASK(NFRAMES_,NEWBUFF1,BG_,FG_)
                BG=BG_
                FG=FG_
                CALL HISTOGRAM(NEWBUFF1)
                CALL SUBLOOK(.TRUE.,NEWBUFF1,.FALSE.)
                CALL DRAWSEPB9(0,0,256)
                IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
                  CALL SUBLOOK(.FALSE.,NEWBUFF1,.TRUE.)
                    CALL DRAWSEPB9(0,0,256)
                  CALL PGCLOS(IDNEW)
                  CALL PGSLCT(IDOLD)
                END IF
              END DO
              IF(CINTERACTIVE.EQ.'y')THEN
                CMORE(1:1)=READC('More fits (y/n)',CMORE,'yn')
                LASK=(CMORE.EQ.'y')
              ELSE
                LASK=.FALSE.
              END IF
            ELSE
              LASK=.FALSE.
            END IF
          END DO
        END IF
C------------------------------------------------------------------------------
C Estudio estadistico del ruido frente a la sen~al en cada subcuadrante.
        IF((CINTERACTIVE.EQ.'y').AND.(CCONT.EQ.'y'))THEN
          CCONT(1:1)=READC('Perform a detailed statistical'//
     +     ' analysis of any subquadrant (y/n)','n','yn')
          IF(CCONT.EQ.'y')THEN
            IF(CCONT.EQ.'y')THEN
              NIMSTAT=10
              NSUBQUA=0
              DO WHILE(NIMSTAT.NE.0)
                IF(NIMSTAT.LE.9)THEN
                  WRITE(CDUMMY,*) NIMSTAT
                ELSE
                  CDUMMY='@'
                END IF
                WRITE(*,101) ' 6  9  4'
                WRITE(*,101) ' 3  1  7'
                WRITE(*,101) ' 8  5  2'
                NIMSTAT=READILIM('Image number (0=exit)',CDUMMY,
     +           0,NFRAMES_)
                IF(NIMSTAT.NE.0) CALL DETAILB9(NEWBUFF1,NCBUFF,
     +           NIMSTAT)
              END DO
            END IF
          END IF
        END IF
C------------------------------------------------------------------------------
        IF((CINTERACTIVE.EQ.'y').OR.(CINTERACTIVE.EQ.'m'))THEN
          IF(CINTERACTIVE.EQ.'y')THEN
            CCONT(1:1)=READC('Do you want to coadd images (y/m/x)',
     +       'y','ymx')
          ELSE
            CCONT='m'
          END IF
          IF(CCONT.EQ.'x')THEN
            CPOST='n'
            RETURN
          END IF
          IF(CCONT.EQ.'m')THEN
            CALL MIDEOFFSET(NFRAMES_,NEWBUFF1,NEWBUFF2,NEWBUFF3,
     +       NEWBUFF4,FOFFSETX,FOFFSETY)
            CPOST='n'
            RETURN
          END IF
C decidimos si eliminamos pixels extran~os
          WRITE(*,*)
          CREJECT(1:1)=
     +     READC('Reject hot/cool pixels (y/n)',CREJECT,'yn')
          IF(CREJECT.EQ.'y')THEN
C numero de pixels umbral para deteccion de pixels calientes
            WRITE(CDUMMY,*) FPIXELTHRESHOLDCR
            FPIXELTHRESHOLDCR=READF('Threshold in no. of pixels >'//
     +       ' Times_Sigma for hot/cool pixels',CDUMMY)
          END IF
        ELSE
          IF(CREJECT.EQ.'y')THEN
            WRITE(*,101) 
     +       'Detecting hot/cool pixels and coadding images...'
          ELSE
            WRITE(*,101) 'Coadding images...'
          END IF
        END IF
C------------------------------------------------------------------------------
        IF(CREJECT.EQ.'y')THEN !redibujamos imagen para superponer pixels
          IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
            CALL SUBLOOK(.FALSE.,NEWBUFF1,.TRUE.)
              CALL DRAWSEPB9(0,0,256)
            CALL PGSLCT(IDOLD)
          END IF
        END IF
C it's time to coadd the 9 frames into a single image
        NAXIS(1,NEWBUFF2)=NXMAXB9_
        NAXIS(2,NEWBUFF2)=NYMAXB9_
        NFRAMES(NEWBUFF2)=1
        NAXIS(1,NEWBUFF2+NMAXBUFF/2)=NXMAXB9_
        NAXIS(2,NEWBUFF2+NMAXBUFF/2)=NYMAXB9_
        NFRAMES(NEWBUFF2+NMAXBUFF/2)=1
        DO I=1,NAXIS(2,NEWBUFF2)
          DO J=1,NAXIS(1,NEWBUFF2)
            IMAGEN(J,I,NEWBUFF2)=0. !inicializamos a cero el nuevo frame
            FPIXUSED(J,I)=0. !aqui almacenamos el numero de pixels usados en
                             !el calculo de cada nuevo pixel en la imagen final
            IMAGEN(J,I,NEWBUFF2+NMAXBUFF/2)=0. !imagen de errores
          END DO
        END DO
        IF(NEWBUFF3.GT.0)THEN
          NAXIS(1,NEWBUFF3)=NXMAXB9_
          NAXIS(2,NEWBUFF3)=NYMAXB9_
          NFRAMES(NEWBUFF2)=1
          DO I=1,NAXIS(2,NEWBUFF3)
            DO J=1,NAXIS(1,NEWBUFF3)
              IF(MASKBOX9(J,I))THEN
                IMAGEN(J,I,NEWBUFF3)=1.
              ELSE
                IMAGEN(J,I,NEWBUFF3)=0.
              END IF
            END DO
          END DO
        END IF
C recorremos cada pixel en la nueva imagen y determinamos cuantos pixels de 
C los frames iniciales se pueden usar en cada caso
        LANYCR=.FALSE.
        NEXTINFO=0
        DO I=1,NAXIS(2,NEWBUFF2)
          DO J=1,NAXIS(1,NEWBUFF2)
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
                IF(INSIDE(JJ,1,256,II,1,256))THEN
                  KK=KK+1
                  PIXEL(KK)=IMAGEN(JJ+DJ(K),II+DI(K),NEWBUFF1)
                  EPIXEL(KK)=
     +             IMAGEN(JJ+DJ(K),II+DI(K),NEWBUFF1+NMAXBUFF/2)
                  PIXELF(KK)=FRACCION_PIXEL
                  K_LAST(KK)=K !necesario para eliminar hot/cool pixels
                  JJ_LAST(KK)=JJ !idem
                  II_LAST(KK)=II !idem
                END IF
              END DO
            END DO
C..............................................................................
C..............................................................................
C En caso necesario, eliminamos hot/cool pixels.
            IF(CREJECT.EQ.'y')THEN
C----->
              IF(KK.LE.4)THEN !si tenemos menos de cinco pixels
                DO NK=1,KK
                  NPIXELSVECINOS=0
                  DO II=II_LAST(NK)-2,II_LAST(NK)+2
                    IF((II.GE.1).AND.(II.LE.256))THEN
                      DO JJ=JJ_LAST(NK)-2,JJ_LAST(NK)+2
                        IF((JJ.GE.1).AND.(JJ.LE.256))THEN
                          IF((II.NE.II_LAST(NK)).OR.
     +                       (JJ.NE.JJ_LAST(NK)))THEN
                            NPIXELSVECINOS=NPIXELSVECINOS+1
                            PIXELVECINO(NPIXELSVECINOS)=IMAGEN(
     +                       JJ+DJ(K_LAST(NK)),
     +                       II+DI(K_LAST(NK)),
     +                       NEWBUFF1)
                          END IF
                        END IF
                      END DO
                    END IF
                  END DO
                  CALL ORDENA1F(NPIXELSVECINOS,PIXELVECINO)
ccc                  DO NN=1,NPIXELSVECINOS-4
ccc                    PIXELVECINO(NN)=PIXELVECINO(NN+2)
ccc                  END DO
ccc                  NPIXELSVECINOS=NPIXELSVECINOS-4
                  FMEAN_PIX=FMEAN0(NPIXELSVECINOS,PIXELVECINO,
     +             FSIGMA_PIX)
                  TSHOT=FTSTUDENTI(NPIXELSVECINOS-2, !grados de libertad
     +             FPIXELTHRESHOLDCR/REAL(NPIXREGION)) !probabilidad
                  IF(ABS(PIXEL(NK)-FMEAN_PIX).GT.TSHOT*FSIGMA_PIX)THEN
                    PIXEL(NK)=FMEDIAN1(NPIXELSVECINOS,PIXELVECINO)
                    EPIXEL(NK)=FSIGMA_PIX
                    CALL PGSCI(4)
                    CALL PGPOINT(1,REAL(JJ_LAST(NK)+DJ(K_LAST(NK))),
     +               REAL(II_LAST(NK)+DI(K_LAST(NK))),21)
                    CALL PGSCI(1)
                    IF(CPOST.EQ.'y')THEN
                      CALL PGSLCT(IDNEW)
                      CALL PGSCI(4)
                      CALL PGPOINT(1,REAL(JJ_LAST(NK)+DJ(K_LAST(NK))),
     +                 REAL(II_LAST(NK)+DI(K_LAST(NK))),21)
                      CALL PGSCI(1)
                      CALL PGSLCT(IDOLD)
                    END IF
                    LANYCR=.TRUE.
                  END IF
                END DO
C----->
              ELSE !si tenemos mas de 4 pixels
C ordenamos los pixels por sen~al creciente
                CALL ORDENA3F3I(KK,PIXEL,EPIXEL,PIXELF,K_LAST,
     +           II_LAST,JJ_LAST)
C calculamos la media y la dispersion eliminando el pixel mas brillante y
C el mas debil
                FMEAN_PIX=FMEAN0(KK-2,PIXEL(2),FSIGMA_PIX)
                TSHOT=FTSTUDENTI(KK-4, !grados de libertad
     +           FPIXELTHRESHOLDCR/REAL(NPIXREGION)) !/ probabilidad
C eliminamos el pixel mas brillante si excede en TSHOT veces sigma
                IF(PIXEL(KK)-FMEAN_PIX.GT.TSHOT*FSIGMA_PIX)THEN
                  CALL PGSCI(5)
                  CALL PGPOINT(1,REAL(JJ_LAST(KK)+DJ(K_LAST(KK))),
     +             REAL(II_LAST(KK)+DI(K_LAST(KK))),21)
                  CALL PGSCI(1)
                  IF(CPOST.EQ.'y')THEN
                    CALL PGSLCT(IDNEW)
                    CALL PGSCI(5)
                    CALL PGPOINT(1,REAL(JJ_LAST(KK)+DJ(K_LAST(KK))),
     +               REAL(II_LAST(Kk)+DI(K_LAST(KK))),21)
                    CALL PGSCI(1)
                    CALL PGSLCT(IDOLD)
                  END IF
                  KK=KK-1
                  LANYCR=.TRUE.
                END IF
C eliminamos el pixel mas debil si excede en -TSHOT veces sigma
                IF(FMEAN_PIX-PIXEL(1).GT.TSHOT*FSIGMA_PIX)THEN
                  CALL PGSCI(6)
                  CALL PGPOINT(1,REAL(JJ_LAST(1)+DJ(K_LAST(1))),
     +             REAL(II_LAST(1)+DI(K_LAST(1))),21)
                  CALL PGSCI(1)
                  IF(CPOST.EQ.'y')THEN
                    CALL PGSLCT(IDNEW)
                    CALL PGSCI(6)
                    CALL PGPOINT(1,REAL(JJ_LAST(1)+DJ(K_LAST(1))),
     +               REAL(II_LAST(1)+DI(K_LAST(1))),21)
                    CALL PGSCI(1)
                    CALL PGSLCT(IDOLD)
                  END IF
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
C..............................................................................
C..............................................................................
C normalizamos cada pixel por el numero total de frames que hemos utilizado
C para calcular la sen~al total en dicho pixel
            IF(KK.EQ.0)THEN !OJO
              IMAGEN(J,I,NEWBUFF2)=0.
              IMAGEN(J,I,NEWBUFF2+NMAXBUFF/2)=0.
              FPIXUSED(J,I)=0.
            ELSE
              DO K=1,KK
                IMAGEN(J,I,NEWBUFF2)=IMAGEN(J,I,NEWBUFF2)+
     +           PIXEL(K)*PIXELF(K)
                IMAGEN(J,I,NEWBUFF2+NMAXBUFF/2)=
     +           IMAGEN(J,I,NEWBUFF2+NMAXBUFF/2)+
     +           EPIXEL(K)*EPIXEL(K)*PIXELF(K)*PIXELF(K)
                FPIXUSED(J,I)=FPIXUSED(J,I)+PIXELF(K)
              END DO
              IMAGEN(J,I,NEWBUFF2)=IMAGEN(J,I,NEWBUFF2)/FPIXUSED(J,I)
              IMAGEN(J,I,NEWBUFF2+NMAXBUFF/2)=
     +         SQRT(IMAGEN(J,I,NEWBUFF2+NMAXBUFF/2))/FPIXUSED(J,I)
            END IF
          END DO
          CALL SHOWPERC(1,NAXIS(2,NEWBUFF2),1,I,NEXTINFO)
        END DO
C si se ha pedido, salvamos en otro buffer el numero de pixels utilizados para
C el calculo de cada pixel
        NAXIS(1,NEWBUFF4)=NXMAXB9_
        NAXIS(2,NEWBUFF4)=NYMAXB9_
        NFRAMES(NEWBUFF4)=1
        DO I=1,NAXIS(2,NEWBUFF4)
          DO J=1,NAXIS(1,NEWBUFF4)
            IMAGEN(J,I,NEWBUFF4)=FPIXUSED(J,I)
            IF(FPIXUSED(J,I).LE.0.0)THEN
              IMAGEN(J,I,NEWBUFF2+NMAXBUFF/2)=-1.
            END IF
          END DO
        END DO
C
        IF(CREJECT.EQ.'y')THEN !cerramos el fichero postscript
          IF(CPOST.EQ.'y')THEN
            CALL PGSLCT(IDNEW)
            CALL PGCLOS(IDNEW)
            CALL PGSLCT(IDOLD)
          END IF
        END IF
C
        IF((CINTERACTIVE.EQ.'y').AND.LANYCR)THEN
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          IF(LECHO) WRITE(*,*)
        END IF
C
        NX1=1
        NX2=NAXIS(1,NEWBUFF2)
        NY1=1
        NY2=NAXIS(2,NEWBUFF2)
        CALL STATISTIC(NEWBUFF2,NX1,NX2,NY1,NY2,.FALSE.,.TRUE.,.TRUE.,
     +   0.0,.FALSE.)
        IF(FSIGMA.GT.0.0)THEN
          BG=FMEAN-5.*FSIGMA
          FG=FMEAN+5.*FSIGMA
        ELSE
          BG=FMEAN-1.0
          FG=FMEAN+1.0
        END IF
        CALL HISTOGRAM(NEWBUFF2)
        CALL SUBLOOK(.FALSE.,NEWBUFF2,.FALSE.)
        IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
          CALL SUBLOOK(.FALSE.,NEWBUFF2,.TRUE.)
            CALL DRAWSEPB9(0,0,256)
          CALL PGCLOS(IDNEW)
          CALL PGSLCT(IDOLD)
        END IF
C------------------------------------------------------------------------------
C Al iterar:
C (1) Calculamos el histograma en las diferentes regiones segun el numero de 
C     imagenes utilizadas (1,..., 9). El histograma se realiza usando
C     pixels comprendidos entre la media +/- 3 sigma (ambas obtenidas
C     eliminando tambien pixels a +/- 3 sigma). Estos histogramas quedan
C     caracterizados por la media y sigma (NOTAR que lo que aqui llamamos sigma
C     es en realidad el error en la media y no la desviacion tipica ---ya que
C     al calcular la imagen suma, hemos dividido por el numero de pixels
C     utilizados, por lo que en dicha imagen cada pixel es la media de una
C     serie de valores; al comparar diferentes pixels estamos midiendo el
C     error en la media---).
C (2) Con el valor promedio de sigma escalado al numero de imagenes sumadas
C     en cada region, podemos crear una mascara. Para ello calculamos el numero
C     de veces sigma que un pixel debe estar por encima de la media para que
C     la probabilidad de que ello ocurra sea muy pequen~a, teniendo en cuenta
C     el numero de pixels presentes. En otras palabras, considerando el numero
C     total de pixels en cada region, existe un valor Times_Sigma tal que
C     que el numero de pixels que estadisticamente puede superar la media en
C     Times_Sigma es igual a un numero fijado por el usuario (normalmente 1).
C..............................................................................
C compute mean +/- 3 sigma in the regions with different numbers of coadded
C frames (using the last available mask)
        DO NC=1,9
          CALL COMPHISTOFF(NC,NXMAXB9_,NYMAXB9_,NEWBUFF2,
     +     FFMEAN(NC),FFSIGMA(NC),NP(NC))
        END DO
C..............................................................................
C compute weighted FFSIGMA using the number of pixels
        SUMNP=0
        DO NC=1,9
          SUMNP=SUMNP+NP(NC)
        END DO
        WRITE(*,100) '>>> Test in total no. of pixels (? = ?, diff): '
C los siguientes numeros no tienen por que ser identicos: la estadistica
C se realiza usando una tolerancia en numero de pixels TOLNPIXEL, con lo
C que habra pixels que no coincidan con un numero entero dentro del intervalo
C establecido (mientras que la diferencia sea pequen~a, la estadistica seguira
C siendo muy fiable)
        WRITE(*,*) NXMAXB9_*NYMAXB9_,SUMNP+NM2,
     +   NXMAXB9_*NYMAXB9_-SUMNP+NM2
        FFSIGMA1=0.
        DO NC=1,9
          FFSIGMA1=FFSIGMA1+FFSIGMA(NC)*SQRT(REAL(NC))*
     +     REAL(NP(NC))/REAL(SUMNP)
        END DO
        DO NC=1,9
          FFSIGMANOR(NC)=FFSIGMA1/SQRT(REAL(NC))
        END DO
C..............................................................................
C compute histogram limits to be employed
        HMIN=+1.E30
        HMAX=-1.E30
        DO NC=1,9
          IF(NP(NC).GT.0)THEN
            IF(FFMEAN(NC)-3.*FFSIGMA(NC).LT.HMIN) 
     +       HMIN=FFMEAN(NC)-3.*FFSIGMA(NC)
            IF(FFMEAN(NC)+3.*FFSIGMA(NC).GT.HMAX) 
     +       HMAX=FFMEAN(NC)+3.*FFSIGMA(NC)
          END IF
        END DO
C..............................................................................
C compute histograms in regions with different number of coadded frames (note
C that we use HMIN and HMAX as histogram limits)
        DBIN=(HMAX-HMIN)/REAL(NBINMAX)
        DO K=1,NBINMAX
          XPLOT(K)=HMIN+REAL(K-1)*DBIN
        END DO
        XMIN=HMIN
        XMAX=HMAX
        DX=XMAX-XMIN
        XMIN=XMIN-DX/20.
        XMAX=XMAX+DX/20.
C..............................................................................
C compute histograms and fit gaussians to mean +/- 3 sigma
        DO NC=1,9
          CALL COMPHISTO(NC,NXMAXB9_,NYMAXB9_,NEWBUFF2,HMIN,DBIN,
     +     NPIX(1,NC))
        END DO
C..............................................................................
C plot histograms
        IF(CINTERACTIVE.EQ.'y')THEN
          WRITE(*,100) 'Press <CR> to plot histograms...'
          READ(*,*)
          IF(LECHO) WRITE(*,*)
        ELSE
          WRITE(*,101) 'Plotting histograms...'
        END IF
        CALL RPGERASW(0.4,1.0,0.0,0.8,0)
        WRITE(*,200)
        DO NC=1,9
          CALL MINIHISTO(NC,XMIN,XMAX,XPLOT,NPIX(1,NC),DBIN,
     +     FFMEAN(NC),FFSIGMA(NC),.FALSE.)
          IF(NP(NC).GT.0)THEN
            TIMESSIGMA(NC)=FNORTIPI(FPIXELTHRESHOLD/REAL(NPIXREGION))
            WRITE(*,*) NC,NP(NC),FFMEAN(NC),FFSIGMA(NC),
     +       FFSIGMA(NC)*SQRT(REAL(NC)),FFSIGMANOR(NC),
     +       TIMESSIGMA(NC)
          ELSE
            TIMESSIGMA(NC)=0.
          END IF
        END DO
        IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
          DO NC=1,9
            CALL MINIHISTO(NC,XMIN,XMAX,XPLOT,NPIX(1,NC),DBIN,
     +       FFMEAN(NC),FFSIGMA(NC),.TRUE.)
          END DO
          CALL PGCLOS(IDNEW)
          CALL PGSLCT(IDOLD)
        END IF
        WRITE(*,*)
        WRITE(*,101) '#1: no. of coadded frames'
        WRITE(*,101) '#2: number of pixels'
        WRITE(*,101) '#3: mean  (removing +/- 3 sigma)'
        WRITE(*,101) '#4: error on the mean (removing +/- 3 sigma)'
        WRITE(*,101) '#5: error on the mean scaled to 1 frame'
        WRITE(*,101) '#6: weighted-<error on the mean> scaled to '//
     +   'no. of frames'
        WRITE(*,100) '#7: Times_Sigma to have the '
        WRITE(*,101) '(threshold no. of pixels)/image'
        WRITE(*,200)
        DO NC=1,9
          XP(NC)=REAL(NC)
          YP(NC)=FFSIGMA(NC)
        END DO
        DO NC=10,18
          XP(NC)=REAL(NC-9)
          YP(NC)=FFSIGMA(NC-9)*SQRT(XP(NC))
        END DO
        CALL SUBPLOT(18,1,18,XP,YP,XP,YP,
     +               .TRUE.,.TRUE.,.FALSE.,.FALSE.,
     +               'no. of coadded frames','error on the mean',
     +               'blue: mean error, green: scaled mean error',
     +               2,-1,0.0)
        DO NC=1,9
          XP(NC)=REAL(NC)
          YP(NC)=FFSIGMA(NC)
        END DO
        CALL SUBPLOTBIS(9,1,9,XP,YP,XP,YP,.FALSE.,.FALSE.,5,21,1.5)
        CALL SUBPLOTBIS(9,1,9,XP,YP,XP,YP,.FALSE.,.FALSE.,4,101,1.)
        DO NC=1,9
          XP(NC)=REAL(NC)
          YP(NC)=FFSIGMA(NC)*SQRT(XP(NC))
        END DO
        CALL SUBPLOTBIS(18,1,9,XP,YP,XP,YP,.FALSE.,.FALSE.,3,17,1.5)
        XP(2)=XP(9)
        YP(1)=FFSIGMANOR(1)
        YP(2)=YP(1)
        CALL SUBPLOTBIS(2,1,2,XP,YP,XP,YP,.FALSE.,.FALSE.,6,101,1.)
C..............................................................................
C replot image with new BG and FG
        IF(CINTERACTIVE.EQ.'y')THEN
          WRITE(*,100) 'Press <CR> to replot image...'
          READ(*,*)
          IF(LECHO) WRITE(*,*)
        END IF
        BG=HMIN
        FG=HMAX
        CALL HISTOGRAM(NEWBUFF2)
        CALL SUBLOOK(.FALSE.,NEWBUFF2,.FALSE.)
        IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
          CALL SUBLOOK(.FALSE.,NEWBUFF2,.TRUE.)
          CALL PGCLOS(IDNEW)
          CALL PGSLCT(IDOLD)
        END IF
C..............................................................................
C replot image with number of pixels 
        IF(CINTERACTIVE.EQ.'y')THEN
          WRITE(*,100) 'Press <CR> to plot image with no. of pixels...'
          READ(*,*)
          IF(LECHO) WRITE(*,*)
        END IF
        NX1=1
        NX2=NAXIS(1,NEWBUFF2)
        NY1=1
        NY2=NAXIS(2,NEWBUFF2)
        DO I=NY1,NY2
          DO J=NX1,NX2
            IMAGEN_(J,I)=FPIXUSED(J,I)
          END DO
        END DO
        CALL PGIMAG(IMAGEN_,NXMAX,NYMAX,NX1,NX2,NY1,NY2,
     +   REAL(NFRAMES_),0.,TR)
        CALL PGQWIN(XW1,XW2,YW1,YW2)
        IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
          CALL PGENV(XW1,XW2,YW1,YW2,JUST,-2)
          CALL PGBOX('IBCTSN',0.0,0,'IBCTSN',0.0,0)
          CALL PGIMAG(IMAGEN_,NXMAX,NYMAX,NX1,NX2,NY1,NY2,0.,9.,TR)
          CALL PGCLOS(IDNEW)
          CALL PGSLCT(IDOLD)
        END IF
C..............................................................................
C replot image with new BG and FG
        IF(CINTERACTIVE.EQ.'y')THEN
          WRITE(*,100) 'Press <CR> to replot image...'
          READ(*,*)
          IF(LECHO) WRITE(*,*)
        END IF
        BG=HMIN
        FG=HMAX
        CALL HISTOGRAM(NEWBUFF2)
        CALL SUBLOOK(.FALSE.,NEWBUFF2,.FALSE.)
        IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
          CALL SUBLOOK(.FALSE.,NEWBUFF2,.TRUE.)
          CALL PGCLOS(IDNEW)
          CALL PGSLCT(IDOLD)
        END IF
C..............................................................................
C save last mask, compute new mask and compare if they are different
        DO I=1,NYMAXB9_
          DO J=1,NXMAXB9_
            MASKBOX9_(J,I)=MASKBOX9(J,I)
          END DO
        END DO
        CALL COMPMASK(NXMAXB9_,NYMAXB9_,NEWBUFF2,
     +   FFSIGMANOR(1),TIMESSIGMA(1),FTAILS,CINTERACTIVE,.TRUE.)
        IF(LPOST(CPOST,CBASEPOST,NPOST,IDOLD,IDNEW))THEN
          CALL SUBLOOK(.FALSE.,NEWBUFF2,.TRUE.)
          CALL COMPMASK(NXMAXB9_,NYMAXB9_,NEWBUFF2,
     +     FFSIGMANOR(1),TIMESSIGMA(1),FTAILS,CINTERACTIVE,.FALSE.)
          CALL PGSLCT(IDOLD)
        END IF
        LANY=.FALSE. !check if there is any difference between masks
        NM1=0 !number of pixels masked in previous mask
        NM2=0 !number of pixels masked in new mask
        IF(CINTERACTIVE.EQ.'y')THEN
          WRITE(*,100) 'Press <CR> to compare masks...'
          READ(*,*)
          IF(LECHO) WRITE(*,*)
        END IF
        DO I=1,NYMAXB9_
          DO J=1,NXMAXB9_
            IF(.NOT.MASKBOX9_(J,I)) NM1=NM1+1
            IF(.NOT.MASKBOX9(J,I)) NM2=NM2+1
            IF(MASKBOX9_(J,I).NEQV.MASKBOX9(J,I))THEN
              LANY=.TRUE.
              IF(.NOT.MASKBOX9_(J,I))THEN
                CALL PGSCI(4)
              ELSE
                CALL PGSCI(6)
              END IF
              CALL PGRECT(REAL(J)-0.5,REAL(J)+0.5,
     +                    REAL(I)-0.5,REAL(I)+0.5)
            END IF
          END DO
        END DO
        CALL PGSCI(1)
        IF(CPOST.EQ.'y')THEN
          CALL PGSLCT(IDNEW)
          DO I=1,NYMAXB9_
            DO J=1,NXMAXB9_
              IF(.NOT.MASKBOX9_(J,I)) NM1=NM1+1
              IF(.NOT.MASKBOX9(J,I)) NM2=NM2+1
              IF(MASKBOX9_(J,I).NEQV.MASKBOX9(J,I))THEN
                LANY=.TRUE.
                IF(.NOT.MASKBOX9_(J,I))THEN
                  CALL PGSCI(4)
                ELSE
                  CALL PGSCI(6)
                END IF
                CALL PGRECT(REAL(J)-0.5,REAL(J)+0.5,
     +                      REAL(I)-0.5,REAL(I)+0.5)
              END IF
            END DO
          END DO
          CALL PGCLOS(IDNEW)
          CALL PGSLCT(IDOLD)
        END IF
C..............................................................................
C comparison between new and previous mask
        WRITE(*,100) '>>> No. of masked pixels in previous mask: '
        WRITE(*,*) NM1
        WRITE(*,100) '>>> No. of masked pixels in new mask.....: '
        WRITE(*,*) NM2
        IF(LANY)THEN
          WRITE(*,101) '>>> New mask is different.'
          CCONT(1:1)=READC('Do you want to apply this new mask '//
     +     '(y/n/p=previous)','y','ynp')
        ELSE
          WRITE(*,101) '>>> New mask is identical to previous mask'
          CCONT(1:1)=
     +     READC('Do you want to apply this new mask anyway (y/n)',
     +     'y','yn')
        END IF
        IF(CCONT.EQ.'p')THEN
          DO I=1,NYMAXB9_
            DO J=1,NXMAXB9_
              MASKBOX9(J,I)=MASKBOX9_(J,I)
            END DO
          END DO
          CCONT='y'
        END IF
        IF(CCONT.EQ.'y')THEN
          CCONT(1:1)=READC('Change parameters in new iteration (y/n)',
     +     'n','yn')
          IF(CCONT.EQ.'y')THEN
            WRITE(CDUMMY,*) FPIXELTHRESHOLD
            FPIXELTHRESHOLD=READF(
     +       'Threshold in no. of pixels > Times_Sigma',CDUMMY)
            WRITE(CDUMMY,'(F8.6)') FTAILS
            FTAILS=READF('Factor to search for extended mask'//
     +       ' pixels (0 <= factor <= 1)',CDUMMY)
            CINTERACTIVE(1:1)=
     +       READC('Run procedure in interactive mode (y/m/n)',
     +       CINTERACTIVE,'ymn')
            IF(CINTERACTIVE.NE.'y')THEN
              CSCALE1(1:1)=
     +         READC('Fit#1 [c]onstant, [s]urface or [n]one'//
     +         ' (c/s/n)',CSCALE1,'csn')
              IF(CSCALE1.EQ.'s')THEN
                WRITE(CDUMMY,*) GX1
                GX1=READILIM('Polynomial degree in X direction',
     +           CDUMMY,0,INT(SQRT(REAL(MAXNCOEFF)))-1)
                WRITE(CDUMMY,*) GY1
                GY1=READILIM('Polynomial degree in Y direction',
     +           CDUMMY,0,INT(SQRT(REAL(MAXNCOEFF)))-1)
                CMEDIANA1(1:1)=
     +           READC('Refine fit minimizing median (y/n)',
     +           CMEDIANA1,'yn')
                IF(CMEDIANA1.EQ.'y')THEN
                  WRITE(CDUMMY,*) YRMSTOL1
                  YRMSTOL1=READF('YRMSTOL for DOWNHILL',CDUMMY)
                END IF
                CREFINE1(1:1)=
     +           READC('Refine fit removing deviating points'
     +           //' (y/n)',CREFINE1,'yn')
              END IF
              IF(CSCALE1.NE.'n')THEN
                CSCALE1OPER(1:1)=
     +           READC('[s]ubtract fit or [d]ivide by fit (s/d)',
     +           CSCALE1OPER,'sd')
                IF(CSCALE1OPER.EQ.'d')THEN
                  WRITE(*,101) 'WARNING: normalization is '//
     +             'recomended to preserve number of counts'
                  CNOR1(1:1)=READC('Normalize fit (y/n)','y','yn')
                ELSE
                  CNOR1='n'
                END IF
              END IF
c
              CSCALE2(1:1)=
     +         READC('Fit#2 [c]onstant, [s]urface or [n]one'//
     +         ' (c/s/n)',CSCALE2,'csn')
              IF(CSCALE2.EQ.'s')THEN
              WRITE(CDUMMY,*) GX2
                GX2=READILIM('Polynomial degree in X direction',
     +           CDUMMY,0,INT(SQRT(REAL(MAXNCOEFF)))-1)
                WRITE(CDUMMY,*) GY2
                GY2=READILIM('Polynomial degree in Y direction',
     +           CDUMMY,0,INT(SQRT(REAL(MAXNCOEFF)))-1)
                CMEDIANA2(1:1)=
     +           READC('Refine fit minimizing median (y/n)',
     +           CMEDIANA2,'yn')
                IF(CMEDIANA2.EQ.'y')THEN
                  WRITE(CDUMMY,*) YRMSTOL2
                  YRMSTOL2=READF('YRMSTOL for DOWNHILL',CDUMMY)
                END IF
                CREFINE2(1:1)=
     +           READC('Refine fit removing deviating points'
     +           //' (y/n)',CREFINE2,'yn')
              END IF
              IF(CSCALE2.NE.'n')THEN
                CSCALE2OPER(1:1)=
     +           READC('[s]ubtract fit or [d]ivide by fit (s/d)',
     +           CSCALE2OPER,'sd')
                IF(CSCALE2OPER.EQ.'d')THEN
                  WRITE(*,101) 'WARNING: normalization is '//
     +             'recomended to preserve number of counts'
                  CNOR2(1:1)=READC('Normalize fit (y/n)','y','yn')
                ELSE
                  CNOR2='n'
                END IF
              END IF
c
              CFILL(1:1)=
     +         READC('Fill masked regions with 2D polynomial '//
     +         'fits (y/n)',CFILL,'yn')
              IF(CFILL.EQ.'y')THEN
                WRITE(CDUMMY,*) GXFILL
                GXFILL=READILIM('Polynomial degree in X direction',
     +           CDUMMY,0,INT(SQRT(REAL(MAXNCOEFF)))-1)
                WRITE(CDUMMY,*) GYFILL
                GYFILL=READILIM('Polynomial degree in Y direction',
     +           CDUMMY,0,INT(SQRT(REAL(MAXNCOEFF)))-1)
                CMEDIANAFILL(1:1)=
     +           READC('Refine fit minimizing median (y/n)',
     +           CMEDIANAFILL,'yn')
                IF(CMEDIANAFILL.EQ.'y')THEN
                  WRITE(CDUMMY,*) YRMSTOLFILL
                  YRMSTOLFILL=READF('YRMSTOL for DOWNHILL',CDUMMY)
                END IF
                CREFINEFILL(1:1)=
     +           READC('Refine fit removing deviating '//
     +           'points (y/n)',CREFINEFILL,'yn')
              END IF
              CERRSKY(1:1)=
     +         READC('Compute error in median sky (y/n)',CERRSKY,
     +         'yn')
c
              IF(CINTERACTIVE.EQ.'n')THEN
                WRITE(*,*)
                CREJECT(1:1)=
     +           READC('Reject hot pixels (y/n)',CREJECT,'yn')
                IF(CREJECT.EQ.'y')THEN
                  WRITE(CDUMMY,*) FPIXELTHRESHOLDCR
                  FPIXELTHRESHOLDCR=READF('Threshold in no. of pixels'//
     +             ' > Times_Sigma for hot pixels',CDUMMY)
                END IF
              END IF
c
              WRITE(*,*)
              WRITE(*,100) '>>> Total no. of operations in stack: '
              WRITE(*,*) NOPER_STACK
              IF(NOPER_STACK.GT.0)THEN
                CSTACK(1:1)=READC('Apply operations in stack (y/n)',
     +           CSTACK,'yn')
              END IF
            END IF
          END IF
          GOTO 10
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
!!!900     WRITE(*,100) 'Press <CR> to return to initial buffer...'
        READ(*,*)
        IF(LECHO) WRITE(*,*)
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
200     FORMAT(79('-'))
        END
C
C******************************************************************************
C Dibuja lineas separando los diferentes frames en el mosaico de 768x768. Si
C NIMA o NQUA es/son distintos de cero, dibuja con otro color los cuadrantes
C apropiados
        SUBROUTINE DRAWSEPB9(NIMA,NQUA,NSIZE)
        IMPLICIT NONE
        INTEGER NIMA,NQUA,NSIZE
C
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
C
        INTEGER I
        INTEGER K,NQUAD
        INTEGER DI(NBOXMAX),DJ(NBOXMAX)
        INTEGER I1(4),I2(4),J1(4),J2(4)
C------------------------------------------------------------------------------
C Note: the pattern of the frames in box-9 is the following:
C       6 9 4
C       3 1 7
C       8 5 2
C offsets of each frame in the 768x768 composite mosaic
        DATA (DI(K),DJ(K),K=1,NBOXMAX) /
ccc     +   256,256,  !1
ccc     +   000,512,  !2
ccc     +   256,000,  !3
ccc     +   512,512,  !4
ccc     +   000,256,  !5
ccc     +   512,000,  !6
ccc     +   256,512,  !7
ccc     +   000,000,  !8
ccc     +   512,256/  !9
     +   1,1,  !1
     +   0,2,  !2
     +   1,0,  !3
     +   2,2,  !4
     +   0,1,  !5
     +   2,0,  !6
     +   1,2,  !7
     +   0,0,  !8
     +   2,1/  !9
C limites en pixels de cada cuadrante
        I1(1)=1
        I2(1)=NSIZE/2
        J1(1)=1
        J2(1)=NSIZE/2
c
        I1(2)=NSIZE/2+1
        I2(2)=NSIZE
        J1(2)=1
        J2(2)=NSIZE/2
c
        I1(3)=1
        I2(3)=NSIZE/2
        J1(3)=NSIZE/2+1
        J2(3)=NSIZE
c
        I1(4)=NSIZE/2+1
        I2(4)=NSIZE
        J1(4)=NSIZE/2+1
        J2(4)=NSIZE
C------------------------------------------------------------------------------
        CALL PGSCI(3)
        DO I=NSIZE,2*NSIZE,NSIZE
          CALL PGMOVE(REAL(I)+0.5,0.5)
          CALL PGDRAW(REAL(I)+0.5,REAL(3*NSIZE)+0.5)
          CALL PGMOVE(0.5,REAL(I)+0.5)
          CALL PGDRAW(REAL(3*NSIZE)+0.5,REAL(I)+0.5)
        END DO
        CALL PGSCI(1)
C------------------------------------------------------------------------------
        IF((NIMA.EQ.0).AND.(NQUA.EQ.0)) RETURN
C------------------------------------------------------------------------------
        CALL PGSCI(2)
        DO K=1,NBOXMAX
          IF((NIMA.EQ.0).OR.(NIMA.EQ.K))THEN
            DO NQUAD=1,4
              IF((NQUA.EQ.0).OR.(NQUA.EQ.NQUAD))THEN
                CALL PGMOVE(DJ(K)+J1(NQUAD)-0.5,DI(K)+I1(NQUAD)-0.5)
                CALL PGDRAW(DJ(K)+J2(NQUAD)+0.5,DI(K)+I1(NQUAD)-0.5)
                CALL PGDRAW(DJ(K)+J2(NQUAD)+0.5,DI(K)+I2(NQUAD)+0.5)
                CALL PGDRAW(DJ(K)+J1(NQUAD)-0.5,DI(K)+I2(NQUAD)+0.5)
                CALL PGDRAW(DJ(K)+J1(NQUAD)-0.5,DI(K)+I1(NQUAD)-0.5)
              END IF
            END DO
          END IF
        END DO
        CALL PGSCI(1)
C------------------------------------------------------------------------------
        END
C
C******************************************************************************
C La estadistica solo se realiza sobre los pixels no afectados por la mascara
        SUBROUTINE COMPHISTOFF(NHISTO,NXMAXB9_,NYMAXB9_,NEWBUFF,
     +   FMEAN,FSIGMA,NP)
        IMPLICIT NONE
C
        INCLUDE 'dimensions.inc'
C
        INTEGER NHISTO
        INTEGER NXMAXB9_,NYMAXB9_
        INTEGER NEWBUFF
        REAL FMEAN,FSIGMA
        INTEGER NP
C
        REAL TOLNPIXEL
        PARAMETER (TOLNPIXEL=0.1) !tolerance when looking for no. of pixels
C
        REAL FMEAN2
C
        INTEGER I,J
        REAL FPIXUSED(NXMAXB9,NYMAXB9)
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL PIXEL(NXMAX*NYMAX)
        LOGICAL MASKBOX9(NXMAXB9,NYMAXB9)
C
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKIMAGEN1_/PIXEL        !usamos este common para ahorra memoria
        COMMON/BLKFPIXUSED/FPIXUSED
        COMMON/BLKMASKBOX9/MASKBOX9
C------------------------------------------------------------------------------
        NP=0
        DO I=1,NYMAXB9_
          DO J=1,NXMAXB9_
            IF((FPIXUSED(J,I).GE.REAL(NHISTO)-TOLNPIXEL).AND.
     +         (FPIXUSED(J,I).LE.REAL(NHISTO)+TOLNPIXEL))THEN
              IF(MASKBOX9(J,I))THEN
                NP=NP+1
                PIXEL(NP)=IMAGEN(J,I,NEWBUFF)
              END IF
            END IF
          END DO
        END DO
C------------------------------------------------------------------------------
C calculamos la media +/- 3 sigma
        IF(NP.EQ.0)THEN
          FMEAN=0.
          FSIGMA=0.
        ELSEIF(NP.EQ.1)THEN
          FMEAN=PIXEL(1)
          FSIGMA=0.
        ELSE
          FMEAN=FMEAN2(NP,PIXEL,3.0,FSIGMA)
        END IF
C
        END
C
C******************************************************************************
C
        SUBROUTINE COMPHISTO(NHISTO,NXMAXB9_,NYMAXB9_,NEWBUFF,
     +   HMIN,DBIN,NPIX)
        IMPLICIT NONE
C
        INCLUDE 'dimensions.inc'
        INTEGER NBINMAX
        PARAMETER (NBINMAX=100)
C
        INTEGER NHISTO
        INTEGER NXMAXB9_,NYMAXB9_
        INTEGER NEWBUFF
        REAL HMIN
        REAL DBIN
        INTEGER NPIX(NBINMAX)
C
        REAL TOLNPIXEL
        PARAMETER (TOLNPIXEL=0.1) !tolerance when looking for no. of pixels
C
        INTEGER I,J,K
        REAL FPIXUSED(NXMAXB9,NYMAXB9)
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        LOGICAL MASKBOX9(NXMAXB9,NYMAXB9)
C
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKFPIXUSED/FPIXUSED
        COMMON/BLKMASKBOX9/MASKBOX9
C------------------------------------------------------------------------------
        DO K=1,NBINMAX
          NPIX(K)=0
        END DO
        DO I=1,NYMAXB9_
          DO J=1,NXMAXB9_
            IF((FPIXUSED(J,I).GE.REAL(NHISTO)-TOLNPIXEL).AND.
     +         (FPIXUSED(J,I).LE.REAL(NHISTO)+TOLNPIXEL))THEN
              IF(MASKBOX9(J,I))THEN
                K=NINT((IMAGEN(J,I,NEWBUFF)-HMIN)/DBIN)+1
                IF((K.GE.1).AND.(K.LE.NBINMAX)) NPIX(K)=NPIX(K)+1
              END IF
            END IF
          END DO
        END DO
C------------------------------------------------------------------------------
        END
C
C******************************************************************************
C
        SUBROUTINE MINIHISTO(NHISTO,XMIN,XMAX,XPLOT,NPIX,DBIN,
     +   FMEAN,FSIGMA,LPOSTSCRIPT)
        IMPLICIT NONE
C
        INTEGER NBINMAX
        PARAMETER (NBINMAX=100)
C
        INTEGER NHISTO
        REAL XMIN,XMAX
        REAL XPLOT(NBINMAX)
        INTEGER NPIX(NBINMAX)
        REAL DBIN
        REAL FMEAN,FSIGMA
        LOGICAL LPOSTSCRIPT
C
        INTEGER K
        INTEGER ISUM
        REAL YMIN,YMAX,DY
        REAL FPIX(NBINMAX)
        REAL XV1,XV2,YV1,YV2
        REAL XW1,XW2,YW1,YW2
        REAL OLD_CH
        REAL AMP,FDUM
        REAL YGAUSS(NBINMAX)
C------------------------------------------------------------------------------
        YMIN=-1.
        YMAX=YMIN
        DO K=1,NBINMAX
          IF(NPIX(K).EQ.0)THEN
            FPIX(K)=-1.
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
        CALL PGSCH(0.7)
        IF(NHISTO.EQ.1)THEN
          CALL PGSVP(0.45,0.59,0.58,0.78)
        ELSEIF(NHISTO.EQ.2)THEN
          CALL PGSVP(0.65,0.79,0.58,0.78)
        ELSEIF(NHISTO.EQ.3)THEN
          CALL PGSVP(0.85,0.99,0.58,0.78)
        ELSEIF(NHISTO.EQ.4)THEN
          CALL PGSVP(0.45,0.59,0.325,0.525)
        ELSEIF(NHISTO.EQ.5)THEN
          CALL PGSVP(0.65,0.79,0.325,0.525)
        ELSEIF(NHISTO.EQ.6)THEN
          CALL PGSVP(0.85,0.99,0.325,0.525)
        ELSEIF(NHISTO.EQ.7)THEN
          CALL PGSVP(0.45,0.59,0.07,0.27)
        ELSEIF(NHISTO.EQ.8)THEN
          CALL PGSVP(0.65,0.79,0.07,0.27)
        ELSEIF(NHISTO.EQ.9)THEN
          CALL PGSVP(0.85,0.99,0.07,0.27)
        END IF
        CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
        CALL PGBOX('BCTSN',0.0,0,'BCTSN',0.0,0)
        CALL PGMTXT('B',2.0,0.5,0.5,'signal')
        CALL PGMTXT('L',2.0,0.5,0.5,'Log[No.pixels]')
        IF(LPOSTSCRIPT)THEN
          CALL PGSCI(3)
        ELSE
          CALL PGSCI(7)
        END IF
        DO K=1,NBINMAX
          CALL PGMOVE(XPLOT(K)-0.5*DBIN,-1.0)
          CALL PGDRAW(XPLOT(K)-0.5*DBIN,FPIX(K))
          CALL PGDRAW(XPLOT(K)+0.5*DBIN,FPIX(K))
          CALL PGDRAW(XPLOT(K)+0.5*DBIN,-1.0)
        END DO
C------------------------------------------------------------------------------
C dibujamos la gaussiana +/- 3 sigma
        IF(FSIGMA.GT.0.0)THEN
          ISUM=0
          DO K=1,NBINMAX
            ISUM=ISUM+NPIX(K)
          END DO
          AMP=REAL(ISUM)*DBIN/(FSIGMA*SQRT(2.*3.141593))
          DO K=1,NBINMAX
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
          CALL PGLINE(NBINMAX,XPLOT,YGAUSS)
        END IF
C------------------------------------------------------------------------------
C recuperamos region de dibujo inicial
        CALL PGSCI(1)
        CALL PGSVP(XV1,XV2,YV1,YV2)
        CALL PGSWIN(XW1,XW2,YW1,YW2)
        CALL PGSCH(OLD_CH)
C
        END
C
C******************************************************************************
C Calcula la nueva mascara. El proceso se realiza en dos pasos. En el primero
C se calculan que pixels exceden la media FFMEAN en TIMESSIGMA FFSIGMANOR. En
C el segundo paso se refina la mascara, permitiendo la inclusion de pixels
C vecinos usando un critero mas ajustado.
C Si LASK=.FALSE., no pregunta nada (util cuando estamos generando el grafico
C en postscript)
        SUBROUTINE COMPMASK(NXMAXB9_,NYMAXB9_,NEWBUFF,
     +   FFSIGMANOR,TIMESSIGMA,FTAILS,CINTERACTIVE,LASK)
        IMPLICIT NONE
C
        INCLUDE 'dimensions.inc'
C
        INTEGER NXMAXB9_,NYMAXB9_
        INTEGER NEWBUFF
        REAL FFSIGMANOR
        REAL TIMESSIGMA
        REAL FTAILS
        CHARACTER*1 CINTERACTIVE
        LOGICAL LASK
C functions
        LOGICAL INSIDE
C variables
        INTEGER I,J
        INTEGER II,JJ,KK,LL
        INTEGER NPIX
        INTEGER NXYTOTAL
        REAL FPIXUSED(NXMAXB9,NYMAXB9)
        REAL PIXEL
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        LOGICAL MASKBOX9(NXMAXB9,NYMAXB9) !mascara
        LOGICAL MASKBOX9_(NXMAXB9,NYMAXB9) !extension de la mascara
        LOGICAL PFOUND
        LOGICAL LECHO
C
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKMASKBOX9/MASKBOX9
        COMMON/BLKFPIXUSED/FPIXUSED
        COMMON/BLKLECHO/LECHO
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C Primer paso: calculamos mascara inicial.
C------------------------------------------------------------------------------
        DO I=1,NYMAXB9_
          DO J=1,NXMAXB9_
            PIXEL=TIMESSIGMA*FFSIGMANOR/SQRT(FPIXUSED(J,I))
            MASKBOX9(J,I)=(IMAGEN(J,I,NEWBUFF).LE.PIXEL)
          END DO
        END DO
C dibujamos mascara inicial
        IF(LASK)THEN
          IF(CINTERACTIVE.EQ.'y')THEN
            WRITE(*,100) 'Press <CR> to overplot new mask...'
            READ(*,*)
            IF(LECHO) WRITE(*,*)
          END IF
        END IF
        CALL PGSCI(2)
        DO I=1,NYMAXB9_
          DO J=1,NXMAXB9_
            IF(.NOT.MASKBOX9(J,I))THEN
              CALL PGRECT(REAL(J)-0.5,REAL(J)+0.5,
     +                    REAL(I)-0.5,REAL(I)+0.5)
            END IF
          END DO
        END DO
        CALL PGSCI(1)
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C Segundo paso: nos movemos en coronas circulares alrededor de cada pixel ya
C incluido en la mascara para comprobar si los pixels vecinos pueden tambien
C an~adirse a la mascara (estamos usando el mismo algoritmo que la subroutina
C WALKER en cleanest).
C NOTA: recordar que MASK(J,I)=.TRUE. significa que el pixel (J,I) NO forma
C parte de la mascara.
C------------------------------------------------------------------------------
        DO I=1,NYMAXB9_
          DO J=1,NXMAXB9_
            MASKBOX9_(J,I)=.TRUE.
          END DO
        END DO
        NXYTOTAL=NXMAXB9_*NYMAXB9_
        DO I=1,NYMAXB9_
          DO J=1,NXMAXB9_
            IF(.NOT.MASKBOX9(J,I))THEN !si pertenece a la mascara
C inicializamos variables
              NPIX=0
              KK=0
              PFOUND=.TRUE.
              DO WHILE((NPIX.LT.NXYTOTAL).AND.PFOUND)
C de momento no hemos encontrado ningun nuevo pixel para la mascara extendida
                PFOUND=.FALSE.
                KK=KK+2
C posicion de arranque (sobre la diagonal)
                II=I-KK/2
                JJ=J-KK/2
C nos movemos hacia arriba
                DO LL=1,KK
                  II=II+1
                  IF(INSIDE(JJ,1,NXMAXB9_,II,1,NYMAXB9_))THEN !dentro
                    NPIX=NPIX+1
                    !no forma parte de la mascara
                    IF(MASKBOX9(JJ,II).AND.MASKBOX9_(JJ,II))THEN 
                      PIXEL=FTAILS*TIMESSIGMA*FFSIGMANOR/
     +                 SQRT(FPIXUSED(JJ,II))
                      IF(IMAGEN(JJ,II,NEWBUFF).GT.PIXEL)THEN
                        MASKBOX9_(JJ,II)=.FALSE.
                        PFOUND=.TRUE.
                      END IF
                    END IF
                  END IF
                END DO
C nos movemos hacia la derecha
                DO LL=1,KK
                  JJ=JJ+1
                  IF(INSIDE(JJ,1,NXMAXB9_,II,1,NYMAXB9_))THEN !dentro
                    NPIX=NPIX+1
                    !no forma parte de la mascara
                    IF(MASKBOX9(JJ,II).AND.MASKBOX9_(JJ,II))THEN 
                      PIXEL=FTAILS*TIMESSIGMA*FFSIGMANOR/
     +                 SQRT(FPIXUSED(JJ,II))
                      IF(IMAGEN(JJ,II,NEWBUFF).GT.PIXEL)THEN
                        MASKBOX9_(JJ,II)=.FALSE.
                        PFOUND=.TRUE.
                      END IF
                    END IF
                  END IF
                END DO
C nos movemos hacia abajo
                DO LL=1,KK
                  II=II-1
                  IF(INSIDE(JJ,1,NXMAXB9_,II,1,NYMAXB9_))THEN !dentro
                    NPIX=NPIX+1
                    !no forma parte de la mascara
                    IF(MASKBOX9(JJ,II).AND.MASKBOX9_(JJ,II))THEN 
                      PIXEL=FTAILS*TIMESSIGMA*FFSIGMANOR/
     +                 SQRT(FPIXUSED(JJ,II))
                      IF(IMAGEN(JJ,II,NEWBUFF).GT.PIXEL)THEN
                        MASKBOX9_(JJ,II)=.FALSE.
                        PFOUND=.TRUE.
                      END IF
                    END IF
                  END IF
                END DO
C nos movemos hacia la izquierda
                DO LL=1,KK
                  JJ=JJ-1
                  IF(INSIDE(JJ,1,NXMAXB9_,II,1,NYMAXB9_))THEN !dentro
                    NPIX=NPIX+1
                    !no forma parte de la mascara
                    IF(MASKBOX9(JJ,II).AND.MASKBOX9_(JJ,II))THEN 
                      PIXEL=FTAILS*TIMESSIGMA*FFSIGMANOR/
     +                 SQRT(FPIXUSED(JJ,II))
                      IF(IMAGEN(JJ,II,NEWBUFF).GT.PIXEL)THEN
                        MASKBOX9_(JJ,II)=.FALSE.
                        PFOUND=.TRUE.
                      END IF
                    END IF
                  END IF
                END DO
              END DO
            END IF
          END DO
        END DO
C dibujamos mascara adicional
        IF(LASK)THEN
          IF(CINTERACTIVE.EQ.'y')THEN
            WRITE(*,100) 'Press <CR> to overplot extended mask...'
            READ(*,*)
            IF(LECHO) WRITE(*,*)
          END IF
        END IF
        CALL PGSCI(3)
        DO I=1,NYMAXB9_
          DO J=1,NXMAXB9_
            IF(.NOT.MASKBOX9_(J,I))THEN
              CALL PGRECT(REAL(J)-0.5,REAL(J)+0.5,
     +                    REAL(I)-0.5,REAL(I)+0.5)
            END IF
          END DO
        END DO
        CALL PGSCI(1)
C------------------------------------------------------------------------------
C Fusionamos ambas mascaras
        DO I=1,NYMAXB9_
          DO J=1,NXMAXB9_
            IF(.NOT.MASKBOX9_(J,I)) MASKBOX9(J,I)=.FALSE.
          END DO
        END DO
C------------------------------------------------------------------------------
100     FORMAT(A,$)
        END
C
C******************************************************************************
C Si el punto J0,I0 esta dentro del rectangulo definido por J1,J2,I1,I2
C la funcion devuelve .TRUE.
        LOGICAL FUNCTION INSIDE(J0,J1,J2,I0,I1,I2)
        IMPLICIT NONE
        INTEGER J0,J1,J2
        INTEGER I0,I1,I2
C
        INSIDE=.TRUE.
        IF(J0.LT.J1) GOTO 70
        IF(J0.GT.J2) GOTO 70
        IF(I0.LT.I1) GOTO 70
        IF(I0.GT.I2) GOTO 70
        RETURN
70      CONTINUE
        INSIDE=.FALSE.
        RETURN
        END
C
C******************************************************************************
C A partir de la mascara de la imagen Box-9 compuesta (normalmente de 312x312 
C pixels), genera las submascaras correspondientes a cada uno de los frames 
C individuales (cada uno de 256x256 pixels). Al utilizar offsets fraccionarios,
C un pixel de la mascara puede implicar llegar a inutilizar hasta 4 pixels
C (notar que tomamos la opcion de eliminar todo pixel "tocado" por la mascara).
        SUBROUTINE GIVESUBMASK(NFRAMES_,NEXTRA_XMINB9,NEXTRA_YMINB9,
     +   I0,J0,FI0,FJ0,NXMAXB9_,NYMAXB9_)
        IMPLICIT NONE
C
        INCLUDE 'dimensions.inc'
C
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
C
        INTEGER NFRAMES_
        INTEGER NEXTRA_XMINB9,NEXTRA_YMINB9
        INTEGER I0(NBOXMAX),J0(NBOXMAX)
        REAL FI0(NBOXMAX),FJ0(NBOXMAX)
        INTEGER NXMAXB9_,NYMAXB9_
C
        LOGICAL INSIDE
C
        INTEGER I,J,K
        INTEGER DII1,DII2,DJJ1,DJJ2
        LOGICAL LOK1,LOK2,LOK3,LOK4
C
        LOGICAL MASKBOX9(NXMAXB9,NYMAXB9)
        LOGICAL SUBMASKBOX9(256,256,9)
        COMMON/BLKMASKBOX9/MASKBOX9
        COMMON/BLKSUBMASKBOX9/SUBMASKBOX9
C------------------------------------------------------------------------------
        DO K=1,NFRAMES_
          DO I=1,256
            DO J=1,256
              DJJ1=NEXTRA_XMINB9+J0(K) !ojo: NEXTRA_XMINB9 es <=0
              IF(FJ0(K).EQ.0.0)THEN
                DJJ2=0
              ELSEIF(FJ0(K).GT.0.0)THEN
                DJJ2=DJJ1+1
              ELSEIF(FJ0(K).LT.0.0)THEN
                DJJ2=DJJ1-1
              END IF
              DII1=NEXTRA_YMINB9+I0(K) !ojo: NEXTRA_YMINB9 es <=0
              IF(FI0(K).EQ.0.0)THEN
                DII2=0
              ELSEIF(FI0(K).GT.0.0)THEN
                DII2=DII1+1
              ELSEIF(FI0(K).LT.0.0)THEN
                DII2=DII1-1
              END IF
C..............................................................................
C Las siguientes lineas no deberian hacer falta. Si todo se ha hecho bien,
C nunca deberiamos estar fuera de limites.
              IF(.NOT.INSIDE(J-DJJ1,1,NXMAXB9_,I-DII1,1,NYMAXB9_))THEN
                WRITE(*,*) J0(K),FJ0(K),I0(K),FI0(K)
                WRITE(*,*) K,J,I,DJJ1,DJJ2,DII1,DII2,NXMAXB9_,NYMAXB9_
                STOP 'FATAL ERROR: LOK1'
              END IF
              IF(.NOT.INSIDE(J-DJJ1,1,NXMAXB9_,I-DII2,1,NYMAXB9_))THEN
                WRITE(*,*) J0(K),FJ0(K),I0(K),FI0(K)
                WRITE(*,*) K,J,I,DJJ1,DJJ2,DII1,DII2,NXMAXB9_,NYMAXB9_
                STOP 'FATAL ERROR: LOK2'
              END IF
              IF(.NOT.INSIDE(J-DJJ2,1,NXMAXB9_,I-DII1,1,NYMAXB9_))THEN
                WRITE(*,*) J0(K),FJ0(K),I0(K),FI0(K)
                WRITE(*,*) K,J,I,DJJ1,DJJ2,DII1,DII2,NXMAXB9_,NYMAXB9_
                STOP 'FATAL ERROR: LOK3'
              END IF
              IF(.NOT.INSIDE(J-DJJ2,1,NXMAXB9_,I-DII2,1,NYMAXB9_))THEN
                WRITE(*,*) J0(K),FJ0(K),I0(K),FI0(K)
                WRITE(*,*) K,J,I,DJJ1,DJJ2,DII1,DII2,NXMAXB9_,NYMAXB9_
                STOP 'FATAL ERROR: LOK4'
              END IF
C..............................................................................
              LOK1=MASKBOX9(J-DJJ1,I-DII1)
              IF(DII2.NE.0)THEN
                LOK2=MASKBOX9(J-DJJ1,I-DII2)
              ELSE
                LOK2=.TRUE.
              END IF
              IF(DJJ2.NE.0)THEN
                LOK3=MASKBOX9(J-DJJ2,I-DII1)
              ELSE
                LOK3=.TRUE.
              END IF
              IF((DII2.NE.0).AND.(DJJ2.NE.0))THEN
                LOK4=MASKBOX9(J-DJJ2,I-DII2)
              ELSE
                LOK4=.TRUE.
              END IF
              SUBMASKBOX9(J,I,K)=(LOK1.AND.LOK2.AND.LOK3.AND.LOK4)
            END DO
          END DO
        END DO
C
        END
C
C******************************************************************************
C Realiza un analisis estadistico detallado de un cuadrante de una imagen del
C Box9, para asi poder discriminar entre variaciones de sen~al debido a cambios
C en el nivel de BIAS o cambios en la sensibilidad/sen~al residual/variaciones
C del nivel del cielo. NBUFF es el buffer con la imagen actual a estudiar.
C NBUFFORIGINAL es el buffer que contiene la informacion original de la imagen
C (el numero de cuentas real, ligeramente alteradas por las operaciones del
C stack, para poder asi tener una buena estimacion del r.m.s. esperado a partir
C de la ganancia, el ruido de lectura y el numero de coadds).
        SUBROUTINE DETAILB9(NBUFF,NBUFFORIGINAL,NIMSTAT)
        IMPLICIT NONE
        INTEGER NBUFF,NBUFFORIGINAL,NIMSTAT
C
        INCLUDE 'dimensions.inc'
C
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
C
        INTEGER READI
        REAL READF
        REAL FMEAN0,FMEAN2
        CHARACTER*255 READC
C
        INTEGER NI,NJ,NK,NK1,NK2,NK3,NK4
        INTEGER I,J,K
        INTEGER I1,I2,J1,J2
        INTEGER DI(NBOXMAX),DJ(NBOXMAX)
        INTEGER NCOADDS
        INTEGER PGOPEN,IDNEW,IDOLD
        INTEGER NSYMB
        INTEGER LMEAN,LSIGMA
        REAL FMEAN(1024),FSIGMA(1024)
        REAL FMEANQ1(256),FSIGMAQ1(256)
        REAL FMEANQ2(256),FSIGMAQ2(256)
        REAL FMEANQ3(256),FSIGMAQ3(256)
        REAL FMEANQ4(256),FSIGMAQ4(256)
        REAL FMEAN_,FSIGMA_
        REAL PIXEL(64),PIXELORIGINAL(64)
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL GAIN,RNOISE
        REAL XP1(256),YP1(256)
        REAL XP2(256),YP2(256)
        REAL XP3(256),YP3(256)
        REAL XP4(256),YP4(256)
        REAL XMIN,XMAX,YMIN,YMAX,YMIN_,YMAX_
        REAL XV1,XV2,YV1,YV2
        REAL XW1,XW2,YW1,YW2
        REAL OLD_CH
        CHARACTER*1 CPOST
        CHARACTER*50 CMEAN,CSIGMA
        LOGICAL SUBMASKBOX9(256,256,9)
C
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKSUBMASKBOX9/SUBMASKBOX9
C------------------------------------------------------------------------------
C Note: the pattern of the frames in box-9 is the following:
C       6 9 4
C       3 1 7
C       8 5 2
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
C------------------------------------------------------------------------------
C Dibujamos valor experado de rms vs sen~al promedio
        GAIN=READF('Gain (e/ADU)..........','@')
        RNOISE=READF('Readout Noise (ADU)...','@')
        NCOADDS=READI('No. of coadds.........','@')
        CPOST(1:1)=READC('Create a Postscript plot (y/n)','n','yn')
C------------------------------------------------------------------------------
        NK=0
        NK1=0
        NK2=0
        NK3=0
        NK4=0
        DO NI=1,32
          I1=REAL(NI-1)*8+1
          I2=I1+7
          DO NJ=1,32
            J1=REAL(NJ-1)*8+1
            J2=J1+7
            IF(NI.LE.16)THEN
              IF(NJ.LE.16)THEN
                CALL PGSCI(5)
              ELSE
                CALL PGSCI(6)
              END IF
            ELSE
              IF(NJ.LE.16)THEN
                CALL PGSCI(3)
              ELSE
                CALL PGSCI(4)
              END IF
            END IF
            CALL PGMOVE(REAL(J1+DJ(NIMSTAT)),REAL(I1+DI(NIMSTAT)))
            CALL PGDRAW(REAL(J2+DJ(NIMSTAT)),REAL(I1+DI(NIMSTAT)))
            CALL PGDRAW(REAL(J2+DJ(NIMSTAT)),REAL(I2+DI(NIMSTAT)))
            CALL PGDRAW(REAL(J1+DJ(NIMSTAT)),REAL(I2+DI(NIMSTAT)))
            CALL PGDRAW(REAL(J1+DJ(NIMSTAT)),REAL(I1+DI(NIMSTAT)))
            K=0
            DO I=I1,I2
              DO J=J1,J2
                IF(SUBMASKBOX9(J,I,NIMSTAT))THEN
                  K=K+1
                  PIXEL(K)=IMAGEN(J+DJ(NIMSTAT),I+DI(NIMSTAT),NBUFF)
                  PIXELORIGINAL(K)=IMAGEN(J+DJ(NIMSTAT),I+DI(NIMSTAT),
     +             NBUFFORIGINAL)
                END IF
              END DO
            END DO
            IF(K.GT.5)THEN
              NK=NK+1
              FMEAN(NK)=FMEAN2(K,PIXEL,3.0,FSIGMA(NK))
              IF(NI.LE.16)THEN
                IF(NJ.LE.16)THEN
                  NK3=NK3+1
                  FMEANQ3(NK3)=FMEAN(NK)
                  FSIGMAQ3(NK3)=FSIGMA(NK)
                  XP3(NK3)=FMEAN2(K,PIXELORIGINAL,3.0,YP3(NK3))
                  YP3(NK3)=SQRT(XP3(NK3)/GAIN+REAL(NCOADDS)*
     +             RNOISE*RNOISE)
                  XP3(NK3)=FMEAN(NK) !damos el "cambiazo" a la abscisa
                ELSE
                  NK4=NK4+1
                  FMEANQ4(NK4)=FMEAN(NK)
                  FSIGMAQ4(NK4)=FSIGMA(NK)
                  XP4(NK4)=FMEAN2(K,PIXELORIGINAL,3.0,YP4(NK4))
                  YP4(NK4)=SQRT(XP4(NK4)/GAIN+REAL(NCOADDS)*
     +             RNOISE*RNOISE)
                  XP4(NK4)=FMEAN(NK) !damos el "cambiazo" a la abscisa
                END IF
              ELSE
                IF(NJ.LE.16)THEN
                  NK1=NK1+1
                  FMEANQ1(NK1)=FMEAN(NK)
                  FSIGMAQ1(NK1)=FSIGMA(NK)
                  XP1(NK1)=FMEAN2(K,PIXELORIGINAL,3.0,YP1(NK1))
                  YP1(NK1)=SQRT(XP1(NK1)/GAIN+REAL(NCOADDS)*
     +             RNOISE*RNOISE)
                  XP1(NK1)=FMEAN(NK) !damos el "cambiazo" a la abscisa
                ELSE
                  NK2=NK2+1
                  FMEANQ2(NK2)=FMEAN(NK)
                  FSIGMAQ2(NK2)=FSIGMA(NK)
                  XP2(NK2)=FMEAN2(K,PIXELORIGINAL,3.0,YP2(NK2))
                  YP2(NK2)=SQRT(XP2(NK2)/GAIN+REAL(NCOADDS)*
     +             RNOISE*RNOISE)
                  XP2(NK2)=FMEAN(NK) !damos el "cambiazo" a la abscisa
                END IF
              END IF
            END IF
          END DO
        END DO
        CALL PGSCI(1)
C------------------------------------------------------------------------------
C dibujamos resultado
        FMEAN_=FMEAN2(NK,FMEAN,3.0,FSIGMA_)
        XMIN=FMEAN_-4.*FSIGMA_
        XMAX=FMEAN_+4.*FSIGMA_
        FMEAN_=FMEAN2(NK,FSIGMA,3.0,FSIGMA_)
        YMIN=FMEAN_-4.*FSIGMA_
        YMAX=FMEAN_+4.*FSIGMA_
        DO K=1,4
          IF(K.EQ.1)THEN
            FMEAN_=FMEAN2(NK1,YP1,3.0,FSIGMA_)
          ELSEIF(K.EQ.2)THEN
            FMEAN_=FMEAN2(NK2,YP2,3.0,FSIGMA_)
          ELSEIF(K.EQ.3)THEN
            FMEAN_=FMEAN2(NK3,YP3,3.0,FSIGMA_)
          ELSEIF(K.EQ.4)THEN
            FMEAN_=FMEAN2(NK4,YP4,3.0,FSIGMA_)
          END IF
          YMIN_=FMEAN_-4*FSIGMA_
          YMAX_=FMEAN_+4*FSIGMA_
          YMIN=AMIN1(YMIN,YMIN_)
          YMAX=AMAX1(YMAX,YMAX_)
        END DO
C almacenamos region de dibujo actual
        CALL PGQVP(0,XV1,XV2,YV1,YV2)
        CALL PGQWIN(XW1,XW2,YW1,YW2)
        CALL PGQCH(OLD_CH)
        IF(CPOST.EQ.'y')THEN
          NSYMB=17
          CALL PGQID(IDOLD)
          CALL PGSAVE
          IDNEW=PGOPEN('?')
          CALL PGSLW(3)
          CALL PGSCF(2)
          CALL PGSVP(0.10,0.50,0.15,0.85)
        ELSE
          NSYMB=1
C borramos plot previo
          CALL RPGERASW(0.00,0.40,0.38,0.80,0)
          CALL PGSCH(0.6)
C definimos nueva region de dibujo para el plot #1
          CALL PGSVP(0.05,0.21,0.46,0.75)
        END IF
C dibujamos caja y datos
        CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
        CALL PGBOX('BCTSN',0.0,0,'BCTSN',0.0,0)
        CALL PGLABEL('signal','r.m.s. of the mean',' ')
        CALL PGMTXT('T',1.5,0.5,0.5,'(from data)')
        IF(NK1.GT.0)THEN
          CALL PGSCI(3)
          CALL PGPOINT(NK1,FMEANQ1,FSIGMAQ1,NSYMB)
          FMEAN_=FMEAN0(NK1,FSIGMAQ1,FSIGMA_)
          WRITE(CMEAN,'(F8.2)') FMEAN_
          CALL RMBLANK(CMEAN,CMEAN,LMEAN)
          WRITE(CSIGMA,'(F8.2)') FSIGMA_
          CALL RMBLANK(CSIGMA,CSIGMA,LSIGMA)
          CALL PGMTXT('T',-2.0,0.05,0.0,
     +     '<r.m.s.>='//CMEAN(1:LMEAN)//'\\(2233)'//CSIGMA(1:LSIGMA))
        END IF
        IF(NK2.GT.0)THEN
          CALL PGSCI(4)
          CALL PGPOINT(NK2,FMEANQ2,FSIGMAQ2,NSYMB)
          FMEAN_=FMEAN0(NK2,FSIGMAQ2,FSIGMA_)
          WRITE(CMEAN,'(F8.2)') FMEAN_
          CALL RMBLANK(CMEAN,CMEAN,LMEAN)
          WRITE(CSIGMA,'(F8.2)') FSIGMA_
          CALL RMBLANK(CSIGMA,CSIGMA,LSIGMA)
          CALL PGMTXT('T',-3.5,0.95,1.0,
     +     '<r.m.s.>='//CMEAN(1:LMEAN)//'\\(2233)'//CSIGMA(1:LSIGMA))
        END IF
        IF(NK3.GT.0)THEN
          CALL PGSCI(5)
          CALL PGPOINT(NK3,FMEANQ3,FSIGMAQ3,NSYMB)
          FMEAN_=FMEAN0(NK3,FSIGMAQ3,FSIGMA_)
          WRITE(CMEAN,'(F8.2)') FMEAN_
          CALL RMBLANK(CMEAN,CMEAN,LMEAN)
          WRITE(CSIGMA,'(F8.2)') FSIGMA_
          CALL RMBLANK(CSIGMA,CSIGMA,LSIGMA)
          CALL PGMTXT('B',-3.5,0.05,0.0,
     +     '<r.m.s.>='//CMEAN(1:LMEAN)//'\\(2233)'//CSIGMA(1:LSIGMA))
        END IF
        IF(NK4.GT.0)THEN
          CALL PGSCI(6)
          CALL PGPOINT(NK4,FMEANQ4,FSIGMAQ4,NSYMB)
          FMEAN_=FMEAN0(NK4,FSIGMAQ4,FSIGMA_)
          WRITE(CMEAN,'(F8.2)') FMEAN_
          CALL RMBLANK(CMEAN,CMEAN,LMEAN)
          WRITE(CSIGMA,'(F8.2)') FSIGMA_
          CALL RMBLANK(CSIGMA,CSIGMA,LSIGMA)
          CALL PGMTXT('B',-2.0,0.95,1.0,
     +     '<r.m.s.>='//CMEAN(1:LMEAN)//'\\(2233)'//CSIGMA(1:LSIGMA))
        END IF
        CALL PGSCI(1)
C definimos nueva region de dibujo para el plot #2
        IF(CPOST.EQ.'y')THEN
          CALL PGSVP(0.50,0.90,0.15,0.85)
        ELSE
          CALL PGSVP(0.21,0.37,0.46,0.75)
        END IF
C dibujamos caja y datos
        CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
        CALL PGBOX('BCTSN',0.0,0,'BCTSM',0.0,0)
        CALL PGLABEL('signal',' ',' ')
        CALL PGMTXT('T',1.5,0.5,0.5,'(from gain and RN)')
        IF(NK1.GT.0)THEN
          CALL PGSCI(3)
          CALL PGPOINT(NK1,XP1,YP1,NSYMB)
          FMEAN_=FMEAN0(NK1,YP1,FSIGMA_)
          WRITE(CMEAN,'(F8.2)') FMEAN_
          CALL RMBLANK(CMEAN,CMEAN,LMEAN)
          WRITE(CSIGMA,'(F8.2)') FSIGMA_
          CALL RMBLANK(CSIGMA,CSIGMA,LSIGMA)
          CALL PGMTXT('T',-2.0,0.05,0.0,
     +     '<r.m.s.>='//CMEAN(1:LMEAN)//'\\(2233)'//CSIGMA(1:LSIGMA))
        END IF
        IF(NK2.GT.0)THEN
          CALL PGSCI(4)
          CALL PGPOINT(NK2,XP2,YP2,NSYMB)
          FMEAN_=FMEAN0(NK2,YP2,FSIGMA_)
          WRITE(CMEAN,'(F8.2)') FMEAN_
          CALL RMBLANK(CMEAN,CMEAN,LMEAN)
          WRITE(CSIGMA,'(F8.2)') FSIGMA_
          CALL RMBLANK(CSIGMA,CSIGMA,LSIGMA)
          CALL PGMTXT('T',-3.5,0.95,1.0,
     +     '<r.m.s.>='//CMEAN(1:LMEAN)//'\\(2233)'//CSIGMA(1:LSIGMA))
        END IF
        IF(NK3.GT.0)THEN
          CALL PGSCI(5)
          CALL PGPOINT(NK3,XP3,YP3,NSYMB)
          FMEAN_=FMEAN0(NK3,YP3,FSIGMA_)
          WRITE(CMEAN,'(F8.2)') FMEAN_
          CALL RMBLANK(CMEAN,CMEAN,LMEAN)
          WRITE(CSIGMA,'(F8.2)') FSIGMA_
          CALL RMBLANK(CSIGMA,CSIGMA,LSIGMA)
          CALL PGMTXT('B',-3.5,0.05,0.0,
     +     '<r.m.s.>='//CMEAN(1:LMEAN)//'\\(2233)'//CSIGMA(1:LSIGMA))
        END IF
        IF(NK4.GT.0)THEN
          CALL PGSCI(6)
          CALL PGPOINT(NK4,XP4,YP4,NSYMB)
          FMEAN_=FMEAN0(NK4,YP4,FSIGMA_)
          WRITE(CMEAN,'(F8.2)') FMEAN_
          CALL RMBLANK(CMEAN,CMEAN,LMEAN)
          WRITE(CSIGMA,'(F8.2)') FSIGMA_
          CALL RMBLANK(CSIGMA,CSIGMA,LSIGMA)
          CALL PGMTXT('B',-2.0,0.95,1.0,
     +     '<r.m.s.>='//CMEAN(1:LMEAN)//'\\(2233)'//CSIGMA(1:LSIGMA))
        END IF
        CALL PGSCI(1)
        IF(CPOST.EQ.'y')THEN
          CALL PGCLOS(IDNEW)
          WRITE(*,101) '   ...OK! Plot created and closed'
          CALL PGSLCT(IDOLD)
          CALL PGUNSA
        END IF
C recuperamos la region de dibujo inicial
        CALL PGSVP(XV1,XV2,YV1,YV2)
        CALL PGSWIN(XW1,XW2,YW1,YW2)
        CALL PGSCH(OLD_CH)
C------------------------------------------------------------------------------
101     FORMAT(A)
        END
C
C******************************************************************************
C Genera una imagen con el valor de la mediana en cada cuadrante para todos
C los frames del Box-9. La salida se realiza a traves de la variable global
C IMAGEN_(). En QMIN y QMAX retorna el maximo y el minimo en cada cuadrante (en
C este caso QMIN()=QMAX(), porque estamos ajustando una constante).
        SUBROUTINE SCALE_CONSTANT(NBUFF,NFRAMES_,NIMA,NQUA,QMIN,QMAX)
        IMPLICIT NONE
        INTEGER NBUFF
        INTEGER NFRAMES_
        INTEGER NIMA,NQUA
        REAL QMIN(36),QMAX(36)
C
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
C
        INCLUDE 'dimensions.inc'
C
        REAL FMEDIAN1
C
        INTEGER I,J,K,KK,NQUAD
        INTEGER I1(4),I2(4),J1(4),J2(4)
        INTEGER DI(NBOXMAX),DJ(NBOXMAX)
        REAL FMEDIAN
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL IMAGEN_(NXMAX,NYMAX)
        REAL PIXEL(128*128)
        LOGICAL SUBMASKBOX9(256,256,9)
C
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKIMAGEN1_/IMAGEN_
        COMMON/BLKSUBMASKBOX9/SUBMASKBOX9
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
C limites en pixels de cada cuadrante
        DATA (I1(NQUAD),I2(NQUAD),J1(NQUAD),J2(NQUAD),NQUAD=1,4) /
     +   001,128,001,128,
     +   129,256,001,128,
     +   001,128,129,256,
     +   129,256,129,256/
C------------------------------------------------------------------------------
C si no ajustamos toda la imagen (o si no es un verdadero box-9), 
C inicializamos a cero todo el buffer de salida
        IF((NFRAMES_.NE.9).OR.(NIMA.NE.0).OR.(NQUA.NE.0))THEN
          DO K=1,NBOXMAX
            DO NQUAD=1,4
              DO I=I1(NQUAD),I2(NQUAD)
                DO J=J1(NQUAD),J2(NQUAD)
                  IMAGEN_(J+DJ(K),I+DI(K))=0.
                END DO
              END DO
              QMIN((K-1)*4+NQUAD)=0.
              QMAX((K-1)*4+NQUAD)=0.
            END DO
          END DO
        END IF
C
        DO K=1,NFRAMES_
          IF((NIMA.EQ.0).OR.(NIMA.EQ.K))THEN
            DO NQUAD=1,4
              IF((NQUA.EQ.0).OR.(NQUA.EQ.NQUAD))THEN
                KK=0
                DO I=I1(NQUAD),I2(NQUAD)
                  DO J=J1(NQUAD),J2(NQUAD)
                    IF(SUBMASKBOX9(J,I,K))THEN
                      KK=KK+1
                      PIXEL(KK)=IMAGEN(J+DJ(K),I+DI(K),NBUFF)
                    END IF
                  END DO
                END DO
                IF(KK.GT.0)THEN
                  FMEDIAN=FMEDIAN1(KK,PIXEL)
                ELSE
                  FMEDIAN=0.0
                END IF
                DO I=I1(NQUAD),I2(NQUAD)
                  DO J=J1(NQUAD),J2(NQUAD)
                    IMAGEN_(J+DJ(K),I+DI(K))=FMEDIAN
                  END DO
                END DO
                QMIN((K-1)*4+NQUAD)=FMEDIAN
                QMAX((K-1)*4+NQUAD)=FMEDIAN
              END IF
            END DO
          END IF
        END DO
C------------------------------------------------------------------------------
        END
C
C******************************************************************************
C Genera una imagen con el valor de la mediana en cada cuadrante para todos
C los frames del Box-9. La salida se realiza a traves de la variable global
C IMAGEN_(). Si LMEDIANA=.TRUE., la rutina minimiza la mediana de los residuos,
C para lo cual utiliza el ajuste inicial por minimos cuadrados y refina los
C coeficientes mediante DOWNHILL. Si LREFINE=.TRUE., la rutina tambien elimina
C puntos que se desvian demasiado iterando el ajuste realizado cada vez que
C algun punto nuevo ha sido eliminado. En QMIN y QMAX retorna el maximo y el 
C minimo en cada cuadrante. En XSOLOUT retornan los coeficientes ajustados en 
C cada cuadrante.
        SUBROUTINE SCALE_SURFACE(NBUFF,NFRAMES_,NIMA,NQUA,GX,GY,
     +   FFSIGMANOR,TIMESSIGMA,
     +   LMEDIANA,YRMSTOL,LREFINE,QMIN,QMAX,XSOLOUT)
        IMPLICIT NONE
C
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
        INTEGER MAXNCOEFF
        PARAMETER (MAXNCOEFF=16)
C
        INTEGER NBUFF
        INTEGER NFRAMES_
        INTEGER NIMA,NQUA
        INTEGER GX,GY
        REAL FFSIGMANOR,TIMESSIGMA
        LOGICAL LMEDIANA
        REAL YRMSTOL
        LOGICAL LREFINE
        REAL QMIN(36),QMAX(36)
        REAL XSOLOUT(MAXNCOEFF,36)
C
        INCLUDE 'dimensions.inc'
C
        REAL YFUNK_SURFACE
        EXTERNAL YFUNK_SURFACE
C
        INTEGER I,J,K,KK,NQUAD
        INTEGER II,JJ,P,Q,L,M,L_,M_
        INTEGER NCOEFF
        INTEGER I1(4),I2(4),J1(4),J2(4)
        INTEGER DI(NBOXMAX),DJ(NBOXMAX)
        INTEGER ORDER(MAXNCOEFF),IOK,IPAR
        INTEGER NEVAL
        INTEGER NBUFF_,GX_,GY_
        INTEGER NDISPLAY
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL IMAGEN_(NXMAX,NYMAX)
        REAL X(128),Y(128)
        REAL A(MAXNCOEFF,MAXNCOEFF),B(MAXNCOEFF)
        REAL SCALEROW(MAXNCOEFF),XSOL(MAXNCOEFF)
        REAL DXSOL(MAXNCOEFF),XXSOL(MAXNCOEFF),DXXSOL(MAXNCOEFF)
        REAL FFACTOR
        LOGICAL SUBMASKBOX9(256,256,9)
        LOGICAL IFGOOD(128,128)
        LOGICAL LANY
C
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKIMAGEN1_/IMAGEN_
        COMMON/BLKSUBMASKBOX9/SUBMASKBOX9
        COMMON/BLKFUNKSURF1/NBUFF_,NQUAD,GX_,GY_,K
        COMMON/BLKFUNKSURF2/IFGOOD
        COMMON/BLKFUNKSURF3/X,Y
C------------------------------------------------------------------------------
C Note: the pattern of the frames in box-9 is the following:
C       6 9 4
C       3 1 7
C       8 5 2
C offsets of each frame in the 768x768 composite mosaic
        DATA (DI(I),DJ(I),I=1,NBOXMAX) /
     +   256,256,  !1
     +   000,512,  !2
     +   256,000,  !3
     +   512,512,  !4
     +   000,256,  !5
     +   512,000,  !6
     +   256,512,  !7
     +   000,000,  !8
     +   512,256/  !9
C limites en pixels de cada cuadrante
        DATA (I1(I),I2(I),J1(I),J2(I),I=1,4) /
     +   001,128,001,128,
     +   129,256,001,128,
     +   001,128,129,256,
     +   129,256,129,256/
C------------------------------------------------------------------------------
C proteccion
        NCOEFF=(GX+1)*(GY+1)
        IF(NCOEFF.GT.MAXNCOEFF)THEN
          WRITE(*,100) 'FATAL ERROR in subroutine SCALE_SURFACE:'
          WRITE(*,101) 'NCOEFF.GT.MAXNCOEFF'
          WRITE(*,100) 'GX,GY,NCOEFF,MAXNCOEFF: '
          WRITE(*,*) GX,GY,NCOEFF,MAXNCOEFF
          STOP
        END IF
C duplicamos variables para evitar problemas con los COMMONs
        NBUFF_=NBUFF
        GX_=GX
        GY_=GY
C si no ajustamos toda la imagen (o si no es un verdadero box-9), 
C inicializamos a cero todo el buffer de salida
        IF((NFRAMES_.NE.9).OR.(NIMA.NE.0).OR.(NQUA.NE.0))THEN
          DO K=1,NBOXMAX
            DO NQUAD=1,4
              DO I=I1(NQUAD),I2(NQUAD)
                DO J=J1(NQUAD),J2(NQUAD)
                  IMAGEN_(J+DJ(K),I+DI(K))=0.
                END DO
              END DO
              QMIN((K-1)*4+NQUAD)=0.
              QMAX((K-1)*4+NQUAD)=0.
              DO KK=1,NCOEFF
                XSOLOUT(KK,(K-1)*4+NQUAD)=0.
              END DO
            END DO
          END DO
        END IF
C normalizamos el recorrido de las variables X e Y al intervalo [-1,1]
        DO J=1,128
          X(J)=(REAL(J)-64.5)/63.5
        END DO
        DO I=1,128
          Y(I)=(REAL(I)-64.5)/63.5
        END DO
C------------------------------------------------------------------------------
        WRITE(*,101) '---> System:'
        NDISPLAY=0
        DO K=1,NFRAMES_
          IF((NIMA.EQ.0).OR.(NIMA.EQ.K))THEN
            DO NQUAD=1,4
              IF((NQUA.EQ.0).OR.(NQUA.EQ.NQUAD))THEN
C..............................................................................
C inicialmente usamos todos los puntos del cuadrante
                DO I=1,128
                  DO J=1,128
                    IFGOOD(J,I)=.TRUE.
                  END DO
                END DO
C..............................................................................
C superponemos mascara
                CALL PGSCI(2)
                DO I=I1(NQUAD),I2(NQUAD)
                  DO J=J1(NQUAD),J2(NQUAD)
                    IF(.NOT.SUBMASKBOX9(J,I,K))THEN
                      CALL PGPOINT(1,REAL(J+DJ(K)),REAL(I+DI(K)),1)
                    END IF
                  END DO
                END DO
                CALL PGSCI(1)
C..............................................................................
C definimos un sistema de (GX+1)*(GY+1) ecuaciones con (GX+1)*(GY+1) incognitas
C II es el numero de ecuacion, y JJ el numero de incognita, por lo que
C el sistema de ecuaciones es A(II,JJ) * X = B(II)
10                II=0
                DO P=0,GX
                  DO Q=0,GY
                    II=II+1
                    JJ=0
                    DO I=0,GX
                      DO J=0,GY
                        JJ=JJ+1
                        A(II,JJ)=0.
                        DO L=J1(NQUAD),J2(NQUAD)
                          L_=L-J1(NQUAD)+1
                          DO M=I1(NQUAD),I2(NQUAD)
                            M_=M-I1(NQUAD)+1
                            IF(SUBMASKBOX9(L,M,K).AND.IFGOOD(L_,M_))THEN
                              FFACTOR=1.
                              IF(I+P.NE.0)FFACTOR=FFACTOR*(X(L_)**(I+P))
                              IF(J+Q.NE.0)FFACTOR=FFACTOR*(Y(M_)**(J+Q))
                              A(II,JJ)=A(II,JJ)+FFACTOR
                            END IF
                          END DO
                        END DO
                      END DO
                    END DO
                    B(II)=0.
                    DO L=J1(NQUAD),J2(NQUAD)
                      L_=L-J1(NQUAD)+1
                      DO M=I1(NQUAD),I2(NQUAD)
                        M_=M-I1(NQUAD)+1
                        IF(SUBMASKBOX9(L,M,K).AND.IFGOOD(L_,M_))THEN
                          FFACTOR=IMAGEN(L+DJ(K),M+DI(K),NBUFF)
                          IF(P.NE.0) FFACTOR=FFACTOR*(X(L_)**(P))
                          IF(Q.NE.0) FFACTOR=FFACTOR*(Y(M_)**(Q))
                          B(II)=B(II)+FFACTOR
                        END IF
                      END DO
                    END DO
                  END DO
                END DO
C..............................................................................
C resolvemos el sistema de ecuaciones
                WRITE(*,'(1X,I1,A1,I1,$)') K,'-',NQUAD
                NDISPLAY=NDISPLAY+1
                IF(NDISPLAY.EQ.19)THEN
                  WRITE(*,*)
                  NDISPLAY=0
                END IF
                CALL LUDCMP(A,NCOEFF,MAXNCOEFF,ORDER,SCALEROW,IOK,IPAR)
                CALL LUSOLV(A,NCOEFF,MAXNCOEFF,ORDER,SCALEROW,B,XSOL)
C..............................................................................
C si se requiere, podemos refinar el ajuste minimizando la mediana
                IF(LMEDIANA)THEN
                  DO KK=1,NCOEFF
                    DXSOL(KK)=XSOL(KK)/10.
                    IF(DXSOL(KK).EQ.0.0) DXSOL(KK)=0.1
                  END DO
                  CALL DOWNHILL(NCOEFF,XSOL,DXSOL,YFUNK_SURFACE,
     +             1.0,0.5,2.0,YRMSTOL,XXSOL,DXXSOL,NEVAL,100)
                  DO KK=1,NCOEFF
                    XSOL(KK)=XXSOL(KK)
                  END DO
                END IF
C..............................................................................
C calculamos la superficie ajustada
                DO M=I1(NQUAD),I2(NQUAD)
                  M_=M-I1(NQUAD)+1
                  DO L=J1(NQUAD),J2(NQUAD)
                    L_=L-J1(NQUAD)+1
                    IMAGEN_(L+DJ(K),M+DI(K))=0.
                    KK=0
                    DO I=0,GX
                      DO J=0,GY
                        KK=KK+1
                        FFACTOR=XSOL(KK)
                        IF(I.NE.0) FFACTOR=FFACTOR*(X(L_)**(I))
                        IF(J.NE.0) FFACTOR=FFACTOR*(Y(M_)**(J))
                        IMAGEN_(L+DJ(K),M+DI(K))=
     +                   IMAGEN_(L+DJ(K),M+DI(K))+FFACTOR
                      END DO
                        END DO
                  END DO
                END DO
C..............................................................................
C decidimos si eliminamos puntos
                IF(LREFINE)THEN
                  IF(FFSIGMANOR*TIMESSIGMA.GT.0.0)THEN
                    LANY=.FALSE. !no hemos eliminado ningun punto mas
                    DO I=I1(NQUAD),I2(NQUAD)
                      M_=I-I1(NQUAD)+1
                      DO J=J1(NQUAD),J2(NQUAD)
                        L_=J-J1(NQUAD)+1
                        IF(SUBMASKBOX9(J,I,K).AND.IFGOOD(L_,M_))THEN
                          IF(ABS(IMAGEN(J+DJ(K),I+DI(K),NBUFF)-
     +                     IMAGEN_(J+DJ(K),I+DI(K))).GT.
     +                     FFSIGMANOR*TIMESSIGMA)THEN
                            LANY=.TRUE.
                           IF(IMAGEN(J+DJ(K),I+DI(K),NBUFF).GT.
     +                      IMAGEN_(J+DJ(K),I+DI(K)))THEN
                             CALL PGSCI(4)
                           ELSE
                             CALL PGSCI(3)
                           END IF
                           CALL PGPOINT(1,REAL(J+DJ(K)),REAL(I+DI(K)),1)
                           IFGOOD(L_,M_)=.FALSE.
                          END IF
                        END IF
                      END DO
                    END DO
                    CALL PGSCI(1)
                    IF(LANY) GOTO 10
                  END IF
                END IF
C..............................................................................
C calculamos maximo y minimo en cada cuadrante
                QMIN((K-1)*4+NQUAD)=
     +           IMAGEN_(J1(NQUAD)+DJ(K),I1(NQUAD)+DI(K))
                QMAX((K-1)*4+NQUAD)=
     +           IMAGEN_(J1(NQUAD)+DJ(K),I1(NQUAD)+DI(K))
                DO M=I1(NQUAD),I2(NQUAD)
                  DO L=J1(NQUAD),J2(NQUAD)
                    IF(IMAGEN_(L+DJ(K),M+DI(K)).LT.QMIN((K-1)*4+NQUAD))
     +               QMIN((K-1)*4+NQUAD)=IMAGEN_(L+DJ(K),M+DI(K))
                    IF(IMAGEN_(L+DJ(K),M+DI(K)).GT.QMAX((K-1)*4+NQUAD))
     +               QMAX((K-1)*4+NQUAD)=IMAGEN_(L+DJ(K),M+DI(K))
                  END DO
                END DO
C..............................................................................
C almacenamos coeficientes ajustados en matriz de salida
                DO KK=1,NCOEFF
                  XSOLOUT(KK,(K-1)*4+NQUAD)=XSOL(KK)
                END DO
C..............................................................................
              END IF
            END DO
          END IF
        END DO
        IF(NDISPLAY.GT.0) WRITE(*,*)
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
        REAL FUNCTION YFUNK_SURFACE(XSOL)
        IMPLICIT NONE
C
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
        INTEGER MAXNCOEFF
        PARAMETER (MAXNCOEFF=100)
C
        REAL XSOL(MAXNCOEFF)
C
        INCLUDE 'dimensions.inc'
C
        REAL FMEDIAN1
C
        INTEGER I,J,L,M,L_,M_,KK
        INTEGER DI(NBOXMAX),DJ(NBOXMAX)
        INTEGER I1(4),I2(4),J1(4),J2(4)
        INTEGER NBUFF,NQUAD,GX,GY,K
        REAL FFACTOR
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL IMAGEN_(NXMAX,NYMAX)
        REAL X(128),Y(128)
        REAL PIXEL(128*128)
        LOGICAL SUBMASKBOX9(256,256,9)
        LOGICAL IFGOOD(128,128)
C
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKIMAGEN1_/IMAGEN_
        COMMON/BLKSUBMASKBOX9/SUBMASKBOX9
        COMMON/BLKFUNKSURF1/NBUFF,NQUAD,GX,GY,K
        COMMON/BLKFUNKSURF2/IFGOOD
        COMMON/BLKFUNKSURF3/X,Y
C------------------------------------------------------------------------------
C Note: the pattern of the frames in box-9 is the following:
C       6 9 4
C       3 1 7
C       8 5 2
C offsets of each frame in the 768x768 composite mosaic
        DATA (DI(I),DJ(I),I=1,NBOXMAX) /
     +   256,256,  !1
     +   000,512,  !2
     +   256,000,  !3
     +   512,512,  !4
     +   000,256,  !5
     +   512,000,  !6
     +   256,512,  !7
     +   000,000,  !8
     +   512,256/  !9
C limites en pixels de cada cuadrante
        DATA (I1(I),I2(I),J1(I),J2(I),I=1,4) /
     +   001,128,001,128,
     +   129,256,001,128,
     +   001,128,129,256,
     +   129,256,129,256/
C------------------------------------------------------------------------------
C calculamos la superficie ajustada
        DO M=I1(NQUAD),I2(NQUAD)
          M_=M-I1(NQUAD)+1
          DO L=J1(NQUAD),J2(NQUAD)
            L_=L-J1(NQUAD)+1
            IMAGEN_(L+DJ(K),M+DI(K))=0.
            KK=0
            DO I=0,GX
              DO J=0,GY
                KK=KK+1
                FFACTOR=XSOL(KK)
                IF(I.NE.0) FFACTOR=FFACTOR*(X(L_)**(I))
                IF(J.NE.0) FFACTOR=FFACTOR*(Y(M_)**(J))
                IMAGEN_(L+DJ(K),M+DI(K))=IMAGEN_(L+DJ(K),M+DI(K))+
     +           FFACTOR
              END DO
                END DO
          END DO
        END DO
C------------------------------------------------------------------------------
C calculamos residuos
        KK=0
        DO I=I1(NQUAD),I2(NQUAD)
          M_=I-I1(NQUAD)+1
          DO J=J1(NQUAD),J2(NQUAD)
            L_=J-J1(NQUAD)+1
            IF(SUBMASKBOX9(J,I,K).AND.IFGOOD(L_,M_))THEN
              KK=KK+1
              PIXEL(KK)=IMAGEN(J+DJ(K),I+DI(K),NBUFF)-
     +                  IMAGEN_(J+DJ(K),I+DI(K))
            END IF
          END DO
        END DO
C------------------------------------------------------------------------------
C calculamos mediana
        YFUNK_SURFACE=FMEDIAN1(KK,PIXEL)
        YFUNK_SURFACE=YFUNK_SURFACE*YFUNK_SURFACE !tomamos el cuadrado para
                                                  !buscar el minimo
C------------------------------------------------------------------------------
        END
C
C******************************************************************************
C Normaliza la sen~al de cada cuadrante al mismo factor constante para poder 
C asi sustraer luego la imagen al normalizar).
C La imagen de entrada es IMAGEN_(), y la de salida tambien.
        SUBROUTINE NORMSUMA(NFRAMES_,NIMA,NQUA)
        IMPLICIT NONE
        INTEGER NFRAMES_
        INTEGER NIMA,NQUA
C
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
C
        INCLUDE 'dimensions.inc'
C
        INTEGER I,J,K,NQUAD
        INTEGER I1(4),I2(4),J1(4),J2(4)
        INTEGER DI(NBOXMAX),DJ(NBOXMAX)
        INTEGER NSIZE
        REAL IMAGEN_(NXMAX,NYMAX)
        REAL FACTOR
        DOUBLE PRECISION DMEAN
C
        COMMON/BLKIMAGEN1_/IMAGEN_
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
C limites en pixels de cada cuadrante
        DATA (I1(NQUAD),I2(NQUAD),J1(NQUAD),J2(NQUAD),NQUAD=1,4) /
     +   001,128,001,128,
     +   129,256,001,128,
     +   001,128,129,256,
     +   129,256,129,256/
C------------------------------------------------------------------------------
C normalizamos cada cuadrante a cero
        DO K=1,NFRAMES_
          IF((NIMA.EQ.0).OR.(NIMA.EQ.K))THEN
            DO NQUAD=1,4
              IF((NQUA.EQ.0).OR.(NQUA.EQ.NQUAD))THEN
                DMEAN=0.D0
                DO I=I1(NQUAD),I2(NQUAD)
                  DO J=J1(NQUAD),J2(NQUAD)
                    DMEAN=DMEAN+IMAGEN_(J+DJ(K),I+DI(K))
                  END DO
                END DO
                NSIZE=(J2(NQUAD)-J1(NQUAD)+1)*(I2(NQUAD)-I1(NQUAD)+1)
                FACTOR=REAL(DMEAN/DBLE(NSIZE))
                DO I=I1(NQUAD),I2(NQUAD)
                  DO J=J1(NQUAD),J2(NQUAD)
                    IMAGEN_(J+DJ(K),I+DI(K))=IMAGEN_(J+DJ(K),I+DI(K))-
     +               FACTOR
                  END DO
                END DO
              END IF
            END DO
          END IF
        END DO
C------------------------------------------------------------------------------
        END
C
C******************************************************************************
C Normaliza la sen~al de cada cuadrante a uno para poder asi dividir luego 
C la imagen al normalizar).
C La imagen de entrada es IMAGEN_(), y la de salida tambien.
        SUBROUTINE NORMDIVI(NFRAMES_,NIMA,NQUA)
        IMPLICIT NONE
        INTEGER NFRAMES_
        INTEGER NIMA,NQUA
C
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
C
        INCLUDE 'dimensions.inc'
C
        INTEGER I,J,K,NQUAD
        INTEGER I1(4),I2(4),J1(4),J2(4)
        INTEGER DI(NBOXMAX),DJ(NBOXMAX)
        INTEGER NSIZE
        REAL FACTOR
        REAL IMAGEN_(NXMAX,NYMAX)
        DOUBLE PRECISION DMEAN
C
        COMMON/BLKIMAGEN1_/IMAGEN_
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
C limites en pixels de cada cuadrante
        DATA (I1(NQUAD),I2(NQUAD),J1(NQUAD),J2(NQUAD),NQUAD=1,4) /
     +   001,128,001,128,
     +   129,256,001,128,
     +   001,128,129,256,
     +   129,256,129,256/
C------------------------------------------------------------------------------
C normalizamos cada cuadrante a uno
        DO K=1,NFRAMES_
          IF((NIMA.EQ.0).OR.(NIMA.EQ.K))THEN
            DO NQUAD=1,4
              IF((NQUA.EQ.0).OR.(NQUA.EQ.NQUAD))THEN
                DMEAN=0.D0
                DO I=I1(NQUAD),I2(NQUAD)
                  DO J=J1(NQUAD),J2(NQUAD)
                    DMEAN=DMEAN+IMAGEN_(J+DJ(K),I+DI(K))
                  END DO
                END DO
                NSIZE=(J2(NQUAD)-J1(NQUAD)+1)*(I2(NQUAD)-I1(NQUAD)+1)
                FACTOR=DMEAN/DBLE(NSIZE)
                IF(FACTOR.GT.0.0)THEN
                  DO I=I1(NQUAD),I2(NQUAD)
                    DO J=J1(NQUAD),J2(NQUAD)
                      IMAGEN_(J+DJ(K),I+DI(K))=IMAGEN_(J+DJ(K),I+DI(K))-
     +                 FACTOR
                    END DO
                  END DO
                ELSE
                  DO I=I1(NQUAD),I2(NQUAD)
                    DO J=J1(NQUAD),J2(NQUAD)
                      IMAGEN_(J+DJ(K),I+DI(K))=1.0
                    END DO
                  END DO
                END IF
              END IF
            END DO
          END IF
        END DO
C------------------------------------------------------------------------------
        END
C
C******************************************************************************
C Calcula la media y la sigma para determinar BG y FG en plots, teniendo en
C cuenta que quiza no queremos utilizar todos los cuadrantes. El calculo no
C es exacto, porque utiliza un promedio de medias y un promedio de sigmas, pero
C para determinar cortes de la imagen, es suficiente.
        SUBROUTINE GIVEME_BGFG(NFRAMES_,NIMA,NQUA,BG,FG)
        IMPLICIT NONE
        INTEGER NFRAMES_
        INTEGER NIMA
        INTEGER NQUA
        REAL BG,FG
C
        INCLUDE 'dimensions.inc'
C
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
C
        INTEGER NQUAD
        INTEGER K,KK
        INTEGER DI(NBOXMAX),DJ(NBOXMAX)
        INTEGER I1(4),I2(4),J1(4),J2(4)
        REAL FMEAN,FSIGMA,FMEDIAN,FMIN,FMAX
        REAL IMAGEN_(NXMAX,NYMAX)
C
        COMMON/BLKIMAGEN1_/IMAGEN_
        COMMON/BLKESTADISTICA/FMEAN,FSIGMA,FMEDIAN,FMIN,FMAX
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
C limites en pixels de cada cuadrante
        DATA (I1(NQUAD),I2(NQUAD),J1(NQUAD),J2(NQUAD),NQUAD=1,4) /
     +   001,128,001,128,
     +   129,256,001,128,
     +   001,128,129,256,
     +   129,256,129,256/
C------------------------------------------------------------------------------
C si tenemos toda la imagen (y es un verdadero box-9), usamos el metodo tipico
        IF((NFRAMES_.EQ.9).AND.(NIMA.EQ.0).AND.(NQUA.EQ.0))THEN
          CALL STATISTIC(0,1,768,1,768,.FALSE.,.FALSE.,.TRUE.,
     +     0.0,.FALSE.)
          BG=FMEAN-5.*FSIGMA
          FG=FMEAN+5.*FSIGMA
        ELSE
C------------------------------------------------------------------------------
C si no tenemos toda la imagen (o no es un verdadero box-9), hacemos uso de 
C los extremos en los cuadrantes a utilizar
          KK=0
          DO K=1,NFRAMES_
            IF((NIMA.EQ.0).OR.(NIMA.EQ.K))THEN
              DO NQUAD=1,4
                IF((NQUA.EQ.0).OR.(NQUA.EQ.NQUAD))THEN
                  KK=KK+1
                  CALL STATISTIC(0,J1(NQUAD)+DJ(K),J2(NQUAD)+DJ(K),
     +             I1(NQUAD)+DI(K),I2(NQUAD)+DI(K),
     +             .FALSE.,.FALSE.,.TRUE.,0.0,.FALSE.)
                  IF(KK.EQ.1)THEN
                    BG=FMEAN-5.*FSIGMA
                    FG=FMEAN+5.*FSIGMA
                  ELSE
                    BG=AMIN1(BG,FMEAN-5.*FSIGMA)
                    FG=AMAX1(FG,FMEAN+5.*FSIGMA)
                  END IF
                END IF
              END DO
            END IF
          END DO
        END IF
C------------------------------------------------------------------------------
        IF(BG.EQ.FG)THEN
          BG=BG-1.0
          FG=FG+1.0
        END IF
C
        END
C
C******************************************************************************
C Calcula el BG y el FG en la imagen del buffer NBUFF, eliminando de la
C estadistica los pixels de la mascara
        SUBROUTINE BGFG_MASK(NFRAMES_,NBUFF,BG,FG)
        IMPLICIT NONE
        INTEGER NFRAMES_
        INTEGER NBUFF
        REAL BG,FG
C
        INCLUDE 'dimensions.inc'
C
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
C
        INTEGER I,J,K,KK
        INTEGER DI(NBOXMAX),DJ(NBOXMAX)
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        DOUBLE PRECISION DMEAN,DSIGMA
        LOGICAL SUBMASKBOX9(256,256,9)
C
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKSUBMASKBOX9/SUBMASKBOX9
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
C------------------------------------------------------------------------------
        DMEAN=0.D0
        DSIGMA=0.D0
C
        KK=0
        DO K=1,NFRAMES_
          DO I=1,256
            DO J=1,256
              IF(SUBMASKBOX9(J,I,K))THEN
                KK=KK+1
                DMEAN=DMEAN+DBLE(IMAGEN(J+DJ(K),I+DI(K),NBUFF))
              END IF
            END DO
          END DO
        END DO
C
        IF(KK.GT.0.)THEN
          DMEAN=DMEAN/DBLE(KK)
          IF(KK.GT.1)THEN
            DO K=1,NFRAMES_
              DO I=1,256
                DO J=1,256
                  IF(SUBMASKBOX9(J,I,K))THEN
                    DSIGMA=DSIGMA+
     +               (DBLE(IMAGEN(J+DJ(K),I+DI(K),NBUFF))-DMEAN)*
     +               (DBLE(IMAGEN(J+DJ(K),I+DI(K),NBUFF))-DMEAN)
                  END IF
                END DO
              END DO
            END DO
          END IF
          DSIGMA=SQRT(DSIGMA/DBLE(KK-1))
        END IF
C
        BG=REAL(DMEAN-5.*DSIGMA)
        FG=REAL(DMEAN+5.*DSIGMA)
        IF(BG.EQ.FG)THEN
          BG=BG-1.0
          FG=FG+1.0
        END IF
C
        END
C
C******************************************************************************
C
        LOGICAL FUNCTION LPOST(CPOST,CBASEPOST,NPOST,
     +   IDOLD,IDNEW)
        IMPLICIT NONE
        CHARACTER*1 CPOST
        CHARACTER*(*) CBASEPOST
        INTEGER NPOST
        INTEGER IDOLD,IDNEW
C
        INTEGER PGOPEN
        INTEGER TRUEBEG,TRUELEN
C
        INTEGER L1,L2
        CHARACTER*4 CNUMBER
        CHARACTER*80 FILEPOSTSCRIPT
        LOGICAL LECHO
C
        COMMON/BLKLECHO/LECHO
C------------------------------------------------------------------------------
        IF(CPOST.EQ.'n')THEN
          LPOST=.FALSE.
          RETURN
        END IF
C
        L1=TRUEBEG(CBASEPOST)
        L2=TRUELEN(CBASEPOST)
C
        NPOST=NPOST+1
        WRITE(CNUMBER,'(I4.4)') NPOST
        FILEPOSTSCRIPT=CBASEPOST(L1:L2)//CNUMBER(1:4)//'.ps/cps'
        IDNEW=PGOPEN(FILEPOSTSCRIPT)
        IF(IDNEW.LE.0)THEN
          WRITE(*,100) 'ERROR invalid graphics device: '
          WRITE(*,101) FILEPOSTSCRIPT
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          IF(LECHO) WRITE(*,*)
          CALL PGSLCT(IDOLD)
          NPOST=NPOST-1
          LPOST=.FALSE.
        ELSE
          LPOST=.TRUE.
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
