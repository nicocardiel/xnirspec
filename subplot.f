C Dibuja un plot de XP(N) y YP(N), con N entre N1 y N2
C EXP,EYP: errores en XP e YP (solo se dibujan si LEX y/o LEY son .TRUE.)
C LX,LY indica si hay que calcular limites automaticamente en X y/o en Y
C CX,CY,CG: etiquetas en eje X, Y y en Grafica
C ICOLOR: color para lineas/puntos
C NSYMB: numero de simbolo 
C        101,102,...=linea de tipo 1,2,.. con PGLINE
C        201,202,...=linea de tipo 1,2,.. con PGBIN
C CH: current height para simbolos
C
        SUBROUTINE SUBPLOT(N,N1,N2,XP,YP,EXP,EYP,
     +                     LX,LY,LEX,LEY,
     +                     CX,CY,CG,
     +                     ICOLOR,NSYMB,CH)
        IMPLICIT NONE
        INTEGER N
        INTEGER N1,N2
        REAL XP(N),YP(N)
        REAL EXP(N),EYP(N)
        LOGICAL LX,LY,LEX,LEY
        CHARACTER*(*) CX
        CHARACTER*(*) CY
        CHARACTER*(*) CG
        INTEGER ICOLOR
        INTEGER NSYMB
        REAL CH
C
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
C
        INTEGER NLAST,N1LAST,N2LAST
        INTEGER ICOLORLAST,NSYMBLAST
        INTEGER NB
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XPLAST(NXYMAX*NOVERSAMPMAX),YPLAST(NXYMAX*NOVERSAMPMAX)
        REAL CHLAST
        CHARACTER*2 CBUFF
        CHARACTER*255 CLABX,CLABY,CLABG
C
        INTEGER I
        REAL DX,DY
        REAL XV1,XV2,YV1,YV2
        REAL XW1,XW2,YW1,YW2
        REAL OLD_CH
C
        COMMON/BLKPLIMITS/XMIN,XMAX,YMIN,YMAX           !para guardar variables
        COMMON/BLKPNLAST/NLAST,N1LAST,N2LAST
        COMMON/BLKPLAST/XPLAST,YPLAST
        COMMON/BLKPLABELS/CLABX,CLABY,CLABG
        COMMON/BLKPCOLSYMB/ICOLORLAST,NSYMBLAST
        COMMON/BLKPCH/CHLAST
C------------------------------------------------------------------------------
C comenzamos buffering
        CALL PGBBUF
C almacenamos region de dibujo actual
        CALL PGQVP(0,XV1,XV2,YV1,YV2)
        CALL PGQWIN(XW1,XW2,YW1,YW2)
        CALL PGQCH(OLD_CH)
C borramos plot previo
        CALL RPGERASW(0.00,0.40,0.38,0.80,0)
C definimos nueva region de dibujo para el plot
        CALL PGSCH(0.8)
        CALL PGSVP(0.05,0.37,0.46,0.75)
C en caso necesario, calculamos limites
        IF(LX)THEN
          IF(LEX)THEN
            XMIN=XP(N1)-EXP(N1)
            XMAX=XP(N1)+EXP(N1)
            IF(N2.GT.N1)THEN
              DO I=N1+1,N2
                IF(XMIN.GT.XP(I)-EXP(I)) XMIN=XP(I)-EXP(I)
                IF(XMAX.LT.XP(I)+EXP(I)) XMAX=XP(I)+EXP(I)
              END DO
            END IF
          ELSE
            CALL FINDMML(N,N1,N2,XP,XMIN,XMAX)
          END IF
          DX=XMAX-XMIN
          IF(DX.EQ.0.0) DX=1.
          XMIN=XMIN-DX/20.
          XMAX=XMAX+DX/20.
        END IF
        IF(LY)THEN
          IF(LEY)THEN
            YMIN=YP(N1)-EYP(N1)
            YMAX=YP(N1)+EYP(N1)
            IF(N2.GT.N1)THEN
              DO I=N1+1,N2
                IF(YMIN.GT.YP(I)-EYP(I)) YMIN=YP(I)-EYP(I)
                IF(YMAX.LT.YP(I)+EYP(I)) YMAX=YP(I)+EYP(I)
              END DO
            END IF
          ELSE
            CALL FINDMML(N,N1,N2,YP,YMIN,YMAX)
          END IF
          DY=YMAX-YMIN
          IF(DY.EQ.0.0) DY=1.
          YMIN=YMIN-DY/20.
          YMAX=YMAX+DY/20.
        END IF
C dibujamos caja
        CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
        CALL PGBOX('BCTSN',0.0,0,'BCTSN',0.0,0)
        CALL PGLABEL(CX,CY,' ')
        CALL PGMTXT('T',1.0,0.5,0.5,CG)
        DO NB=1,NMAXBUFF/2
          CALL PGSCI(NB)
          WRITE(CBUFF,'(A1,I1)') '#',NB
          CALL PGMTXT('R',1.2,REAL(NB-1)/REAL(NMAXBUFF/2),0.0,CBUFF)
        END DO
        CALL PGSCI(1)
C dibujamos datos
        CALL PGSCI(ICOLOR)
        CALL PGSCH(CH)
        IF(NSYMB.LE.100)THEN
          CALL PGPOINT(N2-N1+1,XP(N1),YP(N1),NSYMB)
        ELSEIF(NSYMB.LE.105)THEN
          CALL PGSLS(NSYMB-100)
          CALL PGLINE(N2-N1+1,XP(N1),YP(N1))
          CALL PGSLS(1)
        ELSEIF(NSYMB.LE.205)THEN
          CALL PGSLS(NSYMB-200)
          CALL PGBIN(N2-N1+1,XP(N1),YP(N1),.TRUE.)
          CALL PGSLS(1)
        END IF
        IF(LEX)THEN
          DO I=N1,N2
            CALL PGERRX(1,XP(I)-EXP(I),XP(I)+EXP(I),YP(I),1.0)
          END DO
        END IF
        IF(LEY)THEN
          DO I=N1,N2
            CALL PGERRY(1,XP(I),YP(I)-EYP(I),YP(I)+EYP(I),1.0)
          END DO
        END IF
        CALL PGSCI(1)
C recuperamos region de dibujo inicial
        CALL PGSVP(XV1,XV2,YV1,YV2)
        CALL PGSWIN(XW1,XW2,YW1,YW2)
        CALL PGSCH(OLD_CH)
C terminamos buffering
        CALL PGEBUF
C
        DO I=1,N
          XPLAST(I)=XP(I)
          YPLAST(I)=YP(I)
        END DO
        NLAST=N
        N1LAST=N1
        N2LAST=N2
        CLABX=CX
        CLABY=CY
        CLABG=CG
        ICOLORLAST=ICOLOR
        NSYMBLAST=NSYMB
        CHLAST=CH
C
        END
