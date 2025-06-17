! Dibuja un plot de XP(N) y YP(N), con N entre N1 y N2
! ERRXP,ERRYP: errores en XP e YP (solo se dibujan si LEX y/o LEY son .TRUE.)
! LX,LY indica si hay que calcular limites automaticamente en X y/o en Y
! CX,CY,CG: etiquetas en eje X, Y y en Grafica
! ICOLOR: color para lineas/puntos
! NSYMB: numero de simbolo 
!        101,102,...=linea de tipo 1,2,.. con PGLINE
!        201,202,...=linea de tipo 1,2,.. con PGBIN
! CH: current height para simbolos
!
        SUBROUTINE SUBPLOT(N,N1,N2,XP,YP,ERRXP,ERRYP,LX,LY,LEX,LEY,CX,CY,CG,ICOLOR,NSYMB,CH)
        IMPLICIT NONE
        INTEGER N
        INTEGER N1,N2
        REAL XP(N),YP(N)
        REAL ERRXP(N),ERRYP(N)
        LOGICAL LX,LY,LEX,LEY
        CHARACTER*(*) CX
        CHARACTER*(*) CY
        CHARACTER*(*) CG
        INTEGER ICOLOR
        INTEGER NSYMB
        REAL CH
!
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
!
        INTEGER NLAST,N1LAST,N2LAST
        INTEGER ICOLORLAST,NSYMBLAST
        INTEGER NB
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XPLAST(NXYMAX*NOVERSAMPMAX),YPLAST(NXYMAX*NOVERSAMPMAX)
        REAL CHLAST
        CHARACTER*2 CBUFF
        CHARACTER*255 CLABX,CLABY,CLABG
!
        INTEGER I
        REAL DX,DY
        REAL XV1,XV2,YV1,YV2
        REAL XW1,XW2,YW1,YW2
        REAL OLD_CH
!
        COMMON/BLKPLIMITS/XMIN,XMAX,YMIN,YMAX           !para guardar variables
        COMMON/BLKPNLAST/NLAST,N1LAST,N2LAST
        COMMON/BLKPLAST/XPLAST,YPLAST
        COMMON/BLKPLABELS/CLABX,CLABY,CLABG
        COMMON/BLKPCOLSYMB/ICOLORLAST,NSYMBLAST
        COMMON/BLKPCH/CHLAST
!------------------------------------------------------------------------------
! comenzamos buffering
        CALL PGBBUF
! almacenamos region de dibujo actual
        CALL PGQVP(0,XV1,XV2,YV1,YV2)
        CALL PGQWIN(XW1,XW2,YW1,YW2)
        CALL PGQCH(OLD_CH)
! borramos plot previo
        CALL RPGERASW(0.00,0.40,0.38,0.80,0)
! definimos nueva region de dibujo para el plot
        CALL PGSCH(0.8)
        CALL PGSVP(0.05,0.37,0.46,0.75)
! en caso necesario, calculamos limites
        IF(LX)THEN
          IF(LEX)THEN
            XMIN=XP(N1)-ERRXP(N1)
            XMAX=XP(N1)+ERRXP(N1)
            IF(N2.GT.N1)THEN
              DO I=N1+1,N2
                IF(XMIN.GT.XP(I)-ERRXP(I)) XMIN=XP(I)-ERRXP(I)
                IF(XMAX.LT.XP(I)+ERRXP(I)) XMAX=XP(I)+ERRXP(I)
              END DO
            END IF
          ELSE
            CALL FINDMML(N,N1,N2,XP,XMIN,XMAX)
          END IF
          DX=XMAX-XMIN
          IF(DX.LE.1.0E-6) DX=1.E-6
          XMIN=XMIN-DX/20.
          XMAX=XMAX+DX/20.
        END IF
        IF(LY)THEN
          IF(LEY)THEN
            YMIN=YP(N1)-ERRYP(N1)
            YMAX=YP(N1)+ERRYP(N1)
            IF(N2.GT.N1)THEN
              DO I=N1+1,N2
                IF(YMIN.GT.YP(I)-ERRYP(I)) YMIN=YP(I)-ERRYP(I)
                IF(YMAX.LT.YP(I)+ERRYP(I)) YMAX=YP(I)+ERRYP(I)
              END DO
            END IF
          ELSE
            CALL FINDMML(N,N1,N2,YP,YMIN,YMAX)
          END IF
          DY=YMAX-YMIN
          IF(DY.LE.1.0E-6) DY=1.0E-6
          YMIN=YMIN-DY/20.
          YMAX=YMAX+DY/20.
        END IF
! dibujamos caja
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
! dibujamos datos
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
            CALL PGERRX(1,XP(I)-ERRXP(I),XP(I)+ERRXP(I),YP(I),1.0)
          END DO
        END IF
        IF(LEY)THEN
          DO I=N1,N2
            CALL PGERRY(1,XP(I),YP(I)-ERRYP(I),YP(I)+ERRYP(I),1.0)
          END DO
        END IF
        CALL PGSCI(1)
! recuperamos region de dibujo inicial
        CALL PGSVP(XV1,XV2,YV1,YV2)
        CALL PGSWIN(XW1,XW2,YW1,YW2)
        CALL PGSCH(OLD_CH)
! terminamos buffering
        CALL PGEBUF
!
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
!
        END
