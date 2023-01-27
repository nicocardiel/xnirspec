! Overplot of XP(N) y YP(N), with N in the range [N1,N2]
!
! Esta subrutina no modifica el ultimo conjunto de datos dibujado, por lo que
! lo representado con esta rutina no vuelve a mostrarse con un [zoom] o [whole]
!
! Dibuja un plot de XP(N) y YP(N), con N entre N1 y N2
! ICOLOR: color para lineas/puntos
! NSYMB: numero de simbolo
!        101,102,...=linea de tipo 1,2,.. con PGLINE
!        201,202,...=linea de tipo 1,2,.. con PGBIN
! CH: current height para simbolos
!
        SUBROUTINE SUBPLOTBIS(N,N1,N2,XP,YP,ERRXP,ERRYP,LEX,LEY,ICOLOR,NSYMB,CH)
        IMPLICIT NONE
        INTEGER N
        INTEGER N1,N2
        REAL XP(N),YP(N)
        REAL ERRXP(N),ERRYP(N)
        LOGICAL LEX,LEY
        INTEGER ICOLOR
        INTEGER NSYMB
        REAL CH
!
        INTEGER I
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XV1,XV2,YV1,YV2
        REAL XW1,XW2,YW1,YW2
        REAL OLD_CH
!
        COMMON/BLKPLIMITS/XMIN,XMAX,YMIN,YMAX               !limites anteriores
!------------------------------------------------------------------------------
! comenzamos buffering
        CALL PGBBUF
! almacenamos region de dibujo actual
        CALL PGQVP(0,XV1,XV2,YV1,YV2)
        CALL PGQWIN(XW1,XW2,YW1,YW2)
        CALL PGQCH(OLD_CH)
! definimos nueva region de dibujo para el plot
        CALL PGSVP(0.05,0.37,0.46,0.75)
! dibujamos datos
        CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
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
        END
