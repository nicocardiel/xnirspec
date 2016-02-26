C Overplot of XP(N) y YP(N), with N in the range [N1,N2]
C
C Esta subrutina no modifica el ultimo conjunto de datos dibujado, por lo que
C lo representado con esta rutina no vuelve a mostrarse con un [zoom] o [whole]
C
C Dibuja un plot de XP(N) y YP(N), con N entre N1 y N2
C ICOLOR: color para lineas/puntos
C NSYMB: numero de simbolo
C        101,102,...=linea de tipo 1,2,.. con PGLINE
C        201,202,...=linea de tipo 1,2,.. con PGBIN
C CH: current height para simbolos
C
        SUBROUTINE SUBPLOTBIS(N,N1,N2,XP,YP,EXP,EYP,LEX,LEY,
     +   ICOLOR,NSYMB,CH)
        IMPLICIT NONE
        INTEGER N
        INTEGER N1,N2
        REAL XP(N),YP(N)
        REAL EXP(N),EYP(N)
        LOGICAL LEX,LEY
        INTEGER ICOLOR
        INTEGER NSYMB
        REAL CH
C
        INTEGER I
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XV1,XV2,YV1,YV2
        REAL XW1,XW2,YW1,YW2
        REAL OLD_CH
C
        COMMON/BLKPLIMITS/XMIN,XMAX,YMIN,YMAX               !limites anteriores
C------------------------------------------------------------------------------
C comenzamos buffering
        CALL PGBBUF
C almacenamos region de dibujo actual
        CALL PGQVP(0,XV1,XV2,YV1,YV2)
        CALL PGQWIN(XW1,XW2,YW1,YW2)
        CALL PGQCH(OLD_CH)
C definimos nueva region de dibujo para el plot
        CALL PGSVP(0.05,0.37,0.46,0.75)
C dibujamos datos
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
        END
