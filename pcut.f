C CAXIS: indica si el corte se va a hacer en el eje X o en el eje Y
C MODECUT: single, region, all
C K1,K2: (output) region sumada para hacer el corte
C XP,YP: corte calculado y dibujado
        SUBROUTINE PCUT(CAXIS,MODECUT,K1,K2,XP,YP)
        IMPLICIT NONE
        CHARACTER*1 CAXIS
        INTEGER MODECUT
C
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
C
        INTEGER I,J,L
        INTEGER K1,K2
        INTEGER NX1,NX2,NY1,NY2
        INTEGER IX1,IX2,IY1,IY2
        INTEGER NCBUFF,NBUFF
        INTEGER NAXIS(2,NMAXBUFF)
        REAL XC,YC
        REAL XP(NXYMAX*NOVERSAMPMAX),YP(NXYMAX*NOVERSAMPMAX)
        REAL YP_OVER(NXYMAX*NOVERSAMPMAX)
        REAL XP_ERR(NXYMAX*NOVERSAMPMAX),YP_ERR(NXYMAX*NOVERSAMPMAX)
        CHARACTER*1 CH
        CHARACTER*50 GLABEL
        LOGICAL LOVERPLOT
C
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKIMAGEN2/NCBUFF
        COMMON/BLKXYLIMPLOT/NX1,NX2,NY1,NY2
C------------------------------------------------------------------------------
        DO I=1,NXYMAX
          XP_ERR(I)=0.0
          YP_ERR(I)=0.0
        END DO
C Segun la opcion elegida, seleccionamos el corte
        IF(CAXIS.EQ.'X')THEN
          IF(MODECUT.EQ.1)THEN
            CALL PGSCI(5)
            CALL RPGBAND(5,0,0.,0.,XC,YC,CH)
            CALL PGSCI(1)
            IY1=NINT(YC)
            IF(IY1.LT.NY1) IY1=NY1
            IF(IY1.GT.NY2) IY1=NY2
            IY2=IY1
          ELSEIF(MODECUT.EQ.2)THEN
            CALL PGSCI(5)
            CALL RPGBAND(5,0,0.,0.,XC,YC,CH)
            CALL PGSCI(1)
            IY1=NINT(YC)
            IF(IY1.LT.NY1) IY1=NY1
            IF(IY1.GT.NY2) IY1=NY2
            CALL PGSCI(5)
            CALL RPGBAND(3,0,0.,REAL(IY1),XC,YC,CH)
            CALL PGSCI(1)
            IY2=NINT(YC)
            IF(IY2.LT.NY1) IY2=NY1
            IF(IY2.GT.NY2) IY2=NY2
          ELSE
            IY1=NY1
            IY2=NY2
          END IF
          CALL XCUT(NCBUFF,IY1,IY2,YP)
          DO J=1,NAXIS(1,NCBUFF)
            XP(J)=REAL(J)
          END DO
        ELSE
          IF(MODECUT.EQ.1)THEN
            CALL PGSCI(5)
            CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
            CALL PGSCI(1)
            IX1=NINT(XC)
            IF(IX1.LT.NX1) IX1=NX1
            IF(IX1.GT.NX2) IX1=NX2
            IX2=IX1
          ELSEIF(MODECUT.EQ.2)THEN
            CALL PGSCI(5)
            CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
            CALL PGSCI(1)
            IX1=NINT(XC)
            IF(IX1.LT.NX1) IX1=NX1
            IF(IX1.GT.NX2) IX1=NX2
            CALL PGSCI(5)
            CALL RPGBAND(4,0,REAL(IX1),0.,XC,YC,CH)
            CALL PGSCI(1)
            IX2=NINT(XC)
            IF(IX2.LT.NX1) IX2=NX1
            IF(IX2.GT.NX2) IX2=NX2
          ELSE
            IX1=NX1
            IX2=NX2
          END IF
          CALL YCUT(NCBUFF,IX1,IX2,YP)
          DO I=1,NAXIS(2,NCBUFF)
            XP(I)=REAL(I)
          END DO
        END IF
C------------------------------------------------------------------------------
        IF(CAXIS.EQ.'X')THEN
          WRITE(GLABEL,'(A1,I5,A1,I5,A1)') '[',IY1,',',IY2,']'
          CALL RMBLANK(GLABEL,GLABEL,L)
          CALL SUBPLOT(NAXIS(1,NCBUFF),NX1,NX2,XP,YP,XP_ERR,YP_ERR,
     +     .TRUE.,.TRUE.,.FALSE.,.FALSE.,
     +     'x axis','signal',GLABEL(1:L),NCBUFF+1,201,1.0)
          K1=IY1
          K2=IY2
        ELSE
          WRITE(GLABEL,'(A1,I5,A1,I5,A1)') '[',IX1,',',IX2,']'
          CALL RMBLANK(GLABEL,GLABEL,L)
          CALL SUBPLOT(NAXIS(2,NCBUFF),NY1,NY2,XP,YP,XP_ERR,YP_ERR,
     +     .TRUE.,.TRUE.,.FALSE.,.FALSE.,
     +     'y axis','signal',GLABEL(1:L),NCBUFF+1,201,1.0)
          K1=IX1
          K2=IX2
        END IF
C------------------------------------------------------------------------------
C si hay mas buffers con las mismas dimensiones, dibujamos cortes superpuestos
C (solo lo hacemos con los buffers de datos)
        LOVERPLOT=.FALSE.
        DO NBUFF=1,NMAXBUFF/2
          IF(NBUFF.NE.NCBUFF)THEN
            IF((NAXIS(1,NBUFF).EQ.NAXIS(1,NCBUFF)).AND.
     +         (NAXIS(2,NBUFF).EQ.NAXIS(2,NCBUFF)))THEN
              LOVERPLOT=.TRUE.
              IF(CAXIS.EQ.'X')THEN
                CALL XCUT(NBUFF,IY1,IY2,YP_OVER)
                CALL SUBPLOTBIS(NAXIS(1,NBUFF),NX1,NX2,XP,YP_OVER,
     +           XP_ERR,YP_ERR,
     +           .FALSE.,.FALSE.,NBUFF+1,201,1.0)
              ELSE
                CALL YCUT(NBUFF,IX1,IX2,YP_OVER)
                CALL SUBPLOTBIS(NAXIS(2,NBUFF),NY1,NY2,XP,YP_OVER,
     +           XP_ERR,YP_ERR,
     +           .FALSE.,.FALSE.,NBUFF+1,201,1.0)
              END IF
            END IF
          END IF
        END DO
C dibujamos encima el corte del buffer activo
        IF(LOVERPLOT)THEN
          IF(CAXIS.EQ.'X')THEN
            CALL XCUT(NBUFF,IY1,IY2,YP_OVER)
            CALL SUBPLOTBIS(NAXIS(1,NCBUFF),NX1,NX2,XP,YP,
     +       XP_ERR,YP_ERR,
     +       .FALSE.,.FALSE.,NCBUFF+1,201,1.0)
          ELSE
            CALL SUBPLOTBIS(NAXIS(2,NCBUFF),NY1,NY2,XP,YP,
     +       XP_ERR,YP_ERR,
     +       .FALSE.,.FALSE.,NCBUFF+1,201,1.0)
          END IF
        END IF
C------------------------------------------------------------------------------
        END
