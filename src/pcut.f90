! CAXIS: indica si el corte se va a hacer en el eje X o en el eje Y
! MODECUT: single, region, all
! K1,K2: (output) region sumada para hacer el corte
! XP,YP: corte calculado y dibujado
        SUBROUTINE PCUT(CAXIS,MODECUT,K1,K2,XP,YP)
        IMPLICIT NONE
        CHARACTER*1 CAXIS
        INTEGER MODECUT
!
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
!
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
        REAL XPB(NXYMAX*NOVERSAMPMAX,NMAXBUFF)
        REAL YPB(NXYMAX*NOVERSAMPMAX,NMAXBUFF)
        CHARACTER*1 CH
        CHARACTER*50 GLABEL
        LOGICAL LOVERCUTS,LOVERPLOT_USED
        LOGICAL LPB(NMAXBUFF)
!
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKIMAGEN2/NCBUFF
        COMMON/BLKXYLIMPLOT/NX1,NX2,NY1,NY2
        COMMON/BLKP_ALLBUFF_L/LPB
        COMMON/BLKP_ALLBUFF_XY/XPB,YPB
        COMMON/BLKLOVERCUTS/LOVERCUTS
!------------------------------------------------------------------------------
        DO I=1,NXYMAX
          XP_ERR(I)=0.0
          YP_ERR(I)=0.0
        END DO
        DO I=1,NMAXBUFF
          LPB(I)=.FALSE.
        END DO
! Segun la opcion elegida, seleccionamos el corte
        IF(CAXIS.EQ.'X')THEN
          IF(MODECUT.EQ.1)THEN
            CALL PGSCI(5)
            CALL RPGBAND(5,0,0.,0.,XC,YC,CH)
            CALL PGSCI(1)
            IY1=NINT(YC)
            IF(IY1.LT.NY1) IY1=NY1
            IF(IY1.GT.NY2) IY1=NY2
            IY2=IY1
            WRITE(*,*)
            WRITE(*,100) 'Selecting Y='
            WRITE(*,*) IY1
          ELSEIF(MODECUT.EQ.2)THEN
            CALL PGSCI(5)
            CALL RPGBAND(5,0,0.,0.,XC,YC,CH)
            CALL PGSCI(1)
            IY1=NINT(YC)
            IF(IY1.LT.NY1) IY1=NY1
            IF(IY1.GT.NY2) IY1=NY2
            WRITE(*,*)
            WRITE(*,100) 'Selecting Y='
            WRITE(*,*) IY1
            CALL PGSCI(5)
            CALL RPGBAND(3,0,0.,REAL(IY1),XC,YC,CH)
            CALL PGSCI(1)
            IY2=NINT(YC)
            IF(IY2.LT.NY1) IY2=NY1
            IF(IY2.GT.NY2) IY2=NY2
            WRITE(*,100) 'Selecting Y='
            WRITE(*,*) IY2
          ELSE
            IY1=NY1
            IY2=NY2
            WRITE(*,*)
            WRITE(*,100) 'Selecting Y='
            WRITE(*,*) IY1
            WRITE(*,100) 'Selecting Y='
            WRITE(*,*) IY2
          END IF
          CALL XCUT(NCBUFF,IY1,IY2,YP)
          DO J=1,NAXIS(1,NCBUFF)
            XP(J)=REAL(J)
            XPB(J,NCBUFF)=XP(J)
            YPB(J,NCBUFF)=YP(J)
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
            WRITE(*,*)
            WRITE(*,100) 'Selecting X='
            WRITE(*,*) IX1
          ELSEIF(MODECUT.EQ.2)THEN
            CALL PGSCI(5)
            CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
            CALL PGSCI(1)
            IX1=NINT(XC)
            IF(IX1.LT.NX1) IX1=NX1
            IF(IX1.GT.NX2) IX1=NX2
            WRITE(*,*)
            WRITE(*,100) 'Selecting X='
            WRITE(*,*) IX1
            CALL PGSCI(5)
            CALL RPGBAND(4,0,REAL(IX1),0.,XC,YC,CH)
            CALL PGSCI(1)
            IX2=NINT(XC)
            IF(IX2.LT.NX1) IX2=NX1
            IF(IX2.GT.NX2) IX2=NX2
            WRITE(*,100) 'Selecting X='
            WRITE(*,*) IX2
          ELSE
            IX1=NX1
            IX2=NX2
            WRITE(*,*)
            WRITE(*,100) 'Selecting X='
            WRITE(*,*) IX1
            WRITE(*,100) 'Selecting X='
            WRITE(*,*) IX2
          END IF
          CALL YCUT(NCBUFF,IX1,IX2,YP)
          DO I=1,NAXIS(2,NCBUFF)
            XP(I)=REAL(I)
            XPB(I,NCBUFF)=XP(I)
            YPB(I,NCBUFF)=YP(I)
          END DO
        END IF
!
        LPB(NCBUFF)=.TRUE.
!------------------------------------------------------------------------------
        IF(CAXIS.EQ.'X')THEN
          WRITE(GLABEL,'(A1,I5,A1,I5,A1)') '[',IY1,',',IY2,']'
          CALL RMBLANK(GLABEL,GLABEL,L)
          CALL SUBPLOT(NAXIS(1,NCBUFF),NX1,NX2,XP,YP,XP_ERR,YP_ERR, &
           .TRUE.,.TRUE.,.FALSE.,.FALSE.,'x axis','signal',GLABEL(1:L),NCBUFF,201,1.0)
          K1=IY1
          K2=IY2
        ELSE
          WRITE(GLABEL,'(A1,I5,A1,I5,A1)') '[',IX1,',',IX2,']'
          CALL RMBLANK(GLABEL,GLABEL,L)
          CALL SUBPLOT(NAXIS(2,NCBUFF),NY1,NY2,XP,YP,XP_ERR,YP_ERR, &
           .TRUE.,.TRUE.,.FALSE.,.FALSE.,'y axis','signal',GLABEL(1:L),NCBUFF,201,1.0)
          K1=IX1
          K2=IX2
        END IF
!------------------------------------------------------------------------------
! si hay mas buffers con las mismas dimensiones, dibujamos cortes superpuestos
! (solo lo hacemos con los buffers de datos)
        IF(.NOT.LOVERCUTS) RETURN
        LOVERPLOT_USED=.FALSE.
        DO NBUFF=1,NMAXBUFF/2
          IF(NBUFF.NE.NCBUFF)THEN
            IF((NAXIS(1,NBUFF).EQ.NAXIS(1,NCBUFF)).AND.(NAXIS(2,NBUFF).EQ.NAXIS(2,NCBUFF)))THEN
              LOVERPLOT_USED=.TRUE.
              LPB(NBUFF)=.TRUE.
              IF(CAXIS.EQ.'X')THEN
                CALL XCUT(NBUFF,IY1,IY2,YP_OVER)
                CALL SUBPLOTBIS(NAXIS(1,NBUFF),NX1,NX2,XP,YP_OVER,XP_ERR,YP_ERR,.FALSE.,.FALSE.,NBUFF,201,1.0)
                DO J=1,NAXIS(1,NBUFF)
                  XPB(J,NBUFF)=XP(J)
                  YPB(J,NBUFF)=YP_OVER(J)
                END DO
              ELSE
                CALL YCUT(NBUFF,IX1,IX2,YP_OVER)
                CALL SUBPLOTBIS(NAXIS(2,NBUFF),NY1,NY2,XP,YP_OVER,XP_ERR,YP_ERR,.FALSE.,.FALSE.,NBUFF,201,1.0)
                DO I=1,NAXIS(2,NBUFF)
                  XPB(I,NBUFF)=XP(I)
                  YPB(I,NBUFF)=YP_OVER(I)
                END DO
              END IF
            END IF
          END IF
        END DO
! dibujamos encima el corte del buffer activo
        IF(LOVERPLOT_USED)THEN
          IF(CAXIS.EQ.'X')THEN
            CALL XCUT(NBUFF,IY1,IY2,YP_OVER)
            CALL SUBPLOTBIS(NAXIS(1,NCBUFF),NX1,NX2,XP,YP,XP_ERR,YP_ERR,.FALSE.,.FALSE.,NCBUFF,201,1.0)
          ELSE
            CALL SUBPLOTBIS(NAXIS(2,NCBUFF),NY1,NY2,XP,YP,XP_ERR,YP_ERR,.FALSE.,.FALSE.,NCBUFF,201,1.0)
          END IF
        END IF
!------------------------------------------------------------------------------
100     FORMAT(A,$)
        END
