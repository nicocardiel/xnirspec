!
!******************************************************************************
! Hace un zoom en el histograma
        SUBROUTINE ZOOMHISTOG(BG,FG)
        IMPLICIT NONE
        REAL BG,FG
!
        REAL XV1,XV2,YV1,YV2
        REAL XW1,XW2,YW1,YW2
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XC,YC
        REAL XC1,XC2
        CHARACTER*1 CH
!
        COMMON/BLKHLIMITS/XMIN,XMAX,YMIN,YMAX
!------------------------------------------------------------------------------
! almacenamos region de dibujo actual
        CALL PGQVP(0,XV1,XV2,YV1,YV2)
        CALL PGQWIN(XW1,XW2,YW1,YW2)
!------------------------------------------------------------------------------
! definimos la region de dibujo para el histograma
        CALL PGSVP(0.05,0.37,0.08,0.30)
        CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
!
        CALL PGSCI(5)
        CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
        CALL PGSCI(1)
        XC1=XC
        IF(XMIN.LE.XMAX)THEN
          IF(XC1.LT.XMIN) XC1=XMIN
          IF(XC1.GT.XMAX) XC1=XMAX
        ELSE
          IF(XC1.GT.XMIN) XC1=XMIN
          IF(XC1.LT.XMAX) XC1=XMAX
        END IF
        CALL PGSCI(5)
        CALL RPGBAND(4,0,XC1,0.,XC,YC,CH)
        CALL PGSCI(1)
        XC2=XC
        IF(XMIN.LE.XMAX)THEN
          IF(XC2.LT.XMIN) XC2=XMIN
          IF(XC2.GT.XMAX) XC2=XMAX
        ELSE
          IF(XC2.GT.XMIN) XC2=XMIN
          IF(XC2.LT.XMAX) XC2=XMAX
        END IF
        IF(XC1.NE.XC2)THEN
          BG=AMIN1(XC1,XC2)
          FG=AMAX1(XC1,XC2)
        ELSE
          WRITE(*,101) '***ERROR***'
          WRITE(*,101) '=> BG and FG must be different.'
        END IF
!------------------------------------------------------------------------------
! recuperamos region de dibujo inicial
        CALL PGSVP(XV1,XV2,YV1,YV2)
        CALL PGSWIN(XW1,XW2,YW1,YW2)
!
101     FORMAT(A)
        END
