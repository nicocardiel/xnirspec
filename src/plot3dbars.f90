!
!******************************************************************************
! Realiza un dibujo tridimensional de la region alrededor del pixel XC,YC
! IMODE=0: lineas ocultas
! IMODE=1: lineas transparentes
! MM3D=.FALSE. limites automaticos (se devuelven como output en BG3D,FG3D)
! MM3D=.TRUE. limites fijos (se toman BG3D y FG3D como input)
        SUBROUTINE PLOT3DBARS(NC1,NC2,NS1,NS2,IMODE,MM3D,BG3D,FG3D)
        USE Dynamic_Array_IMAGEN_
        IMPLICIT NONE
        INCLUDE 'interface_imagen_.inc'
! subroutine arguments
        INTEGER NC1,NC2,NS1,NS2
        INTEGER IMODE
        LOGICAL MM3D                           !si .FALSE., limites automaticos
        REAL BG3D,FG3D
! functions
!
        INTEGER I,J
        INTEGER K,L
        INTEGER NCOLOR
        INTEGER ICILO,ICIHI   !lowest and highest color index to use for images
        INTEGER OLD_FS
!delete REAL IMAGEN_(NXMAX,NYMAX)
        REAL XMIN,XMAX,YMIN,YMAX,YMIN0,YMAX0,DX
        REAL XOFFSET,YOFFSET
        REAL X(5),Y(5),YB(4)
        REAL XX(7),YY(7)
        REAL F
        REAL XV1,XV2,YV1,YV2
        CHARACTER*50 CDUMMY
!
!delete COMMON/BLKIMAGEN1_/IMAGEN_              !es global para ahorrar memoria
!------------------------------------------------------------------------------
        CALL PGQCIR(ICILO, ICIHI)
!------------------------------------------------------------------------------
        YMIN0=IMAGEN_(NC1,NS1)
        YMAX0=YMIN0
        DO I=NS1,NS2
          DO J=NC1,NC2
            IF(IMAGEN_(J,I).LT.YMIN0) YMIN0=IMAGEN_(J,I)
            IF(IMAGEN_(J,I).GT.YMAX0) YMAX0=IMAGEN_(J,I)
          END DO
        END DO
        IF(.NOT.MM3D)THEN
          BG3D=YMIN0
          FG3D=YMAX0
          IF(BG3D.EQ.FG3D)THEN
            IF(BG3D.EQ.0.)THEN
              BG3D=-1.
              FG3D=1.
            ELSE
              BG3D=0.9*BG3D
              FG3D=1.1*FG3D
            END IF
          END IF
        END IF
!
        XOFFSET=REAL(NC2-NC1+1)/(4*REAL(NS2-NS1+1))
        YOFFSET=0.95*(FG3D-BG3D)/REAL(NS2-NS1+1)
        XMIN=REAL(NC1)-1.
        XMAX=REAL(NC2)+REAL(NS2-NS1)*XOFFSET+1.
        YMAX=FG3D+REAL(NS2-NS1-1)*YOFFSET
        YMAX=YMAX+YOFFSET
        YMIN=BG3D-YOFFSET
        DX=XMAX-XMIN
!------------------------------------------------------------------------------
! obtenemos la region de dibujo actual
        CALL PGQVP(0,XV1,XV2,YV1,YV2)
        CALL PGQFS(OLD_FS)
        CALL PGSFS(1)
! definimos region de dibujo para imagen 3D
        CALL PGSVP(XV1,XV2,YV1+(YV2-YV1)*0.2,YV2)
        CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
!------------------------------------------------------------------------------
        DO I=NS2,NS1,-1
          DO J=NC1,NC2
            X(1)=REAL(J)-.5+(REAL(I-NS1)-.5)*XOFFSET
            X(2)=X(1)+1.
            X(3)=REAL(J)+.5+(REAL(I-NS1)+.5)*XOFFSET
            X(4)=X(3)-1.
            Y(1)=IMAGEN_(J,I)+(REAL(I-NS1)-.5)*YOFFSET
            Y(2)=Y(1)
            Y(3)=IMAGEN_(J,I)+(REAL(I-NS1)+.5)*YOFFSET
            Y(4)=Y(3)
            X(5)=X(1)
            Y(5)=Y(1)
            YB(1)=BG3D+(REAL(I-NS1)-.5)*YOFFSET
            YB(2)=YB(1)
            YB(3)=BG3D+(REAL(I-NS1)+.5)*YOFFSET
            YB(4)=YB(3)
            IF(IMAGEN_(J,I).LT.BG3D)THEN
              NCOLOR=ICIHI
            ELSEIF(IMAGEN_(J,I).GT.FG3D)THEN
              NCOLOR=ICILO
            ELSE
              F=(IMAGEN_(J,I)-BG3D)/(FG3D-BG3D)
              NCOLOR=ICIHI-NINT(F*REAL(ICIHI-ICILO))
            END IF
            CALL PGSCI(NCOLOR)
!..............................................................................
              IF(IMODE.EQ.0)THEN
              XX(1)=X(1)
              XX(2)=X(1)
              XX(3)=X(4)
              XX(4)=X(3)
              XX(5)=X(3)
              XX(6)=X(2)
              XX(7)=XX(1)
              YY(1)=YB(1)
              YY(2)=Y(1)
              YY(3)=Y(4)
              YY(4)=Y(3)
              YY(5)=YB(3)
              YY(6)=YB(2)
              YY(7)=YY(1)
              CALL PGPOLY(7,XX,YY)
!!!              CALL PGSCI(0)
!!!              CALL PGSCI((ICIHI+ICILO)/2)
              CALL PGSCI(8)
              CALL PGLINE(5,X,Y)
              DO K=1,3
                CALL PGMOVE(X(K),Y(K))
                CALL PGDRAW(X(K),YB(K))
              END DO
!..............................................................................
            ELSEIF(IMODE.EQ.1)THEN
              CALL PGPOLY(5,X,Y)
              DO K=1,4
                CALL PGMOVE(X(K),Y(K))
                CALL PGDRAW(X(K),YB(K))
              END DO
!..............................................................................
            ELSE
            END IF
          END DO
        END DO
!------------------------------------------------------------------------------
        CALL PGSVP(XV1,XV2,YV1+(YV2-YV1)*0.10,YV1+(YV2-YV1)*0.15)
        CALL PGWINDOW(REAL(ICILO),REAL(ICIHI+1),0.,1.)
        DO K=ICILO,ICIHI
          CALL PGSCI(ICIHI-K+ICILO)
          CALL PGRECT(REAL(K),REAL(K+1),0.,1.)
        END DO
        CALL PGSCI(8)
        CALL PGBOX('BC',0.,0,'BC',0.,0)
        CALL PGSCI(1)
        WRITE(CDUMMY,*) BG3D
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
          DO WHILE(CDUMMY(L:L).EQ.'0')
            L=L-1
          END DO
        END IF
        CALL PGMTEXT('B',1.2,0.,0.,CDUMMY(1:L))
        WRITE(CDUMMY,*) FG3D
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
          DO WHILE(CDUMMY(L:L).EQ.'0')
            L=L-1
          END DO
        END IF
        CALL PGMTEXT('B',1.2,1.,1.,CDUMMY(1:L))
        CALL PGSCI(1)
!------------------------------------------------------------------------------
        CALL PGSFS(OLD_FS)
!------------------------------------------------------------------------------
        END
