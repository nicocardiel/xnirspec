!------------------------------------------------------------------------------
! Version 21-January-1998                                        File: button.f
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This routine is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!omment
!
! SUBROUTINE BUTTON(N,TEXT,MODE)
!
! Input: N,TEXT,MODE
! Input (COMMON): global variables in button.inc
!
! Plot buttons and button text in different modes.
!
! INTEGER       N -> button number in the range of available buttons (which
!               runs from 1 to MAX_XBUTT x MAX_YBUTT)
! CHARACTER*(*) TEXT -> the text that will appear in the button
! INTEGER       MODE -> determine the button mode:
!               MODE=-2,-3,...: only text is plotted with PGPLOT color=-NMODE-1
!               (i.e. 1,2,3...)
!               MODE=-1 erase the button
!               MODE=0 whole button is plotted (text in black)
!               MODE=1 only text is plotted (white)
!               MODE=2 only text is plotted (black)
!               MODE=3 only text is plotted (gray, button disabled)
!               MODE=4 whole button with reversed colors (text in black)
!               MODE=5 whole button with reversed colors (text in white)
!
!omment
!------------------------------------------------------------------------------
! Ancho del boton = DX/10; separacion entre botones = DX/40
        SUBROUTINE BUTTON(N,TEXT,MODE)
        IMPLICIT NONE
        INTEGER N
        CHARACTER*(*) TEXT
        INTEGER MODE,LMODE
!------------------------------------------------------------------------------
        INCLUDE 'button.inc'
        INTEGER TRUELEN
!------------------------------------------------------------------------------
! Variables locales
        INTEGER L
        INTEGER NLIN,NCOL
        INTEGER PGSCF_OLD
        INTEGER NCOLORTEXT
        REAL PGSCH_OLD
        REAL X1,X2,Y1,Y2
        REAL XP(4),YP(4)
        REAL DX,DY
        REAL DDX,DDY                               !anchura del canto del boton
        REAL XGAP,YGAP                                     !hueco entre botones
        REAL OFFY
        REAL X1W,X2W,Y1W,Y2W
        REAL X1V,X2V,Y1V,Y2V
        CHARACTER*80 CWHERE,CDUMMY
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
        LMODE=MODE
        NCOLORTEXT=0
        IF((N.GT.MAX_XBUTT*MAX_YBUTT).OR.(N.LT.1))THEN
          WRITE(*,110)'BUTTON NUMBER = ',N
          INCLUDE 'deallocate_arrays.inc'
          STOP 'ERROR: invalid button number in SUBROUTINE BUTTON'
        END IF
        IF(LMODE.GT.5)THEN
          WRITE(*,110)'TYPE = ',LMODE
          INCLUDE 'deallocate_arrays.inc'
          STOP 'ERROR: invalid button type in SUBROUTINE BUTTON'
        END IF
        IF(LMODE.LT.-1)THEN
          NCOLORTEXT=-LMODE-1
          LMODE=2
        END IF
!
        NLIN=(N-1)/MAX_XBUTT+1                              !posicion del boton
        NCOL=N-(NLIN-1)*MAX_XBUTT                           !posicion del boton
!
        IF(LMODE.EQ.-1)THEN
          EXIST_BUTT(N)=.FALSE.
        ELSEIF(LMODE.EQ.3)THEN
          EXIST_BUTT(N)=.FALSE.
        ELSE
          EXIST_BUTT(N)=.TRUE.
        END IF
!
        IF(MODOTEXT_BUTT)THEN
          IF(.NOT.MODOTEXT_PLOTBUTT) GOTO 90    !modo texto sin dibujar botones
        END IF
!------------------------------------------------------------------------------
        DX=X4VPORT-X3VPORT
        DY=Y4VPORT-Y3VPORT
!
        X1=X3VPORT+REAL(NCOL-1)*DX/REAL(MAX_XBUTT)
        X2=X1+DX/REAL(MAX_XBUTT)
        Y2=Y4VPORT-REAL(NLIN-1)*DY/REAL(MAX_YBUTT)
        Y1=Y2-DY/REAL(MAX_YBUTT)
!
        DDX=(X2-X1)/20.
        DDY=(Y2-Y1)/7.
        XGAP=DDX/10.
        YGAP=DDY/10.
        OFFY=(Y2-Y1)*YTEXT_BUTT
!
        X1=X1+XGAP
        X2=X2-XGAP
        Y1=Y1+YGAP
        Y2=Y2-YGAP
!
        CALL PGQWIN(X1W,X2W,Y1W,Y2W)
        CALL PGQVP(0,X1V,X2V,Y1V,Y2V)
        CALL PGVPORT(X3VPORT,X4VPORT,Y3VPORT,Y4VPORT)
        CALL PGWINDOW(X3VPORT,X4VPORT,Y3VPORT,Y4VPORT)
!
        IF(LMODE.EQ.-1)THEN
          CALL PGSCI(0)
          CALL PGRECT(X1-XGAP,X2+XGAP,Y1-YGAP,Y2+YGAP)
          CALL PGSCI(1)
          GOTO 80
        END IF
!
        CALL PGQCF(PGSCF_OLD)
        CALL PGSCF(PGSCF_BUTT)
        CALL PGQCH(PGSCH_OLD)
        CALL PGSCH(PGSCH_BUTT)
!------------------------------------------------------------------------------
        IF((LMODE.GE.1).AND.(LMODE.LE.3))THEN
          IF(LMODE.EQ.1)THEN
            CALL PGSCI(1)
          ELSEIF(LMODE.EQ.2)THEN
            CALL PGSCI(NCOLORTEXT)
          ELSEIF(LMODE.EQ.3)THEN
            EXIST_BUTT(N)=.FALSE.
            CALL PGSCI(15)
          END IF
          CALL PGPTEXT((X1+X2)/2.,Y1+OFFY,0.,0.5,TEXT)
          IF(MODOOVER_BUTT)THEN
            WRITE(CDUMMY,*) N
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            CALL PGSCI(1)
            CALL PGPTEXT((X1+X2)/2.,Y1+OFFY,0.,0.5,CDUMMY(1:L))
          END IF
        ELSE
          CALL PGSCI(12)
          CALL PGRECT(X1+DDX,X2-DDX,Y1+DDY,Y2-DDY)
          IF(LMODE.EQ.0)THEN
            CALL PGSCI(13)
          ELSE
            CALL PGSCI(14)
          END IF
          CALL SUBDATA(XP,X1,X1+DDX,X1+DDX,X1)
          CALL SUBDATA(YP,Y1,Y1+DDY,Y2-DDY,Y2)
          CALL PGPOLY(4,XP,YP)
          CALL SUBDATA(XP,X1,X1+DDX,X2-DDX,X2)
          CALL SUBDATA(YP,Y2,Y2-DDY,Y2-DDY,Y2)
          CALL PGPOLY(4,XP,YP)
          IF(LMODE.EQ.0)THEN
            CALL PGSCI(14)
          ELSE
            CALL PGSCI(13)
          END IF
          CALL SUBDATA(XP,X1,X1+DDX,X2-DDX,X2)
          CALL SUBDATA(YP,Y1,Y1+DDY,Y1+DDY,Y1)
          CALL PGPOLY(4,XP,YP)
          CALL SUBDATA(XP,X2,X2-DDX,X2-DDX,X2)
          CALL SUBDATA(YP,Y1,Y1+DDY,Y2-DDY,Y2)
          CALL PGPOLY(4,XP,YP)
          IF(LMODE.EQ.0)THEN
            CALL PGSCI(1)
          ELSE
            CALL PGSCI(0)
          END IF
          CALL PGMOVE(X1,Y2)
          CALL PGDRAW(X1+DDX,Y2-DDY)
          IF(LMODE.EQ.0)THEN
            CALL PGSCI(0)
          ELSE
            CALL PGSCI(1)
          END IF
          CALL PGMOVE(X2,Y1)
          CALL PGDRAW(X2-DDX,Y1+DDY)
          IF(LMODE.EQ.5)THEN
            CALL PGSCI(1)
          ELSE
            CALL PGSCI(0)
          END IF
          CALL PGPTEXT((X1+X2)/2.,Y1+OFFY,0.,0.5,TEXT)
          IF(MODOOVER_BUTT)THEN
            WRITE(CDUMMY,*) N
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            CALL PGSCI(1)
            CALL PGPTEXT((X1+X2)/2.,Y1+OFFY,0.,0.5,CDUMMY(1:L))
          END IF
          CALL PGSCI(0)
          CALL PGMOVE(X1,Y1)
          CALL PGDRAW(X2,Y1)
          CALL PGDRAW(X2,Y2)
          CALL PGDRAW(X1,Y2)
          CALL PGDRAW(X1,Y1)
          CALL PGSCI(1)
        END IF
!------------------------------------------------------------------------------
        CALL PGSCF(PGSCF_OLD)
        CALL PGSCH(PGSCH_OLD)
        IF(LMODE.EQ.3) CALL PGSCI(1)
        IF(LMODE.EQ.2) CALL PGSCI(1)
!
80      CALL PGVPORT(X1V,X2V,Y1V,Y2V)
        CALL PGWINDOW(X1W,X2W,Y1W,Y2W)
        RETURN
!------------------------------------------------------------------------------
! Botones en modo texto
90      WRITE(CWHERE,'(A,I2,A1,I3,A1)')'[',NLIN,';',(NCOL-1)*20,'f'
        CALL RMBLANK(CWHERE,CWHERE,L)
        WRITE(*,'(A)')CWHERE(1:L)//'                    '
        IF(LMODE.EQ.-1)THEN
        ELSEIF(LMODE.EQ.3)THEN
          WRITE(*,'(A,I3,A,$)')CWHERE(1:L),N,'-'
          WRITE(*,101)TEXT(1:TRUELEN(TEXT))
        ELSEIF(LMODE.EQ.5)THEN
          WRITE(*,100)'[5m'
          WRITE(*,'(A,I3,A,$)')CWHERE(1:L),N,'-'
          WRITE(*,100)TEXT(1:TRUELEN(TEXT))
          WRITE(*,101)'[0m'
        ELSE
          WRITE(*,100)'[5m'
          WRITE(*,'(A,I3,A,$)')'[5m'//CWHERE(1:L),N,'[0m'
          WRITE(*,100)'-'
          WRITE(*,101)TEXT(1:TRUELEN(TEXT))
        END IF
        WRITE(CWHERE,'(A,I2,A3)')'[',MAX_YBUTT+1,';1f'
        CALL RMBLANK(CWHERE,CWHERE,L)
        WRITE(*,'(A)')CWHERE(1:L)//'[J'
        WRITE(*,'(A)')CWHERE(1:L)
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I4)
        END
!
!******************************************************************************
!
        SUBROUTINE SUBDATA(XX,X1,X2,X3,X4)
        IMPLICIT NONE
!
! Introduce en la matriz XX(4) los valores X1,X2,X3 y X4 (ahorra codigo)
!------------------------------------------------------------------------------
! Parametros de la subrutina
        REAL XX(4)
        REAL X1,X2,X3,X4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
        XX(1)=X1
        XX(2)=X2
        XX(3)=X3
        XX(4)=X4
        END
