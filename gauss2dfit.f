C------------------------------------------------------------------------------
C Version 23-July-1998                                         File: gaussfit.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This routine is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C SUBROUTINE GAUSS2DFIT(NBUFF,IX1,IX2,IY1,IY2,LGUESS,
C                       X0,Y0,SIGMAX,SIGMAY,BETA,AMP,CTE,
C                       EX0,EY0,ESIGMAX,ESIGMAY,EBETA,EAMP,ECTE,
C                       EEX0,EEY0,EESIGMAX,EESIGMAY,EEBETA,EEAMP,EECTE,
C                       X0_,Y0_,SIGMAX_,SIGMAY_,AMPX_,AMPY_,CTEX_,CTEY_,
C                       AREA1,AREA2,AREA3,EAREA1,EAREA2,EAREA3,
C                       YRMSTOL,NSIMUL,IMODE,LSHOW)
C
C Input: NBUFF,IX1,IX2,IY1,IY2,YRMSTOL,NSIMUL,IMODE,LGUESS,LSHOW
C Input (COMMON): IMAGEN()
C Input/Output: X0,Y0,SIGMAX,SIGMAY,AMP,CTE,
C Output: EX0,EY0,ESIGMAX,ESIGMAY,EAMP,ECTE
C Output: EEX0,EEY0,EESIGMAX,EESIGMAY,EEAMP,EECTE
C Output: X0_,Y0_,SIGMAX_,SIGMAY_,AMPX_,AMPY_,CTEX_,CTEY_
C Output: AREA1,AREA2,AREA3,EAREA1,EAREA2,EAREA3
C
C Numerical fit of a 2D elliptical gaussian (using DOWNHILL):
C Z(X,Y)=CTE+AMP*
C  EXP[
C       -(X-X0)^2/(2*SIGMAX^2)
C       -(Y-Y0)^2/(2*SIGMAY^2)
C       -BETA*(X-X0)*(Y-Y0)/(SIGMAX*SIGMAY)
C     ]
C
C INTEGER NBUFF -> buffer to be fitted
C INTEGER IX1,IX2,IY1,IY2 -> limits of the box to be fitted (pixels)
C LOGICAL LGUESS -> if .TRUE., X0, Y0, SIGMAX, SIGMAY, BETA, AMP and CTE are
C                   an initial guess of the final fit
C REAL X0,Y0 -> center of the fitted gaussian
C REAL SIGMAX,SIGMAY -> rms lengths of the major and minor axes
C REAL BETA -> measure of the position-angle difference between the principal
C              axes of the ellipse and the coordinate axes (x,y)
C REAL AMP -> peak amplitude
C REAL CTE -> additive constant
C REAL EX0,EY0,ESIGMAX,ESIGMAY,EBETA,EAMP,ECTE -> errors due to the propagation
C                                                 of pixel errors (simulations)
C REAL EEX0,EEY0,EESIGMAX,EESIGMAY,EEBETA,EEAMP,EECTE -> errors due to the use
C                                     of a numerical method (rms from DOWNHILL)
C REAL YRMSTOL -> stopping criterion for DOWNHILL
C INTEGER NSIMUL -> number of simulations to compute errors
C INTEGER IMODE -> IMODE=0 minimize distance; IMODE=1 minimize area
C REAL X0_,Y0_,SIGMAX_,SIGMAY_,AMPX_,AMPY_,CTEX_,CTEY_ -> parameters of 1D fits
C AREA1,EAREA1 -> integrated signal in the data (& error through simulations)
C                 (note that the cte factor is subtracted)
C AREA2,EAREA2 -> integrated signal in the fit (& error through simulations)
C                 (note that the cte factor is subtracted)
C AREA3,EAREA3 -> integrated signal in the fit (& error through simulations)
C                 (numerical integration)
C LOGICAL LSHOW -> if .TRUE. show intermediate results
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE GAUSS2DFIT(NBUFF,IX1,IX2,IY1,IY2,LGUESS,
     +   X0,Y0,SIGMAX,SIGMAY,BETA,AMP,CTE,
     +   EX0,EY0,ESIGMAX,ESIGMAY,EBETA,EAMP,ECTE,
     +   EEX0,EEY0,EESIGMAX,EESIGMAY,EEBETA,EEAMP,EECTE,
     +   X0_,Y0_,SIGMAX_,SIGMAY_,AMPX_,AMPY_,CTEX_,CTEY_,
     +   AREA1,AREA2,AREA3,EAREA1,EAREA2,EAREA3,
     +   YRMSTOL,NSIMUL,IMODE,LSHOW) 
C
        IMPLICIT NONE
C
        INTEGER NBUFF
        INTEGER IX1,IX2,IY1,IY2
        LOGICAL LGUESS
        REAL X0,Y0,SIGMAX,SIGMAY,BETA,AMP,CTE
        REAL EX0,EY0,ESIGMAX,ESIGMAY,EBETA,EAMP,ECTE
        REAL EEX0,EEY0,EESIGMAX,EESIGMAY,EEBETA,EEAMP,EECTE
        REAL X0_,Y0_,SIGMAX_,SIGMAY_,AMPX_,AMPY_,CTEX_,CTEY_
        REAL AREA1,AREA2,AREA3,EAREA1,EAREA2,EAREA3
        REAL YRMSTOL
        INTEGER NSIMUL
        INTEGER IMODE
        LOGICAL LSHOW
C
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
C
        REAL FMEAN0
C
        INTEGER NSIMULMAX
        PARAMETER(NSIMULMAX=1000)                !numero maximo de simulaciones
        REAL PI2
        PARAMETER(PI2=6.283185307)               !2 x pi
C
        INTEGER IX1_,IX2_,IY1_,IY2_
        INTEGER I,J,ISIMUL
        INTEGER NP
        INTEGER NSEED
        INTEGER NPFITX,NPFITY
        INTEGER NEXTINFO
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL IMAGEN_(NXMAX,NYMAX)
        REAL EEX0_,EEY0_
        REAL EESIGMAX_,EESIGMAY_
        REAL EEAMPX_,EEAMPY_
        REAL EECTEX_,EECTEY_
        REAL XF(NXYMAX),YF(NXYMAX)
        REAL X0_SIMUL(NSIMULMAX),Y0_SIMUL(NSIMULMAX)
        REAL SIGMAX_SIMUL(NSIMULMAX),SIGMAY_SIMUL(NSIMULMAX)
        REAL AMP_SIMUL(NSIMULMAX),CTE_SIMUL(NSIMULMAX)
        REAL BETA_SIMUL(NSIMULMAX)
        REAL AREA1_SIMUL(NSIMULMAX)
        REAL AREA2_SIMUL(NSIMULMAX)
C esta linea esta comentada porque en las simulaciones hay problemas al
C calcular la integral numericamente
ccc        REAL AREA3_SIMUL(NSIMULMAX)
        REAL RANDOMNUMBER,R1,R2,ERRYF
        REAL FMEAN
        REAL CC(7)                 !parametros de la gaussiana para integral 2D
        LOGICAL LSIMUL
C
        INTEGER NDIM,NEVAL
        REAL XX0(7),DXX0(7)      !valores iniciales y desplazamientos de prueba
        REAL XX(7),DXX(7)
        EXTERNAL FUNKGSS2D
        REAL FUNKGSS2D
        EXTERNAL FUNKGSS2DA
        REAL FUNKGSS2DA
C
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKIMAGEN1_/IMAGEN_
        COMMON/BLKFIT2DG1/IX1_,IX2_,IY1_,IY2_
        COMMON/BLKFUNK2DY/CC
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        IF(.NOT.LGUESS)THEN
          CTE=0.
          AMP=0.
          X0=0.
          Y0=0.
          SIGMAX=0.
          SIGMAY=0.
          BETA=0.
        END IF
        CTEX_=0.
        CTEY_=0.
        AMPX_=0.
        AMPY_=0.
        X0_=0.
        Y0_=0.
        SIGMAX_=0.
        SIGMAY_=0.
        ECTE=0.
        EECTE=0.
        EAMP=0.
        EEAMP=0.
        EX0=0.
        EEX0=0.
        EY0=0.
        EEY0=0.
        ESIGMAX=0.
        EESIGMAX=0.
        ESIGMAY=0.
        EESIGMAY=0.
        EBETA=0.
        EEBETA=0.
        AREA1=0.
        EAREA1=0.
        AREA2=0.
        EAREA2=0.
        AREA3=0.
        EAREA3=0.
C
        NPFITX=IX2-IX1+1
        NPFITY=IY2-IY1+1
        NP=NPFITX*NPFITY
        IF(NP.LT.8)THEN
          WRITE(*,100) 'NP:'
          WRITE(*,*) NP
          WRITE(*,101) 'ERROR: in subroutine GAUSS2DFIT:'
          WRITE(*,101) ' No. of points for fit < 8'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
C
        IF(NSIMUL.GT.NSIMULMAX)THEN
          WRITE(*,100) 'NSIMULMAX:'
          WRITE(*,*) NSIMUL
          WRITE(*,101) 'ERROR: in subroutine GAUSS2DFIT:'
          WRITE(*,101) ' No. of simulations > NSIMULMAX'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
C
        IX1_=IX1
        IX2_=IX2
        IY1_=IY1
        IY2_=IY2
C------------------------------------------------------------------------------
C Valores de partida
        IF(LGUESS)THEN
          XX0(1)=CTE
          XX0(2)=AMP
          XX0(3)=X0
          XX0(4)=Y0
          XX0(5)=SIGMAX
          XX0(6)=SIGMAY
          XX0(7)=BETA
        ELSE
C si no tenemos una buena estimacion, obtenemos una haciendo ajustes a los 
C cortes promedio de la imagen
C 1) corte en X
          DO J=IX1,IX2
            XF(J-IX1+1)=REAL(J-IX1+1)
            YF(J-IX1+1)=0.
            IF(NBUFF.EQ.0)THEN
              DO I=IY1,IY2
                YF(J-IX1+1)=YF(J-IX1+1)+IMAGEN_(J,I)
              END DO
            ELSE
              DO I=IY1,IY2
                YF(J-IX1+1)=YF(J-IX1+1)+IMAGEN(J,I,NBUFF)
              END DO
            END IF
            YF(J-IX1+1)=YF(J-IX1+1)/REAL(NPFITX)
          END DO
          CALL GAUSCFIT(NPFITX,XF,YF,X0_,SIGMAX_,AMPX_,CTEX_,
     +     EEX0_,EESIGMAX_,EEAMPX_,EECTEX_,YRMSTOL)
          X0_=X0_+REAL(IX1-1)
C 2) corte en Y
          DO I=IY1,IY2
            XF(I-IY1+1)=REAL(I-IY1+1)
            YF(I-IY1+1)=0.
            IF(NBUFF.EQ.0)THEN
              DO J=IX1,IX2
                YF(I-IY1+1)=YF(I-IY1+1)+IMAGEN_(J,I)
              END DO
            ELSE
              DO J=IX1,IX2
                YF(I-IY1+1)=YF(I-IY1+1)+IMAGEN(J,I,NBUFF)
              END DO
            END IF
            YF(I-IY1+1)=YF(I-IY1+1)/REAL(NPFITY)
          END DO
          CALL GAUSCFIT(NPFITY,XF,YF,Y0_,SIGMAY_,AMPY_,CTEY_,
     +     EEY0_,EESIGMAY_,EEAMPY_,EECTEY_,YRMSTOL)
          Y0_=Y0_+REAL(IY1-1)
C Con estas estimaciones ya podemos calcular los parametros de entrada de
C la subrutina DOWNHILL (nota: recordar que los ajustes 1D corresponden a los
C perfiles promedio en X y en Y, y no al perfil que pasa por el centro; por
C eso hacemos una estimacion de AMP usando AMPX y AMPY multiplicadas por un
C factor=SQRT(2*PI)*SIGMA)
          XX0(1)=(CTEX_+CTEY_)/2.
          XX0(2)=(AMPX_*SIGMAY_*2.5+AMPY_*SIGMAX_*2.5)/2.
          XX0(3)=X0_
          XX0(4)=Y0_
          XX0(5)=SIGMAX_
          XX0(6)=SIGMAY_
          XX0(7)=0.0
        END IF
C------------------------------------------------------------------------------
        NDIM=7
        DXX0(1)=XX0(1)/100.
        IF(DXX0(1).EQ.0.0) DXX0(1)=0.1
        DXX0(2)=XX0(2)/100.
        IF(DXX0(2).EQ.0.0) DXX0(2)=0.1
        DXX0(3)=0.1
        DXX0(4)=0.1
        DXX0(5)=0.1
        DXX0(6)=0.1
        DXX0(7)=0.01
C
        IF(NBUFF.GT.0)THEN !si es cero, ya tenemos los datos a ajustar
          DO I=IY1,IY2
            DO J=IX1,IX2
              IMAGEN_(J,I)=IMAGEN(J,I,NBUFF)
            END DO
          END DO
        ELSEIF(NBUFF.LT.0)THEN
          STOP 'FATAL ERROR: NBUFF.LT.0!'
        END IF
        IF(LSHOW)THEN
          print*,'>Before downhill:'
          do i=1,7
            print*,i,xx0(i),dxx0(i)
          end do
        END IF
        CALL DOWNHILL(NDIM,XX0,DXX0,FUNKGSS2D,1.0,0.5,2.0,YRMSTOL,
     +   XX,DXX,NEVAL,500)
        IF(LSHOW)THEN
          print*,'>After downhill:'
          do i=1,7
            print*,i,xx(i),dxx(i)
          end do
        END IF
        IF(IMODE.EQ.1)THEN !como valores iniciales tomo el ajuste anterior
          DO I=1,7
            XX0(I)=XX(I)
            DXX0(I)=DXX(I)
          END DO
          CALL DOWNHILL(NDIM,XX0,DXX0,FUNKGSS2DA,1.0,0.5,2.0,YRMSTOL,
     +     XX,DXX,NEVAL,500)
          IF(LSHOW)THEN
            do i=1,7
              print*,i,xx(i),dxx(i)
            end do
          END IF
        END IF
C
        CTE=XX(1)
        AMP=XX(2)
        X0=XX(3)
        Y0=XX(4)
        SIGMAX=XX(5)
        SIGMAY=XX(6)
        BETA=XX(7)
C
        EECTE=DXX(1)
        EEAMP=DXX(2)
        EEX0=DXX(3)
        EEY0=DXX(4)
        EESIGMAX=DXX(5)
        EESIGMAY=DXX(6)
        EEBETA=DXX(7)
C
        AREA1=0.
        IF(NBUFF.EQ.0)THEN
          DO I=IY1,IY2
            DO J=IX1,IX2
              AREA1=AREA1+IMAGEN_(J,I)-CTE
            END DO
          END DO
        ELSE
          DO I=IY1,IY2
            DO J=IX1,IX2
              AREA1=AREA1+IMAGEN(J,I,NBUFF)-CTE
            END DO
          END DO
        END IF
        AREA2=PI2*SIGMAX*SIGMAY*AMP
        DO I=1,7
          CC(I)=XX(I)
        END DO
ccc        CALL QSIMP2D(REAL(IX1)-0.5,REAL(IX2)+0.5,
ccc     +               REAL(IY1)-0.5,REAL(IY2)+0.5,AREA3)
ccc        AREA3=AREA3-XX(1)*REAL(IX2-IX1+1)*REAL(IY2-IY1+1)
C------------------------------------------------------------------------------
        IF(NSIMUL.LT.2) RETURN      !si no hay que calcular errores, regresamos
        IF(NBUFF.GT.NMAXBUFF/2) RETURN        !idem si no hay imagen de errores
        IF(NBUFF.EQ.0) RETURN !idem si estamos con IMAGEN_()
        LSIMUL=.FALSE.
        DO I=IY1,IY2
          DO J=IX1,IX2
            IF(IMAGEN(J,I,NBUFF+NMAXBUFF/2).GT.0.) LSIMUL=.TRUE.
          END DO
        END DO
        IF(.NOT.LSIMUL) RETURN       !si todos los errores son cero, regresamos
C------------------------------------------------------------------------------
C guardamos los valores finales como valores de prueba para el calculo de
C errores
        DO J=1,7
          XX0(J)=XX(J)
          IF(DXX(J).GT.0.0)DXX0(J)=DXX(J) !de lo contrario usa el valor inicial
        END DO
C------------------------------------------------------------------------------
        NSEED=-1
        NEXTINFO=0
        DO ISIMUL=1,NSIMUL
          DO I=IY1,IY2
            DO J=IX1,IX2
              R1=RANDOMNUMBER(NSEED)
              R2=RANDOMNUMBER(NSEED)
              ERRYF=1.41421356*IMAGEN(J,I,NBUFF+NMAXBUFF/2)*
     +         SQRT(-1.*LOG(1.-R1))*COS(PI2*R2)
              IMAGEN_(J,I)=IMAGEN(J,I,NBUFF)+ERRYF
            END DO
          END DO
          IF(IMODE.EQ.0)THEN
            CALL DOWNHILL(NDIM,XX0,DXX0,FUNKGSS2D,1.0,0.5,2.0,YRMSTOL,
     +       XX,DXX,NEVAL,500)
          ELSE
            CALL DOWNHILL(NDIM,XX0,DXX0,FUNKGSS2DA,1.0,0.5,2.0,YRMSTOL,
     +       XX,DXX,NEVAL,500)
          END IF
          CTE_SIMUL(ISIMUL)=XX(1)
          AMP_SIMUL(ISIMUL)=XX(2)
          X0_SIMUL(ISIMUL)=XX(3)
          Y0_SIMUL(ISIMUL)=XX(4)
          SIGMAX_SIMUL(ISIMUL)=XX(5)
          SIGMAY_SIMUL(ISIMUL)=XX(6)
          BETA_SIMUL(ISIMUL)=XX(7)
          CALL SHOWPERC(1,NSIMUL,1,ISIMUL,NEXTINFO)
          AREA1_SIMUL(ISIMUL)=0.
          DO I=IY1,IY2
            DO J=IX1,IX2
              AREA1_SIMUL(ISIMUL)=AREA1_SIMUL(ISIMUL)+IMAGEN_(J,I)-XX(1)
            END DO
          END DO
          AREA2_SIMUL(ISIMUL)=PI2*XX(5)*XX(6)*XX(2)
          DO I=1,7
            CC(I)=XX(I)
          END DO
C las siguientes lineas estan comentadas porque en las simulaciones la integral
C numerica puede introducirse en un bucle cuasi infinito
ccc        print*,'entrando en QSIMP2D...'
ccc          CALL QSIMP2D(REAL(IX1)-0.5,REAL(IX2)+0.5,
ccc     +                 REAL(IY1)-0.5,REAL(IY2)+0.5,AREA3_SIMUL(ISIMUL))
ccc        print*,'...saliendo de QSIMP2D'
ccc          AREA3_SIMUL(ISIMUL)=AREA3_SIMUL(ISIMUL)-
ccc     +     XX(1)*REAL(IX2-IX1+1)*REAL(IY2-IY1+1)
        END DO
C calculamos rms
        FMEAN=FMEAN0(NSIMUL,CTE_SIMUL,ECTE)
        FMEAN=FMEAN0(NSIMUL,AMP_SIMUL,EAMP)
        FMEAN=FMEAN0(NSIMUL,X0_SIMUL,EX0)
        FMEAN=FMEAN0(NSIMUL,Y0_SIMUL,EY0)
        FMEAN=FMEAN0(NSIMUL,SIGMAX_SIMUL,ESIGMAX)
        FMEAN=FMEAN0(NSIMUL,SIGMAY_SIMUL,ESIGMAY)
        FMEAN=FMEAN0(NSIMUL,BETA_SIMUL,EBETA)
        FMEAN=FMEAN0(NSIMUL,AREA1_SIMUL,EAREA1)
        FMEAN=FMEAN0(NSIMUL,AREA2_SIMUL,EAREA2)
C la siguiente linea tambien esta comentada por lo mismo de antes
ccc        FMEAN=FMEAN0(NSIMUL,AREA3_SIMUL,EAREA3)
        EAREA3=0.
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C CTE=XX(1) AMP=XX(2) X0=XX(3) Y0=XX(4) SIGMAX=XX(5) SIGMAY=XX(6) BETA=XX(7)
        REAL FUNCTION FUNKGSS2D(XX)
        IMPLICIT NONE
        REAL XX(7)
C
        INCLUDE 'dimensions.inc'
C
        INTEGER IX1,IX2,IY1,IY2
        INTEGER I,J
        REAL IMAGEN_(NXMAX,NYMAX)
        REAL FACTOR,FACTOR1,FACTOR2,FACTOR3
        REAL FF
C
        COMMON/BLKIMAGEN1_/IMAGEN_
        COMMON/BLKFIT2DG1/IX1,IX2,IY1,IY2
C------------------------------------------------------------------------------
        FUNKGSS2D=0.
        DO I=IY1,IY2
          DO J=IX1,IX2
            FACTOR1=(REAL(J)-XX(3))*(REAL(J)-XX(3))/(2.*XX(5)*XX(5))
            FACTOR2=(REAL(I)-XX(4))*(REAL(I)-XX(4))/(2.*XX(6)*XX(6))
            FACTOR3=XX(7)*(REAL(J)-XX(3))*(REAL(I)-XX(4))/(XX(5)*XX(6))
            FACTOR=FACTOR1+FACTOR2+FACTOR3
            IF(FACTOR.GT.60.)THEN
              FF=XX(1)
            ELSE
              FF=XX(1)+XX(2)*EXP(-FACTOR)
            END IF
            IF(ABS(IMAGEN_(J,I)-FF).LT.1.E-15)THEN
            ELSE
              FUNKGSS2D=FUNKGSS2D+
     +         (IMAGEN_(J,I)-FF)*(IMAGEN_(J,I)-FF)
            END IF
          END DO
        END DO
        END
C
C******************************************************************************
C CTE=XX(1) AMP=XX(2) X0=XX(3) Y0=XX(4) SIGMAX=XX(5) SIGMAY=XX(6) BETA=XX(7)
        REAL FUNCTION FUNKGSS2DA(XX)
        IMPLICIT NONE
        REAL XX(7)
C
        INCLUDE 'dimensions.inc'
C
        INTEGER IX1,IX2,IY1,IY2
        INTEGER I,J
        REAL IMAGEN_(NXMAX,NYMAX)
        REAL CC(7)
        REAL FF
C
        COMMON/BLKIMAGEN1_/IMAGEN_
        COMMON/BLKFIT2DG1/IX1,IX2,IY1,IY2
        COMMON/BLKFUNK2DY/CC
C------------------------------------------------------------------------------
        DO I=1,7
          CC(I)=XX(I)
        END DO
        FUNKGSS2DA=0.
        DO I=IY1,IY2
          DO J=IX1,IX2
            CALL QSIMP2D(REAL(J)-0.5,REAL(J)+0.5,
     +                   REAL(I)-0.5,REAL(I)+0.5,FF)
            IF(ABS(IMAGEN_(J,I)-FF).LT.1.E-15)THEN
            ELSE
              FUNKGSS2DA=FUNKGSS2DA+
     +         (IMAGEN_(J,I)-FF)*(IMAGEN_(J,I)-FF)
            END IF
          END DO
        END DO
        END
C
C******************************************************************************
C******************************************************************************
C Subrutina del Numerical Recipes modificada para hacer integrales en 2D
      SUBROUTINE QSIMP2D(X1,X2,Y1,Y2,S)
      IMPLICIT NONE
      REAL X1,X2,Y1,Y2,S
C
      INTEGER J,JMAX
      REAL EPS
      REAL OST,OS
      REAL ST
C------------------------------------------------------------------------------
      PARAMETER (EPS=1.E-6, JMAX=20)
      OST=-1.E30
      OS= -1.E30
      DO 11 J=1,JMAX
        CALL TRAPZD2DX(X1,X2,Y1,Y2,ST,J)
        S=(4.*ST-OST)/3.
        IF (ABS(S-OS).LT.EPS*ABS(OS)) RETURN
        OS=S
        OST=ST
11    CONTINUE
ccc      PAUSE 'Too many steps.'
      END
C******************************************************************************
C Subrutina del Numerical Recipes modificada para hacer integrales en 2D
      SUBROUTINE TRAPZD2DX(X1,X2,Y1,Y2,S,N)
      IMPLICIT NONE
      REAL X1,X2,Y1,Y2,S
      INTEGER N
C
      REAL FUNK2DX
      EXTERNAL FUNK2DX
C
      INTEGER J,IT
      REAL DEL,TNM,X,SUM
C------------------------------------------------------------------------------
      IT=1 !evita un WARNING de compilacion
      IF (N.EQ.1) THEN
        S=0.5*(X2-X1)*(FUNK2DX(X1,Y1,Y2)+FUNK2DX(X2,Y1,Y2))
        IT=1
      ELSE
        TNM=IT
        DEL=(X2-X1)/TNM
        X=X1+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNK2DX(X,Y1,Y2)
          X=X+DEL
11      CONTINUE
        S=0.5*(S+(X2-X1)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END
C
C******************************************************************************
C integral entre Y1 y Y2 de f(x,y)dy
      REAL FUNCTION FUNK2DX(X,Y1,Y2)
      IMPLICIT NONE
      REAL X,Y1,Y2
C
      REAL S
C------------------------------------------------------------------------------
      CALL QSIMP2DY(X,Y1,Y2,S)
      FUNK2DX=S
      END
C
C******************************************************************************
C Subrutina del Numerical Recipes modificada para hacer integrales en 2D
      SUBROUTINE QSIMP2DY(X,Y1,Y2,S)
      IMPLICIT NONE
      REAL X,Y1,Y2,S
C
      INTEGER J,JMAX
      REAL EPS
      REAL OST,OS
      REAL ST
C------------------------------------------------------------------------------
      PARAMETER (EPS=1.E-6, JMAX=20)
      OST=-1.E30
      OS= -1.E30
      DO 11 J=1,JMAX
        CALL TRAPZD2DY(X,Y1,Y2,ST,J)
        S=(4.*ST-OST)/3.
        IF (ABS(S-OS).LT.EPS*ABS(OS)) RETURN
        OS=S
        OST=ST
11    CONTINUE
ccc      PAUSE 'Too many steps.'
      END
C******************************************************************************
C Subrutina del Numerical Recipes modificada para hacer integrales en 2D
      SUBROUTINE TRAPZD2DY(X,Y1,Y2,S,N)
      IMPLICIT NONE
      REAL X,Y1,Y2,S
      INTEGER N
C
      REAL FUNK2DY
      EXTERNAL FUNK2DY
C
      INTEGER J,IT
      REAL DEL,TNM,Y,SUM
C------------------------------------------------------------------------------
      IT=1 !evita WARNING de compilacion
      IF (N.EQ.1) THEN
        S=0.5*(Y2-Y1)*(FUNK2DY(X,Y1)+FUNK2DY(X,Y2))
        IT=1
      ELSE
        TNM=IT
        DEL=(Y2-Y1)/TNM
        Y=Y1+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNK2DY(X,Y)
          Y=Y+DEL
11      CONTINUE
        S=0.5*(S+(Y2-Y1)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END
C
C******************************************************************************
C f(x,y) 
      REAL FUNCTION FUNK2DY(X,Y)
      IMPLICIT NONE
      REAL X,Y
C
      REAL FACTOR,FACTOR1,FACTOR2,FACTOR3
      REAL CC(7)
      COMMON/BLKFUNK2DY/CC
C------------------------------------------------------------------------------
      FACTOR1=(X-CC(3))*(X-CC(3))/(2.*CC(5)*CC(5))
      FACTOR2=(Y-CC(4))*(Y-CC(4))/(2.*CC(6)*CC(6))
      FACTOR3=CC(7)*(X-CC(3))*(Y-CC(4))/(CC(5)*CC(6))
      FACTOR=FACTOR1+FACTOR2+FACTOR3
      IF(FACTOR.GT.60.)THEN
        FUNK2DY=CC(1)
      ELSE
        FUNK2DY=CC(1)+CC(2)*EXP(-FACTOR)
      END IF
C
      END
