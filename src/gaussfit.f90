!------------------------------------------------------------------------------
! Version 23-July-1998                                         File: gaussfit.f
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
! SUBROUTINE GAUSSFIT(NPFIT,XFIT,YFIT,EYFIT,X0,SIGMA,AMP,
!                     EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP,YRMSTOL,NSIMUL)
!
! Input: NPFIT,XFIT,YFIT,EYFIT,YRMSTOL,NSIMUL
! Output: X0,SIGMA,AMP,EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP
!
! Numerical fit of a gaussian (using DOWNHILL):
! Y=AMP*EXP[-((X-X0)^2/(2*SIGMA^2))]
!
! INTEGER NPFIT -> number of points to be fitted
! REAL XFIT,YFIT,EYFIT -> x, y and error
! REAL X0 -> center of the fitted gaussian
! REAL SIGMA -> sigma value of the fitted gaussian
! REAL AMP -> maximum of the fitted gaussian
! REAL EX0,ESIGMA,EAMP -> errors in X0,SIGMA,AMP (due to EYFIT --simulations--)
! REAL EEX0,EESIGMA,EEAMP -> errors in X0,SIGMA,AMP (rms from DOWNHILL)
! REAL YRMSTOL -> stopping criterion for DOWNHILL
! INTEGER NSIMUL -> number of simulations to compute errors
!
!omment
!------------------------------------------------------------------------------
        SUBROUTINE GAUSSFIT(NPFIT,XFIT,YFIT,EYFIT,X0,SIGMA,AMP, &
         EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP,YRMSTOL,NSIMUL) 
        IMPLICIT NONE
        INTEGER NPFIT
        REAL XFIT(NPFIT),YFIT(NPFIT),EYFIT(NPFIT)
        REAL X0,SIGMA,AMP,EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP
        REAL YRMSTOL
        INTEGER NSIMUL
!
        INCLUDE 'largest.inc'
!
        INTEGER NSIMULMAX
        PARAMETER(NSIMULMAX=1000)                !numero maximo de simulaciones
        REAL PI2
        PARAMETER(PI2=6.283185307)               !2 x pi
!
        INTEGER I,J,ISIMUL
        INTEGER NP
        INTEGER NSEED
        REAL X0INI
        REAL XMIN,XMAX
        REAL XF(NXYMAX),YF(NXYMAX),ERRYF
        REAL A(3),CHISQR
        REAL X0_SIMUL(NSIMULMAX)
        REAL SIGMA_SIMUL(NSIMULMAX)
        REAL AMP_SIMUL(NSIMULMAX)
        REAL RANDOMNUMBER,R1,R2
        DOUBLE PRECISION MEAN,DISPER
!
        INTEGER NDIM,NEVAL
        REAL XX0(3),DXX0(3)      !valores iniciales y desplazamientos de prueba
        REAL XX(3),DXX(3)
        EXTERNAL FUNKGSS
        REAL FUNKGSS
!
        COMMON/BLKFITG1/NP
        COMMON/BLKFITG2/XF,YF
        COMMON/BLKX0INI/X0INI
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
        AMP=0.
        EAMP=0.
        EEAMP=0.
        X0=0.
        EX0=0.
        EEX0=0.
        SIGMA=0.
        ESIGMA=0.
        EESIGMA=0.
!
        IF(NPFIT.LT.3)THEN
          WRITE(*,101)'ERROR: in subroutine GAUSSFIT:'
          WRITE(*,101)' No. of points for fit < 3'
          WRITE(*,100)'Press <CR> to continue...'
          RETURN
        ELSEIF(NPFIT.GT.NXYMAX)THEN
          WRITE(*,101)'ERROR: in subroutine GAUSSFIT:'
          WRITE(*,101)' No. of points for fit > NXYMAX'
          WRITE(*,100)'Press <CR> to continue...'
          RETURN
        END IF
!
        IF(NSIMUL.GT.NSIMULMAX)THEN
          WRITE(*,101)'ERROR: in subroutine GAUSSFIT:'
          WRITE(*,101)' No. of simulations > NSIMULMAX'
          WRITE(*,100)'Press <CR> to continue...'
          RETURN
        END IF
!
        NP=NPFIT
!------------------------------------------------------------------------------
        DO I=1,NP
          XF(I)=XFIT(I)
          YF(I)=YFIT(I)
        END DO
!------------------------------------------------------------------------------
! como primera estimacion del centro de la gaussiana calculamos el ajuste
! a una parabola y determinamos el maximo/minimo
        CALL POLFIT(XF,YF,YF,NP,3,0,A,CHISQR)
        X0INI=-REAL(A(2)/(2*A(3)))
        X0=0.0           !DOWNHILL supondra que el maximo esta alrededor de x=0
! como estimacion de la amplitud tomamos 3 veces el valor de la dispersion
! de los datos alrededor de la media
        MEAN=0.D0
        DO I=1,NP
          MEAN=MEAN+DBLE(YF(I))
        END DO
        MEAN=MEAN/DBLE(NP)
        DISPER=0.D0
        DO I=1,NP
          DISPER=DISPER+(DBLE(YF(I))-MEAN)*(DBLE(YF(I))-MEAN)
        END DO
        DISPER=DSQRT(DISPER/DBLE(NP-1))
        AMP=3.*REAL(DISPER)
! si la parabola tiene un minimo (derivada segunda positiva), la 
! amplitud es negativa
        IF(A(3).GT.0.D0) AMP=-AMP
!
! Como estimacion del valor de SIGMA tomamos un cuarto del recorrido en la
! variable X
        XMIN=XF(1)
        XMAX=XMIN
        DO I=2,NP
          IF(XF(I).LT.XMIN) XMIN=XF(I)
          IF(XF(I).GT.XMAX) XMAX=XF(I)
        END DO
        SIGMA=(XMAX-XMIN)/4.
!------------------------------------------------------------------------------
! Con estas estimaciones ya podemos calcular los parametros de entrada de
! la subrutina DOWNHILL
        NDIM=3
        XX0(1)=AMP
        XX0(2)=X0
        XX0(3)=SIGMA
        DXX0(1)=0.1*AMP
        DXX0(2)=5.
        DXX0(3)=0.1*SIGMA
!
        CALL DOWNHILL(NDIM,XX0,DXX0,FUNKGSS,1.0,0.5,2.0,YRMSTOL,XX,DXX,NEVAL,500)
!
        AMP=XX(1)
        X0=XX(2)+X0INI
        SIGMA=XX(3)
        EEAMP=DXX(1)
        EEX0=DXX(2)
        EESIGMA=DXX(3)
!------------------------------------------------------------------------------
        IF(NSIMUL.LT.2)RETURN       !si no hay que calcular errores, regresamos
!------------------------------------------------------------------------------
! guardamos los valores finales como valores de prueba para el calculo de
! errores
        DO J=1,3
          XX0(J)=XX(J)
          IF(DXX(J).GT.0.0)DXX0(J)=DXX(J) !de lo contrario usa el valor inicial
        END DO
!------------------------------------------------------------------------------
        NSEED=-1
        DO ISIMUL=1,NSIMUL
          DO I=1,NP
            R1=RANDOMNUMBER(NSEED)
            R2=RANDOMNUMBER(NSEED)
            ERRYF=1.41421356*EYFIT(I)*SQRT(-1.*LOG(1.-R1))*COS(PI2*R2)
            YF(I)=YFIT(I)+ERRYF
          END DO
          CALL DOWNHILL(NDIM,XX0,DXX0,FUNKGSS,1.0,0.5,2.0,YRMSTOL,XX,DXX,NEVAL,500)
          AMP_SIMUL(ISIMUL)=XX(1)
          X0_SIMUL(ISIMUL)=XX(2)+X0INI
          SIGMA_SIMUL(ISIMUL)=XX(3)
        END DO
!
! EAMP: error en AMP
        MEAN=0.D0
        DO ISIMUL=1,NSIMUL
          MEAN=MEAN+DBLE(AMP_SIMUL(ISIMUL))
        END DO
        MEAN=MEAN/DBLE(NSIMUL)
        DISPER=0.D0
        DO ISIMUL=1,NSIMUL
          DISPER=DISPER+(DBLE(AMP_SIMUL(ISIMUL))-MEAN)*(DBLE(AMP_SIMUL(ISIMUL))-MEAN)
        END DO
        DISPER=DSQRT(DISPER/DBLE(NSIMUL-1))
        EAMP=REAL(DISPER)
!
! EX0: error en X0
        MEAN=0.D0
        DO ISIMUL=1,NSIMUL
          MEAN=MEAN+DBLE(X0_SIMUL(ISIMUL))
        END DO
        MEAN=MEAN/DBLE(NSIMUL)
        DISPER=0.D0
        DO ISIMUL=1,NSIMUL
          DISPER=DISPER+(DBLE(X0_SIMUL(ISIMUL))-MEAN)*(DBLE(X0_SIMUL(ISIMUL))-MEAN)
        END DO
        DISPER=DSQRT(DISPER/DBLE(NSIMUL-1))
        EX0=REAL(DISPER)
!
! ESIGMA: error en SIGMA
        MEAN=0.D0
        DO ISIMUL=1,NSIMUL
          MEAN=MEAN+DBLE(SIGMA_SIMUL(ISIMUL))
        END DO
        MEAN=MEAN/DBLE(NSIMUL)
        DISPER=0.D0
        DO ISIMUL=1,NSIMUL
          DISPER=DISPER+(DBLE(SIGMA_SIMUL(ISIMUL))-MEAN)*(DBLE(SIGMA_SIMUL(ISIMUL))-MEAN)
        END DO
        DISPER=DSQRT(DISPER/DBLE(NSIMUL-1))
        ESIGMA=REAL(DISPER)
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
!
!******************************************************************************
!
! XX(1)=AMP, XX(2)=X0, XX(3)=SIGMA
        REAL FUNCTION FUNKGSS(XX)
        IMPLICIT NONE
        REAL XX(3)
!
        INCLUDE 'largest.inc'
!
        INTEGER I
        INTEGER NF
        REAL XF(NXYMAX),YF(NXYMAX)
        COMMON/BLKFITG1/NF
        COMMON/BLKFITG2/XF,YF
        COMMON/BLKX0INI/X0INI
        REAL FF,X0INI,FACTOR
!
! Introducimos X0INI para que el ajuste se realice alrededor de x=0. De esta
! forma el ajuste es mejor
        FUNKGSS=0.
        DO I=1,NF
          FACTOR=((XF(I)-XX(2)-X0INI)*(XF(I)-XX(2)-X0INI))/(2.*XX(3)*XX(3))
          IF(FACTOR.GT.60)THEN !evitamos: Floating point exception 5, underflow
            FF=0.
          ELSE
            FF=XX(1)*EXP(-FACTOR)
          END IF
          IF(ABS(YF(I)-FF).LT.1.E-15)THEN
          ELSE
            FUNKGSS=FUNKGSS+(YF(I)-FF)*(YF(I)-FF)
          END IF
        END DO
        END
