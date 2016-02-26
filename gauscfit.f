C------------------------------------------------------------------------------
C Version 28-February-1997                                     File: gauscfit.f
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
C SUBROUTINE GAUSCFIT(NPFIT,XFIT,YFIT,X0,SIGMA,AMP,Y0,
C                     EEX0,EESIGMA,EEAMP,EEY0,YRMSTOL)
C
C Input: NPFIT,XFIT,YFIT
C Input: YRMSTOL
C Output: X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAM,EEY0
C
C Fit numerically a gaussian + constant (using DOWNHILL):
C Y=Y0+AMP*EXP[-((X-X0)^2/(2*SIGMA^2))]
C
C INTEGER NPFIT -> number of points to be fitted
C REAL XFIT,YFIT -> x, y
C REAL X0 -> center of the fitted gaussian
C REAL SIGMA -> sigma value of the fitted gaussian
C REAL AMP -> maximum of the fitted gaussian
C REAL Y0 -> constant
C REAL EEX0,EESIGMA,EEAMP,EEY0 -> errors in X0,SIGMA,AMP,Y0 (rms from DOWNHILL)
C REAL YRMSTOL -> stopping criterion for DOWNHILL
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE GAUSCFIT(NPFIT,XFIT,YFIT,X0,SIGMA,AMP,Y0,
     +   EEX0,EESIGMA,EEAMP,EEY0,YRMSTOL)
        IMPLICIT NONE
        INTEGER NPFIT
        REAL XFIT(NPFIT),YFIT(NPFIT)
        REAL X0,SIGMA,AMP,Y0
        REAL EEX0,EESIGMA,EEAMP,EEY0
        REAL YRMSTOL
        REAL X0INI
C
        INCLUDE 'largest.inc'
C
        INTEGER NP
        INTEGER I
        REAL XF(NXYMAX),YF(NXYMAX)
        REAL XMIN,XMAX
        REAL A(3),CHISQR
        DOUBLE PRECISION MEAN,DISPER
C
        INTEGER NDIM,NEVAL
        REAL XX0(4),DXX0(4)
        REAL XX(4),DXX(4)
        EXTERNAL FUNKGC
        REAL FUNKGC
        COMMON/BLKFITG1/NP
        COMMON/BLKFITG2/XF,YF
        COMMON/BLKX0INI/X0INI
C
C------------------------------------------------------------------------------
        AMP=0.
        EEAMP=0.
        X0=0.
        EEX0=0.
        SIGMA=0.
        EESIGMA=0.
        Y0=0.
        EEY0=0.
C
        NP=NPFIT
        DO I=1,NP
          XF(I)=XFIT(I)
          YF(I)=YFIT(I)
        END DO
C
        IF(NP.LT.4)THEN
          WRITE(*,101)'FATAL ERROR: in subroutine GAUSCFIT.'
          WRITE(*,101)' No. of points for fit < 4'
          STOP
        END IF
C
C como primera estimacion del centro de la gaussiana calculamos el ajuste
C a una parabola y determinamos el maximo/minimo
        CALL POLFIT(XF,YF,YF,NP,3,0,A,CHISQR)
        X0INI=-REAL(A(2)/(2*A(3)))
        X0=0.1       !DOWNHILL supondra que la gaussiana esta alrededor de cero
C
C Como estimacion de la amplitud tomamos 3 veces el valor de la dispersion
C de los datos alrededor de la media
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
C si la parabola tiene un minimo (derivada segunda positiva), la 
C amplitud es negativa
        IF(A(3).GT.0.D0) AMP=-AMP
C
C Como estimacion del valor de SIGMA tomamos un cuarto del recorrido en la
C variable X
        XMIN=XF(1)
        XMAX=XMIN
        DO I=2,NP
          IF(XF(I).LT.XMIN) XMIN=XF(I)
          IF(XF(I).GT.XMAX) XMAX=XF(I)
        END DO
        SIGMA=(XMAX-XMIN)/4.
C
C Como estimacion de la cte. tomamos el valor medio de los datos
        Y0=REAL(MEAN)
ccc        type*,amp,x0,sigma,y0
C------------------------------------------------------------------------------
C Con estas estimaciones ya podemos calcular los parametros de entrada de
C la subrutina DOWNHILL
        NDIM=4
        XX0(1)=AMP
        XX0(2)=X0
        XX0(3)=SIGMA
        XX0(4)=Y0
        DXX0(1)=0.1*AMP
        DXX0(2)=5.
        DXX0(3)=0.1*SIGMA
        DXX0(4)=Y0+0.1*AMP
C
        CALL DOWNHILL(NDIM,XX0,DXX0,FUNKGC,1.0,0.5,2.0,YRMSTOL,
     +   XX,DXX,NEVAL,500)
C
        AMP=XX(1)
        X0=XX(2)+X0INI
        SIGMA=XX(3)
        Y0=XX(4)
        EEAMP=DXX(1)
        EEX0=DXX(2)
        EESIGMA=DXX(3)
        EEY0=DXX(4)
C------------------------------------------------------------------------------
101     FORMAT(A)
        END
C
C******************************************************************************
C
C XX(1)=AMP, XX(2)=X0, XX(3)=SIGMA, XX(4)=Y0
        REAL FUNCTION FUNKGC(XX)
        IMPLICIT NONE
        INCLUDE 'largest.inc'
        REAL XX(4)
C
        INTEGER I
        INTEGER NP
        REAL XF(NXYMAX),YF(NXYMAX)
        COMMON/BLKFITG1/NP
        COMMON/BLKFITG2/XF,YF
        COMMON/BLKX0INI/X0INI
        DOUBLE PRECISION DFUNK
        REAL FF,X0INI,FACTOR
C
C Introducimos un offset de X0INI para que el ajuste se realiza alrededor
C de x=0 (funciona mucho mejor)
        DFUNK=0.D0
        DO I=1,NP
          FACTOR=((XF(I)-XX(2)-X0INI)*(XF(I)-XX(2)-X0INI))/
     +     (2.*XX(3)*XX(3))
          IF(FACTOR.GT.60)THEN !evitamos: Floating point exception 5, underflow
            FF=XX(4)
          ELSE
            FF=XX(4)+XX(1)*EXP(-FACTOR)
          END IF
          DFUNK=DFUNK+DBLE((YF(I)-FF)*(YF(I)-FF))
        END DO
        FUNKGC=REAL(DFUNK)
        END
