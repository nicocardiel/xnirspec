! Calcula la distancia de un punto a un polinomio, devolviendo asimismo las
! coordenadas del punto mas cercano sobre dicho polinomio.
! XC,YC: coordenadas de un punto generico
! NDEG: grado del polinomio
! COEFF: coeficientes del polinomio
! X0: valor inicial para iniciar la busqueda por DOWNHILL
! D: distancia minima
! XFIN,YFIN: valor final del punto
!------------------------------------------------------------------------------
        SUBROUTINE DISTPP(XC,YC,NDEG,COEFF,X0,D,XFIN,YFIN)
        IMPLICIT NONE
        REAL XC,YC
        INTEGER NDEG
        REAL COEFF(NDEG+1)
        REAL X0
        REAL D
        REAL XFIN,YFIN
!
        REAL FPOLY
        EXTERNAL YFUNKPMIN
        REAL YFUNKPMIN
! duplicamos las variables que tienen que pasar a YFUNKPMIN en COMMONs
        INTEGER NDEG_
        REAL COEFF_(20)
        REAL XC_,YC_
! variables locales
        INTEGER NEVAL
        INTEGER K
        REAL Y0
        REAL X1,Y1,X2,Y2
        REAL DIST2
        REAL ERRXFIN
!
        COMMON/BLKYFUNKPMIN1/NDEG_
        COMMON/BLKYFUNKPMIN2/COEFF_
        COMMON/BLKYFUNKPMIN3/XC_,YC_
!------------------------------------------------------------------------------
! caso trivial (el punto esta sobre el polinomio)
        Y0=FPOLY(NDEG,COEFF,X0)
        DIST2=(X0-XC)*(X0-XC)+(Y0-YC)*(Y0-YC)
        IF(DIST2.EQ.0.0)THEN
          XFIN=X0
          YFIN=Y0
          D=0.
          RETURN
        END IF
!------------------------------------------------------------------------------
! otro caso trivial: el polinomio es una recta
        IF(NDEG.EQ.1)THEN
          X1=X0
          Y1=FPOLY(NDEG,COEFF,X1)
          X2=X0+10.
          Y2=FPOLY(NDEG,COEFF,X2)
          CALL DISTPR(XC,YC,X1,Y1,X2,Y2,D,XFIN,YFIN)
          RETURN
        END IF
!------------------------------------------------------------------------------
! vamos con los casos no triviales
        NDEG_=NDEG
        DO K=1,NDEG+1
          COEFF_(K)=COEFF(K)
        END DO
        XC_=XC
        YC_=YC
! usamos DOWNHILL
        CALL DOWNHILL(1,X0,.1,YFUNKPMIN,1.0,0.5,2.0,1.E-6,XFIN,ERRXFIN,NEVAL,500)
        YFIN=FPOLY(NDEG,COEFF,XFIN)
        D=SQRT((XFIN-XC)*(XFIN-XC)+(YFIN-YC)*(YFIN-YC))
!!!        print*,'YFUNKPMIN> NEVAL: ',NEVAL
!!!        print*,'YFUNKPMIN> XFIN,ERRXFIN: ',XFIN,ERRXFIN,YFIN,D
!------------------------------------------------------------------------------
        END
!
!******************************************************************************
!
        REAL FUNCTION YFUNKPMIN(X)
        IMPLICIT NONE
        REAL X
!
        INTEGER NDEG
        REAL COEFF(20)
        REAL XC,YC
!
        REAL FPOLY
!
        REAL Y
!
        COMMON/BLKYFUNKPMIN1/NDEG
        COMMON/BLKYFUNKPMIN2/COEFF
        COMMON/BLKYFUNKPMIN3/XC,YC
!------------------------------------------------------------------------------
        Y=FPOLY(NDEG,COEFF,X)
        YFUNKPMIN=ABS(X-XC)+ABS(Y-YC)
!!!        print*,'FUNKPMIN>',ABS(X-XC)+ABS(Y-YC)
!
        END
