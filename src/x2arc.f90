! Sea un polinomio definido entre dos puntos A y B. Si S es la fraccion de arco
! que va desde el punto A hasta el punto B, S varia entre 0 y 1. Dado un valor
! inicial de S, la subrutina determina el valor X0 en el cual se obtiene dicha
! fraccion de arco.
! NOTA: S puede en realidad ser mayor que 0 o menor que 1.
! 
        SUBROUTINE X2ARC(NDEG,COEFF,XA,XB,S,X0)
        IMPLICIT NONE
        INTEGER NDEG                                       !grado del polinomio
        REAL COEFF(NDEG+1)                          !coeficientes del polinomio
        REAL XA,XB              !limites en X entre los que S varia entre 0 y 1
        REAL S                          !valor para el cual queremos conocer X0
        REAL X0                                          !solucion del problema
!
        EXTERNAL YFUNKDARC
        REAL YFUNKDARC
! duplicamos las variables que tienen que pasar a YFUNKDARC en COMMONs
        INTEGER NDEG_
        REAL COEFF_(20)
        REAL XA_,XB_,S_
!
        INTEGER K
        INTEGER NEVAL
        REAL XINI,XFIN,ERRXFIN
!
        COMMON/BLKYFUNKDARC1/NDEG_
        COMMON/BLKYFUNKDARC2/COEFF_
        COMMON/BLKYFUNKDARC3/XA_,XB_,S_
!------------------------------------------------------------------------------
! casos triviales
        IF(NDEG.LT.0)THEN
          WRITE(*,101) '***FATAL ERROR***'
          WRITE(*,101) '=> NDEG<0 in X2ARC'
          INCLUDE 'deallocate_arrays.inc'
          STOP
        ELSEIF((NDEG.EQ.0).OR.(NDEG.EQ.1))THEN
          X0=XA+S*(XB-XA)
          RETURN
        END IF
! mas casos triviales (independientes del grado del polinomio)
        IF(S.EQ.0.0)THEN
          X0=XA
          RETURN
        ELSEIF(S.EQ.1.0)THEN
          X0=XB
          RETURN
        END IF
!------------------------------------------------------------------------------
! Y ahora lidiemos con los casos no triviales (que seran los normales)
!------------------------------------------------------------------------------
! duplicamos las variables que pasan a YFUNKDARC en COMMONs
        NDEG_=NDEG
        DO K=1,NDEG+1
          COEFF_(K)=COEFF(K)
        END DO
        XA_=XA
        XB_=XB
        S_=S
! usamos DOWNHILL para resolver nuestro problema
        XINI=XA+S*(XB-XA)       !valor inicial: suponemos una recta entre A y B
        CALL DOWNHILL(1,XINI,0.1,YFUNKDARC,1.0,0.5,2.0,1.E-7,XFIN,ERRXFIN,NEVAL,500)
        X0=XFIN
!!!        print*,'YFUNKDARC> NEVAL: ',NEVAL
!!!        print*,'YFUNKDARC> XFIN,ERRXFIN: ',XFIN,ERRXFIN
!
101     FORMAT(A)
        END
!
!******************************************************************************
!
        REAL FUNCTION YFUNKDARC(X)
        IMPLICIT NONE
        REAL X
!
        REAL ARCLENGTH
!
        INTEGER NDEG                                       !grado del polinomio
        REAL COEFF(20)                              !coeficientes del polinomio
        REAL XA,XB              !limites en X entre los que S varia entre 0 y 1
        REAL S                          !valor para el cual queremos conocer X0
!
        REAL S0,S1
        COMMON/BLKYFUNKDARC1/NDEG
        COMMON/BLKYFUNKDARC2/COEFF
        COMMON/BLKYFUNKDARC3/XA,XB,S
!------------------------------------------------------------------------------
        S0=ARCLENGTH(NDEG,COEFF,XA,XB,100)
        S1=ARCLENGTH(NDEG,COEFF,XA,X,100)
        YFUNKDARC=ABS(S1/S0-S)
!!!        print*,'YFUNKDARC>',s,s1/s0,ABS(S1/S0-S)
!
        END
