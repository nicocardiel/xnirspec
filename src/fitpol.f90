! Version May 8, 2000: 
!       Modifico un poco la subrutina para que los parametros
!       de entrada esten en formato real. El calculo sigue haciendose en doble
!       precision.
!******************************************************************************
! Ajusta un polinomio generico por minimos cuadrados (sin errores)
        SUBROUTINE FITPOL(N,X_,Y_,NDEG,IFCOEF,B_,VARB_,CHISQR_,SR2_,COVAR_)
        IMPLICIT NONE
!
        INCLUDE 'largest.inc'
        INTEGER NDEGMAX         !grado maximo del polinomio que puede ajustarse
        PARAMETER (NDEGMAX=19)
!
        INTEGER N
        REAL X_(N),Y_(N)
        INTEGER NDEG
        LOGICAL IFCOEF(NDEGMAX+1)
        REAL B_(NDEGMAX+1),VARB_(NDEGMAX+1)
        REAL CHISQR_,SR2_
        REAL COVAR_(NDEGMAX+1,NDEGMAX+1)
!
        DOUBLE PRECISION X(NXYMAX),Y(NXYMAX)
        DOUBLE PRECISION B(NDEGMAX+1),VARB(NDEGMAX+1)
        DOUBLE PRECISION CHISQR,SR2
        DOUBLE PRECISION COVAR(NDEGMAX+1,NDEGMAX+1)
!
        INTEGER K
        INTEGER I,J,L
        INTEGER II,JJ
        INTEGER NCOEF
        INTEGER ORDER(NDEGMAX+1)
        INTEGER IOK,IPAR
        DOUBLE PRECISION A(NDEGMAX+1,NDEGMAX+1)
        DOUBLE PRECISION G(NDEGMAX+1)
        DOUBLE PRECISION C1(NDEGMAX+1),C2(NDEGMAX+1)
        DOUBLE PRECISION SCALEROW(NDEGMAX+1)
        DOUBLE PRECISION DSUM
        DOUBLE PRECISION BB(NDEGMAX+1)
        DOUBLE PRECISION X0,Y0
!------------------------------------------------------------------------------
        IF(N.GT.NXYMAX)THEN
          INCLUDE 'deallocate_arrays.inc'
          STOP 'FATAL ERROR: N.GT.NXYMAX in FITPOL'
        END IF
        IF(NDEG.GT.NDEGMAX)THEN
          INCLUDE 'deallocate_arrays.inc'
          STOP 'FATAL ERROR: NDEG.GT.NDEGMAX in FITPOL'
        END IF
!------------------------------------------------------------------------------
! pasamos los datos a doble precision
        DO K=1,N
          X(K)=DBLE(X_(K))
          Y(K)=DBLE(Y_(K))
        END DO
! pasamos los coeficientes a doble precision (esta subrutina puede hacer ajuste
! fijando el valor de algunos parametros que en ese caso no se ajustan, tomandose
! el valor de los coeficientes de entrada)
        DO I=1,NDEG+1
          B(I)=DBLE(B_(I))
        END DO
!------------------------------------------------------------------------------
! numero de coeficientes a ajustar
        NCOEF=0
        DO I=1,NDEG+1
          IF(IFCOEF(I)) NCOEF=NCOEF+1
        END DO
!------------------------------------------------------------------------------
! creamos matrices para plantear el sistema de ecuaciones a resolver
        II=0
        DO I=1,NDEG+1
          IF(IFCOEF(I))THEN
            II=II+1
            JJ=0
            DO J=1,NDEG+1
              IF(IFCOEF(J))THEN
                JJ=JJ+1
                IF((I.EQ.1).AND.(J.EQ.1))THEN
                  A(1,1)=DBLE(N)
                ELSE
                  A(JJ,II)=0.D0
                  DO K=1,N
                    A(JJ,II)=A(JJ,II)+(X(K)**(I+J-2))
                  END DO
                END IF
              END IF
            END DO
            G(II)=0.D0
            DO K=1,N
              DSUM=Y(K)
              DO L=1,NDEG+1
                IF(.NOT.IFCOEF(L))THEN
                  IF(L.EQ.1)THEN
                    DSUM=DSUM-B(L)
                  ELSE
                    DSUM=DSUM-(B(L)*X(K)**(L-1))
                  END IF
                END IF
              END DO
              IF(I.EQ.1)THEN
                G(II)=G(II)+DSUM
              ELSE
                G(II)=G(II)+DSUM*(X(K)**(I-1))
              END IF
            END DO
          END IF
        END DO
!------------------------------------------------------------------------------
! resolvemos el sistema de ecuaciones
        CALL LUDCMPD(A,NCOEF,NDEGMAX+1,ORDER,SCALEROW,IOK,IPAR)
        CALL LUSOLVD(A,NCOEF,NDEGMAX+1,ORDER,SCALEROW,G,BB)
!------------------------------------------------------------------------------
! introducimos los coeficientes calculados en la matriz de salida, sin
! modificar los coeficientes fijos
        II=0
        DO I=1,NDEG+1
          IF(IFCOEF(I))THEN
            II=II+1
            B(I)=BB(II)
          END IF
        END DO
!------------------------------------------------------------------------------
! varianza residual y Chi-cuadrado
        SR2=0.D0
        IF(N.GT.NDEG+1)THEN
          DO K=1,N
            X0=X(K)
            Y0=B(NDEG+1)
            DO L=NDEG,1,-1
              Y0=Y0*X0+B(L)
            END DO
            SR2=SR2+(Y(K)-Y0)*(Y(K)-Y0)
          END DO
          CHISQR=SR2
          SR2=SR2/DBLE(N-(NDEG+1))
        ELSE
          CHISQR=0.D0
        END IF
!------------------------------------------------------------------------------
! calculamos la matriz de varianza-covarianza (=matriz inversa de A)
        DO I=1,NCOEF
          DO J=1,NCOEF
            IF(I.EQ.J)THEN
              C1(J)=1.D0
            ELSE
              C1(J)=0.D0
            END IF
          END DO
          CALL LUSOLVD(A,NCOEF,NDEGMAX+1,ORDER,SCALEROW,C1,C2)
          DO J=1,NCOEF
            COVAR(J,I)=C2(J)
          END DO
        END DO
!------------------------------------------------------------------------------
! introducimos las varianzas calculadas en la matriz de salida
        II=0
        DO I=1,NDEG+1
          IF(IFCOEF(I))THEN
            II=II+1
            VARB(I)=COVAR(II,II)*SR2
          ELSE
            VARB(I)=0.D0
          END IF
        END DO
!------------------------------------------------------------------------------
! las variables de salida no estan en doble precision
        DO I=1,NDEG+1
          IF(IFCOEF(I))THEN
            B_(I)=REAL(B(I))
            IF(VARB(I).LE.1.E-30)THEN
              VARB_(I)=0.
!!! las 2 lineas anteriores evitan: Floating point exception 5, underflow
            ELSE
              VARB_(I)=REAL(VARB(I))
            END IF
          ELSE
            VARB_(I)=0.
          END IF
        END DO
        CHISQR_=REAL(CHISQR)
        SR2_=REAL(SR2)
        DO I=1,NCOEF
          DO J=1,NCOEF
            COVAR_(J,I)=REAL(COVAR(J,I))
          END DO
        END DO
!------------------------------------------------------------------------------
        END
