! Ajusta y dibuja un polinomio de la forma x=x(y)
        SUBROUTINE DRAWPOLXY(N,X,Y,NDEG,COEFF)
        IMPLICIT NONE
        INTEGER N
        REAL X(N),Y(N)
        INTEGER NDEG
!
        INCLUDE 'largest.inc'
!
        REAL ARCLENGTH
        REAL FPOLY
!
        INTEGER I,K
        REAL XP(NXYMAX),YP(NXYMAX)
        REAL COEFF(20),CHISQR
        REAL YMIN,YMAX
!------------------------------------------------------------------------------
        IF(N.LT.1)THEN
          WRITE(*,101) '***ERROR***'
          WRITE(*,101) '=> number of points < 1 in DRAWPOLYX.'
          WRITE(*,100) '(press <CR> to continue...)'
          READ(*,*)
          RETURN
        END IF
!
        CALL POLFIT(Y,X,X,N,NDEG+1,0,COEFF,CHISQR)
!
        WRITE(*,101) '=> fit result: '
        DO K=1,NDEG+1
          WRITE(*,'(A2,I2,A3,$)') 'a(',K-1,'): '
          WRITE(*,*) COEFF(K)
        END DO
!
        CALL FINDMML(N,1,N,Y,YMIN,YMAX)
        DO I=1,NXYMAX
          YP(I)=REAL(I-1)/REAL(NXYMAX-1)*(YMAX-YMIN)+YMIN
          XP(I)=FPOLY(NDEG,COEFF,YP(I))
        END DO
        CALL PGSCI(8)
        CALL PGLINE(NXYMAX,XP,YP)
        CALL PGSCI(1)
!
        WRITE(*,100) '=> Arc length: '
        WRITE(*,*) ARCLENGTH(NDEG,COEFF,YMIN,YMAX,100)
!
100     FORMAT(A,$)
101     FORMAT(A)
        END
