! Ajusta y dibuja un polinomio de la forma y=y(x)
        SUBROUTINE DRAWPOLYX(N,X,Y,NDEG,COEFF)
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
        REAL XMIN,XMAX
!------------------------------------------------------------------------------
        IF(N.LT.1)THEN
          WRITE(*,101) '***ERROR***'
          WRITE(*,101) '=> number of points < 1 in DRAWPOLYX.'
          WRITE(*,100) '(press <CR> to continue...)'
          READ(*,*)
          RETURN
        END IF
!
        CALL POLFIT(X,Y,Y,N,NDEG+1,0,COEFF,CHISQR)
!
        WRITE(*,101) '=> fit result: '
        DO K=1,NDEG+1
          WRITE(*,'(A2,I2,A3,$)') 'a(',K-1,'): '
          WRITE(*,*) COEFF(K)
        END DO
!
        CALL FINDMML(N,1,N,X,XMIN,XMAX)
        DO I=1,NXYMAX
          XP(I)=REAL(I-1)/REAL(NXYMAX-1)*(XMAX-XMIN)+XMIN
          YP(I)=FPOLY(NDEG,COEFF,XP(I))
        END DO
        CALL PGSCI(4)
        CALL PGLINE(NXYMAX,XP,YP)
        CALL PGSCI(1)
!
        WRITE(*,100) '=> Arc length: '
        WRITE(*,*) ARCLENGTH(NDEG,COEFF,XMIN,XMAX,100)
!
100     FORMAT(A,$)
101     FORMAT(A)
        END
