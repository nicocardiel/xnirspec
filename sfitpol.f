C Dados una serie de N puntos X(N), Y(N), esta subrutina realiza ajuste
C polinomicos de forma interactiva
        SUBROUTINE SFITPOL(N,X,Y)
        IMPLICIT NONE
        INTEGER N
        REAL X(N),Y(N)
C
        INCLUDE 'largest.inc'
        INTEGER NSIMULMAX
        PARAMETER (NSIMULMAX=100)
C
        INTEGER READILIM
        REAL FPOLY
        REAL RANDOMNUMBER
C
        INTEGER I,IRAN
        INTEGER NDEG
        INTEGER NCOEFF
        INTEGER K
        INTEGER NSIMUL
        INTEGER NSEED
        INTEGER NUSED
        REAL COEFF(20)
        REAL COEFF_SIMU(20,NSIMULMAX)
        REAL VARCOEFF(20)
        REAL CHISQR,SR2
        REAL XP(NXYMAX),YP(NXYMAX)
        REAL COVAR(20,20)
        REAL XMIN,XMAX
        LOGICAL IFCOEFF(20)
        LOGICAL LEXAMINE
        LOGICAL LUSED(NXYMAX)
C------------------------------------------------------------------------------
        LEXAMINE=.TRUE.
        DO WHILE(LEXAMINE)
          NDEG=READILIM('Polynomial degree (0=exit)',
     +     '0',0,MIN0(19,N-1))
          IF(NDEG.EQ.0)THEN
            LEXAMINE=.FALSE.
          ELSE
            DO K=1,NDEG+1
              IFCOEFF(K)=.TRUE.
            END DO
            NCOEFF=0
            DO K=1,NDEG+1
              IF(IFCOEFF(K)) NCOEFF=NCOEFF+1
            END DO
c..............................................................................
C OLS(Y|X) sin pesos
            CALL FITPOL(N,X,Y,NDEG,IFCOEFF,
     +       COEFF,VARCOEFF,CHISQR,SR2,COVAR)
            WRITE(*,101) '=> OLS(Y|X):'
            CALL SHOW_FITPOL(COEFF,VARCOEFF,IFCOEFF,20,NDEG,N)
            CALL FINDMML(N,1,N,X,XMIN,XMAX)
            DO I=1,NXYMAX
              XP(I)=XMIN+REAL(I-1)/REAL(NXYMAX-1)*(XMAX-XMIN)
              YP(I)=FPOLY(NDEG,COEFF,XP(I))
            END DO
            CALL SUBPLOTBIS(NXYMAX,1,NXYMAX,XP,YP,XP,YP,
     +       .FALSE.,.FALSE.,NDEG,101,1.0)
c..............................................................................
C OLS(Y|X) con bootstrap
            NSEED=-1
            IF(NCOEFF.LT.N)THEN
              DO NSIMUL=1,NSIMULMAX
                NUSED=0
                DO WHILE(NUSED.LT.NCOEFF)
                  DO I=1,N
                    LUSED(I)=.FALSE.
                  END DO
                  DO I=1,N
                    IRAN=INT(REAL(N)*RANDOMNUMBER(NSEED))+1
                    XP(I)=X(IRAN)
                    YP(I)=Y(IRAN)
                    LUSED(IRAN)=.TRUE.
                  END DO
                  NUSED=0
                  DO I=1,N
                    IF(LUSED(I)) NUSED=NUSED+1
                  END DO
                END DO
                CALL FITPOL(N,XP,YP,NDEG,IFCOEFF,
     +           COEFF,VARCOEFF,CHISQR,SR2,COVAR)
                DO K=1,NDEG+1
                  COEFF_SIMU(K,NSIMUL)=COEFF(K)
                END DO
              END DO
              DO K=1,NDEG+1
                COEFF(K)=0.
                DO NSIMUL=1,NSIMULMAX
                  COEFF(K)=COEFF(K)+COEFF_SIMU(K,NSIMUL)
                END DO
                COEFF(K)=COEFF(K)/REAL(NSIMULMAX)
                VARCOEFF(K)=0
                DO NSIMUL=1,NSIMULMAX
                  VARCOEFF(K)=VARCOEFF(K)+
     >             (COEFF_SIMU(K,NSIMUL)-COEFF(K))*
     >             (COEFF_SIMU(K,NSIMUL)-COEFF(K))
                END DO
                VARCOEFF(K)=VARCOEFF(K)/REAL(NSIMULMAX-1)
              END DO
              WRITE(*,101) '=> OLS(Y|X) with bootstrap:'
              CALL SHOW_FITPOL(COEFF,VARCOEFF,IFCOEFF,20,NDEG,N)
            END IF
c..............................................................................
C OLS(Y|X) con jacknife
            IF(NCOEFF.LT.N)THEN
              DO NSIMUL=1,N
                DO I=1,N
                  IF(I.LT.NSIMUL)THEN
                    XP(I)=X(I)
                    YP(I)=Y(I)
                  ELSEIF(I.GT.NSIMUL)THEN
                    XP(I-1)=X(I)
                    YP(I-1)=Y(I)
                  END IF
                END DO
                CALL FITPOL(N-1,XP,YP,NDEG,IFCOEFF,
     +           COEFF,VARCOEFF,CHISQR,SR2,COVAR)
                DO K=1,NDEG+1
                  COEFF_SIMU(K,NSIMUL)=COEFF(K)
                END DO
              END DO
              DO K=1,NDEG+1
                COEFF(K)=0.
                DO NSIMUL=1,NSIMULMAX
                  COEFF(K)=COEFF(K)+COEFF_SIMU(K,NSIMUL)
                END DO
                COEFF(K)=COEFF(K)/REAL(NSIMULMAX)
                VARCOEFF(K)=0
                DO NSIMUL=1,NSIMULMAX
                  VARCOEFF(K)=VARCOEFF(K)+
     >             (COEFF_SIMU(K,NSIMUL)-COEFF(K))*
     >             (COEFF_SIMU(K,NSIMUL)-COEFF(K))
                END DO
                VARCOEFF(K)=VARCOEFF(K)/REAL(NSIMULMAX-1)
              END DO
              WRITE(*,101) '=> OLS(Y|X) with jacknife:'
              CALL SHOW_FITPOL(COEFF,VARCOEFF,IFCOEFF,20,NDEG,N)
            END IF
c..............................................................................
          END IF
        END DO
C
101     FORMAT(A)
        END
