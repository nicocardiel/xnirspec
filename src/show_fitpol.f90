! Los coeficientes son reales.
!******************************************************************************
! Muestra los coeficientes del ajuste polinomico en pantalla, realizando
! un analisis estadistico de la significacion de dichos coeficientes.
        SUBROUTINE SHOW_FITPOL(A,VARA,IFCOEF,NDEGMAX,NDEG,NF)
        IMPLICIT NONE
!
        INTEGER NDEGMAX
        REAL A(NDEGMAX+1)
        REAL VARA(NDEGMAX+1)
        LOGICAL IFCOEF(NDEGMAX+1)
        INTEGER NDEG,NF
!
        INTEGER I
        REAL T,FTSTUDENT
!------------------------------------------------------------------------------
        DO I=1,NDEG+1
          IF(IFCOEF(I))THEN
            IF(I.LE.10)THEN
              WRITE(*,'(A4,I1,A3,$)') '  a(',I-1,'): '
            ELSE
              WRITE(*,'(A3,I2,A3,$)') ' a(',I-1,'): '
            END IF
            IF(VARA(I).EQ.0.0)THEN
              WRITE(*,*) A(I),SQRT(VARA(I))
            ELSE
              T=ABS(A(I))/SQRT(VARA(I))
              WRITE(*,*) A(I),SQRT(VARA(I)),T,2.0*FTSTUDENT(NF-(NDEG+1),T)
            END IF
          ELSE
            IF(I.LE.10)THEN
              WRITE(*,'(A4,I1,A3,$)') ' *a(',I-1,'): '
            ELSE
              WRITE(*,'(A3,I2,A3,$)') '*a(',I-1,'): '
            END IF
            WRITE(*,*) A(I),SQRT(VARA(I))
          END IF
        END DO
!
        END
