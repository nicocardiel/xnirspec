! Busca los picos en un corte de la imagen. Da igual si el corte ha sido
! en el eje X o en el eje Y (corregido o sin corregir de distorsion).
! N: dimension de XP,YP
! XP(N),YP(N): corte en el cual se va a buscar los picos
! N1,N2: region en la que se van a buscar los picos
! NPEAKS: numero de picos encontrados
! XPEAKS,YPEAKS: posicion de los picos y sen~al en los mismos
!
        SUBROUTINE FINDMAX(N,N1,N2,XP,YP,NPEAKS,XPEAKS,YPEAKS)
        IMPLICIT NONE
!
        INCLUDE 'largest.inc'
!
        INTEGER N,N1,N2
        REAL XP(N),YP(N)
        INTEGER NPEAKS
        REAL XPEAKS(NXYMAX),YPEAKS(NXYMAX)
!
        INTEGER NWIDTHMAX
        PARAMETER (NWIDTHMAX=101)
!
        INTEGER I,J,I1,I2
        INTEGER NWIDTH
        INTEGER NMED
        REAL XF(NWIDTHMAX),YF(NWIDTHMAX)
        REAL COEFFPAR(3),CHISQR
        LOGICAL LGOOD
!
        COMMON/BLKDEFAULTS2/NWIDTH
!------------------------------------------------------------------------------
! usamos el algoritmo de findmax.f (REDUCEME) para la primera estimacion
! de los picos
        IF(MOD(NWIDTH,2).EQ.0)THEN
          NWIDTH=NWIDTH+1
          WRITE(*,101) '***WARNING***'
          WRITE(*,100) '=> Effective NWIDTH to be employed: '
          WRITE(*,*) NWIDTH
        END IF
        NMED=NWIDTH/2
!------------------------------------------------------------------------------
        IF(N2-N1+1.LT.NWIDTH)THEN
          WRITE(*,101) '***ERROR***'
          WRITE(*,100) '=> N1, N2, NWIDTH: '
          WRITE(*,*) N1,N2,NWIDTH
          WRITE(*,101) '=> region to search peaks shorter than NWITDH'
          WRITE(*,100) '=> (press <CR>...)'
          READ(*,*)
          RETURN
        END IF
        IF(NWIDTH.LT.3)THEN
          WRITE(*,101) '***WARNING***'
          WRITE(*,100) '=> NWIDTH: '
          WRITE(*,*) NWIDTH
          WRITE(*,101) '=> NWITDH < 3'
          NWIDTH=5
          WRITE(*,100) '=> Effective NWIDTH to be employed: '
          WRITE(*,*) NWIDTH
        END IF
        NPEAKS=0
        I1=N1+NMED
        I2=N2-NMED
        DO I=I1,I2
          LGOOD=.TRUE.
          J=0
          DO WHILE((J.LT.NMED).AND.(LGOOD))
            J=J+1
            IF(YP(I-NMED+J-1).GT.YP(I-NMED+J)) LGOOD=.FALSE.
          END DO
          J=NMED
          DO WHILE((J.LT.NWIDTH-1).AND.(LGOOD))
            J=J+1
            IF(YP(I-NMED+J-1).LT.YP(I-NMED+J)) LGOOD=.FALSE.
          END DO
          IF(LGOOD)THEN
            NPEAKS=NPEAKS+1
            XPEAKS(NPEAKS)=XP(I)
            YPEAKS(NPEAKS)=YP(I)
          END IF
        END DO
!------------------------------------------------------------------------------
! Ahora que ya tenemos una posicion inicial, la refinamos ajustado una
! parabola
        IF(NPEAKS.GT.0)THEN
          DO J=1,NPEAKS
            I1=NINT(XPEAKS(J))-NMED
            I2=NINT(XPEAKS(J))+NMED
            DO I=I1,I2
              XF(I-I1+1)=XP(I)
              YF(I-I1+1)=YP(I)
            END DO
            CALL POLFIT(XF,YF,YF,NWIDTH,3,0,COEFFPAR,CHISQR)
            XPEAKS(J)=-COEFFPAR(2)/2./COEFFPAR(3)     !el maximo de la parabola
          END DO
        END IF
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
