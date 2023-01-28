! Refina la posicion de una linea BA del boundary.
! Nota: el polinomio de entrada en la subrutina es modificados a la salida.
! X0 es un valor de X apropiado para buscar las intersecciones con los
! espectros del borde.
        SUBROUTINE REFINEBA(NDEG,COEFF,X0,FMEAN,FSIGMA)
        USE Dynamic_Array_IMAGEN
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
! subroutine arguments
        INTEGER NDEG
        REAL COEFF(NDEG+1)
        REAL X0,FMEAN,FSIGMA
!
        INCLUDE 'largest.inc'
        INCLUDE 'dimensions.inc'
!
        REAL FPOLY
        REAL FMEAN0,FMEAN2
!
        INTEGER I,J,I1,I2,J1,J2
        INTEGER NWIDTH,NMED
        INTEGER NLINBL,NDEGBL00
        INTEGER NDEGBL(NXYMAX)
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER NPEAKS
        INTEGER NCBUFF
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL COEFFBL(20,NXYMAX)
        REAL TSIGMA
        REAL XPEAKS(NXYMAX),YPEAKS(NXYMAX),RPEAKS(NXYMAX)
        REAL XF(NXYMAX),YF(NXYMAX)
        REAL XMIN,YMIN,XMAX,YMAX
        REAL XXMIN(1),YYMIN(1),XXMAX(1),YYMAX(1)
        REAL XX,YY
!!!        REAL COEFFPAR(3),CHISQR
        REAL X0_,Y0_,AMP,SIGMA
        REAL EEX0_,EEY0_,EEAMP,EESIGMA
        REAL XMIN_,XMAX_,YMIN_,YMAX_
        LOGICAL LDEBUGLOCAL
        LOGICAL LFIT(NXYMAX)
!
!delete COMMON/BLKIMAGEN1/IMAGEN             !imagen FITS leida en formato REAL
        COMMON/BLKIMAGEN2/NCBUFF
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKBOUND1/NLINBL,NDEGBL,NDEGBL00
        COMMON/BLKBOUND2/COEFFBL
        COMMON/BLKDEFAULTS2/NWIDTH
        COMMON/BLKDEFAULTS3/TSIGMA
        COMMON/BLKPLIMITS/XMIN_,XMAX_,YMIN_,YMAX_
!------------------------------------------------------------------------------
        LDEBUGLOCAL=.TRUE.
! intersecciones del polinomio temporal con los espectros del borde
        CALL INTERSEC(NDEG,COEFF,NDEGBL(1),COEFFBL(1,1),X0,XMIN,YMIN)
        CALL INTERSEC(NDEG,COEFF,NDEGBL(NLINBL),COEFFBL(1,NLINBL),X0,XMAX,YMAX)
        IF(LDEBUGLOCAL)THEN
          CALL PGSCI(5)
          !usamos un array unidimensional porque el compilador
          !gfortran-mp-10 da error al usar un escalar en lugar
          !de una matriz
          XXMIN(1)=XMIN
          XXMAX(1)=XMAX
          YYMIN(1)=YMIN
          YYMAX(1)=YMAX
          CALL PGPOINT(1,XXMIN,YYMIN,17)
          CALL PGPOINT(1,XXMAX,YYMAX,17)
          CALL PGSCI(1)
        END IF
!------------------------------------------------------------------------------
        NMED=NWIDTH/2
        I1=NINT(YMIN)
        I2=NINT(YMAX)
        DO I=I1,I2
          YY=REAL(I)
          XX=FPOLY(NDEG,COEFF,YY)
          J1=NINT(XX)-NMED
          J2=NINT(XX)+NMED
          IF(J1.LT.1)THEN
            J1=1
            J2=J1+NWIDTH-1
          ELSEIF(J2.GT.NAXIS(1,NCBUFF))THEN
            J2=NAXIS(1,NCBUFF)
            J1=J2-NWIDTH+1
          END IF
          DO J=J1,J2
            XF(J-J1+1)=REAL(J)
            YF(J-J1+1)=IMAGEN(J,I,NCBUFF)
          END DO
          CALL GAUSCFIT(NWIDTH,XF,YF,X0_,SIGMA,AMP,Y0_,EEX0_,EESIGMA,EEAMP,EEY0_,1.E-6)
!!!          CALL POLFIT(XF,YF,YF,NWIDTH,3,0,COEFFPAR,CHISQR)
!!!          XPEAKS(I-I1+1)=-COEFFPAR(2)/2./COEFFPAR(3)
          XPEAKS(I-I1+1)=X0_
          YPEAKS(I-I1+1)=REAL(I)
        END DO
        NPEAKS=I2-I1+1
        CALL POLFITSIG(NPEAKS,YPEAKS,XPEAKS,TSIGMA,NDEG,COEFF,LFIT)
! en lugar de dibujar el ajuste directamente es mejor mostrar los residuos
        DO J=1,NPEAKS
          RPEAKS(J)=XPEAKS(J)-FPOLY(NDEG,COEFF,YPEAKS(J))
        END DO
        IF(LDEBUGLOCAL)THEN
          CALL PGSCI(2)
          CALL PGPOINT(NPEAKS,XPEAKS,YPEAKS,1)
          CALL PGSCI(1)
          YMIN_=-1.
          YMAX_=+1.
          CALL SUBPLOT(NPEAKS,1,NPEAKS,YPEAKS,RPEAKS,YPEAKS,RPEAKS,.TRUE.,.FALSE.,.FALSE.,.FALSE., &
           'spatial direction','residuals in wavelength direction',' ',3,1,1.5)
          J=0
          DO I=1,NPEAKS
            IF(.NOT.LFIT(I))THEN
              J=J+1
              YF(J)=YPEAKS(I)
              XF(J)=RPEAKS(I)
            END IF
          END DO
          IF(J.GT.0) CALL SUBPLOTBIS(J,1,J,YF,XF,YF,XF,.FALSE.,.FALSE.,2,5,1.0)
          DO I=1,NXYMAX
            YF(I)=YMIN+REAL(I-1)/REAL(NXYMAX-1)*(YMAX-YMIN)
            XF(I)=FPOLY(NDEG,COEFF,YF(I))
          END DO
          CALL PGSCI(3)
          CALL PGLINE(NXYMAX,XF,YF)
          CALL PGSCI(1)
          CALL SUBPLOTBIS(NXYMAX,1,NXYMAX,YF,XF,YF,XF,.FALSE.,.FALSE.,4,101,1.0)
          FMEAN=FMEAN0(NPEAKS,RPEAKS,FSIGMA)
          WRITE(*,100) '> Mean, r.m.s.............: '
          WRITE(*,*) FMEAN,FSIGMA
          FMEAN=FMEAN2(NPEAKS,RPEAKS,3.0,FSIGMA)
          WRITE(*,100) '> Mean, r.m.s. (< 3 sigma): '
          WRITE(*,*) FMEAN,FSIGMA
        END IF
!------------------------------------------------------------------------------
100     FORMAT(A,$)
        END
