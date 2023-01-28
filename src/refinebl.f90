! Refina la posicion de una linea BL del boundary.
! Nota: el polinomio de entrada en la subrutina es modificados a la salida.
! XMIN,XMAX: es el recorrido de la variable X que va a ser examinado.
        SUBROUTINE REFINEBL(NDEG,COEFF,XMIN,XMAX)
        USE Dynamic_Array_IMAGEN
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
! subroutine arguments
        INTEGER NDEG
        REAL COEFF(NDEG+1)
        REAL XMIN,XMAX
!
        INCLUDE 'largest.inc'
        INCLUDE 'dimensions.inc'
!
        REAL FPOLY
!
        INTEGER I,J,I1,I2,J1,J2
        INTEGER NWIDTH,NMED
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER NPEAKS
        INTEGER NCBUFF
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL TSIGMA
        REAL XPEAKS(NXYMAX),YPEAKS(NXYMAX),RPEAKS(NXYMAX)
        REAL XF(NXYMAX),YF(NXYMAX)
        REAL XX,YY
        REAL COEFFPAR(3),CHISQR
        REAL XMIN_,XMAX_,YMIN_,YMAX_
        LOGICAL LDEBUGLOCAL
        LOGICAL LFIT(NXYMAX)
!
!delete COMMON/BLKIMAGEN1/IMAGEN             !imagen FITS leida en formato REAL
        COMMON/BLKIMAGEN2/NCBUFF
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKDEFAULTS2/NWIDTH
        COMMON/BLKDEFAULTS3/TSIGMA
        COMMON/BLKPLIMITS/XMIN_,XMAX_,YMIN_,YMAX_
!------------------------------------------------------------------------------
        LDEBUGLOCAL=.TRUE.
!------------------------------------------------------------------------------
        NMED=NWIDTH/2
        J1=NINT(XMIN)
        J2=NINT(XMAX)
        DO J=J1,J2
          XX=REAL(J)
          YY=FPOLY(NDEG,COEFF,XX)
          I1=NINT(YY)-NMED
          I2=NINT(YY)+NMED
          IF(I1.LT.1)THEN
            I1=1
            I2=I1+NWIDTH-1
          ELSEIF(I2.GT.NAXIS(2,NCBUFF))THEN
            I2=NAXIS(2,NCBUFF)
            I1=I2-NWIDTH+1
          END IF
          DO I=I1,I2
            XF(I-I1+1)=REAL(I)
            YF(I-I1+1)=IMAGEN(J,I,NCBUFF)
          END DO
          CALL POLFIT(XF,YF,YF,NWIDTH,3,0,COEFFPAR,CHISQR)
          XPEAKS(J-J1+1)=REAL(J)
          YPEAKS(J-J1+1)=-COEFFPAR(2)/2./COEFFPAR(3)
        END DO
        NPEAKS=J2-J1+1
        CALL POLFITSIG(NPEAKS,XPEAKS,YPEAKS,TSIGMA,NDEG,COEFF,LFIT)
! en lugar de dibujar el ajuste directamente es mejor mostrar los residuos
        DO J=1,NPEAKS
          RPEAKS(J)=YPEAKS(J)-FPOLY(NDEG,COEFF,XPEAKS(J))
        END DO
        IF(LDEBUGLOCAL)THEN
          CALL PGSCI(2)
          CALL PGPOINT(NPEAKS,XPEAKS,YPEAKS,1)
          CALL PGSCI(1)
          YMIN_=-1.
          YMAX_=+1.
          CALL SUBPLOT(NPEAKS,1,NPEAKS,XPEAKS,RPEAKS,XPEAKS,RPEAKS,.TRUE.,.FALSE.,.FALSE.,.FALSE., &
           'wavelength direction','residuals in spatial direction',' ',3,1,1.5)
          J=0
          DO I=1,NPEAKS
            IF(.NOT.LFIT(I))THEN
              J=J+1
              XF(J)=XPEAKS(I)
              YF(J)=RPEAKS(I)
            END IF
          END DO
          IF(J.GT.0) CALL SUBPLOTBIS(J,1,J,XF,YF,XF,YF,.FALSE.,.FALSE.,2,5,1.0)
          DO I=1,NXYMAX
            XF(I)=XMIN+REAL(I-1)/REAL(NXYMAX-1)*(XMAX-XMIN)
            YF(I)=FPOLY(NDEG,COEFF,XF(I))
          END DO
          CALL PGSCI(3)
          CALL PGLINE(NXYMAX,XF,YF)
          CALL PGSCI(1)
          CALL SUBPLOTBIS(NXYMAX,1,NXYMAX,XF,YF,XF,YF,.FALSE.,.FALSE.,3,101,1.0)
        END IF
!------------------------------------------------------------------------------
        END
