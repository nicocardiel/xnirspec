!------------------------------------------------------------------------------
! SUBROUTINE WVSCALEBOX(CRPIX1,CRVAL1,CDELT1)
!
! Draw a wavelength scale along the top edge of the current viewport, whose
! world coordinates must be pixels along X (as in SUBLOOK). The wavelength
! of pixel X is CRVAL1+(X-CRPIX1)*CDELT1 (CDELT1 can be negative). Tick
! marks are drawn inwards, and the numeric labels above the top edge, with
! character height 0.8; labels that would extend beyond the left or right
! edge of the viewport are not written (they would overlap with the labels
! of the Y axis). The world coordinates, character height and color are
! restored on exit.
!
! Parameters
! ----------
! CRPIX1, CRVAL1, CDELT1 : REAL
!     Wavelength calibration (FITS convention, 1-based pixels).
!------------------------------------------------------------------------------
        SUBROUTINE WVSCALEBOX(CRPIX1,CRVAL1,CDELT1)
        IMPLICIT NONE
        REAL CRPIX1,CRVAL1,CDELT1
!
        REAL CHLABEL
        PARAMETER (CHLABEL=0.8)        !character height of the labels
        REAL DISPLABEL
        PARAMETER (DISPLABEL=0.7)     !displacement of the labels (as PGBOX)
!
        REAL PGRND
!
        INTEGER K,K1,K2,NSUB,NC,MM,IPOW
        REAL X1,X2,Y1,Y2,W1,W2
        REAL XV1,XV2,YV1,YV2
        REAL WTICK,W,FRAC,XL,YL,HALFW
        REAL OLD_CH
        INTEGER OLD_CI
        CHARACTER*20 CLABEL
!------------------------------------------------------------------------------
        IF(CDELT1.EQ.0.0) RETURN
        CALL PGQWIN(X1,X2,Y1,Y2)
        CALL PGQVP(0,XV1,XV2,YV1,YV2)                     !NDC coordinates
        CALL PGQCH(OLD_CH)
        CALL PGQCI(OLD_CI)
        W1=CRVAL1+(X1-CRPIX1)*CDELT1
        W2=CRVAL1+(X2-CRPIX1)*CDELT1
! major tick interval: a "round" number (1, 2 or 5 x 10**n) close to 1/5 of
! the wavelength range, as chosen by PGBOX by default
        WTICK=PGRND(0.2*ABS(W2-W1),NSUB)
        IF(WTICK.LE.0.0) RETURN
        CALL PGSWIN(W1,W2,Y1,Y2)
        CALL PGBOX('CTS',WTICK,NSUB,' ',0.0,0)
! numeric labels, written as MM*10**IPOW with PGNUMB (IPOW is the power of
! ten of WTICK, so MM is an integer)
        CALL PGSCH(CHLABEL)
        IPOW=FLOOR(LOG10(WTICK)+1.E-5)
        K1=CEILING(MIN(W1,W2)/WTICK-1.E-5)
        K2=FLOOR(MAX(W1,W2)/WTICK+1.E-5)
        DO K=K1,K2
          W=REAL(K)*WTICK
          MM=NINT(W/10.0**IPOW)
          CALL PGNUMB(MM,IPOW,1,CLABEL,NC)
          FRAC=(W-W1)/(W2-W1)                !fraction along the top edge
          CALL PGLEN(0,CLABEL(1:NC),XL,YL)                 !length in NDC
          HALFW=0.5*XL/(XV2-XV1)
          IF((FRAC-HALFW.GE.0.0).AND.(FRAC+HALFW.LE.1.0))THEN
            CALL PGMTXT('T',DISPLABEL,FRAC,0.5,CLABEL(1:NC))
          END IF
        END DO
! restore previous settings
        CALL PGSWIN(X1,X2,Y1,Y2)
        CALL PGSCH(OLD_CH)
        CALL PGSCI(OLD_CI)
        END
