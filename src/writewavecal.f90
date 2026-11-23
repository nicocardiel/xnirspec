!------------------------------------------------------------------------------
! SUBROUTINE WRITEWAVECAL(IUNIT,CRPIX1,CRVAL1,CDELT1,CTYPE1,CUNIT1,ISTATUS)
!
! Write the linear wavelength calibration along NAXIS1 in the header of the
! current HDU of a FITS file opened with CFITSIO: CTYPE1 ('WAVE' if blank),
! CUNIT1 (only if not blank), CRPIX1, CRVAL1 and CDELT1. Keywords already
! present in the header (e.g. copied from a header template) are updated
! instead of duplicated. If the header contains CD1_1, it is also set to
! CDELT1, since FITS readers give CD1_1 precedence over CDELT1.
!
! Parameters
! ----------
! IUNIT : INTEGER
!     CFITSIO unit number of the open file.
! CRPIX1, CRVAL1, CDELT1 : REAL
!     Calibration: the wavelength of pixel X (1-based) is
!     CRVAL1+(X-CRPIX1)*CDELT1.
! CTYPE1, CUNIT1 : CHARACTER*(*)
!     Type and units of the wavelength axis.
! ISTATUS : INTEGER
!     CFITSIO status (input and output; nothing is done if it is not zero
!     on input).
!------------------------------------------------------------------------------
        SUBROUTINE WRITEWAVECAL(IUNIT,CRPIX1,CRVAL1,CDELT1,CTYPE1,CUNIT1, &
         ISTATUS)
        IMPLICIT NONE
        INTEGER IUNIT
        REAL CRPIX1,CRVAL1,CDELT1
        CHARACTER*(*) CTYPE1,CUNIT1
        INTEGER ISTATUS
!
        INTEGER NDECIMALS
        PARAMETER (NDECIMALS=8) !E format with 9 significant digits (REAL*4)
!
        INTEGER ISTATUS2
        REAL RDUMMY
        CHARACTER*72 COMMENT
!------------------------------------------------------------------------------
        IF(ISTATUS.NE.0) RETURN
        IF(LEN_TRIM(CTYPE1).GT.0)THEN
          CALL FTUKYS(IUNIT,'CTYPE1',TRIM(ADJUSTL(CTYPE1)), &
           'wavelength axis',ISTATUS)
        ELSE
          CALL FTUKYS(IUNIT,'CTYPE1','WAVE','wavelength axis',ISTATUS)
        END IF
        IF(LEN_TRIM(CUNIT1).GT.0)THEN
          CALL FTUKYS(IUNIT,'CUNIT1',TRIM(ADJUSTL(CUNIT1)), &
           'wavelength units',ISTATUS)
        END IF
        CALL FTUKYE(IUNIT,'CRPIX1',CRPIX1,NDECIMALS,'reference pixel', &
         ISTATUS)
        CALL FTUKYE(IUNIT,'CRVAL1',CRVAL1,NDECIMALS, &
         'wavelength at the reference pixel',ISTATUS)
        CALL FTUKYE(IUNIT,'CDELT1',CDELT1,NDECIMALS, &
         'wavelength increment per pixel',ISTATUS)
! CD1_1 is only updated if it is already present (checked with a separate
! status, since a missing keyword is not an error here)
        ISTATUS2=0
        CALL FTGKYE(IUNIT,'CD1_1',RDUMMY,COMMENT,ISTATUS2)
        IF(ISTATUS2.EQ.0)THEN
          CALL FTUKYE(IUNIT,'CD1_1',CDELT1,NDECIMALS, &
           'wavelength increment per pixel',ISTATUS)
        END IF
        END
