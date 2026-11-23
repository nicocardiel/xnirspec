!------------------------------------------------------------------------------
! SUBROUTINE FITSWAVECAL(IUNIT,LWAVE,CRPIX1,CRVAL1,CDELT1,CTYPE1,CUNIT1)
!
! Read the linear wavelength calibration along NAXIS1 from the header of the
! current HDU of a FITS file already opened with CFITSIO. The image is
! considered to be wavelength calibrated (LWAVE=.TRUE.) only if:
!
! - CRVAL1 is present;
! - CDELT1 is present (or, if not, CD1_1) and is not zero;
! - CTYPE1 is a linear wavelength axis (see WAVECTYPE).
!
! CRPIX1 defaults to 1.0 if absent. If CRVAL1 and CDELT1 (or CD1_1) are
! present but the calibration is rejected, a warning is displayed. The
! CFITSIO status of the caller is not modified (missing keywords are not
! errors).
!
! Parameters
! ----------
! IUNIT : INTEGER
!     CFITSIO unit number of the open file.
!
! Returns
! -------
! LWAVE : LOGICAL
!     .TRUE. if the image is wavelength calibrated.
! CRPIX1, CRVAL1, CDELT1 : REAL
!     Calibration: the wavelength of pixel X (1-based) is
!     CRVAL1+(X-CRPIX1)*CDELT1. If LWAVE is .FALSE., CRPIX1=1 and
!     CRVAL1=CDELT1=0.
! CTYPE1, CUNIT1 : CHARACTER*(*)
!     Values of CTYPE1 and CUNIT1 (blank if absent).
!------------------------------------------------------------------------------
        SUBROUTINE FITSWAVECAL(IUNIT,LWAVE,CRPIX1,CRVAL1,CDELT1,CTYPE1, &
         CUNIT1)
        IMPLICIT NONE
        INTEGER IUNIT
        LOGICAL LWAVE
        REAL CRPIX1,CRVAL1,CDELT1
        CHARACTER*(*) CTYPE1,CUNIT1
!
        LOGICAL WAVECTYPE
!
        INTEGER ISTATUS
        CHARACTER*72 COMMENT
        CHARACTER*8 CKEYDELT
        LOGICAL LCRVAL,LCDELT
!------------------------------------------------------------------------------
        LWAVE=.FALSE.
! each keyword is read with its own status (a missing keyword sets it to
! KEY_NO_EXIST, which must not affect the following calls)
        ISTATUS=0
        CALL FTGKYS(IUNIT,'CTYPE1',CTYPE1,COMMENT,ISTATUS)
        IF(ISTATUS.NE.0) CTYPE1=' '
        ISTATUS=0
        CALL FTGKYS(IUNIT,'CUNIT1',CUNIT1,COMMENT,ISTATUS)
        IF(ISTATUS.NE.0) CUNIT1=' '
        ISTATUS=0
        CALL FTGKYE(IUNIT,'CRPIX1',CRPIX1,COMMENT,ISTATUS)
        IF(ISTATUS.NE.0) CRPIX1=1.0
        ISTATUS=0
        CALL FTGKYE(IUNIT,'CRVAL1',CRVAL1,COMMENT,ISTATUS)
        LCRVAL=(ISTATUS.EQ.0)
        ISTATUS=0
        CKEYDELT='CDELT1'
        CALL FTGKYE(IUNIT,'CDELT1',CDELT1,COMMENT,ISTATUS)
        IF(ISTATUS.NE.0)THEN
          ISTATUS=0
          CKEYDELT='CD1_1'
          CALL FTGKYE(IUNIT,'CD1_1',CDELT1,COMMENT,ISTATUS)
        END IF
        LCDELT=(ISTATUS.EQ.0)
!
        IF(LCRVAL.AND.LCDELT)THEN
          IF(CDELT1.EQ.0.0)THEN
            WRITE(*,'(A)') 'WARNING: wavelength calibration ignored: '// &
             TRIM(CKEYDELT)//' = 0'
          ELSEIF(LEN_TRIM(CTYPE1).EQ.0)THEN
            WRITE(*,'(A)') 'WARNING: wavelength calibration ignored: '// &
             'CTYPE1 not found (expected WAVE, AWAV, WAVELENGTH or LAMBDA)'
          ELSEIF(.NOT.WAVECTYPE(CTYPE1))THEN
            WRITE(*,'(A)') 'WARNING: wavelength calibration ignored: '// &
             'CTYPE1 = '''//TRIM(CTYPE1)//''' is not a linear '// &
             'wavelength axis (WAVE, AWAV, WAVELENGTH, LAMBDA)'
          ELSE
            LWAVE=.TRUE.
          END IF
        END IF
        IF(.NOT.LWAVE)THEN
          CRPIX1=1.0
          CRVAL1=0.0
          CDELT1=0.0
        END IF
        END
