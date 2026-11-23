!------------------------------------------------------------------------------
! LOGICAL FUNCTION WAVECTYPE(CTYPE)
!
! Return .TRUE. if CTYPE (value of the FITS keyword CTYPE1) corresponds to a
! linear wavelength axis. Leading and trailing blanks and the case are
! ignored. Accepted values:
!
! - 'WAVE' and 'AWAV': vacuum and air wavelength (FITS standard, WCS Paper
!   III), without algorithm code (i.e. linear);
! - 'WAVELENGTH' and 'LAMBDA': non-standard but common in older data and in
!   IRAF.
!
! Non-linear spectral axes (e.g. 'WAVE-LOG', 'WAVE-TAB') and 'LINEAR' (used
! by IRAF both for wavelength-calibrated spectra and for pixel axes) are not
! accepted.
!
! Parameters
! ----------
! CTYPE : CHARACTER*(*)
!     Value of CTYPE1.
!
! Returns
! -------
! WAVECTYPE : LOGICAL
!     .TRUE. if CTYPE is one of the accepted values.
!------------------------------------------------------------------------------
        LOGICAL FUNCTION WAVECTYPE(CTYPE)
        IMPLICIT NONE
        CHARACTER*(*) CTYPE
!
        CHARACTER*80 CLOCAL
!------------------------------------------------------------------------------
        CLOCAL=ADJUSTL(CTYPE)
        CALL CHUPPER(CLOCAL)
        SELECT CASE(TRIM(CLOCAL))
        CASE('WAVE','AWAV','WAVELENGTH','LAMBDA')
          WAVECTYPE=.TRUE.
        CASE DEFAULT
          WAVECTYPE=.FALSE.
        END SELECT
        END
