!------------------------------------------------------------------------------
! SUBROUTINE WRITECUTFITS(FITSFILE,N,Y,LWAVE,CRPIX1,CRVAL1,CDELT1,CTYPE1,
!                         CUNIT1,CHIST,ISTATUS)
!
! Write a 1-D profile Y(1:N) as a one-dimensional FITS image (NAXIS=1,
! NAXIS1=N, BITPIX=-32). If LWAVE is .TRUE., the wavelength calibration
! keywords CTYPE1, CUNIT1 (if not blank), CRPIX1, CRVAL1 and CDELT1 are
! also written, so that the wavelength of pixel I (1-based) is
! CRVAL1+(I-CRPIX1)*CDELT1 (the same convention read by SLEEFITS, which
! requires CTYPE1 to be a linear wavelength axis; 'WAVE' is written if
! CTYPE1 is blank). The string CHIST is written as a HISTORY
! record (CFITSIO splits it into several records if needed). The file must
! not exist.
!
! Parameters
! ----------
! FITSFILE : CHARACTER*(*)
!     Name of the output FITS file.
! N : INTEGER
!     Number of pixels.
! Y : REAL array, shape (N)
!     Profile to be saved.
! LWAVE : LOGICAL
!     If .TRUE., write the wavelength calibration keywords.
! CRPIX1, CRVAL1, CDELT1 : REAL
!     Wavelength calibration (ignored if LWAVE is .FALSE.).
! CTYPE1, CUNIT1 : CHARACTER*(*)
!     Type and units of the wavelength axis (ignored if LWAVE is .FALSE.).
! CHIST : CHARACTER*(*)
!     Text of the HISTORY record (not written if blank).
!
! Returns
! -------
! ISTATUS : INTEGER
!     CFITSIO status (0 if the file has been written; otherwise the CFITSIO
!     error message is displayed).
!------------------------------------------------------------------------------
        SUBROUTINE WRITECUTFITS(FITSFILE,N,Y,LWAVE,CRPIX1,CRVAL1,CDELT1, &
         CTYPE1,CUNIT1,CHIST,ISTATUS)
        IMPLICIT NONE
        CHARACTER*(*) FITSFILE
        INTEGER N
        REAL Y(N)
        LOGICAL LWAVE
        REAL CRPIX1,CRVAL1,CDELT1
        CHARACTER*(*) CTYPE1,CUNIT1
        CHARACTER*(*) CHIST
        INTEGER ISTATUS
!
        INTEGER IUNIT
        PARAMETER (IUNIT=84)               !SLEEFITS/SESCRFITS use 80, 82, 83
        INTEGER NDECIMALS
        PARAMETER (NDECIMALS=8) !E format with 9 significant digits (REAL*4)
!
        INTEGER NAXES(1)
        INTEGER ISTATUS2
        CHARACTER*30 ERRTEXT
!------------------------------------------------------------------------------
        ISTATUS=0
        NAXES(1)=N
        CALL FTINIT(IUNIT,FITSFILE,1,ISTATUS)
        CALL FTPHPR(IUNIT,.TRUE.,-32,1,NAXES,0,1,.FALSE.,ISTATUS)
        IF(LWAVE)THEN
          IF(LEN_TRIM(CTYPE1).GT.0)THEN
            CALL FTPKYS(IUNIT,'CTYPE1',TRIM(ADJUSTL(CTYPE1)), &
             'wavelength axis',ISTATUS)
          ELSE
            CALL FTPKYS(IUNIT,'CTYPE1','WAVE','wavelength axis',ISTATUS)
          END IF
          IF(LEN_TRIM(CUNIT1).GT.0)THEN
            CALL FTPKYS(IUNIT,'CUNIT1',TRIM(ADJUSTL(CUNIT1)), &
             'wavelength units',ISTATUS)
          END IF
          CALL FTPKYE(IUNIT,'CRPIX1',CRPIX1,NDECIMALS, &
           'reference pixel',ISTATUS)
          CALL FTPKYE(IUNIT,'CRVAL1',CRVAL1,NDECIMALS, &
           'wavelength at the reference pixel',ISTATUS)
          CALL FTPKYE(IUNIT,'CDELT1',CDELT1,NDECIMALS, &
           'wavelength increment per pixel',ISTATUS)
        END IF
        IF(LEN_TRIM(CHIST).GT.0) CALL FTPHIS(IUNIT,TRIM(CHIST),ISTATUS)
        CALL FTPPRE(IUNIT,1,1,N,Y,ISTATUS)
! close the file even after an error (with a separate status, so that the
! first error is the one reported)
        ISTATUS2=0
        CALL FTCLOS(IUNIT,ISTATUS2)
        IF(ISTATUS.EQ.0) ISTATUS=ISTATUS2
        IF(ISTATUS.NE.0)THEN
          CALL FTGERR(ISTATUS,ERRTEXT)
          WRITE(*,'(A,I0,A)') 'CFITSIO ERROR (status=',ISTATUS,'): '// &
           TRIM(ERRTEXT)
        END IF
        END
