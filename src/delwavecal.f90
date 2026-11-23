!------------------------------------------------------------------------------
! SUBROUTINE DELWAVECAL(IUNIT,LDELETED,ISTATUS)
!
! Remove the wavelength calibration keywords of axis 1 (CTYPE1, CUNIT1,
! CRPIX1, CRVAL1, CDELT1 and CD1_1) from the header of the current HDU of a
! FITS file opened with CFITSIO, but only when CTYPE1 is a linear wavelength
! axis (see WAVECTYPE). To be used when an image without wavelength
! calibration is saved with a header template that contains one (e.g. the
! header of the original spectrum before rectification), so that the file
! is not read back as wavelength calibrated. Other axis types (e.g.
! celestial coordinates, 'RA---TAN') are not modified.
!
! Parameters
! ----------
! IUNIT : INTEGER
!     CFITSIO unit number of the open file.
! ISTATUS : INTEGER
!     CFITSIO status (input and output; nothing is done if it is not zero
!     on input). Missing keywords are not errors.
!
! Returns
! -------
! LDELETED : LOGICAL
!     .TRUE. if the keywords have been removed.
!------------------------------------------------------------------------------
        SUBROUTINE DELWAVECAL(IUNIT,LDELETED,ISTATUS)
        IMPLICIT NONE
        INTEGER IUNIT
        LOGICAL LDELETED
        INTEGER ISTATUS
!
        LOGICAL WAVECTYPE
!
        INTEGER K,ISTATUS2
        CHARACTER*8 CKEY(6)
        CHARACTER*72 CTYPE1,COMMENT
        DATA CKEY /'CTYPE1','CUNIT1','CRPIX1','CRVAL1','CDELT1','CD1_1'/
!------------------------------------------------------------------------------
        LDELETED=.FALSE.
        IF(ISTATUS.NE.0) RETURN
! keywords are read and deleted with a separate status, since a missing
! keyword (KEY_NO_EXIST) is not an error here
        ISTATUS2=0
        CALL FTGKYS(IUNIT,'CTYPE1',CTYPE1,COMMENT,ISTATUS2)
        IF(ISTATUS2.NE.0) RETURN
        IF(.NOT.WAVECTYPE(CTYPE1)) RETURN
        DO K=1,6
          ISTATUS2=0
          CALL FTDKEY(IUNIT,TRIM(CKEY(K)),ISTATUS2)
        END DO
        LDELETED=.TRUE.
        END
