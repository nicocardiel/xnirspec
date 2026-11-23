!------------------------------------------------------------------------------
! SUBROUTINE WAVECAL_CLEAR(NBUFF)
!
! Remove the wavelength calibration of buffer NBUFF (LWAVECAL=.FALSE.,
! CRPIX1=1, CRVAL1=CDELT1=0, blank CTYPE1 and CUNIT1; COMMON blocks
! BLKWAVECAL1/2/3). To be used when a new image whose X axis does not
! correspond pixel by pixel to that of a calibrated image is stored in the
! buffer. With NBUFF=0, all the buffers are cleared.
!
! Parameters
! ----------
! NBUFF : INTEGER
!     Buffer (1 to NMAXBUFF), or 0 for all the buffers.
!------------------------------------------------------------------------------
        SUBROUTINE WAVECAL_CLEAR(NBUFF)
        IMPLICIT NONE
        INTEGER NBUFF
!
        INCLUDE 'dimensions.inc'
!
        INTEGER N1,N2
        REAL CRPIX1(NMAXBUFF),CRVAL1(NMAXBUFF),CDELT1(NMAXBUFF)
        LOGICAL LWAVECAL(NMAXBUFF)
        CHARACTER*20 CTYPE1(NMAXBUFF),CUNIT1(NMAXBUFF)
        COMMON/BLKWAVECAL1/CRPIX1,CRVAL1,CDELT1
        COMMON/BLKWAVECAL2/LWAVECAL
        COMMON/BLKWAVECAL3/CTYPE1,CUNIT1
!------------------------------------------------------------------------------
        IF(NBUFF.EQ.0)THEN
          N1=1
          N2=NMAXBUFF
        ELSE
          N1=NBUFF
          N2=NBUFF
        END IF
        LWAVECAL(N1:N2)=.FALSE.
        CRPIX1(N1:N2)=1.0
        CRVAL1(N1:N2)=0.0
        CDELT1(N1:N2)=0.0
        CTYPE1(N1:N2)=' '
        CUNIT1(N1:N2)=' '
        END
