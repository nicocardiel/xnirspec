!------------------------------------------------------------------------------
! SUBROUTINE WAVECAL_COPY(NFROM,NTO)
!
! Copy the wavelength calibration of buffer NFROM (LWAVECAL, CRPIX1, CRVAL1,
! CDELT1, CTYPE1 and CUNIT1, in the COMMON blocks BLKWAVECAL1/2/3) to buffer
! NTO. To be used when a new image is stored in buffer NTO with the same X
! axis (pixel by pixel) as the image in buffer NFROM.
!
! Parameters
! ----------
! NFROM, NTO : INTEGER
!     Source and destination buffers (1 to NMAXBUFF).
!------------------------------------------------------------------------------
        SUBROUTINE WAVECAL_COPY(NFROM,NTO)
        IMPLICIT NONE
        INTEGER NFROM,NTO
!
        INCLUDE 'dimensions.inc'
!
        REAL CRPIX1(NMAXBUFF),CRVAL1(NMAXBUFF),CDELT1(NMAXBUFF)
        LOGICAL LWAVECAL(NMAXBUFF)
        CHARACTER*20 CTYPE1(NMAXBUFF),CUNIT1(NMAXBUFF)
        COMMON/BLKWAVECAL1/CRPIX1,CRVAL1,CDELT1
        COMMON/BLKWAVECAL2/LWAVECAL
        COMMON/BLKWAVECAL3/CTYPE1,CUNIT1
!------------------------------------------------------------------------------
        IF(NFROM.EQ.NTO) RETURN
        LWAVECAL(NTO)=LWAVECAL(NFROM)
        CRPIX1(NTO)=CRPIX1(NFROM)
        CRVAL1(NTO)=CRVAL1(NFROM)
        CDELT1(NTO)=CDELT1(NFROM)
        CTYPE1(NTO)=CTYPE1(NFROM)
        CUNIT1(NTO)=CUNIT1(NFROM)
        END
