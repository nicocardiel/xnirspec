!------------------------------------------------------------------------------
! Unit tests for SLEEFITS and SESCRFITS (paths that do not read from stdin).
!
! Test image: F(J,I) = 100*I + J + 0.25 (exactly representable as REAL*4),
! with J = 1..NX1 along NAXIS1 and I = 1..NX2 along NAXIS2 (non-square, to
! detect swapped axes). In the file with null values, pixel (JNULL,INULL)
! is NaN.
!
! The expected orientation for each IROTATE value is the convention
! implemented in SLEEFITS (regression test):
!     IROTATE=0          : IMAGEN(J,I)           = F(J,I)
!     IROTATE=90,-270    : IMAGEN(NX2-I+1,J)     = F(J,I)
!     IROTATE=180,-180   : IMAGEN(NX1-J+1,NX2-I+1) = F(J,I)
!     IROTATE=270,-90    : IMAGEN(I,NX1-J+1)     = F(J,I)
! and the null mask LNULL must follow the same transformation.
!------------------------------------------------------------------------------
        PROGRAM TEST_FITSIO
        USE TEST_UTILS
        USE Dynamic_Array_IMAGEN
        USE Dynamic_Array_IMAGEN_
        USE Dynamic_Array_LNULL
        USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE,IEEE_QUIET_NAN
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
        INCLUDE 'interface_imagen_.inc'
        INCLUDE 'interface_lnull.inc'
        INCLUDE 'dimensions.inc'
!
        INTEGER NX1,NX2
        PARAMETER (NX1=7,NX2=5)
        INTEGER JNULL,INULL
        PARAMETER (JNULL=2,INULL=1)
        INTEGER IUNIT
        PARAMETER (IUNIT=90)          !SLEEFITS and SESCRFITS use 80, 82, 83
        CHARACTER*(*) FILE_RT,FILE_NULL,FILE_I16,FILE_1D,FILE_WCS
        PARAMETER (FILE_RT='tmp_fitsio_roundtrip.fits')
        PARAMETER (FILE_NULL='tmp_fitsio_null.fits')
        PARAMETER (FILE_I16='tmp_fitsio_int16.fits')
        PARAMETER (FILE_1D='tmp_fitsio_1d.fits')
        PARAMETER (FILE_WCS='tmp_fitsio_wcs.fits')
!
        INTEGER I,J,K,NBAD
        INTEGER IROT(7)
        INTEGER ISTATUS
        INTEGER NAXES(2)
        INTEGER JDATA(4*3)
        REAL F(NX1,NX2),FNULL(NX1,NX2)
        REAL FDATA(6)
        CHARACTER*60 CNAME
! global variables (COMMONs) shared with SLEEFITS and SESCRFITS
        INTEGER NAXIS(2,NMAXBUFF)
        LOGICAL ANYNULL
        REAL CRPIX1(NMAXBUFF),CRVAL1(NMAXBUFF),CDELT1(NMAXBUFF)
        LOGICAL LWAVECAL(NMAXBUFF)
        REAL PHOTFLAM(NMAXBUFF),PHOTZPT(NMAXBUFF),EXPTIME(NMAXBUFF)
        LOGICAL L_PHOTFLAM(NMAXBUFF),L_PHOTZPT(NMAXBUFF)
        LOGICAL L_EXPTIME(NMAXBUFF)
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKANYNULL/ANYNULL
        COMMON/BLKWAVECAL1/CRPIX1,CRVAL1,CDELT1
        COMMON/BLKWAVECAL2/LWAVECAL
        COMMON/BLK_HSTFLUX1/L_PHOTFLAM,L_PHOTZPT,L_EXPTIME
        COMMON/BLK_HSTFLUX2/PHOTFLAM,PHOTZPT,EXPTIME
!------------------------------------------------------------------------------
        CALL Initialize_Dynamic_Array_IMAGEN
        CALL Initialize_Dynamic_Array_IMAGEN_
        CALL Initialize_Dynamic_Array_LNULL
!
        DO I=1,NX2
          DO J=1,NX1
            F(J,I)=100.0*REAL(I)+REAL(J)+0.25
          END DO
        END DO
        FNULL=F
        FNULL(JNULL,INULL)=IEEE_VALUE(1.0,IEEE_QUIET_NAN)
!------------------------------------------------------------------------------
! round trip: SESCRFITS (buffer 1) -> SLEEFITS (buffer 2)
        CALL DELETE_FILE(FILE_RT)
        IMAGEN(1:NX1,1:NX2,1)=F
        NAXIS(1,1)=NX1
        NAXIS(2,1)=NX2
        CALL SESCRFITS(FILE_RT,1,'none')
        CALL SLEEFITS(FILE_RT,.FALSE.,0,2,.FALSE.,0)
        CALL CHECK_INT('round trip NAXIS1',NAXIS(1,2),NX1)
        CALL CHECK_INT('round trip NAXIS2',NAXIS(2,2),NX2)
        CALL CHECK_TRUE('round trip pixel values', &
         ALL(IMAGEN(1:NX1,1:NX2,2).EQ.F))
        CALL CHECK_TRUE('round trip ANYNULL=.FALSE.',.NOT.ANYNULL)
        CALL CHECK_TRUE('round trip LNULL all .FALSE.', &
         .NOT.ANY(LNULL(1:NX1,1:NX2,2)))
        CALL CHECK_TRUE('round trip LWAVECAL=.FALSE. (no CRVAL1/CDELT1)', &
         .NOT.LWAVECAL(2))
!------------------------------------------------------------------------------
! null value (NaN) and image orientation
        CALL DELETE_FILE(FILE_NULL)
        IMAGEN(1:NX1,1:NX2,1)=FNULL
        CALL SESCRFITS(FILE_NULL,1,'none')
        IROT=(/0,90,180,270,-90,-180,-270/)
        DO K=1,7
          CALL SLEEFITS(FILE_NULL,.FALSE.,IROT(K),3,.FALSE.,0)
          WRITE(CNAME,'(A,I0)') 'IROTATE=',IROT(K)
          CALL CHECK_TRUE(TRIM(CNAME)//' ANYNULL=.TRUE.',ANYNULL)
          CALL CHECK_ORIENTATION(TRIM(CNAME),IROT(K),3)
        END DO
!------------------------------------------------------------------------------
! LBOX9=.TRUE.: the image is returned in IMAGEN_ and NAXIS is not modified
        NAXIS(1,4)=-1
        NAXIS(2,4)=-1
        IMAGEN_(1:NX1,1:NX2)=-1.0
        CALL SLEEFITS(FILE_RT,.FALSE.,0,4,.TRUE.,0)
        CALL CHECK_TRUE('LBOX9 image returned in IMAGEN_', &
         ALL(IMAGEN_(1:NX1,1:NX2).EQ.F))
        CALL CHECK_TRUE('LBOX9 NAXIS of buffer unchanged', &
         (NAXIS(1,4).EQ.-1).AND.(NAXIS(2,4).EQ.-1))
!------------------------------------------------------------------------------
! BITPIX=16 (integer reading path), 4x3 image with values 10*I+J
        DO I=1,3
          DO J=1,4
            JDATA((I-1)*4+J)=10*I+J
          END DO
        END DO
        ISTATUS=0
        NAXES=(/4,3/)
        CALL FTINIT(IUNIT,'!'//FILE_I16,1,ISTATUS)
        CALL FTPHPR(IUNIT,.TRUE.,16,2,NAXES,0,1,.FALSE.,ISTATUS)
        CALL FTPPRJ(IUNIT,1,1,12,JDATA,ISTATUS)
        CALL FTCLOS(IUNIT,ISTATUS)
        CALL CHECK_INT('CFITSIO status writing BITPIX=16 file',ISTATUS,0)
        CALL SLEEFITS(FILE_I16,.FALSE.,0,5,.FALSE.,0)
        CALL CHECK_INT('BITPIX=16 NAXIS1',NAXIS(1,5),4)
        CALL CHECK_INT('BITPIX=16 NAXIS2',NAXIS(2,5),3)
        NBAD=0
        DO I=1,3
          DO J=1,4
            IF(IMAGEN(J,I,5).NE.REAL(10*I+J)) NBAD=NBAD+1
          END DO
        END DO
        CALL CHECK_INT('BITPIX=16 mismatched pixels',NBAD,0)
!------------------------------------------------------------------------------
! NAXIS=1 (single spectrum)
        FDATA=(/1.5,2.5,3.5,4.5,5.5,6.5/)
        ISTATUS=0
        NAXES(1)=6
        CALL FTINIT(IUNIT,'!'//FILE_1D,1,ISTATUS)
        CALL FTPHPR(IUNIT,.TRUE.,-32,1,NAXES,0,1,.FALSE.,ISTATUS)
        CALL FTPPRE(IUNIT,1,1,6,FDATA,ISTATUS)
        CALL FTCLOS(IUNIT,ISTATUS)
        CALL CHECK_INT('CFITSIO status writing NAXIS=1 file',ISTATUS,0)
        CALL SLEEFITS(FILE_1D,.FALSE.,0,6,.FALSE.,0)
        CALL CHECK_INT('NAXIS=1 NAXIS1',NAXIS(1,6),6)
        CALL CHECK_INT('NAXIS=1 NAXIS2',NAXIS(2,6),1)
        CALL CHECK_TRUE('NAXIS=1 pixel values',ALL(IMAGEN(1:6,1,6).EQ.FDATA))
!------------------------------------------------------------------------------
! wavelength calibration keywords and HST photometry keywords (IBUFF<>0):
! EXPTIME present, PHOTFLAM and PHOTZPT absent
        ISTATUS=0
        NAXES=(/NX1,NX2/)
        CALL FTINIT(IUNIT,'!'//FILE_WCS,1,ISTATUS)
        CALL FTPHPR(IUNIT,.TRUE.,-32,2,NAXES,0,1,.FALSE.,ISTATUS)
        CALL FTPKYE(IUNIT,'CRPIX1',1.5,7,'reference pixel',ISTATUS)
        CALL FTPKYE(IUNIT,'CRVAL1',5000.0,7,'wavelength at CRPIX1',ISTATUS)
        CALL FTPKYE(IUNIT,'CDELT1',2.5,7,'dispersion',ISTATUS)
        CALL FTPKYE(IUNIT,'EXPTIME',300.0,7,'exposure time',ISTATUS)
        CALL FTPPRE(IUNIT,1,1,NX1*NX2,F,ISTATUS)
        CALL FTCLOS(IUNIT,ISTATUS)
        CALL CHECK_INT('CFITSIO status writing WCS file',ISTATUS,0)
        CALL SLEEFITS(FILE_WCS,.FALSE.,0,7,.FALSE.,7)
        CALL CHECK_TRUE('LWAVECAL=.TRUE.',LWAVECAL(7))
        CALL CHECK_REAL('CRPIX1',CRPIX1(7),1.5,1.E-6,0.0)
        CALL CHECK_REAL('CRVAL1',CRVAL1(7),5000.0,1.E-6,0.0)
        CALL CHECK_REAL('CDELT1',CDELT1(7),2.5,1.E-6,0.0)
        CALL CHECK_TRUE('L_EXPTIME=.TRUE.',L_EXPTIME(7))
        CALL CHECK_REAL('EXPTIME',EXPTIME(7),300.0,1.E-6,0.0)
        CALL CHECK_TRUE('L_PHOTFLAM=.FALSE.',.NOT.L_PHOTFLAM(7))
        CALL CHECK_TRUE('L_PHOTZPT=.FALSE.',.NOT.L_PHOTZPT(7))
        CALL CHECK_TRUE('WCS file pixel values', &
         ALL(IMAGEN(1:NX1,1:NX2,7).EQ.F))
!------------------------------------------------------------------------------
        CALL DELETE_FILE(FILE_RT)
        CALL DELETE_FILE(FILE_NULL)
        CALL DELETE_FILE(FILE_I16)
        CALL DELETE_FILE(FILE_1D)
        CALL DELETE_FILE(FILE_WCS)
        CALL Deallocate_Array_IMAGEN
        CALL Deallocate_Array_IMAGEN_
        CALL Deallocate_Array_LNULL
        CALL TEST_SUMMARY()
!
        CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE CHECK_ORIENTATION(LABEL,IROTATE,NBUFF)
!
! Check dimensions, pixel values and null mask of buffer NBUFF after reading
! FILE_NULL with SLEEFITS and the given IROTATE (see the convention in the
! header of this file). The null pixel must have value 0.0.
!
! Parameters
! ----------
! LABEL : CHARACTER*(*)
!     Prefix for the names of the checks.
! IROTATE : INTEGER
!     Rotation angle employed in SLEEFITS.
! NBUFF : INTEGER
!     Buffer where the image was read.
!------------------------------------------------------------------------------
        SUBROUTINE CHECK_ORIENTATION(LABEL,IROTATE,NBUFF)
        CHARACTER*(*), INTENT(IN) :: LABEL
        INTEGER, INTENT(IN) :: IROTATE,NBUFF
!
        INTEGER II,JJ,JD,ID,NBADV,NBADM,NEXP1,NEXP2
        REAL FEXP
        LOGICAL LEXP
!
        SELECT CASE(IROTATE)
        CASE(90,-270,270,-90)
          NEXP1=NX2
          NEXP2=NX1
        CASE DEFAULT
          NEXP1=NX1
          NEXP2=NX2
        END SELECT
        CALL CHECK_INT(LABEL//' NAXIS1',NAXIS(1,NBUFF),NEXP1)
        CALL CHECK_INT(LABEL//' NAXIS2',NAXIS(2,NBUFF),NEXP2)
        IF((NAXIS(1,NBUFF).NE.NEXP1).OR.(NAXIS(2,NBUFF).NE.NEXP2)) RETURN
!
        NBADV=0
        NBADM=0
        DO II=1,NX2
          DO JJ=1,NX1
            SELECT CASE(IROTATE)
            CASE(0)
              JD=JJ
              ID=II
            CASE(90,-270)
              JD=NX2-II+1
              ID=JJ
            CASE(180,-180)
              JD=NX1-JJ+1
              ID=NX2-II+1
            CASE DEFAULT                                          !270,-90
              JD=II
              ID=NX1-JJ+1
            END SELECT
            LEXP=((JJ.EQ.JNULL).AND.(II.EQ.INULL))
            IF(LEXP)THEN
              FEXP=0.0
            ELSE
              FEXP=F(JJ,II)
            END IF
            IF(IMAGEN(JD,ID,NBUFF).NE.FEXP) NBADV=NBADV+1
            IF(LNULL(JD,ID,NBUFF).NEQV.LEXP) NBADM=NBADM+1
          END DO
        END DO
        CALL CHECK_INT(LABEL//' mismatched pixel values',NBADV,0)
        CALL CHECK_INT(LABEL//' mismatched null mask pixels',NBADM,0)
        END SUBROUTINE CHECK_ORIENTATION
!
        END PROGRAM TEST_FITSIO
