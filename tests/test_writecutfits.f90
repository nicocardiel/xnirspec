!------------------------------------------------------------------------------
! Unit tests for WRITECUTFITS: a 1-D profile is written as a one-dimensional
! FITS image, with and without wavelength calibration, and read back with
! CFITSIO.
!------------------------------------------------------------------------------
        PROGRAM TEST_WRITECUTFITS
        USE TEST_UTILS
        IMPLICIT NONE
!
        INTEGER N
        PARAMETER (N=50)
        INTEGER IUNIT
        PARAMETER (IUNIT=90)                     !WRITECUTFITS uses unit 84
        CHARACTER*(*) FILE_WAVE,FILE_NOWAVE
        PARAMETER (FILE_WAVE='tmp_cut_wave.fits')
        PARAMETER (FILE_NOWAVE='tmp_cut_nowave.fits')
!
        INTEGER I,ISTATUS,NAXIS,NAXIS1,BITPIX,NKEYS,NMORE
        REAL Y(N),YREAD(N)
        REAL CRPIX1,CRVAL1,CDELT1
        CHARACTER*80 COMMENT,CARD
        CHARACTER*20 CTYPE1,CUNIT1
        LOGICAL ANYNULL,LHIST
!------------------------------------------------------------------------------
        DO I=1,N
          Y(I)=100.0+REAL(I)**1.5
        END DO
        CALL DELETE_FILE(FILE_WAVE)
        CALL DELETE_FILE(FILE_NOWAVE)
!------------------------------------------------------------------------------
! with wavelength calibration
        CALL WRITECUTFITS(FILE_WAVE,N,Y,.TRUE.,3.5,15234.567,-2.125, &
         'AWAV','Angstrom', &
         'xnirspec: cut along X of buffer #1, rows [10,12]',ISTATUS)
        CALL CHECK_INT('wave: ISTATUS',ISTATUS,0)
        ISTATUS=0
        CALL FTOPEN(IUNIT,FILE_WAVE,0,I,ISTATUS)
        CALL FTGKYJ(IUNIT,'NAXIS',NAXIS,COMMENT,ISTATUS)
        CALL FTGKYJ(IUNIT,'NAXIS1',NAXIS1,COMMENT,ISTATUS)
        CALL FTGKYJ(IUNIT,'BITPIX',BITPIX,COMMENT,ISTATUS)
        CALL FTGKYE(IUNIT,'CRPIX1',CRPIX1,COMMENT,ISTATUS)
        CALL FTGKYE(IUNIT,'CRVAL1',CRVAL1,COMMENT,ISTATUS)
        CALL FTGKYE(IUNIT,'CDELT1',CDELT1,COMMENT,ISTATUS)
        CALL FTGKYS(IUNIT,'CTYPE1',CTYPE1,COMMENT,ISTATUS)
        CALL FTGKYS(IUNIT,'CUNIT1',CUNIT1,COMMENT,ISTATUS)
        CALL FTGPVE(IUNIT,1,1,N,0.0,YREAD,ANYNULL,ISTATUS)
        CALL CHECK_INT('wave: CFITSIO status reading',ISTATUS,0)
        CALL CHECK_STR('wave: CTYPE1',CTYPE1,'AWAV')
        CALL CHECK_STR('wave: CUNIT1',CUNIT1,'Angstrom')
        CALL CHECK_INT('wave: NAXIS',NAXIS,1)
        CALL CHECK_INT('wave: NAXIS1',NAXIS1,N)
        CALL CHECK_INT('wave: BITPIX',BITPIX,-32)
        CALL CHECK_REAL('wave: CRPIX1',CRPIX1,3.5,0.0,0.0)
        CALL CHECK_REAL('wave: CRVAL1',CRVAL1,15234.567,0.0,0.0)
        CALL CHECK_REAL('wave: CDELT1',CDELT1,-2.125,0.0,0.0)
        CALL CHECK_TRUE('wave: data',ALL(YREAD.EQ.Y))
! HISTORY record
        CALL FTGHSP(IUNIT,NKEYS,NMORE,ISTATUS)
        LHIST=.FALSE.
        DO I=1,NKEYS
          CALL FTGREC(IUNIT,I,CARD,ISTATUS)
          IF(CARD(1:8).EQ.'HISTORY ')THEN
            LHIST=(INDEX(CARD,'cut along X of buffer #1, rows [10,12]') &
             .GT.0)
          END IF
        END DO
        CALL CHECK_TRUE('wave: HISTORY',LHIST)
        CALL FTCLOS(IUNIT,ISTATUS)
!------------------------------------------------------------------------------
! without wavelength calibration: no CRVAL1/CDELT1
        CALL WRITECUTFITS(FILE_NOWAVE,N,Y,.FALSE.,1.0,1.0,1.0,'WAVE',' ', &
         ' ',ISTATUS)
        CALL CHECK_INT('nowave: ISTATUS',ISTATUS,0)
        ISTATUS=0
        CALL FTOPEN(IUNIT,FILE_NOWAVE,0,I,ISTATUS)
        CALL FTGKYJ(IUNIT,'NAXIS1',NAXIS1,COMMENT,ISTATUS)
        CALL FTGPVE(IUNIT,1,1,N,0.0,YREAD,ANYNULL,ISTATUS)
        CALL CHECK_INT('nowave: CFITSIO status reading',ISTATUS,0)
        CALL CHECK_INT('nowave: NAXIS1',NAXIS1,N)
        CALL CHECK_TRUE('nowave: data',ALL(YREAD.EQ.Y))
        CALL FTGKYE(IUNIT,'CRVAL1',CRVAL1,COMMENT,ISTATUS)
        CALL CHECK_TRUE('nowave: no CRVAL1',ISTATUS.NE.0)
        ISTATUS=0
        CALL FTGKYE(IUNIT,'CDELT1',CDELT1,COMMENT,ISTATUS)
        CALL CHECK_TRUE('nowave: no CDELT1',ISTATUS.NE.0)
        ISTATUS=0
        CALL FTGKYS(IUNIT,'CTYPE1',CTYPE1,COMMENT,ISTATUS)
        CALL CHECK_TRUE('nowave: no CTYPE1',ISTATUS.NE.0)
        ISTATUS=0
        CALL FTCLOS(IUNIT,ISTATUS)
!------------------------------------------------------------------------------
! the file already exists: CFITSIO error, the file is not modified
        WRITE(*,'(A)') '(a CFITSIO error message is expected next)'
        CALL WRITECUTFITS(FILE_NOWAVE,N/2,Y,.FALSE.,1.0,1.0,1.0,' ',' ', &
         ' ',ISTATUS)
        CALL CHECK_TRUE('existing file: ISTATUS /= 0',ISTATUS.NE.0)
        ISTATUS=0
        CALL FTOPEN(IUNIT,FILE_NOWAVE,0,I,ISTATUS)
        CALL FTGKYJ(IUNIT,'NAXIS1',NAXIS1,COMMENT,ISTATUS)
        CALL FTCLOS(IUNIT,ISTATUS)
        CALL CHECK_INT('existing file: unchanged NAXIS1',NAXIS1,N)
!
        CALL DELETE_FILE(FILE_WAVE)
        CALL DELETE_FILE(FILE_NOWAVE)
        CALL TEST_SUMMARY()
        END PROGRAM TEST_WRITECUTFITS
