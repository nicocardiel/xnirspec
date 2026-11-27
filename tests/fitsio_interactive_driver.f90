!------------------------------------------------------------------------------
! Driver for test_fitsio_interactive.sh: tests the paths of SLEEFITS and
! SESCRFITS that ask questions through stdin. The answers are provided by
! the script, in this order:
!   0        -> SLEEFITS, 3D cube: collapse the whole cube
!   2        -> SLEEFITS, 3D cube: read slice 2
!   -1, 2, 3 -> SLEEFITS, 3D cube: collapse slices 2 to 3
!   (none)   -> SLEEFITS, single HDU with EXTEND=T: no question is asked
!   -1, 2, 0, 1 -> SLEEFITS, file with extensions (numbered from 0): list
!               the HDUs, table HDU 2 rejected, empty primary HDU 0
!               rejected, read extension 1
!   <empty>  -> idem: the default is extension 1 (first image HDU)
!   -1, <empty> -> SLEEFITS, tile-compressed image: list the HDUs (name
!               COMPRESSED_IMAGE, since EXTNAME is missing), read
!               extension 1 (default)
!   <empty>  -> SESCRFITS with header template: skip 5 keywords (default)
!   <empty>  -> idem, template with wavelength calibration (removed, since
!               the buffer is not calibrated)
!   <empty>  -> idem, template with celestial WCS (kept)
!
! Test image: F(J,I) = 100*I + J + 0.25; slice K of the cube is
! F(J,I) + 1000*K (the collapsed cube is the sum of the slices).
!------------------------------------------------------------------------------
        PROGRAM FITSIO_INTERACTIVE_DRIVER
        USE TEST_UTILS
        USE Dynamic_Array_IMAGEN
        USE Dynamic_Array_IMAGEN_
        USE Dynamic_Array_LNULL
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
        INCLUDE 'interface_imagen_.inc'
        INCLUDE 'interface_lnull.inc'
        INCLUDE 'dimensions.inc'
!
        INTEGER NX1,NX2,NX3
        PARAMETER (NX1=4,NX2=3,NX3=3)
        INTEGER IUNIT
        PARAMETER (IUNIT=90)          !SLEEFITS and SESCRFITS use 80, 82, 83
        CHARACTER*(*) FILE_CUBE,FILE_MEF,FILE_HDR,FILE_OUT,FILE_EXT1
        CHARACTER*(*) FILE_CMP
        CHARACTER*(*) FILE_HDRW,FILE_OUTW
        PARAMETER (FILE_CUBE='tmp_fitsio_cube.fits')
        PARAMETER (FILE_MEF='tmp_fitsio_mef.fits')
        PARAMETER (FILE_EXT1='tmp_fitsio_extend1.fits')
        PARAMETER (FILE_CMP='tmp_fitsio_compressed.fits')
        PARAMETER (FILE_HDR='tmp_fitsio_hdrtemplate.fits')
        PARAMETER (FILE_OUT='tmp_fitsio_withhdr.fits')
        PARAMETER (FILE_HDRW='tmp_fitsio_hdrtemplate_wcs.fits')
        PARAMETER (FILE_OUTW='tmp_fitsio_withhdr_wcs.fits')
!
        INTEGER I,J,K
        INTEGER ISTATUS,NKEYS,NSPACE,IREADWRITE,IBLOCK
        INTEGER NKEYS_TEMPLATE,NKEYS_OWN,NNAXIS1,NCOMMENT
        INTEGER NAXES(3)
        INTEGER NROWS
        PARAMETER (NROWS=5)
        CHARACTER*8 TTYPE(2),TFORM(2),TUNIT(2)
        REAL F(NX1,NX2),CUBE(NX1,NX2,NX3)
        CHARACTER*20 COBJECT
        CHARACTER*50 COMMENT
        CHARACTER*80 RECORD
        CHARACTER*20 CVAL
        REAL RVAL
! global variables (COMMONs) shared with SLEEFITS and SESCRFITS
        INTEGER NAXIS(2,NMAXBUFF)
        LOGICAL LWAVECAL(NMAXBUFF)
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKWAVECAL2/LWAVECAL
!------------------------------------------------------------------------------
        CALL Initialize_Dynamic_Array_IMAGEN
        CALL Initialize_Dynamic_Array_IMAGEN_
        CALL Initialize_Dynamic_Array_LNULL
        LWAVECAL=.FALSE.
!
        DO I=1,NX2
          DO J=1,NX1
            F(J,I)=100.0*REAL(I)+REAL(J)+0.25
          END DO
        END DO
        DO K=1,NX3
          CUBE(:,:,K)=F+1000.0*REAL(K)
        END DO
!------------------------------------------------------------------------------
! 3D data cube
        ISTATUS=0
        NAXES=(/NX1,NX2,NX3/)
        CALL FTINIT(IUNIT,'!'//FILE_CUBE,1,ISTATUS)
        CALL FTPHPR(IUNIT,.TRUE.,-32,3,NAXES,0,1,.FALSE.,ISTATUS)
        CALL FTPPRE(IUNIT,1,1,NX1*NX2*NX3,CUBE,ISTATUS)
        CALL FTCLOS(IUNIT,ISTATUS)
        CALL CHECK_INT('CFITSIO status writing 3D cube',ISTATUS,0)
!
        CALL SLEEFITS(FILE_CUBE,.FALSE.,0,1,.FALSE.,0)            !answer: 0
        CALL CHECK_INT('cube NAXIS1',NAXIS(1,1),NX1)
        CALL CHECK_INT('cube NAXIS2',NAXIS(2,1),NX2)
        CALL CHECK_TRUE('cube collapsed (sum of all slices)', &
         ALL(IMAGEN(1:NX1,1:NX2,1).EQ.3.0*F+6000.0))
        CALL SLEEFITS(FILE_CUBE,.FALSE.,0,1,.FALSE.,0)            !answer: 2
        CALL CHECK_TRUE('cube single slice 2', &
         ALL(IMAGEN(1:NX1,1:NX2,1).EQ.F+2000.0))
        CALL SLEEFITS(FILE_CUBE,.FALSE.,0,1,.FALSE.,0)    !answers: -1, 2, 3
        CALL CHECK_TRUE('cube slices 2 to 3 collapsed (sum)', &
         ALL(IMAGEN(1:NX1,1:NX2,1).EQ.2.0*F+5000.0))
!------------------------------------------------------------------------------
! single HDU with EXTEND=T: the extension number is not asked
        ISTATUS=0
        NAXES(1:2)=(/NX1,NX2/)
        CALL FTINIT(IUNIT,'!'//FILE_EXT1,1,ISTATUS)
        CALL FTPHPR(IUNIT,.TRUE.,-32,2,NAXES,0,1,.TRUE.,ISTATUS)
        CALL FTPPRE(IUNIT,1,1,NX1*NX2,F,ISTATUS)
        CALL FTCLOS(IUNIT,ISTATUS)
        CALL CHECK_INT('CFITSIO status writing single HDU with EXTEND', &
         ISTATUS,0)
        CALL SLEEFITS(FILE_EXT1,.FALSE.,0,6,.FALSE.,0)           !no answer
        CALL CHECK_TRUE('single HDU with EXTEND=T: pixel values', &
         ALL(IMAGEN(1:NX1,1:NX2,6).EQ.F))
!------------------------------------------------------------------------------
! file with an empty primary HDU (0), one image extension (1) and one
! binary table (2)
        ISTATUS=0
        CALL FTINIT(IUNIT,'!'//FILE_MEF,1,ISTATUS)
        CALL FTPHPR(IUNIT,.TRUE.,-32,0,NAXES,0,1,.TRUE.,ISTATUS)
        CALL FTCRHD(IUNIT,ISTATUS)
        NAXES(1:2)=(/NX1,NX2/)
        CALL FTPHPR(IUNIT,.TRUE.,-32,2,NAXES,0,1,.FALSE.,ISTATUS)
        CALL FTPKYS(IUNIT,'EXTNAME','SCI','extension name',ISTATUS)
        CALL FTPPRE(IUNIT,1,1,NX1*NX2,F,ISTATUS)
        CALL FTCRHD(IUNIT,ISTATUS)
        TTYPE=(/'A','B'/)
        TFORM=(/'E','J'/)
        TUNIT=(/' ',' '/)
        CALL FTPHBN(IUNIT,NROWS,2,TTYPE,TFORM,TUNIT,'TAB',0,ISTATUS)
        CALL FTPCLE(IUNIT,1,1,1,NROWS,(/1.0,2.0,3.0,4.0,5.0/),ISTATUS)
        CALL FTPCLJ(IUNIT,2,1,1,NROWS,(/1,2,3,4,5/),ISTATUS)
        CALL FTCLOS(IUNIT,ISTATUS)
        CALL CHECK_INT('CFITSIO status writing file with extensions', &
         ISTATUS,0)
!
        CALL SLEEFITS(FILE_MEF,.FALSE.,0,2,.FALSE.,0)  !answers: -1, 2, 0, 1
        CALL CHECK_INT('extension NAXIS1',NAXIS(1,2),NX1)
        CALL CHECK_INT('extension NAXIS2',NAXIS(2,2),NX2)
        CALL CHECK_TRUE('extension pixel values', &
         ALL(IMAGEN(1:NX1,1:NX2,2).EQ.F))
        CALL SLEEFITS(FILE_MEF,.FALSE.,0,7,.FALSE.,0)      !answer: <empty>
        CALL CHECK_INT('default extension NAXIS1',NAXIS(1,7),NX1)
        CALL CHECK_TRUE('default extension pixel values', &
         ALL(IMAGEN(1:NX1,1:NX2,7).EQ.F))
!------------------------------------------------------------------------------
! tile-compressed float image (lossless: GZIP without quantization), written
! by CFITSIO as extension 1 after an empty primary HDU. The keywords BITPIX,
! NAXIS1, ... of this HDU describe the binary table (BITPIX=8, NAXIS1=8),
! not the image (BITPIX=-32, 4x3)
        ISTATUS=0
        NAXES(1:2)=(/NX1,NX2/)
        CALL FTINIT(IUNIT,'!'//FILE_CMP//'[compress GZIP; q 0]',1,ISTATUS)
        CALL FTPHPR(IUNIT,.TRUE.,-32,2,NAXES,0,1,.TRUE.,ISTATUS)
        CALL FTPPRE(IUNIT,1,1,NX1*NX2,F,ISTATUS)
        CALL FTCLOS(IUNIT,ISTATUS)
        CALL CHECK_INT('CFITSIO status writing compressed image', &
         ISTATUS,0)
        CALL SLEEFITS(FILE_CMP,.FALSE.,0,8,.FALSE.,0) !answers: -1, <empty>
        CALL CHECK_INT('compressed image NAXIS1',NAXIS(1,8),NX1)
        CALL CHECK_INT('compressed image NAXIS2',NAXIS(2,8),NX2)
        CALL CHECK_TRUE('compressed image pixel values', &
         ALL(IMAGEN(1:NX1,1:NX2,8).EQ.F))
!------------------------------------------------------------------------------
! SESCRFITS copying the header of a template file. The template starts
! with SIMPLE, BITPIX, NAXIS, NAXIS1 and NAXIS2, followed by the two COMMENT
! lines that CFITSIO adds after NAXIS2, and OBJECT. With the default NSKIP=5
! the remaining records are copied after the NKEYS_OWN records written by
! SESCRFITS itself (the same five keywords and two COMMENT lines), except
! those identical to records already present: only OBJECT must be copied
        ISTATUS=0
        NAXES(1:2)=(/2,2/)
        CALL FTINIT(IUNIT,'!'//FILE_HDR,1,ISTATUS)
        CALL FTPHPR(IUNIT,.TRUE.,-32,2,NAXES,0,1,.FALSE.,ISTATUS)
        CALL FTPKYS(IUNIT,'OBJECT','TESTOBJ','object name',ISTATUS)
        CALL FTPPRE(IUNIT,1,1,4,(/1.0,2.0,3.0,4.0/),ISTATUS)
        CALL FTCLOS(IUNIT,ISTATUS)
        CALL FTOPEN(IUNIT,FILE_HDR,0,IBLOCK,ISTATUS)
        CALL FTGHSP(IUNIT,NKEYS_TEMPLATE,NSPACE,ISTATUS)
        CALL FTCLOS(IUNIT,ISTATUS)
        CALL CHECK_INT('CFITSIO status writing header template',ISTATUS,0)
        NKEYS_OWN=NKEYS_TEMPLATE-1                    !all except OBJECT
!
        CALL DELETE_FILE(FILE_OUT)
        IMAGEN(1:NX1,1:NX2,3)=F
        NAXIS(1,3)=NX1
        NAXIS(2,3)=NX2
        CALL SESCRFITS(FILE_OUT,3,FILE_HDR)              !answer: <empty>
!
        ISTATUS=0
        IREADWRITE=0
        CALL FTOPEN(IUNIT,FILE_OUT,IREADWRITE,IBLOCK,ISTATUS)
        CALL FTGHSP(IUNIT,NKEYS,NSPACE,ISTATUS)
        NNAXIS1=0
        NCOMMENT=0
        DO I=1,NKEYS
          CALL FTGREC(IUNIT,I,RECORD,ISTATUS)
          IF(RECORD(1:8).EQ.'NAXIS1  ') NNAXIS1=NNAXIS1+1
          IF(RECORD(1:8).EQ.'COMMENT ') NCOMMENT=NCOMMENT+1
        END DO
        CALL FTGKYS(IUNIT,'OBJECT',COBJECT,COMMENT,ISTATUS)
        CALL FTCLOS(IUNIT,ISTATUS)
        CALL CHECK_INT('CFITSIO status reading output header',ISTATUS,0)
        CALL CHECK_STR('header template: OBJECT copied',COBJECT,'TESTOBJ')
        CALL CHECK_INT('header template: number of keywords',NKEYS, &
         NKEYS_OWN+1)
        CALL CHECK_INT('header template: NAXIS1 not duplicated',NNAXIS1,1)
        CALL CHECK_INT('header template: COMMENT lines not duplicated', &
         NCOMMENT,2)
        CALL SLEEFITS(FILE_OUT,.FALSE.,0,4,.FALSE.,0)
        CALL CHECK_TRUE('header template: pixel values', &
         ALL(IMAGEN(1:NX1,1:NX2,4).EQ.F))
!------------------------------------------------------------------------------
! header template with wavelength calibration and an uncalibrated buffer
! (#3): the calibration keywords copied from the template are removed, and
! the file is not read back as calibrated
        CALL WRITE_TEMPLATE('WAVE',21000.0,2.5)
        CALL DELETE_FILE(FILE_OUTW)
        CALL SESCRFITS(FILE_OUTW,3,FILE_HDRW)             !answer: <empty>
        ISTATUS=0
        CALL FTOPEN(IUNIT,FILE_OUTW,0,IBLOCK,ISTATUS)
        CALL FTGKYS(IUNIT,'OBJECT',COBJECT,COMMENT,ISTATUS)
        CALL CHECK_INT('wavelength template: CFITSIO status',ISTATUS,0)
        CALL CHECK_STR('wavelength template: OBJECT copied',COBJECT,'WAVEOBJ')
        CALL FTGKYS(IUNIT,'CTYPE1',CVAL,COMMENT,ISTATUS)
        CALL CHECK_TRUE('wavelength template: CTYPE1 removed',ISTATUS.NE.0)
        ISTATUS=0
        CALL FTGKYE(IUNIT,'CRVAL1',RVAL,COMMENT,ISTATUS)
        CALL CHECK_TRUE('wavelength template: CRVAL1 removed',ISTATUS.NE.0)
        ISTATUS=0
        CALL FTGKYE(IUNIT,'CD1_1',RVAL,COMMENT,ISTATUS)
        CALL CHECK_TRUE('wavelength template: CD1_1 removed',ISTATUS.NE.0)
        ISTATUS=0
        CALL FTCLOS(IUNIT,ISTATUS)
        CALL SLEEFITS(FILE_OUTW,.FALSE.,0,5,.FALSE.,0)
        CALL CHECK_TRUE('wavelength template: read back not calibrated', &
         .NOT.LWAVECAL(5))
! header template with celestial WCS: kept
        CALL WRITE_TEMPLATE('RA---TAN',150.0,-1.E-4)
        CALL DELETE_FILE(FILE_OUTW)
        CALL SESCRFITS(FILE_OUTW,3,FILE_HDRW)             !answer: <empty>
        ISTATUS=0
        CALL FTOPEN(IUNIT,FILE_OUTW,0,IBLOCK,ISTATUS)
        CALL FTGKYS(IUNIT,'CTYPE1',CVAL,COMMENT,ISTATUS)
        CALL FTGKYE(IUNIT,'CRVAL1',RVAL,COMMENT,ISTATUS)
        CALL FTCLOS(IUNIT,ISTATUS)
        CALL CHECK_INT('celestial template: CFITSIO status',ISTATUS,0)
        CALL CHECK_STR('celestial template: CTYPE1 kept',CVAL,'RA---TAN')
        CALL CHECK_REAL('celestial template: CRVAL1 kept',RVAL,150.0, &
         0.0,0.0)
!------------------------------------------------------------------------------
        CALL DELETE_FILE(FILE_CUBE)
        CALL DELETE_FILE(FILE_MEF)
        CALL DELETE_FILE(FILE_EXT1)
        CALL DELETE_FILE(FILE_CMP)
        CALL DELETE_FILE(FILE_HDR)
        CALL DELETE_FILE(FILE_OUT)
        CALL DELETE_FILE(FILE_HDRW)
        CALL DELETE_FILE(FILE_OUTW)
        CALL Deallocate_Array_IMAGEN
        CALL Deallocate_Array_IMAGEN_
        CALL Deallocate_Array_LNULL
        CALL TEST_SUMMARY()
!
        CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE WRITE_TEMPLATE(CTYPE,VAL,DELT)
!
! Write the header template FILE_HDRW (2x2 image) with OBJECT='WAVEOBJ' and
! the axis-1 keywords CTYPE1=CTYPE, CUNIT1, CRPIX1=1, CRVAL1=VAL,
! CDELT1=DELT and CD1_1=DELT.
!------------------------------------------------------------------------------
        SUBROUTINE WRITE_TEMPLATE(CTYPE,VAL,DELT)
        CHARACTER*(*) CTYPE
        REAL VAL,DELT
!
        ISTATUS=0
        NAXES(1:2)=(/2,2/)
        CALL FTINIT(IUNIT,'!'//FILE_HDRW,1,ISTATUS)
        CALL FTPHPR(IUNIT,.TRUE.,-32,2,NAXES,0,1,.FALSE.,ISTATUS)
        CALL FTPKYS(IUNIT,'OBJECT','WAVEOBJ','object name',ISTATUS)
        CALL FTPKYS(IUNIT,'CTYPE1',CTYPE,'axis type',ISTATUS)
        CALL FTPKYS(IUNIT,'CUNIT1','Angstrom','units',ISTATUS)
        CALL FTPKYE(IUNIT,'CRPIX1',1.0,7,'reference pixel',ISTATUS)
        CALL FTPKYE(IUNIT,'CRVAL1',VAL,7,'value at CRPIX1',ISTATUS)
        CALL FTPKYE(IUNIT,'CDELT1',DELT,7,'increment',ISTATUS)
        CALL FTPKYE(IUNIT,'CD1_1',DELT,7,'increment',ISTATUS)
        CALL FTPPRE(IUNIT,1,1,4,(/1.0,2.0,3.0,4.0/),ISTATUS)
        CALL FTCLOS(IUNIT,ISTATUS)
        CALL CHECK_INT('CFITSIO status writing template '//TRIM(CTYPE), &
         ISTATUS,0)
        END SUBROUTINE WRITE_TEMPLATE
        END PROGRAM FITSIO_INTERACTIVE_DRIVER
