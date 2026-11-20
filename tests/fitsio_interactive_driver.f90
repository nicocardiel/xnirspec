!------------------------------------------------------------------------------
! Driver for test_fitsio_interactive.sh: tests the paths of SLEEFITS and
! SESCRFITS that ask questions through stdin. The answers are provided by
! the script, in this order:
!   0        -> SLEEFITS, 3D cube: collapse the whole cube
!   2        -> SLEEFITS, 3D cube: read slice 2
!   -1, 2, 3 -> SLEEFITS, 3D cube: collapse slices 2 to 3
!   2        -> SLEEFITS, file with extensions: read extension number 2
!   <empty>  -> SESCRFITS with header template: skip 5 keywords (default)
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
        CHARACTER*(*) FILE_CUBE,FILE_MEF,FILE_HDR,FILE_OUT
        PARAMETER (FILE_CUBE='tmp_fitsio_cube.fits')
        PARAMETER (FILE_MEF='tmp_fitsio_mef.fits')
        PARAMETER (FILE_HDR='tmp_fitsio_hdrtemplate.fits')
        PARAMETER (FILE_OUT='tmp_fitsio_withhdr.fits')
!
        INTEGER I,J,K
        INTEGER ISTATUS,NKEYS,NSPACE,IREADWRITE,IBLOCK
        INTEGER NKEYS_TEMPLATE,NKEYS_OWN,NNAXIS1,NCOMMENT
        INTEGER NAXES(3)
        REAL F(NX1,NX2),CUBE(NX1,NX2,NX3)
        CHARACTER*20 COBJECT
        CHARACTER*50 COMMENT
        CHARACTER*80 RECORD
! global variables (COMMONs) shared with SLEEFITS and SESCRFITS
        INTEGER NAXIS(2,NMAXBUFF)
        COMMON/BLKNAXIS/NAXIS
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
! file with an empty primary HDU and one image extension
        ISTATUS=0
        CALL FTINIT(IUNIT,'!'//FILE_MEF,1,ISTATUS)
        CALL FTPHPR(IUNIT,.TRUE.,-32,0,NAXES,0,1,.TRUE.,ISTATUS)
        CALL FTCRHD(IUNIT,ISTATUS)
        NAXES(1:2)=(/NX1,NX2/)
        CALL FTPHPR(IUNIT,.TRUE.,-32,2,NAXES,0,1,.FALSE.,ISTATUS)
        CALL FTPKYS(IUNIT,'EXTNAME','SCI','extension name',ISTATUS)
        CALL FTPPRE(IUNIT,1,1,NX1*NX2,F,ISTATUS)
        CALL FTCLOS(IUNIT,ISTATUS)
        CALL CHECK_INT('CFITSIO status writing file with extension', &
         ISTATUS,0)
!
        CALL SLEEFITS(FILE_MEF,.FALSE.,0,2,.FALSE.,0)             !answer: 2
        CALL CHECK_INT('extension NAXIS1',NAXIS(1,2),NX1)
        CALL CHECK_INT('extension NAXIS2',NAXIS(2,2),NX2)
        CALL CHECK_TRUE('extension pixel values', &
         ALL(IMAGEN(1:NX1,1:NX2,2).EQ.F))
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
        CALL DELETE_FILE(FILE_CUBE)
        CALL DELETE_FILE(FILE_MEF)
        CALL DELETE_FILE(FILE_HDR)
        CALL DELETE_FILE(FILE_OUT)
        CALL Deallocate_Array_IMAGEN
        CALL Deallocate_Array_IMAGEN_
        CALL Deallocate_Array_LNULL
        CALL TEST_SUMMARY()
        END PROGRAM FITSIO_INTERACTIVE_DRIVER
