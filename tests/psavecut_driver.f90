!------------------------------------------------------------------------------
! Driver for test_psavecut.sh: draws a synthetic 1-D profile with SUBPLOT
! (as the cuts of PCUT, with PGPLOT opened with the /NULL device) and saves
! it with PSAVECUT, which reads the output file name from stdin (fallback
! path of READC_FILE). Only part of the profile (pixels 5 to N-5) is
! plotted, to check that the complete profile is saved.
!
! Test cases (first command-line argument):
!   xwave  : cut along X of buffer #2, which has wavelength calibration:
!            CRPIX1, CRVAL1 and CDELT1 must be written
!   ycut   : cut along Y of buffer #1 (which also has wavelength
!            calibration): no calibration keywords
!   none   : the answer is "none": nothing is saved
!   exists : the first answer is an existing file (created by the driver),
!            the second one a new file, which must be written
!
! The output file is tmp_psavecut.fits (deleted at the end).
!------------------------------------------------------------------------------
        PROGRAM PSAVECUT_DRIVER
        USE TEST_UTILS
        IMPLICIT NONE
        INCLUDE 'dimensions.inc'
        INCLUDE 'button.inc'
!
        INTEGER N
        PARAMETER (N=40)
        INTEGER IUNIT
        PARAMETER (IUNIT=90)
        CHARACTER*(*) OUTFILE,OLDFILE
        PARAMETER (OUTFILE='tmp_psavecut.fits')
        PARAMETER (OLDFILE='tmp_psavecut_old.fits')
!
        INTEGER I,ID,PGOPEN,ISTATUS,NAXIS1
        REAL X(N),Y(N),YREAD(N)
        REAL RVAL
        CHARACTER*20 CCASE
        CHARACTER*80 COMMENT
        LOGICAL LSAVED,ANYNULL,LEXIST
! global variables (COMMONs) shared with PSAVECUT and the input routines
        INTEGER NAXIS(2,NMAXBUFF)
        REAL CRPIX1(NMAXBUFF),CRVAL1(NMAXBUFF),CDELT1(NMAXBUFF)
        LOGICAL LWAVECAL(NMAXBUFF)
        LOGICAL LECHO
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKWAVECAL1/CRPIX1,CRVAL1,CDELT1
        COMMON/BLKWAVECAL2/LWAVECAL
        COMMON/BLKLECHO/LECHO
!------------------------------------------------------------------------------
        CALL GET_COMMAND_ARGUMENT(1,CCASE)
        WRITE(*,'(A)') 'Test case: '//TRIM(CCASE)
        LECHO=.FALSE.
        ID=PGOPEN('/NULL')
        IF(ID.LE.0)THEN
          WRITE(*,'(A)') 'FAIL: PGOPEN(''/NULL'')'
          ERROR STOP 1
        END IF
        MODOTEXT_BUTT=.TRUE.
        MODOTEXT_PLOTBUTT=.FALSE.
        CALL DELETE_FILE(OUTFILE)
        CALL DELETE_FILE(OLDFILE)
! buffers #1 (N rows) and #2 (N columns), both wavelength calibrated
        NAXIS=0
        NAXIS(:,1)=(/100,N/)
        NAXIS(:,2)=(/N,30/)
        LWAVECAL=.FALSE.
        LWAVECAL(1:2)=.TRUE.
        CRPIX1(1:2)=(/1.0,2.0/)
        CRVAL1(1:2)=(/4000.0,21500.25/)
        CDELT1(1:2)=(/1.0,-0.75/)
        DO I=1,N
          X(I)=REAL(I)
          Y(I)=50.0+REAL(I)+10.0*SIN(REAL(I))
        END DO
!------------------------------------------------------------------------------
        SELECT CASE(TRIM(CCASE))
        CASE('xwave','none','exists')
          CALL SUBPLOT(N,5,N-5,X,Y,X,Y,.TRUE.,.TRUE.,.FALSE.,.FALSE., &
           'x axis','signal','[10,12]',2,201,1.0)
        CASE('ycut')
          CALL SUBPLOT(N,5,N-5,X,Y,X,Y,.TRUE.,.TRUE.,.FALSE.,.FALSE., &
           'y axis','signal','[30,31]',1,201,1.0)
        CASE DEFAULT
          WRITE(*,'(A)') 'FAIL: unknown test case '//TRIM(CCASE)
          ERROR STOP 1
        END SELECT
        IF(TRIM(CCASE).EQ.'exists')THEN
          CALL WRITECUTFITS(OLDFILE,3,Y,.FALSE.,1.0,1.0,1.0,' ',ISTATUS)
        END IF
!------------------------------------------------------------------------------
        CALL PSAVECUT(LSAVED)
!------------------------------------------------------------------------------
        IF(TRIM(CCASE).EQ.'none')THEN
          CALL CHECK_TRUE('nothing saved',.NOT.LSAVED)
          INQUIRE(FILE=OUTFILE,EXIST=LEXIST)
          CALL CHECK_TRUE('no output file',.NOT.LEXIST)
        ELSE
          CALL CHECK_TRUE('profile saved',LSAVED)
          ISTATUS=0
          CALL FTOPEN(IUNIT,OUTFILE,0,I,ISTATUS)
          CALL FTGKYJ(IUNIT,'NAXIS1',NAXIS1,COMMENT,ISTATUS)
          CALL FTGPVE(IUNIT,1,1,N,0.0,YREAD,ANYNULL,ISTATUS)
          CALL CHECK_INT('CFITSIO status reading',ISTATUS,0)
          CALL CHECK_INT('complete profile saved (NAXIS1)',NAXIS1,N)
          CALL CHECK_TRUE('data',ALL(YREAD.EQ.Y))
          IF(TRIM(CCASE).EQ.'ycut')THEN
            CALL FTGKYE(IUNIT,'CRVAL1',RVAL,COMMENT,ISTATUS)
            CALL CHECK_TRUE('cut along Y: no CRVAL1',ISTATUS.NE.0)
            ISTATUS=0
          ELSE
            CALL FTGKYE(IUNIT,'CRPIX1',RVAL,COMMENT,ISTATUS)
            CALL CHECK_REAL('CRPIX1',RVAL,2.0,0.0,0.0)
            CALL FTGKYE(IUNIT,'CRVAL1',RVAL,COMMENT,ISTATUS)
            CALL CHECK_REAL('CRVAL1',RVAL,21500.25,0.0,0.0)
            CALL FTGKYE(IUNIT,'CDELT1',RVAL,COMMENT,ISTATUS)
            CALL CHECK_REAL('CDELT1',RVAL,-0.75,0.0,0.0)
            CALL CHECK_INT('CFITSIO status reading keywords',ISTATUS,0)
          END IF
          CALL FTCLOS(IUNIT,ISTATUS)
        END IF
        IF(TRIM(CCASE).EQ.'exists')THEN
          ISTATUS=0
          CALL FTOPEN(IUNIT,OLDFILE,0,I,ISTATUS)
          CALL FTGKYJ(IUNIT,'NAXIS1',NAXIS1,COMMENT,ISTATUS)
          CALL FTCLOS(IUNIT,ISTATUS)
          CALL CHECK_INT('existing file not overwritten',NAXIS1,3)
        END IF
        CALL DELETE_FILE(OUTFILE)
        CALL DELETE_FILE(OLDFILE)
        CALL PGCLOS
        CALL TEST_SUMMARY()
        END PROGRAM PSAVECUT_DRIVER
