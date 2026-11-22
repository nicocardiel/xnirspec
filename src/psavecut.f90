!------------------------------------------------------------------------------
! SUBROUTINE PSAVECUT(LSAVED)
!
! Save the complete 1-D profile displayed in the upper-left plot (the last
! profile drawn by SUBPLOT, all its NLAST points, not only the plotted
! range) as a one-dimensional FITS image (see WRITECUTFITS). The output file
! name is asked with READC_FILE (wildcards list the matching files, "none"
! cancels, and existing files are not overwritten).
!
! The profile is identified as a cut of the buffer NB along X (or Y) when it
! has been drawn as the cuts of PCUT and of the main program: X label
! 'x axis' ('y axis'), color number = buffer number NB, NAXIS1 (NAXIS2)
! points and X(I)=I. If it is a cut along X and the buffer has wavelength
! calibration (LWAVECAL), CRPIX1, CRVAL1 and CDELT1 of the buffer are
! written in the header. A HISTORY record describes the origin of the
! profile. For other profiles only the Y values are saved.
!
! Returns
! -------
! LSAVED : LOGICAL
!     .TRUE. if the file has been written.
!------------------------------------------------------------------------------
        SUBROUTINE PSAVECUT(LSAVED)
        IMPLICIT NONE
        LOGICAL LSAVED
!
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
!
        INTEGER SYSTEMFUNCTION
        CHARACTER*255 READC_FILE
!
        INTEGER I,NB,ISTATUS,ISYSTEM
        CHARACTER*255 OUTFILE
        CHARACTER*255 CHIST
        LOGICAL LPIXEL,LXCUT,LYCUT,LWAVE,LEXIST,LASK
! global variables (COMMONs) shared with SUBPLOT and the main program
        INTEGER NLAST,N1LAST,N2LAST
        INTEGER ICOLORLAST,NSYMBLAST
        INTEGER NAXIS(2,NMAXBUFF)
        REAL XPLAST(NXYMAX*NOVERSAMPMAX),YPLAST(NXYMAX*NOVERSAMPMAX)
        REAL CRPIX1(NMAXBUFF),CRVAL1(NMAXBUFF),CDELT1(NMAXBUFF)
        CHARACTER*255 CLABX,CLABY,CLABG
        LOGICAL LWAVECAL(NMAXBUFF)
        COMMON/BLKPNLAST/NLAST,N1LAST,N2LAST
        COMMON/BLKPLAST/XPLAST,YPLAST
        COMMON/BLKPLABELS/CLABX,CLABY,CLABG
        COMMON/BLKPCOLSYMB/ICOLORLAST,NSYMBLAST
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKWAVECAL1/CRPIX1,CRVAL1,CDELT1
        COMMON/BLKWAVECAL2/LWAVECAL
!------------------------------------------------------------------------------
        LSAVED=.FALSE.
        IF(NLAST.LE.0)THEN
          WRITE(*,101) 'ERROR: there is no 1-D profile to be saved'
          RETURN
        END IF
!------------------------------------------------------------------------------
! identify the profile: X coordinates must be the pixel numbers 1,...,NLAST
        LPIXEL=.TRUE.
        DO I=1,NLAST
          IF(XPLAST(I).NE.REAL(I))THEN
            LPIXEL=.FALSE.
            EXIT
          END IF
        END DO
        NB=ICOLORLAST
        LXCUT=.FALSE.
        LYCUT=.FALSE.
        IF(LPIXEL.AND.(NB.GE.1).AND.(NB.LE.NMAXBUFF))THEN
          LXCUT=((CLABX.EQ.'x axis').AND.(NLAST.EQ.NAXIS(1,NB)))
          LYCUT=((CLABX.EQ.'y axis').AND.(NLAST.EQ.NAXIS(2,NB)))
        END IF
        LWAVE=.FALSE.
        IF(LXCUT) LWAVE=LWAVECAL(NB)
!
        WRITE(*,101) '* Saving the 1-D profile in the upper-left plot '// &
         '(all its points)'
        IF(LXCUT)THEN
          WRITE(CHIST,'(A,I0,A)') 'xnirspec: cut along X of buffer #', &
           NB,', rows '//TRIM(CLABG)
        ELSEIF(LYCUT)THEN
          WRITE(CHIST,'(A,I0,A)') 'xnirspec: cut along Y of buffer #', &
           NB,', columns '//TRIM(CLABG)
        ELSE
          CHIST='xnirspec: 1-D profile ('//TRIM(CLABY)//' versus '// &
           TRIM(CLABX)//')'
        END IF
        WRITE(*,101) '  '//TRIM(CHIST(11:))
        WRITE(*,'(A,I0)') '  number of pixels: ',NLAST
        IF(LWAVE)THEN
          WRITE(*,'(A,3(1X,ES15.7))') '  wavelength calibration '// &
           '(CRPIX1, CRVAL1, CDELT1):',CRPIX1(NB),CRVAL1(NB),CDELT1(NB)
        ELSE
          WRITE(*,101) '  no wavelength calibration'
        END IF
        IF(.NOT.LPIXEL)THEN
          WRITE(*,101) 'WARNING: the X coordinates of this profile are '// &
           'not the pixel numbers and they are not saved'
        END IF
!------------------------------------------------------------------------------
! output file name
        LASK=.TRUE.
        DO WHILE(LASK)
          OUTFILE=READC_FILE('Output FITS file name (none=EXIT)', &
           '*.fits','@')
          IF((INDEX(OUTFILE,'*').NE.0).OR.(INDEX(OUTFILE,'?').NE.0))THEN
            ISYSTEM=SYSTEMFUNCTION('ls '//TRIM(ADJUSTL(OUTFILE)))
          ELSEIF(OUTFILE.EQ.'none')THEN
            WRITE(*,101) '[CANCEL]'
            RETURN
          ELSE
            INQUIRE(FILE=OUTFILE,EXIST=LEXIST)
            IF(LEXIST)THEN
              WRITE(*,101) 'ERROR: this file already exists. Try again.'
            ELSE
              LASK=.FALSE.
            END IF
          END IF
        END DO
!------------------------------------------------------------------------------
! (the calibration is ignored when LWAVE is .FALSE.; NB may then be out of
! range, so buffer #1 is passed instead)
        IF(.NOT.LWAVE) NB=1
        CALL WRITECUTFITS(OUTFILE,NLAST,YPLAST,LWAVE,CRPIX1(NB), &
         CRVAL1(NB),CDELT1(NB),CHIST,ISTATUS)
        IF(ISTATUS.EQ.0)THEN
          LSAVED=.TRUE.
          WRITE(*,101) 'File '//TRIM(OUTFILE)//' written'
        ELSE
          WRITE(*,101) 'ERROR: the file could not be written'
        END IF
101     FORMAT(A)
        END
