!
!******************************************************************************
! SUBROUTINE FITSINFO(IUNIT,FITSFILE,NHDU)
!
! Print a summary of the HDUs of a FITS file already opened with CFITSIO,
! in the same layout as astropy (fitsinfo, HDUList.info()):
!
!   Filename: file.fits
!   No.    Name      Ver    Type      Cards   Dimensions   Format
!     0  PRIMARY       1 PrimaryHDU       6   ()
!     1  SCI           1 ImageHDU         9   (4, 3)   float32
!     2  TAB           1 BinTableHDU     13   5R x 2C   [E, J]
!
! Parameters
! ----------
! IUNIT : INTEGER
!     CFITSIO unit number of the opened file.
! FITSFILE : CHARACTER*(*)
!     File name (only used to print it).
! NHDU : INTEGER
!     Number of HDUs in the file (as returned by FTTHDU).
!
! Notes
! -----
! HDUs are numbered from 0 (primary HDU), as in astropy, whereas CFITSIO
! numbers them from 1. Image dimensions are given in the astropy order
! (NAXIS1, NAXIS2, ...), and the image format is the stored data type
! (BITPIX), followed by the equivalent type returned by FTGIET when BSCALE
! and BZERO change it. For tile-compressed images, Cards is the number of
! records of the underlying binary table (astropy gives that of the image
! header), since this is what CFITSIO (FTGHSP) returns. The file is
! left positioned at an arbitrary HDU: the caller must move to the desired
! HDU (FTMAHD) afterwards. CFITSIO errors are printed (PRINTERROR) but do
! not modify the status of the caller.
!------------------------------------------------------------------------------
        SUBROUTINE FITSINFO(IUNIT,FITSFILE,NHDU)
        IMPLICIT NONE
! subroutine arguments
        INTEGER IUNIT
        CHARACTER*(*) FITSFILE
        INTEGER NHDU
! local parameters
        INTEGER MAXDIM                      !maximum number of axes displayed
        PARAMETER (MAXDIM=99)
! local variables
        INTEGER IHDU,K,L
        INTEGER ISTATUS,HDUTYPE
        INTEGER NKEYS,NSPACE
        INTEGER EXTVER
        INTEGER NDIM,NAXES(MAXDIM)
        INTEGER BITPIX,IMGTYPE
        INTEGER NROWS,NCOLS
        CHARACTER*3 CHDU,CVER
        CHARACTER*5 CCARDS
        CHARACTER*8 KEYNAME
        CHARACTER*12 CTYPE
        CHARACTER*12 CNUM,CNUM2
        CHARACTER*20 EXTNAME
        CHARACTER*20 CTFORM
        CHARACTER*50 COMMENT
        CHARACTER*255 CDIM,CFORMAT,CLINE
        LOGICAL ZIMAGE
!------------------------------------------------------------------------------
        WRITE(*,'(A)') 'Filename: '//TRIM(FITSFILE)
        WRITE(*,'(A)') 'No.    Name      Ver    Type      Cards   '// &
         'Dimensions   Format'
        DO IHDU=1,NHDU
          ISTATUS=0
          NKEYS=0
          NDIM=0
          NROWS=0
          NCOLS=0
          CALL FTMAHD(IUNIT,IHDU,HDUTYPE,ISTATUS)
          IF(ISTATUS.GT.0)THEN
            CALL PRINTERROR(ISTATUS)
            RETURN
          END IF
          CALL FTGHSP(IUNIT,NKEYS,NSPACE,ISTATUS)
! EXTNAME and EXTVER are optional (status 202: keyword not found)
          CALL FTGKYS(IUNIT,'EXTNAME',EXTNAME,COMMENT,ISTATUS)
          IF(ISTATUS.EQ.202)THEN
            ISTATUS=0
            IF(IHDU.EQ.1)THEN
              EXTNAME='PRIMARY'
            ELSE
              EXTNAME=' '
            END IF
          END IF
          CALL FTGKYJ(IUNIT,'EXTVER',EXTVER,COMMENT,ISTATUS)
          IF(ISTATUS.EQ.202)THEN
            ISTATUS=0
            EXTVER=1
          END IF
          CDIM=' '
          CFORMAT=' '
          IF(HDUTYPE.EQ.0)THEN
! image HDU (CFITSIO also presents tile-compressed images as images; they
! are identified by ZIMAGE=T in the header of the binary table)
            CALL FTGKYL(IUNIT,'ZIMAGE',ZIMAGE,COMMENT,ISTATUS)
            IF(ISTATUS.EQ.202)THEN
              ISTATUS=0
              ZIMAGE=.FALSE.
            END IF
            IF(ZIMAGE)THEN
              CTYPE='CompImageHDU'
            ELSEIF(IHDU.EQ.1)THEN
              CTYPE='PrimaryHDU'
            ELSE
              CTYPE='ImageHDU'
            END IF
            CALL FTGIDM(IUNIT,NDIM,ISTATUS)
            NDIM=MIN(NDIM,MAXDIM)
            IF(NDIM.GT.0) CALL FTGISZ(IUNIT,NDIM,NAXES,ISTATUS)
! dimensions as a Python tuple: (), (n,), (n1, n2), ... (L is the length of
! the string built so far: TRIM would remove the blank after each comma)
            CDIM='('
            L=1
            DO K=1,NDIM
              WRITE(CNUM,'(I12)') NAXES(K)
              IF(K.GT.1)THEN
                CDIM(L+1:L+2)=', '
                L=L+2
              END IF
              CDIM(L+1:)=ADJUSTL(CNUM)
              L=LEN_TRIM(CDIM)
            END DO
            IF(NDIM.EQ.1) CDIM=TRIM(CDIM)//','
            CDIM=TRIM(CDIM)//')'
! stored type (BITPIX) and, as astropy does, the equivalent type when
! BSCALE/BZERO change it, e.g. "int16 (rescales to uint16)"
            IF(NDIM.GT.0)THEN
              CALL FTGIDT(IUNIT,BITPIX,ISTATUS)
              CALL FTGIET(IUNIT,IMGTYPE,ISTATUS)
              CFORMAT=CTYPENAME(BITPIX)
              IF(IMGTYPE.NE.BITPIX)THEN
                CFORMAT=TRIM(CFORMAT)//' (rescales to '// &
                 TRIM(CTYPENAME(IMGTYPE))//')'
              END IF
            END IF
          ELSE
! ASCII (HDUTYPE=1) or binary (HDUTYPE=2) table: rows x columns and the
! list of TFORMn values
            IF(HDUTYPE.EQ.1)THEN
              CTYPE='TableHDU'
            ELSE
              CTYPE='BinTableHDU'
            END IF
            CALL FTGNRW(IUNIT,NROWS,ISTATUS)
            CALL FTGNCL(IUNIT,NCOLS,ISTATUS)
            WRITE(CNUM,'(I12)') NROWS
            WRITE(CNUM2,'(I12)') NCOLS
            CDIM=TRIM(ADJUSTL(CNUM))//'R x '//TRIM(ADJUSTL(CNUM2))//'C'
            CFORMAT='['
            L=1
            DO K=1,NCOLS
              WRITE(CNUM,'(I12)') K
              KEYNAME='TFORM'//TRIM(ADJUSTL(CNUM))
              CALL FTGKYS(IUNIT,KEYNAME,CTFORM,COMMENT,ISTATUS)
              IF(K.GT.1)THEN
                CFORMAT(L+1:L+2)=', '
                L=L+2
              END IF
! truncate long lists (the line would be too long anyway)
              IF(L+LEN_TRIM(CTFORM).GT.200)THEN
                CFORMAT(L+1:)='...'
                EXIT
              END IF
              CFORMAT(L+1:)=CTFORM
              L=LEN_TRIM(CFORMAT)
            END DO
            CFORMAT=TRIM(CFORMAT)//']'
          END IF
! same layout as astropy: "{:3d}  {:10}  {:3} {:11}  {:5d}   {}   {}"
! (the name and type fields are widened, not truncated, when necessary)
          WRITE(CHDU,'(I3)') IHDU-1
          WRITE(CVER,'(I3)') EXTVER
          WRITE(CCARDS,'(I5)') NKEYS
          CLINE=CHDU//'  '//EXTNAME(1:MAX(10,LEN_TRIM(EXTNAME)))//'  '// &
           CVER//' '//CTYPE(1:MAX(11,LEN_TRIM(CTYPE)))//'  '//CCARDS// &
           '   '//TRIM(CDIM)//'   '//TRIM(CFORMAT)
          WRITE(*,'(A)') TRIM(CLINE)
          IF(ISTATUS.GT.0)THEN
            CALL PRINTERROR(ISTATUS)
          END IF
        END DO
!
        CONTAINS
!------------------------------------------------------------------------------
! Name of the numpy data type corresponding to a CFITSIO image type code
! (BITPIX values, plus 10, 20, 40 and 80 for the equivalent types int8,
! uint16, uint32 and uint64 returned by FTGIET).
!------------------------------------------------------------------------------
        CHARACTER*20 FUNCTION CTYPENAME(ITYPE)
        INTEGER ITYPE
        CHARACTER*12 CNUMT
!
        SELECT CASE(ITYPE)
        CASE(8)
          CTYPENAME='uint8'
        CASE(10)
          CTYPENAME='int8'
        CASE(16)
          CTYPENAME='int16'
        CASE(20)
          CTYPENAME='uint16'
        CASE(32)
          CTYPENAME='int32'
        CASE(40)
          CTYPENAME='uint32'
        CASE(64)
          CTYPENAME='int64'
        CASE(80)
          CTYPENAME='uint64'
        CASE(-32)
          CTYPENAME='float32'
        CASE(-64)
          CTYPENAME='float64'
        CASE DEFAULT
          WRITE(CNUMT,'(I12)') ITYPE
          CTYPENAME='BITPIX='//TRIM(ADJUSTL(CNUMT))
        END SELECT
        END FUNCTION CTYPENAME
        END
