!------------------------------------------------------------------------------
! Driver for test_hplot.sh: draws a synthetic 1-D profile with SUBPLOT (as
! the cuts of the main program) and runs HPLOT('manual'), which reads the
! new limits of the upper-left plot (Xmin, Xmax, Ymin, Ymax) from stdin, or
! HPLOT('min,max'), which computes the Y limits from the displayed cuts
! within the current X limits. PGPLOT is opened with the /NULL device.
!
! Test cases (first command-line argument):
!   default   : four empty answers; the limits computed by SUBPLOT must
!               not change
!   set       : Xmin=10, Xmax=50, Ymin=-1, Ymax=100
!   retry     : Xmin=50, Xmax=10 (rejected), then Xmin=10, Xmax=50;
!               Ymin=5, Ymax=5 (rejected), then Ymin=0, Ymax=1
!   noprofile : no profile drawn (NLAST=0): nothing is asked (stdin is
!               empty) and the limits do not change
!   minmax    : 'min,max' with X limits [10.5,20.5]: points x=11,...,20,
!               Y limits 11 and 20 enlarged by 5% of the range on each side
!   minmaxover: as minmax, with the cut of another buffer (y = 2x) also
!               displayed: Y limits 11 and 40 (enlarged)
!   minmaxnone: 'min,max' with X limits [200,300] (no points): Y limits
!               unchanged
!
! Profile: y = x, x = 1,...,100.
!------------------------------------------------------------------------------
        PROGRAM HPLOT_DRIVER
        USE TEST_UTILS
        USE Dynamic_Arrays_XYPB
        IMPLICIT NONE
        INCLUDE 'interface_xypb.inc'
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
!
        INTEGER N
        PARAMETER (N=100)
!
        INTEGER I
        INTEGER ID,PGOPEN
        REAL X(N),Y(N)
        REAL XMIN0,XMAX0,YMIN0,YMAX0
        CHARACTER*20 CCASE
! global variables (COMMONs) shared with SUBPLOT, HPLOT and READF
        INTEGER NLAST,N1LAST,N2LAST
        REAL XMIN,XMAX,YMIN,YMAX
        LOGICAL LPB(NMAXBUFF)
        LOGICAL LECHO
        COMMON/BLKPLIMITS/XMIN,XMAX,YMIN,YMAX
        COMMON/BLKPNLAST/NLAST,N1LAST,N2LAST
        COMMON/BLKP_ALLBUFF_L/LPB
        COMMON/BLKLECHO/LECHO
!------------------------------------------------------------------------------
        CALL GET_COMMAND_ARGUMENT(1,CCASE)
        WRITE(*,'(A)') 'Test case: '//TRIM(CCASE)
        LECHO=.FALSE.
        CALL Initialize_Dynamic_Arrays_XYPB
        DO I=1,NMAXBUFF
          LPB(I)=.FALSE.                   !no additional cuts to be drawn
        END DO
        ID=PGOPEN('/NULL')
        IF(ID.LE.0)THEN
          WRITE(*,'(A)') 'FAIL: PGOPEN(''/NULL'')'
          ERROR STOP 1
        END IF
        DO I=1,N
          X(I)=REAL(I)
          Y(I)=REAL(I)
        END DO
        IF(TRIM(CCASE).EQ.'noprofile')THEN
          NLAST=0
          XMIN=1.0
          XMAX=2.0
          YMIN=3.0
          YMAX=4.0
        ELSE
          CALL SUBPLOT(N,1,N,X,Y,X,Y,.TRUE.,.TRUE.,.FALSE.,.FALSE., &
           'x axis','signal','[1,1]',1,201,1.0)
        END IF
! cuts of the data buffers as stored by PCUT: buffer #1 is the displayed
! profile and, in minmaxover, buffer #2 is overplotted
        DO I=1,N
          XPB(I,1)=X(I)
          YPB(I,1)=Y(I)
          XPB(I,2)=X(I)
          YPB(I,2)=2.0*Y(I)
        END DO
        IF(CCASE(1:6).EQ.'minmax')THEN
          LPB(1)=.TRUE.
          LPB(2)=(TRIM(CCASE).EQ.'minmaxover')
          IF(TRIM(CCASE).EQ.'minmaxnone')THEN
            XMIN=200.0
            XMAX=300.0
          ELSE
            XMIN=10.5
            XMAX=20.5
          END IF
        END IF
        XMIN0=XMIN
        XMAX0=XMAX
        YMIN0=YMIN
        YMAX0=YMAX
!------------------------------------------------------------------------------
        IF(CCASE(1:6).EQ.'minmax')THEN
          CALL HPLOT('min,max')
        ELSE
          CALL HPLOT('manual')
        END IF
!------------------------------------------------------------------------------
        SELECT CASE(TRIM(CCASE))
        CASE('default','noprofile')
! the defaults are written with list-directed output, hence the tolerance
          CALL CHECK_REAL('xmin',XMIN,XMIN0,1.E-6,0.0)
          CALL CHECK_REAL('xmax',XMAX,XMAX0,1.E-6,0.0)
          CALL CHECK_REAL('ymin',YMIN,YMIN0,1.E-6,0.0)
          CALL CHECK_REAL('ymax',YMAX,YMAX0,1.E-6,0.0)
        CASE('set')
          CALL CHECK_REAL('xmin',XMIN,10.0,0.0,0.0)
          CALL CHECK_REAL('xmax',XMAX,50.0,0.0,0.0)
          CALL CHECK_REAL('ymin',YMIN,-1.0,0.0,0.0)
          CALL CHECK_REAL('ymax',YMAX,100.0,0.0,0.0)
        CASE('retry')
          CALL CHECK_REAL('xmin',XMIN,10.0,0.0,0.0)
          CALL CHECK_REAL('xmax',XMAX,50.0,0.0,0.0)
          CALL CHECK_REAL('ymin',YMIN,0.0,0.0,0.0)
          CALL CHECK_REAL('ymax',YMAX,1.0,0.0,0.0)
        CASE('minmax')
          CALL CHECK_REAL('xmin unchanged',XMIN,XMIN0,0.0,0.0)
          CALL CHECK_REAL('xmax unchanged',XMAX,XMAX0,0.0,0.0)
          CALL CHECK_REAL('ymin',YMIN,11.0-9.0/20.0,1.E-6,0.0)
          CALL CHECK_REAL('ymax',YMAX,20.0+9.0/20.0,1.E-6,0.0)
        CASE('minmaxover')
          CALL CHECK_REAL('ymin',YMIN,11.0-29.0/20.0,1.E-6,0.0)
          CALL CHECK_REAL('ymax',YMAX,40.0+29.0/20.0,1.E-6,0.0)
        CASE('minmaxnone')
          CALL CHECK_REAL('ymin unchanged',YMIN,YMIN0,0.0,0.0)
          CALL CHECK_REAL('ymax unchanged',YMAX,YMAX0,0.0,0.0)
        CASE DEFAULT
          WRITE(*,'(A)') 'FAIL: unknown test case '//TRIM(CCASE)
          ERROR STOP 1
        END SELECT
        IF((TRIM(CCASE).NE.'noprofile').AND.(CCASE(1:6).NE.'minmax'))THEN
! the whole profile is redrawn
          CALL CHECK_INT('n1last',N1LAST,1)
          CALL CHECK_INT('n2last',N2LAST,N)
        END IF
        CALL PGCLOS
        CALL Deallocate_Arrays_XYPB
        CALL TEST_SUMMARY()
        END PROGRAM HPLOT_DRIVER
