!------------------------------------------------------------------------------
! Driver for test_pcutgauss.sh: draws a synthetic 1-D profile with SUBPLOT
! (as the cuts of the main program) and runs the interactive gaussian fit
! PCUTGAUSS, with PGPLOT opened with the /NULL device and buttons in text
! mode, so that every click of RPGBAND is read from stdin as "n", XC, YC
! and the character (A=left button, X=right button: cancel).
!
! Test cases (first command-line argument):
!   fit    : six clicks (left continuum 5-25, right continuum 90-55, line
!            30-50); the fitted gaussian must be x0=40, sigma=3, amp=50
!   cancel : right button at the third click; no fit
!   nocont : continuum regions without points (200-210, 300-310): error
!            after the fourth click (the line region is not asked for)
!
! Profile: y = 10 + 0.5*x + 50*exp[-(x-40)^2/(2*3^2)], x = 1,...,100.
!------------------------------------------------------------------------------
        PROGRAM PCUTGAUSS_DRIVER
        USE TEST_UTILS
        IMPLICIT NONE
        INCLUDE 'button.inc'
!
        INTEGER N
        PARAMETER (N=100)
!
        INTEGER I
        INTEGER ID,PGOPEN
        REAL X(N),Y(N)
        REAL X0,SIGMA,AMP
        CHARACTER*20 CCASE
        LOGICAL LFIT
! global variables (COMMONs) shared with PCUTGAUSS and the input routines
        REAL YRMSTOL
        LOGICAL LECHO
        COMMON/BLKDEFAULTS9/YRMSTOL
        COMMON/BLKLECHO/LECHO
!------------------------------------------------------------------------------
        CALL GET_COMMAND_ARGUMENT(1,CCASE)
        WRITE(*,'(A)') 'Test case: '//TRIM(CCASE)
        LECHO=.FALSE.
        YRMSTOL=1.E-5                          !as in the main program
        ID=PGOPEN('/NULL')
        IF(ID.LE.0)THEN
          WRITE(*,'(A)') 'FAIL: PGOPEN(''/NULL'')'
          ERROR STOP 1
        END IF
        MODOTEXT_BUTT=.TRUE.
        MODOTEXT_PLOTBUTT=.FALSE.
! synthetic profile drawn in the upper-left plot (stored by SUBPLOT)
        DO I=1,N
          X(I)=REAL(I)
          Y(I)=10.0+0.5*X(I)+50.0*EXP(-(X(I)-40.0)**2/(2.0*9.0))
        END DO
        CALL SUBPLOT(N,1,N,X,Y,X,Y,.TRUE.,.TRUE.,.FALSE.,.FALSE., &
         'x axis','signal','[1,1]',1,201,1.0)
!------------------------------------------------------------------------------
        CALL PCUTGAUSS(LFIT,X0,SIGMA,AMP)
!------------------------------------------------------------------------------
        SELECT CASE(TRIM(CCASE))
        CASE('fit')
          CALL CHECK_TRUE('fit performed',LFIT)
          CALL CHECK_REAL('x0',X0,40.0,0.0,1.E-2)
          CALL CHECK_REAL('sigma',SIGMA,3.0,1.E-2,0.0)
          CALL CHECK_REAL('amp',AMP,50.0,1.E-2,0.0)
        CASE('cancel')
          CALL CHECK_TRUE('fit cancelled',.NOT.LFIT)
          CALL CHECK_REAL('amp set to zero',AMP,0.0,0.0,0.0)
        CASE('nocont')
          CALL CHECK_TRUE('no fit without continuum points',.NOT.LFIT)
        CASE DEFAULT
          WRITE(*,'(A)') 'FAIL: unknown test case '//TRIM(CCASE)
          ERROR STOP 1
        END SELECT
        CALL PGCLOS
        CALL TEST_SUMMARY()
        END PROGRAM PCUTGAUSS_DRIVER
