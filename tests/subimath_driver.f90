!------------------------------------------------------------------------------
! Driver for test_subimath.sh: runs SUBIMATH without a graphic display and
! checks the result of the test case given as the first command-line
! argument.
!
! SUBIMATH is driven by buttons. Here the buttons are used in text mode
! (MODOTEXT_BUTT=.TRUE., without plotting them), so that RPGBAND and IFBUTTON
! read from stdin, for each "click", the answer to "Are you pressing a
! button" (y) and the button number. Other questions (constants, filter
! options, ...) are read from stdin too. PGPLOT is opened with the /NULL
! device. The button layout is the same as in the main program (xnirspec.f90).
! The sequence of answers of each test case is given in test_subimath.sh.
!
! Initial buffers (NX1 x NX2 images; buffers #k+NMAXBUFF/2 contain the errors
! of buffer #k):
!   #1: A(J,I) = 100*I + J        error EA = 3
!   #2: B(J,I) = 2 + J + 10*I     error EB = 4    (no zeros)
!   other buffers: undefined (NAXIS=0)
!
! Operations are performed in place on the first selected buffer (#1), using
! the second one (#2) when needed. Filters write the result in the second
! selected buffer (#3 here). Expected values are computed independently in
! this program.
!------------------------------------------------------------------------------
        PROGRAM SUBIMATH_DRIVER
        USE TEST_UTILS
        USE Dynamic_Array_IMAGEN
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
        INCLUDE 'dimensions.inc'
        INCLUDE 'button.inc'
!
        INTEGER NX1,NX2
        PARAMETER (NX1=6,NX2=4)
        REAL EA0,EB0
        PARAMETER (EA0=3.0,EB0=4.0)
        REAL RTOL
        PARAMETER (RTOL=1.E-5)
!
        INTEGER I,J,II,JJ,K
        INTEGER ID,PGOPEN
        INTEGER NE1,NE2,NE3                         !error buffers of #1,#2,#3
        REAL A(NX1,NX2),B(NX1,NX2),EA(NX1,NX2),EB(NX1,NX2)
        REAL FEXP(NX1,NX2),EEXP(NX1,NX2)
        REAL XCUT(NX1),XCUT_ERR(NX1)
        REAL WORK(9)
        CHARACTER*40 CCASE
! global variables (COMMONs) shared with SUBIMATH and the I/O routines
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER NFRAMES(NMAXBUFF)
        REAL CRPIX1(NMAXBUFF),CRVAL1(NMAXBUFF),CDELT1(NMAXBUFF)
        LOGICAL LWAVECAL(NMAXBUFF)
        LOGICAL LECHO
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKNFRAMES/NFRAMES
        COMMON/BLKWAVECAL1/CRPIX1,CRVAL1,CDELT1
        COMMON/BLKWAVECAL2/LWAVECAL
        COMMON/BLKLECHO/LECHO
!------------------------------------------------------------------------------
        CALL GET_COMMAND_ARGUMENT(1,CCASE)
        IF(LEN_TRIM(CCASE).EQ.0)THEN
          WRITE(*,'(A)') 'FAIL: missing test case name'
          ERROR STOP 1
        END IF
        WRITE(*,'(A)') 'Test case: '//TRIM(CCASE)
!
        CALL Initialize_Dynamic_Array_IMAGEN
        LECHO=.FALSE.
        NE1=1+NMAXBUFF/2
        NE2=2+NMAXBUFF/2
        NE3=3+NMAXBUFF/2
! PGPLOT without graphic display and buttons in text mode
        ID=PGOPEN('/NULL')
        IF(ID.LE.0)THEN
          WRITE(*,'(A)') 'FAIL: PGOPEN(''/NULL'')'
          ERROR STOP 1
        END IF
        MODOTEXT_BUTT=.FALSE.
        MODOTEXT_PLOTBUTT=.FALSE.
        EXIST_BUTT=.FALSE.
        CALL BUTTSPR(0.440,0.965,0.050,0.750)
        CALL BUTTSXB(10)
        CALL BUTTSYB(25)
        CALL BUTTSCF(1)
        CALL BUTTSCH(0.8)
        MODOTEXT_BUTT=.TRUE.
! initial buffers
        DO I=1,NX2
          DO J=1,NX1
            A(J,I)=100.0*REAL(I)+REAL(J)
            B(J,I)=2.0+REAL(J)+10.0*REAL(I)
          END DO
        END DO
        EA=EA0
        EB=EB0
        SELECT CASE(TRIM(CCASE))
        CASE('div_frames_zero')
          B(3,2)=0.0                              !division by zero at (3,2)
        CASE('median_filter')
          A(3,2)=1.E4                       !outlier removed by the median
        CASE('xcut_mask','xcut_div_zero')
          B(2,1)=-999.0                       !masked pixels in buffer #2
          B(2,3)=-999.0
          B(5,:)=-999.0                             !fully masked column
        END SELECT
        NAXIS=0
        NFRAMES=0
        LWAVECAL=.FALSE.
        CALL SET_BUFFER(1,A)
        CALL SET_BUFFER(2,B)
        CALL SET_BUFFER(NE1,EA)
        CALL SET_BUFFER(NE2,EB)
        LWAVECAL(2)=.TRUE.
        CRPIX1(2)=1.0
        CRVAL1(2)=5000.0
        CDELT1(2)=2.0
        LWAVECAL(1)=.TRUE.
        CRPIX1(1)=2.0
        CRVAL1(1)=6000.0
        CDELT1(1)=1.5
        LWAVECAL(NE1)=.TRUE.
        CRPIX1(NE1)=2.0
        CRVAL1(NE1)=6000.0
        CDELT1(NE1)=1.5
!------------------------------------------------------------------------------
        CALL SUBIMATH
!------------------------------------------------------------------------------
        SELECT CASE(TRIM(CCASE))
!..............................................................................
        CASE('add_frames')                                     !#1 = #1 + #2
          CALL CHECK_IMAGE('data #1',1,A+B)
          CALL CHECK_IMAGE('error #1',NE1,SQRT(EA**2+EB**2))
          CALL CHECK_UNCHANGED_2()
!..............................................................................
        CASE('sub_frames')                                     !#1 = #1 - #2
          CALL CHECK_IMAGE('data #1',1,A-B)
          CALL CHECK_IMAGE('error #1',NE1,SQRT(EA**2+EB**2))
          CALL CHECK_UNCHANGED_2()
!..............................................................................
        CASE('mul_frames')                                     !#1 = #1 * #2
          CALL CHECK_IMAGE('data #1',1,A*B)
          CALL CHECK_IMAGE('error #1',NE1,SQRT(A**2*EB**2+B**2*EA**2))
          CALL CHECK_UNCHANGED_2()
!..............................................................................
        CASE('div_frames')                                     !#1 = #1 / #2
          CALL CHECK_IMAGE('data #1',1,A/B)
          CALL CHECK_IMAGE('error #1',NE1,SQRT(A**2*EB**2+B**2*EA**2)/B**2)
          CALL CHECK_UNCHANGED_2()
!..............................................................................
        CASE('div_frames_zero')        !#1 = #1 / #2, with B(3,2)=0 -> -1, -2
          DO I=1,NX2
            DO J=1,NX1
              IF(B(J,I).EQ.0.0)THEN
                FEXP(J,I)=-1.0
                EEXP(J,I)=-2.0
              ELSE
                FEXP(J,I)=A(J,I)/B(J,I)
                EEXP(J,I)=SQRT(A(J,I)**2*EB(J,I)**2+B(J,I)**2*EA(J,I)**2)/ &
                 B(J,I)**2
              END IF
            END DO
          END DO
          CALL CHECK_IMAGE('data #1',1,FEXP)
          CALL CHECK_IMAGE('error #1',NE1,EEXP)
!..............................................................................
        CASE('add_constant')                    !#1 = #1 + 7 (error of 7: 0.5)
          CALL CHECK_IMAGE('data #1',1,A+7.0)
          CALL CHECK_IMAGE('error #1',NE1,SQRT(EA**2+0.5**2))
!..............................................................................
        CASE('div_constant_zero')    !#1 = #1 / 0 (rejected), then #1 = #1 / 4
          CALL CHECK_IMAGE('data #1',1,A/4.0)
          CALL CHECK_IMAGE('error #1',NE1,EA/4.0)
!..............................................................................
        CASE('add_constant_region')         !#1 = #1 + 10 in X=2..3, Y=2..4
          FEXP=A
          FEXP(2:3,2:4)=A(2:3,2:4)+10.0
          CALL CHECK_IMAGE('data #1',1,FEXP)
          CALL CHECK_IMAGE('error #1',NE1,EA)
!..............................................................................
        CASE('copy')                                            !#3 <- #2
          CALL CHECK_INT('NAXIS1 #3',NAXIS(1,3),NX1)
          CALL CHECK_INT('NAXIS2 #3',NAXIS(2,3),NX2)
          CALL CHECK_IMAGE('data #3',3,B)
          CALL CHECK_INT('NAXIS1 error #3',NAXIS(1,NE3),NX1)
          CALL CHECK_IMAGE('error #3',NE3,EB)
          CALL CHECK_TRUE('LWAVECAL #3 copied',LWAVECAL(3))
          CALL CHECK_REAL('CRVAL1 #3 copied',CRVAL1(3),5000.0,0.0,0.0)
          CALL CHECK_REAL('CDELT1 #3 copied',CDELT1(3),2.0,0.0,0.0)
          CALL CHECK_UNCHANGED_2()
!..............................................................................
        CASE('reverse_x')                           !#3 = #1 reversed along X
          DO J=1,NX1
            FEXP(J,:)=A(NX1-J+1,:)
          END DO
          CALL CHECK_FILTER_OUTPUT(FEXP)
! reversed X axis: CRPIX1 -> NAXIS1-CRPIX1+1, CDELT1 -> -CDELT1
          CALL CHECK_WAVECAL('#3',3,.TRUE.,REAL(NX1)-2.0+1.0,6000.0,-1.5)
!..............................................................................
        CASE('reverse_y')                           !#3 = #1 reversed along Y
          DO I=1,NX2
            FEXP(:,I)=A(:,NX2-I+1)
          END DO
          CALL CHECK_FILTER_OUTPUT(FEXP)
          CALL CHECK_WAVECAL('#3',3,.TRUE.,2.0,6000.0,1.5)
!..............................................................................
        CASE('median_filter')        !#3 = 3x3 median of #1 (box clipped at
          DO I=1,NX2                 !the borders of the image)
            DO J=1,NX1
              K=0
              DO II=MAX(1,I-1),MIN(NX2,I+1)
                DO JJ=MAX(1,J-1),MIN(NX1,J+1)
                  K=K+1
                  WORK(K)=A(JJ,II)
                END DO
              END DO
              FEXP(J,I)=MEDIAN(K,WORK)
            END DO
          END DO
          CALL CHECK_FILTER_OUTPUT(FEXP)
!..............................................................................
        CASE('bin_expand')            !#3 = #1 binned by 2 along X, expanded
          DO J=1,NX1,2
            FEXP(J,:)=0.5*(A(J,:)+A(J+1,:))
            FEXP(J+1,:)=FEXP(J,:)
          END DO
          CALL CHECK_FILTER_OUTPUT(FEXP)
!..............................................................................
        CASE('xcut_mask')  !#1 = #1 - mean X-cut of #2 (rows 1-4, mask -999)
          CALL COMPUTE_XCUT()
          DO J=1,NX1
            FEXP(J,:)=A(J,:)-XCUT(J)
            EEXP(J,:)=SQRT(EA(J,:)**2+XCUT_ERR(J)**2)
          END DO
          CALL CHECK_IMAGE('data #1',1,FEXP)
          CALL CHECK_IMAGE('error #1',NE1,EEXP)
          CALL CHECK_IMAGE('data #2 unchanged',2,B)
!..............................................................................
        CASE('xcut_div_zero')     !#1 = #1 / mean X-cut of #2 (as xcut_mask);
          CALL COMPUTE_XCUT()     !column 5 is fully masked, so its X-cut is
          DO J=1,NX1              !zero -> data -1, error -2 (FZERO, FZERO_ERR)
            IF(XCUT(J).EQ.0.0)THEN
              FEXP(J,:)=-1.0
              EEXP(J,:)=-2.0
            ELSE
              FEXP(J,:)=A(J,:)/XCUT(J)
              EEXP(J,:)=SQRT(EA(J,:)**2*XCUT(J)**2+A(J,:)**2*XCUT_ERR(J)**2)/ &
               XCUT(J)**2
            END IF
          END DO
          CALL CHECK_IMAGE('data #1',1,FEXP)
          CALL CHECK_IMAGE('error #1',NE1,EEXP)
!..............................................................................
        CASE('shift_wavecal')        !#3 = #1 shifted by +1.5 pixels along X:
          CALL CHECK_INT('NAXIS1 #3',NAXIS(1,3),NX1) !calibration of #1 with
          CALL CHECK_WAVECAL('#3',3,.TRUE.,3.5,6000.0,1.5) !CRPIX1 + 1.5
          CALL CHECK_WAVECAL('error #3',NE3,.TRUE.,3.5,6000.0,1.5)
!..............................................................................
        CASE('quad_4_wavecal')    !#1 <- 4 quadrants of #2: no calibration
          CALL CHECK_WAVECAL('#1',1,.FALSE.,1.0,0.0,0.0)
          CALL CHECK_WAVECAL('error #1',NE1,.FALSE.,1.0,0.0,0.0)
!..............................................................................
        CASE DEFAULT
          WRITE(*,'(A)') 'FAIL: unknown test case '//TRIM(CCASE)
          ERROR STOP 1
        END SELECT
!
        CALL PGCLOS
        CALL Deallocate_Array_IMAGEN
        CALL TEST_SUMMARY()
!
        CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE SET_BUFFER(NBUFF,F)
!
! Store the NX1 x NX2 image F in buffer NBUFF and set its dimensions.
!
! Parameters
! ----------
! NBUFF : INTEGER
!     Buffer number.
! F : REAL(NX1,NX2)
!     Image.
!------------------------------------------------------------------------------
        SUBROUTINE SET_BUFFER(NBUFF,F)
        INTEGER, INTENT(IN) :: NBUFF
        REAL, INTENT(IN) :: F(NX1,NX2)
!
        IMAGEN(1:NX1,1:NX2,NBUFF)=F
        NAXIS(1,NBUFF)=NX1
        NAXIS(2,NBUFF)=NX2
        NFRAMES(NBUFF)=1
        END SUBROUTINE SET_BUFFER
!------------------------------------------------------------------------------
! SUBROUTINE CHECK_IMAGE(LABEL,NBUFF,FEXPECTED)
!
! Check that buffer NBUFF contains the NX1 x NX2 image FEXPECTED, with
! relative tolerance RTOL (and absolute tolerance RTOL for values close to
! zero). The first mismatched pixel, if any, is printed.
!
! Parameters
! ----------
! LABEL : CHARACTER*(*)
!     Name of the check.
! NBUFF : INTEGER
!     Buffer number.
! FEXPECTED : REAL(NX1,NX2)
!     Expected image.
!------------------------------------------------------------------------------
        SUBROUTINE CHECK_IMAGE(LABEL,NBUFF,FEXPECTED)
        CHARACTER*(*), INTENT(IN) :: LABEL
        INTEGER, INTENT(IN) :: NBUFF
        REAL, INTENT(IN) :: FEXPECTED(NX1,NX2)
!
        INTEGER I0,J0,NBAD
        REAL F0
!
        NBAD=0
        DO I0=1,NX2
          DO J0=1,NX1
            F0=IMAGEN(J0,I0,NBUFF)
            IF(ABS(F0-FEXPECTED(J0,I0)).GT. &
             RTOL*(1.0+ABS(FEXPECTED(J0,I0))))THEN
              NBAD=NBAD+1
              IF(NBAD.EQ.1)THEN
                WRITE(*,'(A,I0,A,I0,A,ES16.8,A,ES16.8)') &
                 '      first mismatch at (',J0,',',I0,'): got ',F0, &
                 ', expected ',FEXPECTED(J0,I0)
              END IF
            END IF
          END DO
        END DO
        CALL CHECK_INT(LABEL//' (mismatched pixels)',NBAD,0)
        END SUBROUTINE CHECK_IMAGE
!------------------------------------------------------------------------------
! SUBROUTINE CHECK_WAVECAL(LABEL,NBUFF,LEXP,CRPIX1EXP,CRVAL1EXP,CDELT1EXP)
!
! Check the wavelength calibration of buffer NBUFF (LWAVECAL, CRPIX1, CRVAL1
! and CDELT1).
!------------------------------------------------------------------------------
        SUBROUTINE CHECK_WAVECAL(LABEL,NBUFF,LEXP,CRPIX1EXP,CRVAL1EXP, &
         CDELT1EXP)
        CHARACTER*(*), INTENT(IN) :: LABEL
        INTEGER, INTENT(IN) :: NBUFF
        LOGICAL, INTENT(IN) :: LEXP
        REAL, INTENT(IN) :: CRPIX1EXP,CRVAL1EXP,CDELT1EXP
!
        CALL CHECK_TRUE('LWAVECAL '//LABEL,LWAVECAL(NBUFF).EQV.LEXP)
        CALL CHECK_REAL('CRPIX1 '//LABEL,CRPIX1(NBUFF),CRPIX1EXP,1.E-6,0.0)
        CALL CHECK_REAL('CRVAL1 '//LABEL,CRVAL1(NBUFF),CRVAL1EXP,1.E-6,0.0)
        CALL CHECK_REAL('CDELT1 '//LABEL,CDELT1(NBUFF),CDELT1EXP,1.E-6,0.0)
        END SUBROUTINE CHECK_WAVECAL
!------------------------------------------------------------------------------
! SUBROUTINE COMPUTE_XCUT()
!
! Expected mean X-cut of buffer #2 over rows 1 to NX2, excluding the pixels
! with the mask value -999, and its error (FMEAN0E: SQRT(sum of squared
! errors)/K). Fully masked columns give XCUT=XCUT_ERR=0.
!------------------------------------------------------------------------------
        SUBROUTINE COMPUTE_XCUT()
        INTEGER I0,J0,K0
!
        DO J0=1,NX1
          K0=0
          XCUT(J0)=0.0
          DO I0=1,NX2
            IF(B(J0,I0).NE.-999.0)THEN
              K0=K0+1
              XCUT(J0)=XCUT(J0)+B(J0,I0)
            END IF
          END DO
          IF(K0.GT.0)THEN
            XCUT(J0)=XCUT(J0)/REAL(K0)
            XCUT_ERR(J0)=EB0*SQRT(REAL(K0))/REAL(K0)
          ELSE
            XCUT(J0)=0.0
            XCUT_ERR(J0)=0.0
          END IF
        END DO
        END SUBROUTINE COMPUTE_XCUT
!------------------------------------------------------------------------------
! SUBROUTINE CHECK_UNCHANGED_2()
!
! Check that buffer #2 and its error buffer have not been modified.
!------------------------------------------------------------------------------
        SUBROUTINE CHECK_UNCHANGED_2()
        CALL CHECK_IMAGE('data #2 unchanged',2,B)
        CALL CHECK_IMAGE('error #2 unchanged',NE2,EB)
        END SUBROUTINE CHECK_UNCHANGED_2
!------------------------------------------------------------------------------
! SUBROUTINE CHECK_FILTER_OUTPUT(FEXPECTED)
!
! Check the result of a filter applied to #1 with output in #3: the data of
! #3 must be FEXPECTED, the error buffer of #3 must be zero (filters do not
! propagate errors), and #1 must not be modified.
!
! Parameters
! ----------
! FEXPECTED : REAL(NX1,NX2)
!     Expected image in buffer #3.
!------------------------------------------------------------------------------
        SUBROUTINE CHECK_FILTER_OUTPUT(FEXPECTED)
        REAL, INTENT(IN) :: FEXPECTED(NX1,NX2)
!
        REAL ZERO(NX1,NX2)
!
        ZERO=0.0
        CALL CHECK_INT('NAXIS1 #3',NAXIS(1,3),NX1)
        CALL CHECK_INT('NAXIS2 #3',NAXIS(2,3),NX2)
        CALL CHECK_IMAGE('data #3',3,FEXPECTED)
        CALL CHECK_IMAGE('error #3 set to zero',NE3,ZERO)
        CALL CHECK_IMAGE('data #1 unchanged',1,A)
        END SUBROUTINE CHECK_FILTER_OUTPUT
!------------------------------------------------------------------------------
! REAL FUNCTION MEDIAN(N,X)
!
! Median of X(1:N) (insertion sort of a copy; independent of FMEDIAN1).
!
! Parameters
! ----------
! N : INTEGER
!     Number of values.
! X : REAL(N)
!     Values (not modified).
!
! Returns
! -------
! MEDIAN : REAL
!     Median value (mean of the two central values if N is even).
!------------------------------------------------------------------------------
        REAL FUNCTION MEDIAN(N,X)
        INTEGER, INTENT(IN) :: N
        REAL, INTENT(IN) :: X(N)
!
        INTEGER L1,L2
        REAL Y(N),T
!
        Y=X
        DO L1=2,N
          T=Y(L1)
          L2=L1-1
          DO WHILE(L2.GE.1)
            IF(Y(L2).LE.T) EXIT
            Y(L2+1)=Y(L2)
            L2=L2-1
          END DO
          Y(L2+1)=T
        END DO
        IF(MOD(N,2).EQ.0)THEN
          MEDIAN=0.5*(Y(N/2)+Y(N/2+1))
        ELSE
          MEDIAN=Y(N/2+1)
        END IF
        END FUNCTION MEDIAN
!
        END PROGRAM SUBIMATH_DRIVER
