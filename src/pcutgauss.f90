!------------------------------------------------------------------------------
! SUBROUTINE PCUTGAUSS(LFIT,X0,SIGMA,AMP)
!
! Interactive gaussian fit to the 1-D profile displayed in the upper-left
! plot (the last profile drawn by SUBPLOT). The user selects with the cursor
! (two clicks each) the left continuum region, the right continuum region
! and the line region. The continuum is fitted with a straight line (see
! CUTCONTFIT) and drawn as soon as both continuum regions are defined; then
! a gaussian is fitted to the profile minus the continuum in the line region
! (see CUTGAUSSFIT). The selected regions, the continuum and the fit are
! drawn on the plot, and the fitted parameters are displayed in the
! terminal. Pressing the right mouse button (character 'X') at any click
! cancels the fit.
!
! Returns
! -------
! LFIT : LOGICAL
!     .TRUE. if the fit has been performed.
! X0, SIGMA, AMP : REAL
!     Center, sigma and amplitude of the fitted gaussian (zero if LFIT is
!     .FALSE.).
!------------------------------------------------------------------------------
        SUBROUTINE PCUTGAUSS(LFIT,X0,SIGMA,AMP)
        IMPLICIT NONE
        LOGICAL LFIT
        REAL X0,SIGMA,AMP
!
        INCLUDE 'largest.inc'
!
        INTEGER NSIMUL
        PARAMETER (NSIMUL=100)       !simulations to estimate the errors
        INTEGER NCURVE
        PARAMETER (NCURVE=500)       !points to draw the fitted function
        INTEGER ICOLCONT,ICOLLINE
        PARAMETER (ICOLCONT=8,ICOLLINE=7)  !colors: continuum (orange), line
                                           !(yellow); not used by the buffers
        REAL FWHMFACTOR                                      !2*SQRT(2*LN(2))
        PARAMETER (FWHMFACTOR=2.354820)
        REAL SQRT2PI
        PARAMETER (SQRT2PI=2.506628)
!
        INTEGER I,K
        INTEGER NCONT,NLINE,IERR
        INTEGER OLD_CI,OLD_LS
        REAL XLIM(6)
        REAL XC,YC
        REAL A0,A1,RMSCONT
        REAL EX0,ESIGMA,EAMP
        REAL XV1,XV2,YV1,YV2,XW1,XW2,YW1,YW2
        REAL XA,XB,XP(NCURVE),YP(NCURVE)
        CHARACTER*1 CH
        CHARACTER*50 CPROMPT(3)
        LOGICAL LERRSIM
! global variables (COMMONs) shared with SUBPLOT and the main program
        INTEGER NLAST,N1LAST,N2LAST
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XPLAST(NXYMAX*NOVERSAMPMAX),YPLAST(NXYMAX*NOVERSAMPMAX)
        REAL YRMSTOL
        COMMON/BLKPLIMITS/XMIN,XMAX,YMIN,YMAX
        COMMON/BLKPNLAST/NLAST,N1LAST,N2LAST
        COMMON/BLKPLAST/XPLAST,YPLAST
        COMMON/BLKDEFAULTS9/YRMSTOL
!
        DATA CPROMPT /'left continuum region', &
                      'right continuum region', &
                      'line region'/
!------------------------------------------------------------------------------
        LFIT=.FALSE.
        X0=0.0
        SIGMA=0.0
        AMP=0.0
        IF(NLAST.LE.0)THEN
          WRITE(*,101) 'ERROR: there is no 1-D profile to be fitted'
          RETURN
        END IF
! use the viewport and world coordinates of the upper-left plot (as defined
! in SUBPLOT)
        CALL PGQVP(0,XV1,XV2,YV1,YV2)
        CALL PGQWIN(XW1,XW2,YW1,YW2)
        CALL PGQCI(OLD_CI)
        CALL PGQLS(OLD_LS)
        CALL PGSVP(0.05,0.37,0.46,0.72)
        CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
!------------------------------------------------------------------------------
! select the regions (two clicks each); every click is marked with a
! vertical line
        WRITE(*,101) '* Gaussian fit to the 1-D profile in the upper-left plot'
        WRITE(*,101) '  (right mouse button: cancel)'
        DO K=1,3
          WRITE(*,101) '=> Click twice to define the '//TRIM(CPROMPT(K))
          DO I=1,2
            CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
            IF(CH.EQ.'X')THEN
              WRITE(*,101) '[CANCEL]'
              GOTO 900
            END IF
            XLIM(2*(K-1)+I)=XC
            IF(K.LE.2)THEN
              CALL PGSCI(ICOLCONT)
            ELSE
              CALL PGSCI(ICOLLINE)
            END IF
            CALL PGMOVE(XC,YMIN)
            CALL PGDRAW(XC,YMAX)
          END DO
! both continuum regions defined: fit and draw the continuum (dashed line
! between the outermost continuum limits)
          IF(K.EQ.2)THEN
            CALL CUTCONTFIT(NXYMAX*NOVERSAMPMAX,N1LAST,N2LAST, &
             XPLAST,YPLAST,XLIM(1:4),A0,A1,NCONT,RMSCONT,IERR)
            IF(IERR.NE.0)THEN
              WRITE(*,101) 'ERROR: fewer than 2 points (with different '// &
               'X) in the continuum regions'
              GOTO 900
            END IF
            XA=MINVAL(XLIM(1:4))
            XB=MAXVAL(XLIM(1:4))
            CALL PGSCI(ICOLCONT)
            CALL PGSLS(2)
            CALL PGMOVE(XA,A0+A1*XA)
            CALL PGDRAW(XB,A0+A1*XB)
            CALL PGSLS(1)
          END IF
        END DO
!------------------------------------------------------------------------------
        CALL CUTGAUSSFIT(NXYMAX*NOVERSAMPMAX,N1LAST,N2LAST,XPLAST,YPLAST, &
         XLIM,YRMSTOL,NSIMUL,A0,A1,NCONT,RMSCONT, &
         X0,SIGMA,AMP,EX0,ESIGMA,EAMP,LERRSIM,NLINE,IERR)
        IF(IERR.EQ.1)THEN
          WRITE(*,101) 'ERROR: fewer than 2 points (with different X) '// &
           'in the continuum regions'
          GOTO 900
        ELSEIF(IERR.EQ.2)THEN
          WRITE(*,101) 'ERROR: fewer than 3 points in the line region'
          GOTO 900
        ELSEIF(IERR.EQ.3)THEN
          WRITE(*,101) 'ERROR: too many points in the line region'
          GOTO 900
        END IF
        LFIT=.TRUE.
!------------------------------------------------------------------------------
! draw the fit (gaussian + continuum) in the line region (the continuum has
! already been drawn)
        XA=MIN(XLIM(5),XLIM(6))
        XB=MAX(XLIM(5),XLIM(6))
        DO I=1,NCURVE
          XP(I)=XA+(XB-XA)*REAL(I-1)/REAL(NCURVE-1)
          YP(I)=A0+A1*XP(I)+AMP*EXP(-(XP(I)-X0)**2/(2.0*SIGMA*SIGMA))
        END DO
        CALL PGSCI(ICOLLINE)
        CALL PGLINE(NCURVE,XP,YP)
!------------------------------------------------------------------------------
! fitted parameters
        WRITE(*,101) '------------------------------------------------------'
        WRITE(*,'(A,I0,A)') 'Continuum: y = a0 + a1*x  (',NCONT,' points)'
        WRITE(*,'(A,ES15.7)') '  a0                  : ',A0
        WRITE(*,'(A,ES15.7)') '  a1                  : ',A1
        WRITE(*,'(A,ES15.7)') '  rms of residuals    : ',RMSCONT
        WRITE(*,'(A,I0,A)') 'Gaussian fit to profile - continuum  (',NLINE, &
         ' points)'
        WRITE(*,'(A)') '  y = amp*exp[-(x-x0)^2/(2*sigma^2)]'
        WRITE(*,'(A,ES15.7,A,ES12.4)') '  x0                  : ',X0, &
         ' +/- ',EX0
        WRITE(*,'(A,ES15.7,A,ES12.4)') '  sigma               : ',SIGMA, &
         ' +/- ',ESIGMA
        WRITE(*,'(A,ES15.7,A,ES12.4)') '  FWHM                : ', &
         FWHMFACTOR*SIGMA,' +/- ',FWHMFACTOR*ESIGMA
        WRITE(*,'(A,ES15.7,A,ES12.4)') '  amp                 : ',AMP, &
         ' +/- ',EAMP
        WRITE(*,'(A,ES15.7)') '  area (amp*sigma*sqrt(2*pi)): ', &
         AMP*SIGMA*SQRT2PI
        IF(LERRSIM)THEN
          WRITE(*,'(A,I0,A)') '  (errors from ',NSIMUL, &
           ' simulations with noise = rms of continuum residuals)'
        ELSE
          WRITE(*,'(A)') '  (errors: rms of the final simplex of DOWNHILL)'
        END IF
        WRITE(*,101) '------------------------------------------------------'
!------------------------------------------------------------------------------
! restore the previous plot settings
900     CALL PGSCI(OLD_CI)
        CALL PGSLS(OLD_LS)
        CALL PGSVP(XV1,XV2,YV1,YV2)
        CALL PGSWIN(XW1,XW2,YW1,YW2)
101     FORMAT(A)
        END
