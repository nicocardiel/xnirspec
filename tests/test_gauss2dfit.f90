!------------------------------------------------------------------------------
! Unit tests for GAUSS2DFIT: fit of a 2D gaussian (without initial guess,
! i.e. starting from the 1D fits to the mean profiles) to synthetic,
! noise-free images in a 15x15 box, as in CENTROID:
!
!   Z(X,Y) = CTE + AMP*EXP[-(X-X0)^2/(2*SIGMAX^2) - (Y-Y0)^2/(2*SIGMAY^2)]
!
! The wide case (SIGMA ~ 2.4 pixels) is a regression test: the initial
! amplitude was overestimated by a factor (2.5*SIGMA)**2/15 ~ 2.4, and
! DOWNHILL converged to a spurious solution with SIGMAX ~ 0.4.
!------------------------------------------------------------------------------
        PROGRAM TEST_GAUSS2DFIT
        USE TEST_UTILS
        USE Dynamic_Array_IMAGEN
        USE Dynamic_Array_IMAGEN_
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
        INCLUDE 'interface_imagen_.inc'
!
        INTEGER NBOX,NHALF
        PARAMETER (NBOX=15,NHALF=7)              !box size and (NBOX-1)/2
        REAL YRMSTOL
        PARAMETER (YRMSTOL=1.E-5)                 !as in the main program
!------------------------------------------------------------------------------
        CALL Initialize_Dynamic_Array_IMAGEN
        CALL Initialize_Dynamic_Array_IMAGEN_
! wide profile, close to the star in afztUCM0070.fits that failed
        CALL CHECK_FIT('wide',450.0,1000.0,89.0,636.2,2.35,2.45)
! compact profile
        CALL CHECK_FIT('compact',100.0,5000.0,50.3,40.6,1.2,1.1)
! very wide profile (wings partially outside the box)
        CALL CHECK_FIT('very wide',0.0,200.0,30.4,30.1,3.5,3.2)
!
        CALL Deallocate_Array_IMAGEN
        CALL Deallocate_Array_IMAGEN_
        CALL TEST_SUMMARY()
!
        CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE CHECK_FIT(LABEL,CTE0,AMP0,X00,Y00,SX0,SY0)
!
! Store in buffer 1 a 2D gaussian with the given parameters (BETA=0) in a
! NBOX x NBOX box centered at the nearest pixel to (X00,Y00), fit it with
! GAUSS2DFIT and check the fitted parameters.
!------------------------------------------------------------------------------
        SUBROUTINE CHECK_FIT(LABEL,CTE0,AMP0,X00,Y00,SX0,SY0)
        CHARACTER*(*), INTENT(IN) :: LABEL
        REAL, INTENT(IN) :: CTE0,AMP0,X00,Y00,SX0,SY0
!
        INTEGER I,J,IX1,IX2,IY1,IY2
        REAL X0,Y0,SIGMAX,SIGMAY,BETA,AMP,CTE
        REAL EX0,EY0,ESIGMAX,ESIGMAY,EBETA,EAMP,ECTE
        REAL EEX0,EEY0,EESIGMAX,EESIGMAY,EEBETA,EEAMP,EECTE
        REAL X0_,Y0_,SIGMAX_,SIGMAY_,AMPX_,AMPY_,CTEX_,CTEY_
        REAL AREA1,AREA2,AREA3,EAREA1,EAREA2,EAREA3
!
        IX1=NINT(X00)-NHALF
        IX2=IX1+NBOX-1
        IY1=NINT(Y00)-NHALF
        IY2=IY1+NBOX-1
        DO I=IY1,IY2
          DO J=IX1,IX2
            IMAGEN(J,I,1)=CTE0+AMP0*EXP(-(REAL(J)-X00)**2/(2.0*SX0**2)- &
             (REAL(I)-Y00)**2/(2.0*SY0**2))
          END DO
        END DO
        CALL GAUSS2DFIT(1,IX1,IX2,IY1,IY2,.FALSE., &
         X0,Y0,SIGMAX,SIGMAY,BETA,AMP,CTE, &
         EX0,EY0,ESIGMAX,ESIGMAY,EBETA,EAMP,ECTE, &
         EEX0,EEY0,EESIGMAX,EESIGMAY,EEBETA,EEAMP,EECTE, &
         X0_,Y0_,SIGMAX_,SIGMAY_,AMPX_,AMPY_,CTEX_,CTEY_, &
         AREA1,AREA2,AREA3,EAREA1,EAREA2,EAREA3, &
         YRMSTOL,0,0,.FALSE.)
        CALL CHECK_REAL(LABEL//': X0',X0,X00,0.0,0.01)
        CALL CHECK_REAL(LABEL//': Y0',Y0,Y00,0.0,0.01)
        CALL CHECK_REAL(LABEL//': SIGMAX',ABS(SIGMAX),SX0,0.01,0.0)
        CALL CHECK_REAL(LABEL//': SIGMAY',ABS(SIGMAY),SY0,0.01,0.0)
        CALL CHECK_REAL(LABEL//': AMP',AMP,AMP0,0.01,0.0)
        CALL CHECK_REAL(LABEL//': CTE',CTE,CTE0,0.0,0.01*AMP0)
        CALL CHECK_REAL(LABEL//': BETA',BETA,0.0,0.0,0.01)
        END SUBROUTINE CHECK_FIT
        END PROGRAM TEST_GAUSS2DFIT
