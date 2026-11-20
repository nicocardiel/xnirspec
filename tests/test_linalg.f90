!------------------------------------------------------------------------------
! Unit tests for LUDCMP/LUSOLV and FITPOL.
! Singular matrices are not tested: LUDCMP then waits for <CR> in stdin.
!------------------------------------------------------------------------------
        PROGRAM TEST_LINALG
        USE TEST_UTILS
        IMPLICIT NONE
!
        INTEGER NDIM
        PARAMETER (NDIM=5)                !physical dimension larger than N=3
        INTEGER NDEGMAX
        PARAMETER (NDEGMAX=19)                         !as defined in FITPOL
!
        INTEGER I,K
        INTEGER ORDER(3),IOK,IPAR
        REAL A(NDIM,NDIM),SCALEROW(3),B(3),X(3)
        REAL XSOL1(3),XSOL2(3)
        REAL XF(10),YF(10)
        LOGICAL IFCOEF(NDEGMAX+1)
        REAL BF(NDEGMAX+1),VARBF(NDEGMAX+1)
        REAL CHISQR,SR2
        REAL COVAR(NDEGMAX+1,NDEGMAX+1)
        CHARACTER*40 CNAME
!------------------------------------------------------------------------------
! A(I,J): row I, column J
!   2x +  y -  z =   8
!  -3x -  y + 2z = -11
!  -2x +  y + 2z =  -3   ->  (x,y,z) = (2,3,-1)
        A=0.0
        A(1,1:3)=(/ 2.0, 1.0,-1.0/)
        A(2,1:3)=(/-3.0,-1.0, 2.0/)
        A(3,1:3)=(/-2.0, 1.0, 2.0/)
        XSOL1=(/2.0,3.0,-1.0/)
        XSOL2=(/1.0,-1.0,2.0/)
        CALL LUDCMP(A,3,NDIM,ORDER,SCALEROW,IOK,IPAR)
        CALL CHECK_INT('LUDCMP IOK',IOK,0)
        B=(/8.0,-11.0,-3.0/)
        CALL LUSOLV(A,3,NDIM,ORDER,SCALEROW,B,X)
        DO I=1,3
          WRITE(CNAME,'(A,I0,A)') 'LUSOLV x(',I,')'
          CALL CHECK_REAL(TRIM(CNAME),X(I),XSOL1(I),1.E-5,1.E-5)
        END DO
! the LU decomposition can be reused with another right-hand side:
! solution (1,-1,2) -> B = (-1,2,1)
        B=(/-1.0,2.0,1.0/)
        CALL LUSOLV(A,3,NDIM,ORDER,SCALEROW,B,X)
        DO I=1,3
          WRITE(CNAME,'(A,I0,A)') 'LUSOLV reused, x(',I,')'
          CALL CHECK_REAL(TRIM(CNAME),X(I),XSOL2(I),1.E-5,1.E-5)
        END DO
!------------------------------------------------------------------------------
! FITPOL with noiseless data from y = 1 - 2x + 0.5x^2, x = 0,...,9
        DO K=1,10
          XF(K)=REAL(K-1)
          YF(K)=1.0-2.0*XF(K)+0.5*XF(K)*XF(K)
        END DO
        IFCOEF=.TRUE.
        BF=0.0
        CALL FITPOL(10,XF,YF,2,IFCOEF,BF,VARBF,CHISQR,SR2,COVAR)
        CALL CHECK_REAL('FITPOL b0',BF(1),1.0,0.0,1.E-4)
        CALL CHECK_REAL('FITPOL b1',BF(2),-2.0,0.0,1.E-4)
        CALL CHECK_REAL('FITPOL b2',BF(3),0.5,0.0,1.E-4)
! fixing the constant term (IFCOEF(1)=.FALSE.): it keeps its input value
! and the remaining coefficients are fitted
        IFCOEF=.TRUE.
        IFCOEF(1)=.FALSE.
        BF=0.0
        BF(1)=1.0
        CALL FITPOL(10,XF,YF,2,IFCOEF,BF,VARBF,CHISQR,SR2,COVAR)
        CALL CHECK_REAL('FITPOL fixed b0 unchanged',BF(1),1.0,0.0,0.0)
        CALL CHECK_REAL('FITPOL fixed b0, b1',BF(2),-2.0,0.0,1.E-4)
        CALL CHECK_REAL('FITPOL fixed b0, b2',BF(3),0.5,0.0,1.E-4)
        CALL CHECK_REAL('FITPOL fixed b0, VARB(1)=0',VARBF(1),0.0,0.0,0.0)
!
        CALL TEST_SUMMARY()
        END PROGRAM TEST_LINALG
