!------------------------------------------------------------------------------
! Unit tests for FPOLY, FDERI and FPOLYINV.
! Polynomial convention: COEFF(1) is the constant term, COEFF(K+1) multiplies
! X**K.
!------------------------------------------------------------------------------
        PROGRAM TEST_FPOLY
        USE TEST_UTILS
        IMPLICIT NONE
!
        REAL FPOLY,FDERI
        REAL COEFF(3)
        REAL X0,XEXP
!------------------------------------------------------------------------------
! p(x) = 1 - 2x + 3x^2
        COEFF=(/1.0,-2.0,3.0/)
        CALL CHECK_REAL('FPOLY p(0)',FPOLY(2,COEFF,0.0),1.0,1.E-6,0.0)
        CALL CHECK_REAL('FPOLY p(2)',FPOLY(2,COEFF,2.0),9.0,1.E-6,0.0)
        CALL CHECK_REAL('FPOLY p(-1.5)',FPOLY(2,COEFF,-1.5),10.75,1.E-6,0.0)
! degree 0: constant polynomial
        CALL CHECK_REAL('FPOLY degree 0',FPOLY(0,(/4.5/),123.0),4.5,1.E-6,0.0)
!------------------------------------------------------------------------------
! p'(x) = -2 + 6x
        CALL CHECK_REAL('FDERI p''(2)',FDERI(2,COEFF,2.0),10.0,1.E-6,0.0)
        CALL CHECK_REAL('FDERI p''(-1.5)',FDERI(2,COEFF,-1.5),-11.0,1.E-6,0.0)
        CALL CHECK_REAL('FDERI degree 1',FDERI(1,(/3.0,7.0/),5.0),7.0,1.E-6,0.0)
        CALL CHECK_REAL('FDERI degree 0',FDERI(0,(/3.0/),5.0),0.0,0.0,1.E-6)
!------------------------------------------------------------------------------
! FPOLYINV uses Newton's method and stops when |x_{n+1}-x_n| <= 1E-4
        CALL FPOLYINV(2,(/-2.0,0.0,1.0/),0.0,1.0,X0)
        CALL CHECK_REAL('FPOLYINV x^2-2=0',X0,SQRT(2.0),0.0,1.E-4)
! 1 - 2x + 3x^2 = 7  ->  3x^2 - 2x - 6 = 0 (positive root)
        XEXP=(2.0+SQRT(76.0))/6.0
        CALL FPOLYINV(2,COEFF,7.0,1.0,X0)
        CALL CHECK_REAL('FPOLYINV p(x)=7',X0,XEXP,0.0,1.E-4)
!
        CALL TEST_SUMMARY()
        END PROGRAM TEST_FPOLY
