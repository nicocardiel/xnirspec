!------------------------------------------------------------------------------
! Unit tests for CUBSPL and CUBSPLX.
!------------------------------------------------------------------------------
        PROGRAM TEST_SPLINE
        USE TEST_UTILS
        IMPLICIT NONE
!
        INTEGER N
        PARAMETER (N=5)
!
        INTEGER I,I0,IMODE
        REAL X(N),Y(N),S(N),A(N),B(N),C(N)
        REAL X0,Y0
        REAL XT(5)
        CHARACTER*60 CNAME
!------------------------------------------------------------------------------
        X=(/0.0,1.0,2.5,4.0,5.0/)               !non equally-spaced abscissae
        XT=(/0.5,1.7,3.3,4.9,2.5/)
!------------------------------------------------------------------------------
! a natural spline (IMODE=1) reproduces a straight line exactly, also when
! extrapolating with the first/last cubic
        Y=3.0*X-1.0
        S=0.0
        CALL CUBSPL(X,Y,N,1,S,A,B,C)
        I0=1
        DO I=1,5
          X0=XT(I)
          CALL CUBSPLX(X,Y,A,B,C,N,I0,X0,Y0)
          WRITE(CNAME,'(A,F4.1)') 'CUBSPL IMODE=1 line at x=',X0
          CALL CHECK_REAL(TRIM(CNAME),Y0,3.0*X0-1.0,1.E-5,1.E-5)
        END DO
        CALL CUBSPLX(X,Y,A,B,C,N,I0,-1.0,Y0)
        CALL CHECK_REAL('CUBSPL IMODE=1 line extrapolated x=-1',Y0,-4.0, &
         1.E-5,1.E-5)
        CALL CUBSPLX(X,Y,A,B,C,N,I0,6.0,Y0)
        CALL CHECK_REAL('CUBSPL IMODE=1 line extrapolated x=6',Y0,17.0, &
         1.E-5,1.E-5)
!------------------------------------------------------------------------------
! for every end condition, the spline passes through the data points
        Y=SIN(X)
        DO IMODE=1,4
          S=0.0
          IF(IMODE.EQ.4)THEN                !first derivatives at both ends
            S(1)=COS(X(1))
            S(N)=COS(X(N))
          END IF
          CALL CUBSPL(X,Y,N,IMODE,S,A,B,C)
          I0=1
          DO I=1,N
            CALL CUBSPLX(X,Y,A,B,C,N,I0,X(I),Y0)
            WRITE(CNAME,'(A,I0,A,I0,A)') 'CUBSPL IMODE=',IMODE, &
             ' through node ',I
            CALL CHECK_REAL(TRIM(CNAME),Y0,Y(I),1.E-5,1.E-5)
          END DO
        END DO
!------------------------------------------------------------------------------
! a clamped spline (IMODE=4, exact end slopes) reproduces a cubic exactly:
! y = x^3 - 2x, y' = 3x^2 - 2
        Y=X**3-2.0*X
        S=0.0
        S(1)=3.0*X(1)**2-2.0
        S(N)=3.0*X(N)**2-2.0
        CALL CUBSPL(X,Y,N,4,S,A,B,C)
        I0=1
        DO I=1,4
          X0=XT(I)
          CALL CUBSPLX(X,Y,A,B,C,N,I0,X0,Y0)
          WRITE(CNAME,'(A,F4.1)') 'CUBSPL IMODE=4 cubic at x=',X0
          CALL CHECK_REAL(TRIM(CNAME),Y0,X0**3-2.0*X0,1.E-4,1.E-4)
        END DO
!
        CALL TEST_SUMMARY()
        END PROGRAM TEST_SPLINE
