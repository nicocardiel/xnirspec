!------------------------------------------------------------------------------
! Unit tests for ORDENA1F, ORDENA1F1I and BINSEARCH.
!------------------------------------------------------------------------------
        PROGRAM TEST_SORT
        USE TEST_UTILS
        IMPLICIT NONE
!
        INTEGER I,N0
        INTEGER IB(6)
        REAL A(6),A0(6),A1(1)
        REAL X(5)
        LOGICAL LOK
!------------------------------------------------------------------------------
! ORDENA1F: unsorted (with a tie), already sorted, reversed, N=1
        A=(/5.0,3.0,9.0,1.0,3.0,7.0/)
        CALL ORDENA1F(6,A)
        CALL CHECK_TRUE('ORDENA1F unsorted with tie', &
         ALL(A.EQ.(/1.0,3.0,3.0,5.0,7.0,9.0/)))
        CALL ORDENA1F(6,A)
        CALL CHECK_TRUE('ORDENA1F already sorted', &
         ALL(A.EQ.(/1.0,3.0,3.0,5.0,7.0,9.0/)))
        A=(/6.0,5.0,4.0,3.0,2.0,1.0/)
        CALL ORDENA1F(6,A)
        CALL CHECK_TRUE('ORDENA1F reversed', &
         ALL(A.EQ.(/1.0,2.0,3.0,4.0,5.0,6.0/)))
        A1=(/2.0/)
        CALL ORDENA1F(1,A1)
        CALL CHECK_REAL('ORDENA1F N=1',A1(1),2.0,0.0,0.0)
!------------------------------------------------------------------------------
! ORDENA1F1I: the integer array (here, the original indices) must be
! permuted in parallel, i.e. A(I) = A0(IB(I)) for every I
        A0=(/5.0,3.0,9.0,1.0,3.0,7.0/)
        A=A0
        IB=(/(I,I=1,6)/)
        CALL ORDENA1F1I(6,A,IB)
        CALL CHECK_TRUE('ORDENA1F1I sorts real array', &
         ALL(A.EQ.(/1.0,3.0,3.0,5.0,7.0,9.0/)))
        LOK=.TRUE.
        DO I=1,6
          IF(A(I).NE.A0(IB(I))) LOK=.FALSE.
        END DO
        CALL CHECK_TRUE('ORDENA1F1I permutes integer array in parallel',LOK)
!------------------------------------------------------------------------------
! BINSEARCH: N0 such that X(N0) <= X0 < X(N0+1); 0 below X(1), N above X(N)
        X=(/1.0,2.0,3.0,4.0,5.0/)
        N0=1
        CALL BINSEARCH(X,5,2.5,N0)
        CALL CHECK_INT('BINSEARCH interior (start N0=1)',N0,2)
        N0=5
        CALL BINSEARCH(X,5,2.5,N0)
        CALL CHECK_INT('BINSEARCH interior (start N0=5)',N0,2)
        N0=3
        CALL BINSEARCH(X,5,3.0,N0)
        CALL CHECK_INT('BINSEARCH X0 equal to X(3)',N0,3)
        N0=3
        CALL BINSEARCH(X,5,1.0,N0)
        CALL CHECK_INT('BINSEARCH X0 equal to X(1)',N0,1)
        N0=3
        CALL BINSEARCH(X,5,0.5,N0)
        CALL CHECK_INT('BINSEARCH below X(1)',N0,0)
        N0=3
        CALL BINSEARCH(X,5,5.5,N0)
        CALL CHECK_INT('BINSEARCH above X(N)',N0,5)
!
        CALL TEST_SUMMARY()
        END PROGRAM TEST_SORT
