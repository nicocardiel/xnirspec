!------------------------------------------------------------------------------
! Unit tests for FMEAN0, FMEDIAN1 and FROBUSTRMS0.
!------------------------------------------------------------------------------
        PROGRAM TEST_STATS
        USE TEST_UTILS
        IMPLICIT NONE
!
        REAL FMEAN0,FMEDIAN1,FROBUSTRMS0
        REAL X8(8),X3(3),X4(4),X1(1)
        REAL FMEAN,SIGMA
!------------------------------------------------------------------------------
! mean and r.m.s. (normalised by N-1): sum of squared deviations = 32
        X8=(/2.0,4.0,4.0,4.0,5.0,5.0,7.0,9.0/)
        FMEAN=FMEAN0(8,X8,SIGMA)
        CALL CHECK_REAL('FMEAN0 mean',FMEAN,5.0,1.E-6,0.0)
        CALL CHECK_REAL('FMEAN0 sigma (N-1)',SIGMA,SQRT(32.0/7.0),1.E-6,0.0)
        X1=(/3.5/)
        FMEAN=FMEAN0(1,X1,SIGMA)
        CALL CHECK_REAL('FMEAN0 N=1 mean',FMEAN,3.5,1.E-6,0.0)
        CALL CHECK_REAL('FMEAN0 N=1 sigma',SIGMA,0.0,0.0,1.E-6)
!------------------------------------------------------------------------------
! median (the input array is returned sorted)
        X3=(/3.0,1.0,2.0/)
        CALL CHECK_REAL('FMEDIAN1 odd N',FMEDIAN1(3,X3),2.0,1.E-6,0.0)
        X4=(/4.0,1.0,3.0,2.0/)
        CALL CHECK_REAL('FMEDIAN1 even N',FMEDIAN1(4,X4),2.5,1.E-6,0.0)
        CALL CHECK_TRUE('FMEDIAN1 returns sorted array', &
         ALL(X4.EQ.(/1.0,2.0,3.0,4.0/)))
        X1=(/-7.0/)
        CALL CHECK_REAL('FMEDIAN1 N=1',FMEDIAN1(1,X1),-7.0,1.E-6,0.0)
!------------------------------------------------------------------------------
! robust r.m.s. = 0.7413*(X(I75)-X(I25)), nearest-rank percentiles:
! N=8 -> I25=INT(2.5)=2, I75=INT(6.5)=6 -> 0.7413*(6-2)
        X8=(/5.0,1.0,8.0,3.0,7.0,2.0,6.0,4.0/)
        CALL CHECK_REAL('FROBUSTRMS0 N=8',FROBUSTRMS0(8,X8),0.7413*4.0, &
         1.E-6,0.0)
        X1=(/2.0/)
        CALL CHECK_REAL('FROBUSTRMS0 N=1',FROBUSTRMS0(1,X1),0.0,0.0,1.E-6)
!
        CALL TEST_SUMMARY()
        END PROGRAM TEST_STATS
