!------------------------------------------------------------------------------
! Weighted mean
!------------------------------------------------------------------------------
!omment
!
! REAL FUNCTION FMEAN0W(N,X,EX,SIGMA)
!
! Input: N,X,EX
! Output: FMEAN0W (function),SIGMA
!
! Calculate the mean value of X(N) and its r.m.s.
!
! INTEGER N -> no. of elements
! REAL    X(N) -> input matrix
! REAL    EX(N) -> input errors
! REAL SIGMA -> expected r.m.s. around the mean value
!
!omment
!------------------------------------------------------------------------------
        REAL FUNCTION FMEAN0W(N,X,EX,SIGMA)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL EX(N)
        REAL SIGMA
!
        REAL FMEAN0
!
        INTEGER I
        REAL SUM,SUME,SUMW
        LOGICAL LNULL
!------------------------------------------------------------------------------
        IF(N.LE.0) STOP 'FATAL ERROR: in function FMEAN0: N.LE.0'
!------------------------------------------------------------------------------
! chequeamos que no hay errores nulos
        LNULL=.FALSE.
        DO I=1,N
          IF(EX(I).LE.0.0) LNULL=.TRUE.
        END DO
        IF(LNULL)THEN
!!!          WRITE(*,100) 'WARNING: negative o null weights in FMEAN0W.'
!!!          WRITE(*,101) ' Using FMEAN0 instead.'
          FMEAN0W=FMEAN0(N,X,SIGMA)
          RETURN
        END IF
!------------------------------------------------------------------------------
        SUM=0.
        SUMW=0.
        DO I=1,N
          SUME=EX(I)*EX(I)
          SUM=SUM+X(I)/SUME
          SUMW=SUMW+1./SUME
        END DO
        FMEAN0W=SUM/SUMW
        SIGMA=SQRT(REAL(N)/SUMW)
!------------------------------------------------------------------------------
!!!100     FORMAT(A,$)
!!!101     FORMAT(A)
        END
