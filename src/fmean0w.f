C------------------------------------------------------------------------------
C Weighted mean
C------------------------------------------------------------------------------
Comment
C
C REAL FUNCTION FMEAN0W(N,X,EX,SIGMA)
C
C Input: N,X,EX
C Output: FMEAN0W (function),SIGMA
C
C Calculate the mean value of X(N) and its r.m.s.
C
C INTEGER N -> no. of elements
C REAL    X(N) -> input matrix
C REAL    EX(N) -> input errors
C REAL SIGMA -> expected r.m.s. around the mean value
C
Comment
C------------------------------------------------------------------------------
        REAL FUNCTION FMEAN0W(N,X,EX,SIGMA)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL EX(N)
        REAL SIGMA
C
        REAL FMEAN0
C
        INTEGER I
        REAL SUM,SUME,SUMW
        LOGICAL LNULL
C------------------------------------------------------------------------------
        IF(N.LE.0) STOP 'FATAL ERROR: in function FMEAN0: N.LE.0'
C------------------------------------------------------------------------------
C chequeamos que no hay errores nulos
        LNULL=.FALSE.
        DO I=1,N
          IF(EX(I).LE.0.0) LNULL=.TRUE.
        END DO
        IF(LNULL)THEN
ccc          WRITE(*,100) 'WARNING: negative o null weights in FMEAN0W.'
ccc          WRITE(*,101) ' Using FMEAN0 instead.'
          FMEAN0W=FMEAN0(N,X,SIGMA)
          RETURN
        END IF
C------------------------------------------------------------------------------
        SUM=0.
        SUMW=0.
        DO I=1,N
          SUME=EX(I)*EX(I)
          SUM=SUM+X(I)/SUME
          SUMW=SUMW+1./SUME
        END DO
        FMEAN0W=SUM/SUMW
        SIGMA=SQRT(REAL(N)/SUMW)
C------------------------------------------------------------------------------
ccc100     FORMAT(A,$)
ccc101     FORMAT(A)
        END
