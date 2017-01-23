C------------------------------------------------------------------------------
C Version 16-June-1998                                         File: fmedian1.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This routine is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C REAL FUNCTION FMEDIAN1E(N,X,XE,EFMEDIAN)
C
C Input: N,X,XE
C Output: FMEDIAN1E (function), X (sorted), EFMEDIAN
C
C Calculate the median value of X(N). It is important to note that this 
C subroutine rearranges the matrix X which is returned sorted. This routine
C also computes the error of FMEDIAN1E using Monte Carlo simulations.
C
C INTEGER N -> no. of elements
C REAL    X(N) -> input matrix
C REAL    EX(N) -> input error matrix
C REAL    EFMEDIAN -> error of FMEDIAN1E
C
Comment
C------------------------------------------------------------------------------
        REAL FUNCTION FMEDIAN1E(N,X,EX,EFMEDIAN)
        IMPLICIT NONE
        INTEGER N
        REAL X(N),EX(N)
        REAL EFMEDIAN
C
        INCLUDE 'largest.inc'
C
        INTEGER NMAXSIMUL
        PARAMETER (NMAXSIMUL=100)
C
        REAL FMEAN0
C variables locales
        INTEGER NN
        INTEGER I,K
        INTEGER NSEED
        REAL X_(NXYMAX)
        REAL FMEDIAN_(NMAXSIMUL)
        REAL XDUM
        REAL RANDOMNUMBER,R1,R2
C------------------------------------------------------------------------------
        IF(N.EQ.0) STOP 'FATAL ERROR: in function FMEDIAN: N=0.'
        CALL ORDENA1F(N,X)
        NN=N/2
        IF(MOD(N,2).EQ.0)THEN
          FMEDIAN1E=(X(NN)+X(NN+1))/2.
        ELSE
          FMEDIAN1E=X(NN+1)
        END IF
C
        NSEED=-1
        DO K=1,NMAXSIMUL
          DO I=1,N
            R1=RANDOMNUMBER(NSEED)
            R2=RANDOMNUMBER(NSEED)*2.*3.141592654
            X_(I)=X(I)+1.414213562*EX(I)*SQRT(-ALOG(1.-R1))*COS(R2)
          END DO
          CALL ORDENA1F(N,X_)
          IF(MOD(N,2).EQ.0)THEN
            FMEDIAN_(K)=(X_(NN)+X_(NN+1))/2.
          ELSE
            FMEDIAN_(K)=X_(NN+1)
          END IF
        END DO
        XDUM=FMEAN0(NMAXSIMUL,FMEDIAN_,EFMEDIAN)
C
        END
