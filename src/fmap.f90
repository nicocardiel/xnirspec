! Rutina que calcula la transformacion segun el mapping realizado
!
        SUBROUTINE FMAP(N,AIJ,BIJ,X,Y,U,V)
        INTEGER N
        REAL AIJ((N+1)*N)
        REAL BIJ((N+1)*N)
        REAL X,Y,U,V
!
        INTEGER II,L,M
        REAL FSCALE
        REAL XS,YS
!
        COMMON/BLKFSCALE/FSCALE
!------------------------------------------------------------------------------
        XS=X*FSCALE
        YS=Y*FSCALE
!
        II=0
        U=0.
        V=0.
        DO L=0,N
          IF(L.EQ.0)THEN
            DO M=0,N
              II=II+1
              IF(M.EQ.0)THEN
                U=U+AIJ(II)
                V=V+BIJ(II)
              ELSE
                U=U+AIJ(II)*(YS**M)
                V=V+BIJ(II)*(YS**M)
              END IF
            END DO
          ELSE
            DO M=0,N-L
              II=II+1
              IF(M.EQ.0)THEN
                U=U+AIJ(II)*(XS**L)
                V=V+BIJ(II)*(XS**L)
              ELSE
                U=U+AIJ(II)*(XS**L)*(YS**M)
                V=V+BIJ(II)*(XS**L)*(YS**M)
              END IF
            END DO
          END IF
        END DO
!
        U=U/FSCALE
        V=V/FSCALE
!
        END
