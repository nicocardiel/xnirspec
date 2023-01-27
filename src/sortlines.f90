! Ordena las lineas espectrales que definen el boundary
!
        SUBROUTINE SORTLINES
        IMPLICIT NONE
!
        INCLUDE 'largest.inc'
!
        INTEGER L,K
        INTEGER NLINBL,NLINBA
        INTEGER NDEGBA(NXYMAX),NDEGBL(NXYMAX)
        INTEGER NDEGBA00,NDEGBL00
        INTEGER NDEGB_(NXYMAX)
        INTEGER LSORT(NXYMAX)
        REAL COEFFBL(20,NXYMAX),COEFFBA(20,NXYMAX)
        REAL COEFFB_(20,NXYMAX)
        REAL XMINBL(NXYMAX),XMAXBL(NXYMAX)
        REAL YMINBL(NXYMAX),YMAXBL(NXYMAX)
        REAL XMINBA(NXYMAX),XMAXBA(NXYMAX)
        REAL YMINBA(NXYMAX),YMAXBA(NXYMAX)
        REAL XMINB_(NXYMAX),XMAXB_(NXYMAX)
        REAL YMINB_(NXYMAX),YMAXB_(NXYMAX)
        REAL FSORT(NXYMAX)
!
        COMMON/BLKBOUND1/NLINBL,NDEGBL,NDEGBL00
        COMMON/BLKBOUND1B/NLINBA,NDEGBA,NDEGBA00
        COMMON/BLKBOUND2/COEFFBL
        COMMON/BLKBOUND2B/COEFFBA
        COMMON/BLKBOUND4/XMINBL,YMINBL,XMAXBL,YMAXBL
        COMMON/BLKBOUND5/XMINBA,YMINBA,XMAXBA,YMAXBA
!------------------------------------------------------------------------------
! ordenamos BL usando YMAXBL
        IF(NLINBL.GT.0)THEN
          DO L=1,NLINBL
            FSORT(L)=YMAXBL(L)
            LSORT(L)=L
          END DO
!
          CALL ORDENA1F1I(NLINBL,FSORT,LSORT)
!
          DO L=1,NLINBL
            NDEGB_(L)=NDEGBL(LSORT(L))
            DO K=1,NDEGB_(L)+1
              COEFFB_(K,L)=COEFFBL(K,LSORT(L))
            END DO
            XMINB_(L)=XMINBL(LSORT(L))
            XMAXB_(L)=XMAXBL(LSORT(L))
            YMINB_(L)=YMINBL(LSORT(L))
            YMAXB_(L)=YMAXBL(LSORT(L))
          END DO
!
          DO L=1,NLINBL
            NDEGBL(L)=NDEGB_(L)
            DO K=1,NDEGBL(L)+1
              COEFFBL(K,L)=COEFFB_(K,L)
            END DO
            XMINBL(L)=XMINB_(L)
            XMAXBL(L)=XMAXB_(L)
            YMINBL(L)=YMINB_(L)
            YMAXBL(L)=YMAXB_(L)
          END DO
!
          NDEGBL00=0
          DO L=1,NLINBL
            IF(NDEGBL(L).GT.NDEGBL00) NDEGBL00=NDEGBL(L)
          END DO
        END IF
!------------------------------------------------------------------------------
! si tenemos intersecciones, calculamos limites
        IF((NLINBA.GT.0).AND.(NLINBL.GT.0))THEN
          DO L=1,NLINBA
            CALL INTERSEC(NDEGBA(L),COEFFBA(1,L),NDEGBL(NLINBL),COEFFBL(1,NLINBL), &
             0.5*(XMINBL(NLINBL)+XMAXBL(NLINBL)),XMAXBA(L),YMAXBA(L))
            CALL INTERSEC(NDEGBA(L),COEFFBA(1,L),NDEGBL(1),COEFFBL(1,1), &
             0.5*(XMINBL(1)+XMAXBL(1)),XMINBA(L),YMINBA(L))
          END DO
        END IF
!------------------------------------------------------------------------------
! ordenamos BA usando XMAXBA
        IF(NLINBA.GT.0)THEN
          DO L=1,NLINBA
            FSORT(L)=XMAXBA(L)
            LSORT(L)=L
          END DO
!
          CALL ORDENA1F1I(NLINBA,FSORT,LSORT)
!
          DO L=1,NLINBA
            NDEGB_(L)=NDEGBA(LSORT(L))
            DO K=1,NDEGB_(L)+1
              COEFFB_(K,L)=COEFFBA(K,LSORT(L))
            END DO
            XMINB_(L)=XMINBA(LSORT(L))
            XMAXB_(L)=XMAXBA(LSORT(L))
            YMINB_(L)=YMINBA(LSORT(L))
            YMAXB_(L)=YMAXBA(LSORT(L))
          END DO
!
          DO L=1,NLINBA
            NDEGBA(L)=NDEGB_(L)
            DO K=1,NDEGBA(L)+1
              COEFFBA(K,L)=COEFFB_(K,L)
            END DO
            XMINBA(L)=XMINB_(L)
            XMAXBA(L)=XMAXB_(L)
            YMINBA(L)=YMINB_(L)
            YMAXBA(L)=YMAXB_(L)
          END DO
!
          NDEGBA00=0
          DO L=1,NLINBA
            IF(NDEGBA(L).GT.NDEGBA00) NDEGBA00=NDEGBA(L)
          END DO
        END IF
!------------------------------------------------------------------------------
        END
