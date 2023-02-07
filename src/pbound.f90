! Dibuja la calibracion actual de la distorsion del espectro (boundary)
        SUBROUTINE PBOUND
        USE Dynamic_Array_COEFFBL
        USE Dynamic_Array_COEFFBA
        IMPLICIT NONE
        INCLUDE 'interface_coeffbl.inc'
        INCLUDE 'interface_coeffba.inc'
!
        INCLUDE 'largest.inc'
!
        REAL FPOLY
!
        INTEGER I,L,LL
        INTEGER NDEGBL(NXYMAX),NDEGBA(NXYMAX)  !polynomial degrees for boundary
        INTEGER NDEGBL00,NDEGBA00
        INTEGER NLINBL,NLINBA                  !number of arc lines in boundary
!delete REAL COEFFBL(20,NXYMAX),COEFFBA(20,NXYMAX)      !pol. coef. in boundary
        REAL XP(NXYMAX),YP(NXYMAX)
        REAL XMIN,XMAX,YMIN,YMAX
        REAL X0,Y0
        REAL XMINBL(NXYMAX),XMAXBL(NXYMAX)
        REAL YMINBL(NXYMAX),YMAXBL(NXYMAX)
        REAL XMINBA(NXYMAX),XMAXBA(NXYMAX)
        REAL YMINBA(NXYMAX),YMAXBA(NXYMAX)
        LOGICAL LBOUNDL,LBOUNDA
!
        COMMON/BLKBOUND1/NLINBL,NDEGBL,NDEGBL00
        COMMON/BLKBOUND1B/NLINBA,NDEGBA,NDEGBA00
!delete COMMON/BLKBOUND2/COEFFBL
!delete COMMON/BLKBOUND2B/COEFFBA
        COMMON/BLKBOUND3/LBOUNDL,LBOUNDA
        COMMON/BLKBOUND4/XMINBL,YMINBL,XMAXBL,YMAXBL
        COMMON/BLKBOUND5/XMINBA,YMINBA,XMAXBA,YMAXBA
!------------------------------------------------------------------------------
        CALL PGBBUF
        IF(LBOUNDL)THEN
          CALL PGSCI(4)
          DO L=1,NLINBL
            XMIN=XMINBL(L)
            XMAX=XMAXBL(L)
            DO I=1,NXYMAX
              XP(I)=REAL(I-1)/REAL(NXYMAX-1)*(XMAX-XMIN)+XMIN
              YP(I)=FPOLY(NDEGBL(L),COEFFBL(1,L),XP(I))
            END DO
            CALL PGLINE(NXYMAX,XP,YP)
          END DO
        END IF
!
        IF(LBOUNDA)THEN
          CALL PGSCI(4)
          DO L=1,NLINBA
            YMIN=YMINBA(L)
            YMAX=YMAXBA(L)
            DO I=1,NXYMAX
              YP(I)=REAL(I-1)/REAL(NXYMAX-1)*(YMAX-YMIN)+YMIN
              XP(I)=FPOLY(NDEGBA(L),COEFFBA(1,L),YP(I))
            END DO
            CALL PGLINE(NXYMAX,XP,YP)
          END DO
        END IF
!
        IF(LBOUNDL.AND.LBOUNDA)THEN
          CALL PGSCI(6)
          DO L=1,NLINBA
            DO LL=1,NLINBL
              CALL INTERSEC(NDEGBA(L),COEFFBA(1,L),NDEGBL(LL),COEFFBL(1,LL),0.5*(XMINBA(L)+XMAXBA(L)),X0,Y0)
              CALL PGPOINT(1,X0,Y0,22)
            END DO
          END DO
        END IF
!
        CALL PGSCI(1)
        CALL PGEBUF
!------------------------------------------------------------------------------
        END
